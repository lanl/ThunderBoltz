"""The main ThunderBoltz wrapper object with relevant helper functions."""
import copy
from datetime import datetime
import io
import json
import multiprocessing as mp
import os
from os import listdir as ls
from os.path import join as pjoin
from pathlib import Path
import re
import shutil
import subprocess
from subprocess import PIPE
import sys
import time
from typing import (
    Any,
    Tuple,
)
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
from scipy.stats import linregress
from scipy.integrate import quadrature as quad

import thunderboltz as tb
from thunderboltz.input import CrossSections
from thunderboltz import parameters
from thunderboltz.parameters import tb_parameters
from thunderboltz.parameters import wrap_parameters
from thunderboltz import parsing
from thunderboltz.parallel import MPRunner
from thunderboltz.parallel import visit
from thunderboltz.plotting.figure_gen import plot_timeseries
from thunderboltz.plotting.figure_gen import plot_vdfs

# Physical Constants
ME = 9.1093837e-31 # kg
QE = 1.60217663e-19 # C or J/eV
AMU_TO_KG = 1.6605e-27 # kg / amu

# Elastic angular distribution function names and parameters
EADF = {"He_Park": ("Park", [0.283, 0.667, 0.0307, 16.971, 6.59, 29.02, 0.0258, 0.295, 0.00328, 1.794]),
        "default": (None, None)}

# Electron impact ionization energy sharing distribution function names and parameters
EESD = {"default": (None, None),
        "equal": ("Equal", None),
        "uniform": ("Uniform", None)}

# Parameters specifically intended for cross section set construction
CS_PARAMS = ["n", "analytic_cs", "proper_sample", "nsamples", "egen",
             "mix_thresh", "eadf", "ECS", "fixed_background"]


class ThunderBoltz(MPRunner):
    """ThunderBoltz 0D DSMC simulation wrapper.

    Args:
        directory (str): The path to a directory that will host
            ThunderBoltz compiled, source, input, and output files.
        cs (CrossSections): The set of cross section information
            required for this simulation. Optionally supplied as an alternative
            to ``indeck``, default is an empty CrossSections object.
        out_file (str): Optional file base name for ThunderBoltz stdout buffer,
            default is ``"thunderboltz"``.
        monitor (bool): Runtime flag, when set to ``True`` an empty ``monitor``
            file will be generated in the simulation directory.
            Deleting this file will cause the ThunderBoltz process to exit,
            but allow the wrapper to continue execution. This is useful
            performing several simulation calculations sequentially,
            but manual exit is required for each one, or if post processing
            is required immediately after ThunderBoltz exits.
        live (bool): Run and update time series plotting GUI during simulation.
        ts_plot_params (list[str]): The default output parameters to be plotted
            by :func:`plot_timeseries`.
        **params:
            ThunderBoltz simulation parameters. Any attributes of
            :class:`~.thunderboltz.parameters.TBParameters` or
            :class:`~.thunderboltz.parameters.WrapParameters` can be passed here.
    """
    #: Name of the simulation output file produced by ThunderBoltz
    logfile = "thunderboltz.log"
    #: Files to be read when tabular data is requested
    output_files = ["thunderboltz.out", "Particle_Type", "Counts.dat", logfile]
    #: Time series plot parameters
    ts_plot_params = ["MEe", "mobN", "a_n"]

    def __init__(self, directory=None, cs=None, out_file=None,
                 monitor=None, live=None, live_rate=None, ts_plot_params=None, **params):
        """Initialize the ThunderBoltz object."""
        # Initialize data structures
        #: Particle-specific times series data
        self.particle_tables = []
        self.kinetic_table = None #: Banner output data
        self.timeseries = None #: All tick-by-tick simulation data
        self.vdfs = None #: Particle velocity dump data
        #: Particle velocity data intended for particle initialization
        self.vdf_init_data = None

        #: Time step at which steady-state calculations are considered converged
        self.time_conv = None
        self.counts = None #: pd.DataFrame: Table of collision counts
        # Initialize runtime logging params
        self.elapsed_time = None #: Elapsed wall-clock time of calculation
        self.runtime_start = None #: Date/Time of calculation start
        self.runtime_end = None #: Date/Time of calculation end
        # Initialize runtime params
        self.out_file = out_file #: Name for ThunderBoltz stdout file

        # Plotting objects and settings
        #: (bool): Run and update time series plotting GUI during simulation.
        self.live = live
        #: (bool): Run and update reaction rate plotting GUI during simulation.
        self.live_rate = live_rate
        #: The figure object for the time series plot.
        self.ts_fig = None
        #: The figure object for the rate plot
        self.rate_fig = None
        #: List of functions that are called every time banner output is
        #: updated.
        self.callbacks = []

        self.directory = directory #: Simulation directory
        #: Recorded thunderboltz warnings read in from output files
        self.err_stack = []
        if ts_plot_params: self.ts_plot_params = ts_plot_params

        # The Popen subprocess that will execute ThunderBoltz
        self._proc = None
        #: (bool): Option to create temp file during run that causes safe exit upon deletion
        self.monitor = monitor

        # Grab hyper param defaults
        all_params = copy.deepcopy(wrap_parameters)
        # Update with default thunderboltz params
        all_params.update(copy.deepcopy(tb_parameters))
        # Update with user initialized thunderboltz params
        all_params.update(copy.deepcopy(params))
        # Remember initial init parameters
        self.initial_params = copy.deepcopy(params)
        self.cs = cs
        # Make empty CrossSections object
        if cs is None:
            self.cs = CrossSections()
        # Set input for the first time
        self.hp = {}
        self.tb_params = {}
        self.set_(**all_params)

    def set_(self, **p):
        """Update parameters, call appropriate functions ensuring
        input parameters are self-consistent.

        Args:
            **p: Optional keyword parameters to update the calculator.
                can be any of :class:`~.thunderboltz.parameters.TBParameters` or
                :class:`~.thunderboltz.parameters.WrapParameters`.
        """
        # Set thunderboltz and hyper params appropriately
        self.hp.update({k: v for k, v in p.items()
            if k in wrap_parameters})
        self.tb_params.update({k: v for k, v in p.items()
            if k in tb_parameters})
        # For now, also allow arbitrary parameters for log files
        self.hp.update({k: v for k, v in p.items()
            if k not in tb_parameters})
        # Update directory
        if "directory" in p:
            self.directory = p["directory"]

        # Determine if special parameters changed
        resample = CS_PARAMS + ["indeck"]
        if any(i in resample for i in p):
            self._load_indeck()

        # Adjust differential models if necessary
        if ("eesd" in p or "eadf" in p) and len(self.cs.table):
            self._update_differential_models()

        # No matter what, input files should be considered unwritten
        # and input re-constrained
        self.written_input = False
        self._constrain_input()

        # Interpret velocity init file request (must occur after NP is set)
        if self.hp["vdf_init"]: self._prepare_vdf_init()

    def reset(self):
        """Reset output data for a new run."""
        self.particle_tables = []
        self.kinetic_table = None
        self.timeseries = None
        self.vdfs = None
        self.vdf_init_data = None
        self.time_conv = None
        self.counts = None
        self.elapsed_time = None
        self.runtime_start = None
        self.runtime_end = None
        self.ts_fig = None
        self.rate_fig = None
        self.callbacks = []
        self.err_stack = []
        self._proc = None

    def get_sim_param(self, key):
        """Return the value of a simulation parameter.

        Raises:
            IndexError: if the parameter is not set.
        """
        if key in self.hp:
            return self.hp[key]
        if key in self.tb_params:
            return self.tb_params[key]

        # If not found, raise error
        raise IndexError(f"'{key}' not a set simulation parameter.")

    def read_tb_params(self, fname, ignore=[]):
        """Takes file name of an input deck, updates the simulation
        parameters and returns the simulation parameters which were
        read from the file ``fname``.

        Args:
            fname (str): The name of the indeck file to read.
            ignore (list[str]): Don't read certain ThunderBoltz
                params, e.g. ``["MP", "QP"]`` would ignore the
                mass and charge parameters in an indeck file.
        Returns:
            dict: ``tb_params``.
        """
        # Read the input file
        with open(fname, "r") as f:
            lines = f.readlines()
        # Drop lines in ignore, comments, and the cross section info table
        # and form parameter key, value pairs in tuples. Keys may repeat at
        # this point
        keys = [l.split()[0] for l in lines if l.split()[0] not in ignore+["CC", "CS"]]
        vals = [l.split()[1:] for l in lines if l.split()[0] not in ignore+["CC", "CS"]]

        # Consolidate the values of repeated keys into a list of lists lines
        multi_keys = [k for k in keys if keys.count(k) > 1]
        multi_vals = [v for k, v in zip(keys, vals) if keys.count(k) > 1]
        tb_params = {k: [] for k in set(multi_keys)}
        for mk, mv in zip(multi_keys, multi_vals):
            tb_params[mk].append([infer_param_type(p) for p in mv])

        # Add the rest of the keys
        for k, v in zip(keys, vals):
            if k in multi_keys: continue
            tb_params[k] = [infer_param_type(p) for p in v]
            if len(v) == 1: tb_params[k] = tb_params[k][0]

        self.tb_params.update(tb_params)
        return tb_params

    def _prepare_vdf_init(self):
        """Infer type of vdf_init request and populate vdf_init_data
        accordingly."""
        # Create directory to store simulation veloctity init data.
        os.makedirs(pjoin(self.directory, "velocity_init"), exist_ok=True)
        # Convert to a list of each particle type
        if np.array(self.hp["vdf_init"], dtype=object).ndim < 2:
            self.hp["vdf_init"] = [self.hp["vdf_init"]]

        self.vdf_init_data = []
        self.tb_params["LV"] = []
        for vdf_info, ptype in self.hp["vdf_init"]:
            if isinstance(vdf_info, int):
                # Interpret as step
                if not vdf_info in self._check_available_vdfs():
                    raise RuntimeError("Requested vdf_init step not available")
                self.get_vdfs(steps=vdf_info, particle_type=ptype, sample_cap=None)
                self.vdf_init_data.append([
                    self.vdfs.iloc[self.vdfs.step == vdf_info, -3:].copy(),
                    ptype,
                ])
            elif isinstance(vdf_info, str):
                # Interpret as path to vdf file
                self.vdf_init_data.append([
                    parsing.extract_vdf(vdf_info, sample_cap=None),
                    ptype,
                ])
            elif isinstance(vdf_info, pd.DataFrame):
                self.vdf_init_data.append([vdf_info.reset_index(drop=True), ptype])
            elif isinstance(vdf_info, np.ndarray) or isinstance(vdf_info, list):
                self.vdf_init_data.append([
                    pd.DataFrame(vdf_info, columns=["vx", "vy", "vz"]),
                    ptype,
                ])

            if self.hp["downsample"]:
                # Truncate init data to match NP
                self.vdf_init_data[-1][0] = (
                    self.vdf_init_data[-1][0].iloc[:self.tb_params["NP"][ptype],:].copy()
                )

            # Add requests to tb input deck
            fname = pjoin("velocity_init", f"ParticleType{ptype}.dat")
            self.tb_params["LV"].append([fname, ptype])


    def _load_indeck(self):
        """Sample and data from indeck options if available."""
        indeck = self.hp["indeck"]
        if callable(indeck):
            # Auto generate CS and indeck with function `indeck`
            # and these parameters
            csps = {k: v for k, v in self.hp.items() if k in CS_PARAMS}
            self.cs, tbp = indeck(**csps)
            # Ignore user specified params
            tb_updates = {k: v for k, v in tbp.items() if k not in self.initial_params}
            self.tb_params.update(tb_updates)
        elif indeck:
            # Interpret as path to pre-written input data
            self.cs = CrossSections(input_path=indeck, read_cs_data=True)
            # Ignore init params set by user
            ignored = list(self.initial_params.keys())
            self.read_tb_params(
                    pjoin(indeck, self.cs.input_fname), ignore=ignored)
        # Otherwise, don't configure indeck data

    def _update_differential_models(self):
        """Update cross section process tables with parameters."""
        # Elastic
        name, params = EADF[self.hp["eadf"]]
        if isinstance(params, list):
            params = " ".join(str(p) for p in params)
        msk = self.cs.table.rtype.str.contains("Elastic")
        self.cs.table.loc[msk, ["model_name", "params"]] = (name, params)
        # Ionization
        name, params = EESD[self.hp["eesd"]]
        if isinstance(params, list):
            params = " ".join(str(p) for p in params)
        msk = self.cs.table.rtype.str.contains("Ionization")
        self.cs.table.loc[msk, ["model_name", "params"]] = (name, params)

    def _impose_nprod_crit(self):
        """Minimum requirement for the product of each species particle count.
        Currently assumes first two species are the electrons and background
        gas respectively."""

        # Aliases
        tb = self.tb_params
        hp = self.hp
        # Use some defaults if not defined
        if "DE" not in hp:
            hp["DE"] = 0.1 # eV
        if "EP_0" not in hp:
            hp["EP_0"] = 10.
        if "Nmin" not in hp: # Npairs multiplier (>=1)
            hp["Nmin"] = 1.0
        # Reduced mass (kg)
        m1, m2 = tb["MP"][0], tb["MP"][1]
        mu = m1*m2/(m1+m2) * AMU_TO_KG # amu to kg
        if not "L" in tb:
            raise RuntimeError("no L specified for DT calculation")

        cs_sets = self.cs.data
        cs_df = pd.DataFrame()
        # Collect all cross sections into one data frame.
        for key, df_ in cs_sets.items():
            df = df_.copy()
            e, cs = df.values.T
            df["rtype"] = key
            # Calculate vsig(eps) for each cs
            df["vsig"] = np.sqrt(2*e*QE/mu)*cs
            cs_df = pd.concat((cs_df, df), ignore_index=True)

        vsig_min_max = np.min(cs_df.groupby("rtype").vsig.agg("max"))
        vsig_sum_max = np.sum(cs_df.groupby("rtype").vsig.agg("max"))
        # Convert reduced field to Vm^2
        EoN = 1e-21*hp["Ered"]
        # Condition: only accelerate EP_0 electrons a max of DE each step
        Ne = hp["Nmin"]*np.sqrt(2*QE*hp["EP_0"]/ME)*2*QE*EoN/(QE*hp["DE"]*vsig_min_max)

        # Take maximum of current number of particles with this condition
        Ne = max(Ne, tb["NP"][0]*tb["NP"][1])

        if "pct_ion" in hp:
            tb["NP"][1] = np.ceil(np.sqrt(Ne/hp["pct_ion"]))
            tb["NP"][0] = np.ceil(hp["pct_ion"]*np.sqrt(Ne/hp["pct_ion"]))

        # Compute E with NP and L
        n_gas = tb["NP"][hp["gas_index"]]/tb["L"]**3
        tb["E"] = EoN*n_gas
        # Compute DT with E, EP_0, DE
        tb["DT"] = QE*hp["DE"]/(QE*tb["E"]) * np.sqrt(ME/(2*QE*hp["EP_0"]))

        # Scale individual parameters beyond calibrated values.
        # Note, this will break E/N initial condition
        for k, v in self.hp["pc_scale"].items():
            if k == "NP":
                for i in range(len(tb[k])):
                    tb[k][i] = int(v*tb[k][i])
            else:
                tb[k] = v*tb[k]

    def _constrain_input(self):
        """Use alternative parameters to derive those required
        for the ThunderBoltz input files."""
        tb = self.tb_params
        hp = self.hp
        # Check gas index
        if hp["gas_index"] > 0 and len(tb["NP"]) == 1:
            # Assume there must be a gas for reduced calculations
            hp["gas_index"] = 0

        if self.hp["autostep"]:
            if not self.hp["DE"]:
                raise RuntimeError("Specify max DeltaE per time step for autostep")
            if not self.hp["Ered"]:
                raise RuntimeError("Specify reduced field (E/N) for autostep")
            if not self.hp["pct_ion"]:
                raise RuntimeError("Specify percent ionization (pct_ion) for autostep")
            # Calculate Nprod constraint and time step
            self._impose_nprod_crit()

        elif hp["pct_ion"] and hp["NN"]:
            # Calculate # of particles
            tb["NP"] = [int(hp["NN"]*hp["pct_ion"]), hp["NN"]]

        # Gas density
        n_gas = tb["NP"][hp["gas_index"]]/tb["L"]**3

        # Calculate other derived parameters, if need be
        if hp["duration"]:
            tb["NS"] = hp["duration"] // tb["DT"]

        if hp["Ered"] is not None:
            # Calculate E-field in V/m such that an electron
            # would be accelerated in the +z direction.
            tb["E"] = -1e-21*hp["Ered"]*n_gas

        if hp["Bred"] is not None:
            if not isinstance(hp["Bred"], list):
                raise RuntimeError("Please specify all three components of the "
                                   "reduced magnetic field.")
            # Calculate absolute B-field in Tesla (1 Gauss = 1e23 Hx m^-3)
            tb["B"] = [1e-27*Bred_i*n_gas for Bred_i in hp["Bred"]]

        # Round NP if its not already an integer
        tb["NP"] = [round(n) for n in tb["NP"]]

    def _validate_input(self):
        """Make sure the ThunderBoltz indeck will be valid before running."""
        # Alias
        tb = self.tb_params
        # Make sure QP and MP are correctly specified
        SP = len(tb["QP"])
        if not SP == len(tb["MP"]):
            raise RuntimeError("QP and MP arguments do not match in length")

        # All other multi-argument settings should not exceed these in length match SP parameter
        if any(len(tb[k]) > SP for k in ["NP", "TP"]):
            raise RuntimeError("Some multi-argument parameters do not match in length")

        # Infer number of species from number of charge and mass arguments
        tb["SP"] = SP

        # If multi-argument settings are shorter in length than inferred SP, set to 0
        arg_diff = SP - len(tb["NP"])
        if arg_diff:
            tb["NP"].extend([0]*arg_diff)
        arg_diff = SP - len(tb["TP"])
        if arg_diff:
            tb["TP"].extend([0.]*arg_diff)
        arg_diff = SP - len(tb["VV"])
        if arg_diff:
            tb["VV"].extend([0.]*arg_diff)

        # Ensure that the indices in the CS indeck portion do not exceed SP
        cs_sp_count = self.cs.table[["r1", "r2", "p1", "p2"]].max().max() + 1
        if cs_sp_count > SP:
            raise IndexError("Collision model species indices exceeds "
                             "specified number of species: "
                             f"({cs_sp_count} > {SP})")

    def set_fixed_tracking(self):
        """Change all reaction species indices of differing reactant values
        to be between only particle 0 and 1 (e.g. 0+1->0+2 is changed to 0+1->0+1)."""
        diff_msk = self.cs.table.r1 != self.cs.table.r2
        if len(diff_msk):
            self.cs.table.loc[diff_msk, ["r1", "r2", "p1", "p2"]] = [0, 1, 0, 1]

    def write_input(self, directory):
        """Write all the input files into a directory with
        the current settings.

        Args:
            directory (str): The path to the simulation directory.
        """
        # Write cross section files
        self.cs.write(directory=directory)
        # Write particle velocity init files if necessary.
        if self.vdf_init_data is not None:
            for i, (vdf_dat, ptype) in enumerate(self.vdf_init_data):
                fname = pjoin(directory, self.tb_params["LV"][i][0])
                vdf_dat.to_csv(fname, header=False, index=False)

        # Write simulation params into input.in
        with open(pjoin(directory, "input.in"), "w") as f:
            for k, v in self.tb_params.items():
                if v is None or v == [] or k not in tb_parameters:
                    continue
                if isinstance(v, list):
                    if isinstance(v[0], list):
                        for line in v:
                            args = " ".join(str(s) for s in line)
                            f.write(f"{k} {args}\n")
                    else:
                        v = " ".join(str(s) for s in  v)
                        f.write(f"{k} {v}\n")
                else:
                    f.write(f"{k} {v}\n")
            # Write CS arguments
            cs_deck = self.cs.get_deck()
            f.write(cs_deck)

        self.written_input = True

    def compile_from(self, src_path, debug=False):
        """Copy TB files from src_path and compile
        in simulation directory.

        Args:
            src_path (str): The location of the source files
                to compile from.
            debug (bool): If ``True``, compile C++ with the ``-g``
                debug flag.

        Raises:
            RuntimeError: if there is a compilation issue.
        """
        # Get source file names
        src_files = ls(src_path)
        # Copy each to simulation directory
        for sf in src_files:
            if ".cpp" in sf or ".h" in sf:
                src_file_path = pjoin(src_path, sf)
                dest_file_path = pjoin(self.directory, sf)
                shutil.copy(src_file_path, dest_file_path)

        # Go to simulation directory and compile
        pwd = os.getcwd()
        os.chdir(self.directory)

        # Compile (hard-coded for security reasons).
        cmd = "g++ -std=c++17 DSMC0D.cpp -o thunderboltz.bin -Wall -Werror -Wsign-compare"
        # Debug
        if debug:
            cmd = "g++ -g -std=c++17 DSMC0D.cpp -o thunderboltz.debug -Wall -Werror -Wsign-compare"
        cerr = os.system(cmd)
        if cerr:
            raise RuntimeError("Failed to compile.")

        os.chdir(pwd)

    def compile_debug(self):
        """Prepare all files and compile ThunderBoltz with -g debug flag."""
        self.run(dryrun=True, debug=True)

    def add_callback(self, f):
        """Add a function to the list of functions that will be
        called during banner output.

        Args:
            f (callable[,]): The function that will be called. It
                should accept no arguments and return no arguments.
        """
        self.callbacks.append(f)

    def run(self, src_path=None, bin_path=None, out_file="thunderboltz",
            monitor=False, dryrun=False, debug=False, std_banner=False,
            live=False, live_rate=False):
        """Execute with the current parameters in the simulation directory.
           The internal API ThunderBoltz version will be used in lieu of
           user-provided binary/source files.

        Args:
            src_path (str): Optional path to source files to copy into the
                simulation directory. The source is then compiled there.
            bin_path (str): Optional path to binary executable to copy
                into the simulation directory.
            out_file (str): The file name for stdout buffer of the
                ThunderBoltz process.
            monitor (bool): Runtime flag, when set to `True` will generate
                an empty `monitor` file in the simulation directory. Deleting this
                file will cause the ThunderBoltz process to exit, but allow
                the wrapper to continue execution. This is useful performing
                several simulation calculations sequentially, but a manual exit
                is required for each one.
            dryrun (bool): Setup all the files for the calculation, but do not
                run the calculation.
            debug: (bool): Compile with C++ ``-g`` debug flag.
            std_banner (bool): Toggle banner output streaming to stdout in
                addition to being written to the ``out_file`` buffer.
            live (bool): Run and update time series plotting GUI during simulation.
            live_rate (bool): Run and update rate plotting GUI during simulation.

        Raises:
            RuntimeError: if there is no simulation directory set or if
                the one provided does not exist.
        """
        if not self.directory:
            raise RuntimeError(
                "Simulation directory has not been established")
        elif not os.path.isdir(self.directory):
            raise RuntimeError(
                "Simulation directory does not exist")

        # Overwrite monitor, out_file, and live options
        if self.monitor is None:
            self.monitor = monitor
        if self.out_file is None:
            self.out_file = out_file
        if self.live is None:
            self.live = live
        if self.live_rate is None:
            self.live_rate = live_rate

        # Validate parameters and write indeck files
        self._validate_input()
        if not self.written_input:
            self.write_input(self.directory)

        # Determine binary TB file
        if not bin_path:
            if not src_path:
                # Use internal package src
                src_path = tb.src_path
            # Create binary file from source
            self.compile_from(src_path, debug=debug)
        else:
            # Copy binary file from user path
            bin_file = Path(bin_path).parts[-1]
            shutil.copy(bin_path, pjoin(self.directory, bin_file))

        # Record time and log params
        self.runtime_start = datetime.now().strftime("%d-%m-%y_%H:%M:%S")
        self._log()

        # Do not run ThunderBoltz for dry runs
        if dryrun: return

        ##### Run the program #####

        # Sub-process commands hard-coded for security reasons
        args = ["./thunderboltz.bin", "input.in"]

        # Create callbacks for communication between the
        # ThunderBoltz thread and this parent thread.

        if self.monitor:
            # Create temporary monitor file
            with open(pjoin(self.directory, "monitor"), "w") as f:
                f.write("") # Write empty file
            self.add_callback(self._monitor_callback)

        if self.live:
            # Live time series plotting
            self.add_callback(self._live_ts_callback)
        if self.live_rate:
            # Live rate coefficient plotting
            self.add_callback(self._live_rate_callback)

        # Start precision timer
        t1 = time.perf_counter()

        # Start the ThunderBoltz subprocess
        self._interactive_subprocess(args, tee=std_banner,
            std_buf=pjoin(self.directory, self.out_file+".out"),
            err_buf=pjoin(self.directory, self.out_file+".err"),
        )

        # After subprocess exits, record wall-clock time.
        self.elapsed_time = time.perf_counter() - t1
        self.runtime_end = datetime.now().strftime("%d-%m-%y_%H:%M:%S")

        # Update log files with elapsed times
        self._log()

    def _interactive_subprocess(self, args, tee=True, std_buf=None,
            err_buf=None):
        """Open a subprocess for the ThunderBoltz simulation and
        pipe the output back to the parent. Optionally allow for
        'tee'-like behavior and also print out process
        stdout at the same time."""
        # Start process in simulation directory
        with visit(self.directory):
            self._proc = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)

        p = self._proc # Alias
        # Ensure file reading doesn't wait for subprocess to close when
        # sending updates
        os.set_blocking(p.stdout.fileno(), False)

        if std_buf: f = open(std_buf, "a", buffering=1)

        def run_callbacks():
            lines = p.stdout.readlines()
            txt = "".join(s.decode() for s in lines)
            if std_buf: f.write(txt)
            if self.callbacks: [c() for c in self.callbacks]
            if tee: print(txt, end="")

        while p.poll() is None:
            # Release lock to the subprocess if running
            # on one core.
            time.sleep(2)

            run_callbacks()

        # Run buffering one last time.
        run_callbacks()

        if std_buf: f.close()

        # Log errors into .err file after failure and raise
        # as warning in the parent
        err = "\n".join(s.decode() for s in p.stderr.readlines())
        if err_buf:
            with open(err_buf, "w") as e:
                e.write(err)
        elif err:
            warnings.warn(err, RuntimeWarning)

        return p

    def _monitor_callback(self):
        """Exit the ThunderBoltz process when the
        monitor file is removed by the user."""
        p = self._proc # Alias

        # Check if monitor file is still there, or if process ended.
        monfile = pjoin(self.directory, "monitor")
        if not os.path.isfile(monfile) or p.poll() is not None:
            # At this point, the process is either terminated, or should be
            if p.poll() is None:
                p.kill()
            # And the monitor file is either deleted, or should be
            if os.path.isfile(monfile):
                os.remove(monfile)

    def _live_ts_callback(self):
        """Update the time series GUI."""
        # Make sure output files have data
        if not os.path.isfile(
                pjoin(self.directory, self.out_file+".out")):
            return # Don't try to plot anything

        # Clear plots for data redraw if not initialized yet
        if self.ts_fig: [a.clear() for a in self.ts_fig.axes]
        self.plot_timeseries()
        if self.ts_fig: self.ts_fig.canvas.draw()
        plt.pause(0.2)

    def _live_rate_callback(self):
        """Update the rate plot GUI."""
        if not os.path.isfile(
                pjoin(self.directory, self.out_file+".out")):
            return # Don't try to plot anything

        # Clear plots for data redraw if not initialized yet
        if self.rate_fig: [a.clear() for a in self.rate_fig.axes]
        self.plot_rates()
        if self.rate_fig: self.rate_fig.canvas.draw()
        plt.pause(0.2)

    def read_stdout(self, fname):
        """Read the banner output data.

        Args:
            fname (str): The name of the ``.out`` file to read.

        Returns:
            :class:`pandas.DataFrame`: The banner data.

        Note:
            If ThunderBoltz warnings are found (e.g. particle overload),
            a message will be appended to ``err_stack``.
        """
        fpath = pjoin(self.directory, fname)
        with open(fpath, "r") as f:
            # Get to data lines
            l = f.readline()
            while "0, t" not in l:
                if "Too many particles!" in l:
                    self.err_stack.append("particle overload")
                if not l:
                    self.kinetic_table = pd.DataFrame()
                    return

                l = f.readline()
            # Read the rest of the lines
            lines = [l] + f.readlines()

        if not lines:
            self.kinetic_table = pd.DataFrame()
            return

        # Convert line data to DataFrame
        if "Too many particles!" in lines[-1]:
            lines = lines[:-1]
            self.err_stack.append("particle overload")

        lines = [l.replace(" ", "") for l in lines]
        # Ignore comments
        lines = [l for l in lines if not l.lstrip().startswith("#")]
        txt = "".join(lines)
        cmax = len(lines[0].split(","))
        a = np.genfromtxt(io.StringIO(txt), delimiter=",",
            usecols=range(0, cmax+1, 2), autostrip=True, ndmin=2)
        columns = ["step"] + [lines[0].split(",")[i] for i in range(1, cmax, 2)]
        self.kinetic_table = pd.DataFrame(a, columns=columns)
        return self.kinetic_table

    def read_particle_table(self, i):
        """Read species specific output data, including
        density, velocity, displacement, energy, and temperature.

        Args:
            i (int): The species index.

        Returns:
            :class:`pandas.DataFrame`: The particle data for
            species ``i``.
        """

        fpath = pjoin(self.directory, f"Particle_Type_{i}.dat")
        with open(fpath, "r") as f:
            column_names = f.readline().replace(",", "").split()
            lines = f.readlines()

        lines = [l.replace(" ", "") for l in lines]
        if not lines: return
        txt = "".join(lines)
        a = np.genfromtxt(io.StringIO(txt), delimiter=",",
                autostrip=True, ndmin=2)
        self.particle_tables[i] = pd.DataFrame(a, columns=column_names)
        return self.particle_tables[i]

    def read(self, directory=None, read_input=True, read_cs_data=False, only=None):
        """Read the simulation directory of a ThunderBoltz run, possibly
        all of its input and output files.

        Args:
            directory (str): The location of the simulation directory
                from which to read.
            read_input (bool): Whether or not to read any input data.
            read_cs_data (bool): Whether or not to read cross section data.
                This can be expensive, and often isn't necessary.
            only (list or None): Only read certain types of files. Default is
                ``["thunderboltz.out", "Particle_Type", "Counts.dat", ""thunderboltz.log"]``.
        """
        if directory: self.directory = directory

        if read_input:
            # Read cross sections input first
            self.cs.read(self.directory, read_cs_data=read_cs_data)
            self.read_tb_params(pjoin(self.directory, self.cs.input_fname))
        # Read other files
        logfiles = ls(self.directory)
        if only:
            logfiles = [o for o in only if o in logfiles]
        for filename in logfiles:
            fpath = pjoin(self.directory, filename)
            if filename == "Counts.dat":
                counts_dat = np.genfromtxt(fpath, delimiter=",")
                # If only one row of data, ensure 2D array
                if counts_dat.ndim == 1 and len(counts_dat) == len(self.cs.table)+1:
                    counts_dat = counts_dat.reshape(1, -1)
                self.counts = pd.DataFrame(counts_dat,
                    columns=["step"]+self.cs.table.csfile.to_list())
            elif filename == "thunderboltz.out":
                self.read_stdout(filename)
            elif filename == "thunderboltz.log":
                self.read_log(filename)

        # Read particle type files in order
        if only and "Particle_Type" not in only:
            return
        pfiles = [f for f in ls(self.directory) if "Particle_Type" in f]
        self.particle_tables = [None]*len(pfiles)
        for i in range(len(pfiles)):
            self.read_particle_table(i)

    def get_particle_tables(self):
        """Return the particle table data for each species
        in a list.

        Returns:
            list[ :class:`pandas.DataFrame`]: The list of particle
            table data. Each table with have columns with keywords
            matching the attributes of
            :class:`~.thunderboltz.parameters.ParticleParameters`.
        """
        # Implicitly call the main output parsing code.
        self.get_timeseries()
        return self.particle_tables

    def get_timeseries(self):
        """Collect the relevant time series data from a ThunderBoltz simulation
        directory and add input parameter columns.

        Returns:
            :class:`pandas.DataFrame`: The table of time series data.
        """
        # Read output files
        self.read(only=self.output_files, read_input=False)
        # Make sure actual output files got read
        pt = self.particle_tables
        kt = self.kinetic_table

        if kt is None:
            raise RuntimeError("Missing output logs, make sure stdout and "
                               "particle tables are being written.")
        if pt == [] or pt[0] is None or (len(pt[0]) < 2 or len(kt) < 2):
            # Do not try to join until more data is available
            self.timeseries = pd.DataFrame()
            return self.timeseries

        # Create time series
        ts = self.timeseries = pt[0].join(kt.drop("t", axis=1), how="inner")

        # Add input and meta parameters
        self._tabulate(ts)

        # Calculate gas density
        ts["n_gas"] = self.tb_params["NP"][self.hp["gas_index"]]/ts.L**3
        # Calculate reduced field
        ts.loc[:, "mobN"] = ts.Vzi.abs()/ts.E.abs()*ts.n_gas

        # Compute some additional parameters
        self.DT = ts["t"].values[1] - ts["t"].values[0]
        # For convenience include electron z position in main time series
        ts["Rzi"] = ts.Rzi

        # Reaction rates
        self._calculate_rates()
        # Diffusion rates
        self._calculate_diffusion()
        # Bulk swarm flows / reaction rates
        self._calculate_bulk_rates()

        # Calculate the reduced Townshend ionization coefficient
        ts["a_n"] = ts.k_i / ts.Vzi.abs()

        return ts

    def _calculate_rates(self, ts=None, ddt=None, std=False):
        """Calculate rate coefficients from reaction counts.
        These calculations require synchronization across different
        ThunderBoltz output files, so the maximum amount of data available across
        all files is used.
        """
        if ts is None: ts = self.timeseries
        pt = self.particle_tables
        c = self.counts.set_index("step")

        # Mask counts and particle tables to match timeseries data
        c = c[c.index >= ts["step"].min()].copy()
        for i in range(len(pt)):
            pt[i] = pt[i][pt[i]["t"] >= ts["t"].min()].copy()

        # Find the shortest list (if files are still being written by process)
        npoints = min(len(ts), *[len(p) for p in pt], len(c))
        # Truncate counts to the shortest
        c = c.iloc[:npoints].copy()

        # Default time derivative is forward difference estimator
        if ddt is None: ddt = lambda x, name=None: (x.diff()/self.DT, 0*x)

        # Compute collision rate data
        dcdt, c_std = ddt(c)

        ctab = self.cs.table
        # Get ionization collision files.
        ionz_files = ctab.loc[(ctab.rtype.str.contains("Ionization"))
                             |(ctab.csfile.str.contains("ion")), "csfile"]
        ionz_c = c[ionz_files].copy()
        # Compute total ionization count
        tot_ionz_c = ionz_c.sum(axis=1)
        # Compute total interaction count over all processes
        total_c = c.sum(axis=1)
        ionz_rate, ionz_std = ddt(tot_ionz_c)
        tot_rate, tot_std = ddt(total_c)
        # Generate specific rates for ionization and total processes
        k_scale = 1/(ts.n_gas*ts.Ni*ts.L**3)[:npoints]
        ts.loc[:npoints,"k_i"] = k_scale.values*ionz_rate.values
        ts.loc[:npoints,"k_tot"] = k_scale.values*tot_rate.values

        if std:
            # Compute std data
            ts.loc[:npoints,"k_i_std"] = k_scale.values*ionz_std.values
            ts.loc[:npoints,"k_tot_std"] = k_scale.values*tot_std.values
        
        # Generate rates for all processes individually
        for i, row in ctab.iterrows():
            k_scale = 1/(pt[row.r1].Ni*pt[row.r2].Ni*ts.L**3)[:npoints]
            k_scale[k_scale == np.inf] = 0
            ts.loc[:npoints,f"k_{i+1}"] = k_scale.values*dcdt.iloc[:npoints,i].values
            if std:
                # Compute std data per process
                ts.loc[:npoints,f"k_{i+1}_std"] = k_scale.values*c_std.iloc[:npoints,i].values

    def _calculate_bulk_rates(self, ts=None, ddt=None, std=False):
        """Calculate the bulk flow rates and the associated bulk
        mobility and ionization coefficients."""

        if ts is None: ts = self.timeseries
        # Default time derivative is forward difference estimator
        if ddt is None: ddt = lambda x, name=None: (x.diff()/self.DT, 0*x)

        # Calculate bulk parameters
        zdrift, zdrift_std = ddt(ts.Rzi, r"$\langle r_z \rangle$")
        ts["Wz"] = np.abs(zdrift)
        ts["mobN_bulk"] = zdrift/ts.E.abs()*ts.n_gas
        ts["a_n_bulk"] = ts.k_i / zdrift

        if std:
            ts["Wz_std"] = zdrift_std
            ts["mobN_bulk_std"] = zdrift_std/ts.E.abs()*ts.n_gas
            # Division of two uncertain parameters propagates error
            # as s(A/B) = |A/B| * sqrt((s(A)/A)^2 + (s(B)/B)^2)
            if "k_i_std" in ts:
                ts["a_n_bulk_std"] = np.abs((ts.k_i/zdrift) * (zdrift_std/zdrift))
            else:
                ts["a_n_bulk_std"] = (np.abs(ts.k_i/zdrift)
                  * np.sqrt((ts.k_i_std/ts.k_i)**2 + (zdrift_std/zdrift)**2))

    def _calculate_diffusion(self, ts=None, ddt=None, std=False):
        r"""Calculate the diffusion coefficients from position/velocity
        correlation data.

        Args:
            ts (:class:`pandas.DataFrame`): Specify timeseries data
                to calcaulte the diffusion rates for. Default is all time
                series data.
            ddt (callable): Specify a custom derivative method.
            std (bool): Option to return standard deviation data propagated
                through derivative opererations (i.e. if fits are used to compute slope values).

        Note:
            These diffusion calculation assume that the magnetic field is
            oriented in the :math:`y`-direction. In particular, the Hall diffusion
            is calculated using <xz> correlation data.
        """
        if ts is None: ts = self.timeseries

        # Default to forward difference time derivative operator
        if ddt is None: ddt = lambda x, name=None: (x.diff()/self.DT, 0*x)

        for axis in ["X", "Y", "Z"]:
            # Compute on-axis bulk diffusion coefficients
            R = ts[f"R{axis.lower()}i"]
            dii, dii_std = ddt(1/2 * (ts[2*axis] - R**2),
                parameters.latex[f"D_fit_{axis}"])
            ts[f"D_{axis}{axis}_bulk"] = dii
            # Compute on-axis flux diffusion coefficients
            V = ts[f"V{axis.lower()}i"]
            RV = ts[f"{axis}V{axis}"]
            ts[f"D_{axis}{axis}"] = RV - R*V

            if std:
                ts[f"D_{axis}{axis}_bulk_std"] = dii_std

        # Compute hall diffusion
        # fig, ax = plt.subplots()
        dh, dh_std = ddt(ts["XZ"] - ts["Rxi"]*ts["Rzi"], parameters.latex["D_fit_H"])
        ts["D_H_bulk"] = dh
        ts["D_H"] = ts["XVZ"] + ts["ZVX"] - ts["Rxi"]*ts["Vzi"] - ts["Rzi"]*ts["Vxi"]

        if std:
            ts["D_H_bulk_std"] = dh_std

    def get_ss_params(self, ss_func=None, fit=False):
        """Get steady-state transport parameter values by averaging last section
        of time series. By default, the last fourth of the available data is
        considered to be steady-state. Standard deviations over this interval
        will be computed for each parameter in a new column with a "_std"
        suffix added to the column name.

        Args:
            ss_func (callable[:class:`pandas.DataFrame`,
                :class:`~.thunderboltz.parameters.TBParameters`]->:class:`pandas.DataFrame`):
                A function that takes in numerical time series data and returns a
                new frame with only data that is considered to be at steady state.
                For example, ``ss_func=lambda df: df[df.t > 1e-6]`` would select only
                times in the simulation after one microsecond for steady state
                calculations, or ``ss_func=lambda df: df.iloc[50:,:]`` would select
                the last 50 time steps.
            fit (bool): Option to use a line of best fit over the steady state
                window to calculate time dependent parameters (bulk swarm parameters
                and rate coefficients).
                rather than averaging derivatives with the forward difference formula.
                ``True`` does so for all rate parameters. Default is False.

        Returns:
            :class:`pandas.DataFrame`: The aggregated steady-state data for
            each output parameter along with columns specifying the input
            parameters.

        Raises:
            RuntimeWarning: if not enough steps are available to compute steady
                state statistics.

        Warning:
            Currently, the last quarter of the time series data is assumed to be
            in steady-state by default when calculating these steady-state parameters.
            One can verify that this is true by viewing the figures produced by
            :py:func:`plot_timeseries`. Otherwise, one may run the simulation for longer,
            or provide the appropriate steady state criteria via ``ss_func``.
        """
        # Require timeseries data
        if self.timeseries is None:
            ts = self.get_timeseries()
            if not len(ts):
                # Return empty data
                return pd.DataFrame()

        ts = self.timeseries
        # Get numerical time series data
        tsn = ts.select_dtypes(include=[np.number]).reset_index(drop=True)

        if ss_func: # Take user steady state function is specified
            ts_last = ss_func(tsn, self).copy().reset_index(drop=True)
        else: # Otherwise take last fourth
            ts_last = tsn[tsn.index > 3*len(tsn)//4].copy().reset_index(drop=True)

        # Must have steady state data
        require_ts = 3
        if len(ts_last) < require_ts:
            raise RuntimeError("Not enough steps to generate steady state statistics")
        # Set time step after which stats are generated
        self.time_conv = ts_last.t.values[0]

        if fit:
            # Recalculate using fit method for time derivative calculations
            self._calculate_rates(ts_last, ddt=self.compute_fit, std=True)
            self._calculate_diffusion(ts_last, ddt=self.compute_fit, std=True)
            self._calculate_bulk_rates(ts_last, ddt=self.compute_fit, std=True)

        # Calculate steady-state statistics by aggregating over time
        mean = ts_last.agg("mean")
        std = ts_last[[c for c in ts_last if c+"_std" not in ts_last]].agg("std")
        stats = pd.concat((mean, std.add_suffix("_std")), axis=0)
        stats = stats.drop(["t", "step"])

        # Discard double aggregated _std_std columns 
        stats = stats.drop([c for c in stats.index if "_std_std" in c])

        # Add input / meta parameters
        stats = pd.DataFrame([stats], index=[0])
        stats = self._tabulate(stats)
        return stats

    def get_slope_fit():
        pass

    def get_etrans(self):
        """Return the energy weighted counts of each reaction computed
        for each time step.

        Returns:
            :class:`pandas.DataFrame`: A table of each process and the
            energy transfer through that channel as a proportion
            to the total energy transfer.
        """
        table = self.cs.table.set_index("csfile")
        weighted = self.counts.copy()
        # Loop through columns of reaction count table.
        for reaction in weighted.columns[1:]:
            weighted[reaction] *= table.loc[reaction, "B"]
        weighted.iloc[:, 1:] = weighted.iloc[:, 1:].divide(weighted.iloc[:, 1:].sum(axis=1), axis=0)
        return weighted

    def get_vdfs(self, steps="last", sample_cap=500000, particle_type=0, v=0):
        """Read the electron velocities arrays within a ThunderBoltz calculation.
        Velocity units are in m/s. If velocity dump files are found
        corresponding to ``steps``, update ``vdfs``.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.
            particle_type (str, list[int], or int): Specify which kinds of species data
                should be read from.

                * ``int``: The particle type to read. Default is ``0``.
                * ``list[int]``: A set of particle types to read.
                * ``"all"``: Read all particle types.

            v (int): Verbosity -- 0: silent, 1: print file paths before reading.

        Returns:
            :class:`pandas.DataFrame`: The particle velocity dump data.

        Warning:
            Large files are truncated to the first ``sample_cap`` lines.
            It is assumed that the particle ordering in the dump files
            is not correlated with any velocity statistics, but
            this may not be the case when ``egen`` is on. In that case, ensure
            the entire velocity dump file is being read.
        """
        # Check what's available out of the requested
        available = self._check_available_vdfs(steps)
        if not len(available):
            return # Nothing to plot
        if self.vdfs is None:
            # Initialize empty DataFrame with named columns
            self.vdfs = pd.DataFrame({k: pd.Series([], dtype=dtype)
                for k, dtype in zip(["ptype", "step", "vx", "vy", "vz"],
                                    [np.int64]*2 + [np.float64]*3)})
                # Keep "vx ... vz" as the last three columns

        # Find which steps are missing from current vdf table
        unique = self.vdfs[self.vdfs.ptype==0].step.unique()
        missing = [step for step in available
                   if step not in unique]
        # Extract them from the simulation directory.
        missing_vdfs = parsing.extract_vdfs(
            self.directory, steps=missing, sample_cap=sample_cap,
            particle_type=particle_type, v=v,
        )

        # Only read the missing ones
        self.vdfs = pd.concat((self.vdfs, missing_vdfs), ignore_index=True)
        return self.vdfs

    def get_edfs(self, steps="last", sample_cap=500000):
        """Read the electron velocity distribution functions and return
        the component and total energy distributions within a
        ThunderBoltz calculation. Energy units are in eV. Invokes
        :func:`get_vdfs`.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.

        Returns:
            :class:`pandas.DataFrame`: A table with the signed and unsigned
                energy components of each particle.
        """

        self.get_vdfs(steps, sample_cap)
        edf = self.vdfs.copy().rename(
            columns={k: k.replace("v", "E") for k in self.vdfs.columns[1:]})
        # Transform columns to eV, first signed energies
        m = self.tb_params["MP"][0]*AMU_TO_KG
        for col in ["Ex", "Ey", "Ez"]:
            edf[col+"_signed"] = 1/2 * m/QE * edf[col] * edf[col].map(np.abs)
        # Then the unsigned
        for col in ["Ex", "Ey", "Ez"]:
            edf[col] = edf[col+"_signed"].abs()
        # Total energy is the sum of its orthogonal components.
        edf["E"] = edf.Ex + edf.Ey + edf.Ez
        return edf

    def get_counts(self):
        """Return the cumulative reaction counts for each process and
        each time step.
        
        Returns:
            :class:`pandas.DataFrame`: A table where each column corresponds to a
            collision process and each row corresponds to a time step.
        """
        return self.counts

    def describe_dist(self, steps="last", sample_cap=500000):
        """Generate percentile and count statistics of the electron
        velocity / energy distribution for various time steps.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.

        Returns:
            :class:`pandas.DataFrame`: A table with statistical descriptions of
            the velocity and energy distributions.
        """

        # Get energy distribution, this will implicitly update self.vdfs
        efs = self.get_edfs(steps, sample_cap)
        # Concatenate the two frames
        joined = pd.concat((self.vdfs, efs.iloc[:,-4:]), axis=1)
        stats = (joined
            .groupby(["ptype", "step"])
            .describe(percentiles=[0.01, .1, .25, .5, .75, .9, .99])
            .stack()
        )
        return stats

    def set_ts_plot_params(self, params):
        """Set the default series plotted by
        :meth:`~thunderboltz.ThunderBoltz.plot_timeseries`"""
        self.ts_plot_params = params

    def plot_timeseries(self, series=None, save=None, stamp=[], v=0, update=True):
        """Create a diagnostic plot of ThunderBoltz time series
        data.

        Args:
            series (list[str|tuple]):
                The y-parameters to plot onto the time series figure. If the
                element is a string, the corresponding parameter will be plotted
                if available in :class:`~.thunderboltz.parameters.OutputParameters`.
                If the element is a tuple, the first argument will be interpreted as the
                name of a new user defined parameter, and the second argument a user defined
                function that calculates it. The function must accept timeseries data
                and return a single series to be plotted.
            save (str): Option to save the plot to a file path.
            stamp (list[str]): Option to stamp the figure with the value of
                descriptive parameters, e.g. the field, or initial
                number of particles. See :class:`~.thunderboltz.parameters.TBParameters`
                and :class:`~.thunderboltz.parameters.WrapParameters`.
            v (int): Verbosity -- 0: silent, 1: print file paths before
                plotting.
            update (bool): If set to ``False``, assume required data has
                already been parsed into ThunderBoltz frames.

        Returns:
            :class:`matplotlib.figure.Figure`: The timeseries figure object.
        """
        if series: self.ts_plot_params = series
        # If no data has been read yet, or program still running, read it
        process_running = self._proc and self._proc.poll() is None
        if update and (self.timeseries is None or process_running):
            self.get_timeseries()

        if not len(self.timeseries):
            return # No data

        # Create figure and ax objects, if necessary
        if not self.ts_fig:
            self.ts_fig, _ = plt.subplots(figsize=(14,9))
        fig = self.ts_fig # Alias

        if self.timeseries is None:
            raise RuntimeError("Timeseries data not read.")

        def format_ts(ts):
            # Convert mobility values
            ts["mobN"] = 1e-24*ts.mobN
            if "mobN_bulk" in ts:
                ts["mobN_bulk"] = 1e-24*ts.mobN_bulk
            if "mobN_bulk_fit" in ts:
                ts["mobN_bulk_fit"] = 1e-24*ts.mobN_bulk_fit
            return ts

        ts = format_ts(self.timeseries.copy())
        # Plot timeseries with this figure
        plot_timeseries(fig, ts, series=self.ts_plot_params, save=save)
        # super title
        stamps = copy.deepcopy(self.hp)
        stamps.update(self.tb_params)
        if "directory" not in stamps:
            stamps["directory"] = Path(self.directory).parts[-1]
        stamps = {k: v for k, v in stamps.items() if k in stamp}
        # Use pandas Series to create a pretty string.
        stamp_fmt = pd.DataFrame([pd.Series(stamps, dtype="object")]).T
        # stamp_fmt.loc[:,0]
        sss = stamp_fmt.to_string(float_format=lambda s: f"{s:.3e}" if s=="DT" else str(s))
        stt = "\n".join(str(sss).split("\n")[1:])
        fig.suptitle(stt, fontsize=10, fontproperties={"family": "monospace"})

        if save:
            # Save with the name of the simulation directory
            sname = Path(self.directory).parts[-1]
            if not os.path.isdir(save):
                os.makedirs(save)
            saveas = pjoin(save, sname+".pdf")
            if v: print(f"Saving time series in {saveas}")
            fig.savefig(saveas)
            plt.close(fig)

        return fig

    def plot_rates(self, save=None, stamp=None, v=0, update=True):
        """Create a diagnostic plot of ThunderBoltz time series
        data.

        Args:
            save (str): Option to save the plot to a file path.
            stamp (list[str]): Option to stamp the figure with the value of
                descriptive parameters, e.g. the field, or initial
                number of particles. See :class:`~.thunderboltz.parameters.TBParameters`
                and :class:`~.thunderboltz.parameters.WrapParameters`.
            v (int): Verbosity -- 0: silent, 1: print file paths before
                plotting.
            update (bool): If set to ``False``, assume required data has
                already been parsed into ThunderBoltz frames.

        Returns:
            :class:`matplotlib.figure.Figure`: The plot_rate figure object.
        """
        if self.rate_fig is None:
            self.rate_fig, ax = plt.subplots(figsize=(14, 10))

        ax = self.rate_fig.axes[0]

        if update: self.get_timeseries()
        ts = self.timeseries

        kkeys = [k for k in ts.columns if k.startswith("k_") and "i"
                 not in k]
        for col in ts.columns:
            if not col.startswith("k_") or "i" in col: continue
            ax.plot(ts.t, ts[col], label=col)
            ax.legend(fontsize=10)
            ax.set_yscale("log")
            ax.set_ylabel("Rate Coefficient (m$^3$/s)")
            ax.set_xlabel("Time (s)")

        if save:
            # Save with the name of the simulation directory
            sname = Path(self.directory).parts[-1]
            if not os.path.isdir(save):
                os.makedirs(save)
            saveas = pjoin(save, sname+"_rates.pdf")
            if v: print(f"Saving rate constants in {saveas}")
            self.rate_fig.savefig(saveas)
            plt.close(self.rate_fig)

        return self.rate_fig


    def plot_edf_comps(self, steps="last", sample_cap=500000,
                       bins=100, maxwellian=True, save=None):
        """Plot the directional components of the energy distribution function.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.
            bins (int): Total number of bins to divide the energy space into.
            maxwellian (bool): Option to draw a maxwellian distribution
                with the same temperature for comparison.
            save (str): Optional location of directory to save the figure in.
        """
        # TODO: stamp with temperature for species.
        # TODO: Add variable-width bins
        # TODO: Add symmetrical log plots.
        # TODO: Add option to choose between the above.

        # XXX: For now, it is assumed that only VDFs for electrons exist.

        edfs = self.get_edfs(steps="last", sample_cap=500000)
        # Save pretty component names
        xlabels = [r"$\e" + f"psilon{c}$ (eV)" for c in ["_x", "_y", "_z", r"_{\rm tot}"]]
        figs = []
        steps = []
        for step, edf in edfs.groupby("step"):
            steps.append(step)
            fig, axs = plt.subplots(2, 2)
            steps.append(step)
            figs.append(fig)
            for i, (comp, ax) in enumerate(
                    zip(["Ex", "Ey", "Ez", "E"], axs.flatten())):
                # Get component series
                edf_c = edf[comp+"_signed"] if comp != "E" else edf["E"]
                # Plot the histogram
                edf_c.plot.hist(ax=ax, bins=bins)
                mean_e = edf["E"].mean()
                Ne = len(edf) # Counts of particles

                # Create maxwellian plot
                _, divs = np.histogram(edf_c.values, bins=bins)
                dE = divs[1] - divs[0]
                centers = divs[:-1] + dE/2
                prefactor=1
                # Mass of this species
                m = self.tb_params["MP"][0]
                if comp == "E":
                    # The speed distribution in units of energy
                    maxwellian = (lambda e: Ne*dE * e/(3/2*mean_e)**2 * np.exp(-abs(e)/(3/2*mean_e)))
                else:
                    # The 1D Gaussian velocity distribution in units of energy
                    maxwellian = lambda e: Ne*dE / (3*mean_e) * np.exp(-2/3*abs(e)/mean_e)
                maxw_y = np.vectorize(maxwellian)(centers)
                ax.plot(centers, maxw_y, label="Maxwellian")

                # View on a log scale
                ax.set_yscale("log")
                ax.set_ylabel(None)
                ax.set_xlabel(xlabels[i])
                ax.set_ylim(.9, None)
            fig.suptitle(f"step {step}", fontsize=13)
            fig.subplots_adjust(hspace=.18, wspace=.15)
            # Label axes
            axs[0,0].set_ylabel("Count")
            axs[1,0].set_ylabel("Count")
            fig.tight_layout()

            if save:
                sname = Path(self.directory).parts[-1]
                if not os.path.isdir(save):
                    os.makedirs(save)
                saveas = pjoin(save, sname+f"{step}_edf.pdf")
                print(f"Saving edf plots in {saveas}")
                sys.stdout.flush()
                fig.savefig(saveas)
                plt.close(fig)

        return figs, steps


    def plot_edfs(self, steps="last", sample_cap=500000,
            bins=100, plot_cs=False, save=None):
        """Plot the electron total energy distribution function, optionally
        include the provided cross sections for comparison.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.
            bins (int): Total number of bins to divide the energy space into.
            maxwellian (bool): Option to draw a maxwellian distribution
                with the same temperature for comparison.
            save (bool): Optional location of directory to save the figure in.

        Returns:
            (Tuple[list[matplotlib.figure.Figure], list[int]): The list of
            figures and a list of their corresponding step indices.

        Note:
            It currently assumed that only data for one particle type is to be
            plotted.
        """
        edfs = self.get_edfs(steps, sample_cap)

        fig, ax = plt.subplots(figsize=(8,7))

        if plot_cs:
            # Plot it on the right y-axes
            self.cs.plot_cs(ax=ax.twinx())

        for step, df in edfs.groupby("step"):
            vals, divs = np.histogram(df.E, bins=bins)
            de = divs[1] - divs[0]
            evals = divs[:-1] + de/2
            ax.plot(evals, vals/de, label=f"step {step}")

        ax.legend()
        ax.set_ylabel("Particle number density (eV$^{-1}$)")
        ax.set_xlabel("Energy (eV)")

        fig.tight_layout()


    def plot_vdfs(self, steps="last", save=None, bins=100,
            sample_cap=500000):
        """Plot the joint distribution heat map between the x-y and x-z
        velocities.

        Args:
            steps (str, list[int], or int): Options for which time steps to
                read:

                * ``"last"``: Only read the VDF of the last time step
                * ``"first"``: Only read the VDF of the first time step
                * ``"all"``: Read a separate VDF for each time step.
                * ``list[int]``: Read VDF for each time step included in list.
                * ``int``: read VDF at one specific time step.

            sample_cap (int): Limit the number of samples read from the dump
                file for very large files. Default is 500000. If bool(sample_cap)
                evaluates to ``False``, then no cap will be imposed.
            bins (int): Total number of bins to divide the energy space into.
            save (str): Optional location of directory to save the figure in.

        Returns:
            (Tuple[list[matplotlib.figure.Figure], list[int]): The list of
            figures and a list of their corresponding step indices.
        """
        # Check whats available out of the requested
        self.get_vdfs(steps, sample_cap)
        available = self._check_available_vdfs(steps)
        # Confine plots to only the requested
        vdfs_data = self.vdfs[self.vdfs.step.isin(available)].copy()

        # Loop through included vdf groups
        figs = []
        steps = []
        for step, vdf in vdfs_data.groupby("step"):
            fig_xy = plot_vdfs(self.vdfs, bins=bins)
            fig_xz = plot_vdfs(self.vdfs, pair=["vx", "vz"], bins=bins)
            if save:
                sname = Path(self.directory).parts[-1]
                if not os.path.isdir(save):
                    os.makedirs(save)
                saveas = pjoin(save, sname+f"_xy{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fig_xy.savefig(saveas)
                plt.close(fig_xy)
                saveas = pjoin(save, sname+f"_xz{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fig_xz.savefig(saveas)
                plt.close(fig_xz)
            figs.extend([fig_xy, fig_xz])
            steps.append(step)
        return figs, steps

    def plot_cs(self, ax=None, legend=True, vsig=False,
                thresholds=False, save=None, **plot_args):
        r"""Plot the cross sections models.

        Args:
            ax (:class:`~.matplotlib.axes.Axes` or None): Optional
                axes object to plot on top of, default is None.
                If ax is None, then a new figure and ax object
                will be created.
            legend (bool): Activate axes legend if true,
                default is True.
            vsig (bool): Plot :math:`\sqrt{\frac{2\epsilon}{m_{\rm e}}}\sigma(\epsilon)`
                rather than :math:`\sigma(\epsilon)`.
            thresholds (bool): if True, plot energy units in thresholds.
            save (str): Optional location of the directory to save the plot in.
            **plot_args: Optional arguments passed to Axes.plot().
        Returns:
            :class:`matplotlib.axes.Axes`: The axes object of the plot.
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        ax = self.cs.plot_cs(
            ax, legend, vsig, thresholds, **plot_args)

        if save:
            sname = Path(self.directory).parts[-1]
            if not os.path.isdir(save):
                os.makedirs(save)
            saveas = pjoin(save, sname+"_cs.pdf")
            print(f"Saving cross section plot in {saveas}", flush=True)
            fig.savefig(saveas)
            plt.close(fig)

        ax.grid(which='major', color='gray', linestyle='-', alpha=0.5)
        ax.grid(which='minor', color='gray', linestyle='-', alpha=0.3)

        return ax

    def _check_available_vdfs(self, steps="all", get_types=False):
        """Check what steps are available base on required criteria."""
        available = []
        # First, check directory for file data
        try:
            for efile in ls(self.directory):
                if "Dump" in efile:
                    if "_0_" not in efile: continue
                    available.append(int(efile.split("_")[-1].split(".")[0]))
        except FileNotFoundError:
            pass
        # Then check for already parsed data
        if self.vdfs is not None:
            p0_only = self.vdfs.ptype == 0
            available = sorted(list(set(
                available + self.vdfs[p0_only].step.unique().tolist()
            )))

        if not len(available) or steps == "all":
            return available # No data available

        # Choose which time steps to include based on `steps`
        if steps == "last":
            return [max(available)]
        elif steps == "first":
            return [min(available)]
        elif isinstance(steps, list):
            return [s for s in steps if s in available]
        elif isinstance(steps, int):
            return [s for s in [steps] if s in available]

    def _tabulate(self, df):
        """Add params to DataFrame."""
        # Log inputs into table format
        for p, v in self.tb_params.items():
            if isinstance(v, list):
                df[p] = None # Init column to avoid type error
                df.loc[:, p] = " ".join(str(t) for t in v)
            else:
                df.loc[:, p] = v
        # Log meta params into table format
        for p, v in self.hp.items():
            if isinstance(v, list):
                df.loc[:, p] = " ".join(str(t) for t in v)
            elif isinstance(v, dict):
                df.loc[:, p] = str(v)
            elif callable(v):
                df.loc[:, p] = v.__name__
            else:
                df.loc[:, p] = v
        return df

    def _log(self, directory=None):
        """Write down the hyper parameters into a file."""
        if not directory:
            directory = self.directory
        with open(pjoin(directory, self.logfile), "w") as f:
            dat = {"hp": self.hp, "tb": self.tb_params, "elapsed_time": self.elapsed_time,
                    "runtime_start": self.runtime_start, "runtime_end": self.runtime_end}
            a = copy.deepcopy(self.hp["indeck"])
            # Convert functions to function names
            if callable(a):
                dat["hp"]["indeck"] = a.__name__
            # Convert any numpy/pandas datatypes into python primitives
            to_primitive(dat)
            # Serialize into log file
            json.dump(dat, f)
            self.hp["indeck"] = a

    def to_pickleable(self):
        """Return a picklable version of this object."""
        # It is currently entirely picklable
        return self

    def get_directory(self):
        """Return the path of the current simulation."""
        return self.directory

    def read_log(self, logfile):
        """Read json file from simulation directory. Update
        the corresponding settings in the ThunderBoltz object.

        Args:
            logfile (str): The path name of the logfile to read.
        """
        lpath = pjoin(self.directory, logfile)
        if os.path.isfile(lpath):
            with open(lpath, "r") as f:
                dat = json.load(f)
                self.hp.update(dat["hp"]) # Update in case params added
                self.tb_params.update(dat["tb"])
                self.elapsed_time = dat["elapsed_time"]
                self.runtime_start = dat["runtime_start"]
                self.runtime_end = dat["runtime_end"]

    def compute_fit(self, x_, name=None):
        """Get the slope and associated error of the line of best fit. If
        x is a Dataframe, do so for each column."""
        is_series = isinstance(x_, pd.Series)
        # Rest index to normalize by DT instead of time
        i = x_.index[0]
        x = copy.deepcopy(x_).reset_index(drop=True)
        if is_series: x = x.to_frame()
        e = copy.deepcopy(x) # store errors 

        for col in x:
            t = x.index * self.DT 
            popt, pcov = np.polyfit(t, x[col], 1, cov=True)
            # Overwrite with slopes and errors, rescaled
            x[col] = popt[0]
            # Slope standard deviation is sqrt(Cov[0,0])
            e[col] = np.sqrt(pcov[0][0])

        # Return in same shape
        if is_series: return to_series(x), to_series(e)
        return x, e

def read(directory, read_cs_data=False):
    """Create a ThunderBoltz object by reading from a
    ThunderBoltz simulation directory.

    Args:
        directory (str): The directory from which to
            initialize the ThunderBoltz object.
        read_cs_data (bool):
            When set to true, the reader will look for cs_data,
            default is ``False``

    Returns:
        :class:`~.thunderboltz.ThunderBoltz`: The ThunderBoltz object with
        tabulated data if available.
    """
    tb = ThunderBoltz()
    tb.read(directory, read_cs_data=read_cs_data)
    return tb

def query_tree(directory, name_req=None, param_req=None,
              read_cs_data=False, callback=None, agg=True):
    r"""Walk a directory tree and search for ThunderBoltz
    simulation directories to read. Either return a list of
    :class:`~.thunderboltz.ThunderBoltz` objects, or a custom aggregation of
    the output data.

    Args:
        directory (str): The root path to search for ThunderBoltz
            data in.

        name_req (callable[str,bool]): A requirement on the file path
            names to be included in the query. The callable accepts
            the file path of a thunderboltz simulation directory and
            should return ``True`` if that directory is to be included
            in the query.

            e.g. ``name_req=lambda s: "test_type_1" in s``
            would return only data in a subfolder ``test_type_1``.

        param_req (dict): A requirement on the
            parameter settings of the ThunderBoltz calculations.
            The dictionary corresponding to simulation parameters
            that must be set by the read ThunderBoltz object.

            e.g. ``param_req={"Ered": 100, "L": 1e-6}`` would only
            return data from calculations with a reduced field of
            100 Td and a cell length of 1 :math:`\mu{\rm m}`.

        callback: (callable[:class:`ThunderBoltz`, Any]):
            A function that accepts a ThunderBoltz object and
            returns the desired data.

        agg: If ``callback`` is set, attempt to aggregate the data
            based on the data type:

            ===============================   ============================
            callback Return Type              Behavior
            ===============================   ============================
            :class:`pandas.DataFrame`         Frames will be
                                              concatenated row-wise and
                                              one larger DataFrame will be
                                              returned.

            list[:class:`pandas.DataFrame`]   A list of frames the same
                                              length of the return value
                                              will be returned. The frame
                                              at index ``i`` will contain
                                              the concatenated data
                                              from each simulation returned
                                              by ``callable(tb)[i]``.

            list[Any]                         A list of lists will be
                                              returned. The list at index
                                              ``i`` will contain a list
                                              of items returned from each
                                              call to ``callable(tb)[i]``.

            Any                               Return values will be
                                              returned in a list.
            ===============================   ============================

            If ``agg`` is set to False, always return a list of callback data
            without any concatenation.

    Returns:
        list[:class:`ThunderBoltz`], or :class:`pandas.DataFrame`, or list[:class:`pandas.DataFrame`], or list[list[Any]]:
        See ``agg`` option for behavior. Default return type
        is list[:class:`ThunderBoltz`].

    """

    tbs = [] # For return type 1
    df = pd.DataFrame() # For return type 2
    dfs = [] # For return type 3
    tb_dat = [] # For return type 4
    # Default return type is list of TB objects
    rtype = tbs
    for path, dirs, files in os.walk(directory):
        # Ensure this is a ThunderBoltz directory
        if "cross_sections" in dirs and "thunderboltz.out" in files:
            # Ensure file path is matches requirements
            if name_req:
                if not name_req(path): continue

            # Read in data for this ThunderBoltz obj.
            tb = read(path, read_cs_data=read_cs_data)
            if param_req:
                valid_params = True
                for k, v in param_req.items():
                    set_val = tb.get_sim_param(k)
                    valid_params = (set_val == v) and valid_params
                if not valid_params: continue

            # Append object once determined to be included
            tbs.append(tb)

            # If desired, extract data from ThunderBoltz objects.
            if not callback:
                continue
            dat = callback(tb)
            if isinstance(dat, pd.DataFrame) and agg:
                # Build one big data frame with all the data
                df = pd.concat((df, dat), ignore_index=True)
                rtype = df
            elif isinstance(dat, (list, tuple)):
                # Build a list of either lists or DataFrames
                n = len(dat)
                is_df_list = all(isinstance(item, pd.DataFrame) for item in dat)

                # Initialize if necessary
                if not len(dfs): dfs = n*[pd.DataFrame()]
                if not len(tb_dat): tb_dat = [[] for _ in range(n)]
                for i, item in enumerate(dat):
                    if is_df_list and agg:
                        dfs[i] = pd.concat((dfs[i], item),
                            ignore_index=True)
                        rtype = dfs
                    else:
                        tb_dat[i].append(item)
                        rtype = tb_dat
            # Otherwise just return the list of ThunderBoltz objects

    # Return the correct datatype
    return rtype

def plot_tree(path, series=["MEe", "mobN", "a_n"], name_req=None, param_req=None,
              save=None, stamp=["directory"]):
    """Query a tree of calculations and make a simple time series plot for each
    one.

    Args:
        series (list[str]):
            The y-parameters to plot onto the time series figure. If the
            element is a string, the corresponding parameter will be plotted
            if available in :class:`~.thunderboltz.parameters.OutputParameters`.
            If the element is a tuple, the first argument will be interpreted as the
            name of a new user defined parameter, and the second argument a user defined
            function that calculates it. The function must accept timeseries data
            and return a single series to be plotted.
        path (str): The root of the tree to plot calculations from.
        name_req (callable[str,bool]): A requirement on the file path
            names to be included in the query. The callable accepts
            the file path of a thunderboltz simulation directory and
            should return ``True`` if that directory is to be included
            in the query.

            e.g. ``name_req=lambda s: "test_type_1" in s``
            would return only data in a subfolder ``test_type_1``.

        param_req (dict): A requirement on the
            parameter settings of the ThunderBoltz calculations.
            The dictionary corresponding to simulation parameters
            that must be set by the read ThunderBoltz object.

            e.g. ``param_req={"Ered": 100, "L": 1e-6}`` would only
            return data from calculations with a reduced field of
            100 Td and a cell length of 1 :math:`\mu{\rm m}`.

        save (str): Option to save the plot to a file path.

        stamp (list[str]): Option to stamp the figure with the value of
            descriptive parameters, e.g. the field, or initial
            number of particles. See :class:`~.thunderboltz.parameters.TBParameters`
            and :class:`~.thunderboltz.parameters.WrapParameters`.
    """
    calcs = tb.query_tree(path, name_req=name_req, param_req=param_req)
    for calc in calcs: calc.plot_timeseries(
        series=series, stamp=stamp, save=save)

def to_primitive(d):
    """Recursively convert collections with non primitive types
    back to python natives."""
    if isinstance(d, dict):
        iterable = d.items()
    elif isinstance(d, list):
        iterable = enumerate(d)

    for k, v in iterable:
        if isinstance(v, dict) or isinstance(v, list):
            # Recurse mutable collections in place
            to_primitive(v)
        elif isinstance(v, np.ndarray):
            d[k] = v.tolist()
        elif isinstance(v, np.number):
            d[k] = infer_param_type(v)
        elif isinstance(v, pd.DataFrame):
            d[k] = v.to_numpy().tolist()

def infer_param_type(x):
    """Sequential type checking for parameter io."""
    types = [int, float, str]
    for t in types:
        try:
            return t(x)
        except:
            continue
    raise RuntimeError("No type satisfied by {x}")

def to_series(df):
    """Attempt to convert dataframe to series."""
    if len(df.columns) > 1:
        raise RuntimeError("Cannot convert multi-column DataFrame to Series.")
    else:
        return df.iloc[:,0]
