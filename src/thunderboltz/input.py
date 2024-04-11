import copy
import json
import os
from os.path import join as pjoin
from os import listdir as ls
import re
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import physical_constants as pc

from thunderboltz import parsing
from thunderboltz import data_path
from thunderboltz.kinetic import Process
from thunderboltz.plotting import styles

# Physical Constants
ME = 9.1093837e-31 # electron mass, kg
QE = 1.60217663e-19 # elementary charge C or J/eV
AMU_TO_KG = 1.6605e-27 # kg / amu
a_0 = 5.29177e-11 # Bohr radius in m
eVRyd = 13.6056980659 # eV per Rydberg

# Map LXCat process types to ThunderBoltz process types.
lxcat_map = {
    "ELASTIC": "Elastic",
    "EFFECTIVE": "Elastic",
    "ROTATIONAL": "Inelastic",
    "EXCITATION": "Inelastic",
    "IONIZATION": "Ionization",
    # "ATTACHMENT": Not yet implemented in ThunderBoltz
}

class CrossSections(object):
    """ThunderBoltz cross section set data type. Consists of
    a set of cross sections each with a file reference and a
    reaction table.

    Args:
        directory (str): The path to a ThunderBoltz simulation
            directory in which input files are to be written.
            Default is ``None``.
        input_path (str): The path to a set of ThunderBoltz
            input files from which input data can be read.
            The file structure should be something like:
            ::
               path/to/input_path
                |———indeck_file.in  <— The main ThunderBoltz indeck file.
                |———cross_sections  <— Cross section directory
                |   |———cs1.dat     <— ThunderBoltz-formatted cross section file.
                |   |———cs2.dat     <— ThunderBoltz-formatted cross section file.
                |   ...
                ...

        cs_dir_name (str): The name of the cross section directory.
        input_fname (str): The name of the main ThunderBoltz indeck
            file. Default is None, in which case the indeck will be
            searched for in `input_path` and must end with `.in`.
    """
    # TODO: Update this API so that cross section files and indeck files may be
    # separate if desired, but defaults to finding first compatible .in file and
    # cross_section directories

    def __init__(self, directory=None, input_path=None, read_cs_data=True,
            cs_dir_name="cross_sections", input_fname=None):
        #: Place for cs files in simulation dir
        self.cs_dir_name = cs_dir_name
        #: Input deck filename default
        self.input_fname = input_fname

        # Store data
        cols = ["csfile", "r1", "r2", "rtype", "B", "p1", "p2", "model_name", "params"]
        types = [str]+[np.int64]*2+[str]+[np.float64]+[np.int64]*2+[str]*2
        #: The reaction table, with columns
        self.table = pd.DataFrame({k: pd.Series([], dtype=dtype)
            for k, dtype in zip(cols, types)})

        #: Data tables for each cross section.
        self.data = {} # Store the cross section data.
        # Hyper parameters
        self.input_path = input_path #: Input path to default input data
        self.def_cs_dir = None
        self.def_indeck = None
        if input_path:
            # Read input from default cs data files
            self.read(input_path, read_cs_data=read_cs_data)
        # Otherwise, do not read/store any CS data

        if directory:
            self.find_infile()

            if input_path:
                # Create cs_dir for new calculation
                self.indeck = pjoin(directory, self.input_fname)
                self.cs_dir = pjoin(directory, self.cs_dir_name)
                # Make dir.
                if not os.path.isdir(self.cs_dir):
                    os.mkdir(self.cs_dir)

    def find_infile(self):
        """Look for indeck file in input_path.

        Raises:
            RuntimeError: If input_path is not set, if multiple indecks
            are found, or if no indecks are found.
        """
        if not self.input_path:
            raise RuntimeError("No input directory to search for indeck file in.")
        if self.input_fname is None:
            # Try to find one on type ".in"
            fs = [f for f in ls(self.input_path) if f.endswith(".in")]
            if not len(fs):
                raise RuntimeError("No input file found in indeck directory.")
            if len(fs) > 1:
                raise RuntimeError("Multiple input files found in indeck directory")
            self.input_fname = fs[0]

    def read(self, input_path, read_cs_data=True):
        """Read ThunderBoltz cross section data from a directory with a
        single input file and a set of cross section files.

        Args:
            input_path (str): The path to the directory
                with cross section data. If specified, CrossSections.input_path
                will be updated as well.
            read_cs_data (bool): If ``False``, only read the process header
                information, and not the actual cs data itself, default
                is ``True``.
        """
        # Make this input path as the reference input directory for this object
        self.input_path = input_path
        # Locate a '.in' file.
        self.find_infile()
        self.def_cs_dir = pjoin(input_path, self.cs_dir_name)
        self.def_indeck = pjoin(input_path, self.input_fname)

        # Reinitialize the reaction table
        cols = ["csfile", "r1", "r2", "rtype", "B", "p1", "p2", "model_name", "params"]
        types = [str]+[np.int64]*2+[str]+[np.float64]+[np.int64]*2+[str]*2
        self.table = pd.DataFrame({k: pd.Series([], dtype=dtype)
            for k, dtype in zip(cols, types)})
        # Read the input deck for cross section defaults.
        with open(self.def_indeck, "r") as f:
            for line in f.readlines():
                if "CS" not in line or "CC" in line:
                    continue
                l = line.split()[1:len(cols)]
                # set csfile relative to cs_dir
                l[0] = "/".join(l[0].replace("./", "").split("/")[-1:])
                lps = " ".join(line.split()[len(cols):])
                # Save to a data table
                items = l + [lps]
                row = {k: v for k, v in zip(cols, items)}
                self.table = pd.concat(
                    (self.table, pd.DataFrame([row])), ignore_index=True)
        # Recast types (they were read in as strings)
        self.table = self.table.astype({k: v for k, v in zip(cols, types)})

        if not read_cs_data:
            return

        # Read the input cross section files.
        for csfile in self.table.csfile.unique():
            with open(pjoin(self.def_cs_dir, csfile), "r") as f:
                # Loop through cross section data
                df = pd.DataFrame()
                for i, line in enumerate(f.readlines()):
                    e, cs = list(map(float, line.split()))
                    row = {"Energy (eV)": e, "Cross Section (m^2)": cs}
                    df = pd.concat((df, pd.DataFrame([row])), ignore_index=True)
                # Save data tables into a hashmap
                self.data[csfile] = df


    def write(self, directory=None):
        """Write cross section files into the simulation
        cross section directory.

        Args:
            directory (str): Option to write to a specific directory.
                If provided, ``CrossSections.cs_dir`` will be updated.
                Default is None.
        Raises:
            RuntimeError: If no directory is set or provided.
        """
        if directory:
            self.cs_dir = pjoin(directory, self.cs_dir_name)
        if not self.cs_dir:
            raise RuntimeError(
                "Nowhere to write, specify directory")
        for cs_file, cs_data in self.data.items():
            # cs_file = reaction + ".dat"
            cs_path = pjoin(self.cs_dir, cs_file.replace("./", ""))
            cs_dir = os.path.split(cs_path)[0]
            if not os.path.isdir(cs_dir):
                os.makedirs(cs_dir)
            with open(cs_path, "w") as f:
                f.write("\n".join("\t".join(str(s) for s in pair) for pair in
                    cs_data.values))

    def add_process(self, p):
        """Take a Process object and update the cross
        section data and cross section reaction table.

        Args:
            p (Process): The process object for a single
                type of interaction.
        """
        # Make sure cross section data is available for this process
        p.require_cs()
        cs_row = p.to_df()
        name = cs_row.csfile.values[0]
        cs_dat = copy.deepcopy(p.data)
        # Add descriptors to cross section table
        self.table = pd.concat(
            (self.table, cs_row), ignore_index=True)
        # Specify cross section data
        self.data[name] = cs_dat

    def add_processes(self, ps):
        """Add multiple cross sections to the reaction table.

        Args:
            ps (list[Process]): A list of process objects to add.
        """
        for p in ps: self.add_process(p)

    def add_differential_model(self, rtype, name, params=None):
        """Add a differential model to a certain type of process.

        Args:
            rtype (str): "Elastic", "Inelastic", or "Ionization",
                the broad collision process type.
            name (str): The name of the differential process model.
                Available built-in options for each ``rtype`` are:

                ========== ====================
                ``rtype``  ``name``
                ========== ====================
                Elastic    Park, Murphy
                Ionization Equal, Uniform
                ========== ====================

            params (list[float]): Optional list of parameters required
                by the differential model.
        """
        msk = self.table.rtype.str.contains(rtype)
        self.table.loc[msk, "model_name"] = name
        if params:
            self.table.loc[msk, "params"] = " ".join(str(p) for p in params)

    def set_fixed_background(self, fixed=True):
        """Set all particle conserving processes to have the
        `FixedParticle2` tag or not."""
        c = self.table
        # Make dynamic fist
        c.loc[:,"rtype"] = c.rtype.str.replace("FixedParticle2", "")
        if fixed:
            # Add fixed tags to non ionizing collisions.
            nion = ~c.rtype.str.contains("Ionization")
            c.loc[nion,"rtype"] += "FixedParticle2"

    def plot_cs(self, ax=None, legend=True, vsig=False,
                thresholds=False, **plot_args):
        r"""Plot the cross sections models.

        Args:
            ax (:class:`~.matplotlib.axes.Axes` or None): Optional
                axes object to plot on top of, default is ``None``.
                If ax is ``None``, then a new figure and Axes object
                will be created.
            legend (bool or dict): Activate axes legend if true,
                default is True. If a dictionary is passed, it is
                interpreted as arguments to :meth:`~.matplotlib.axes.Axes.legend`.
            vsig (bool): Plot :math:`\sqrt{\frac{2\epsilon}{m_{\rm e}}}\sigma(\epsilon)`
                rather than :math:`\sigma(\epsilon)` on the y-axis.
            thresholds (bool): if ``True``, plot
                :math:`\frac{\epsilon}{\epsilon_{\rm ion}} - 1`
                rather than :math:`\epsilon` on the x-axis.
            **plot_args: Optional arguments passed to Axes.plot().
        Returns:
            ax (matplotlib.axes.Axes): The axes object.
        """
        # Column name aliases
        en = "Energy (eV)"
        csn = "Cross Section (m^2)"

        # Generate axes object if not provided
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 10))

        # Separate label from plot args if plotting multiple cross sections
        label = ""
        if "label" in plot_args:
            label = "_" + plot_args["label"]
            del plot_args["label"]
        # Loop through the processes
        for (csfile, B), color in zip(
                self.table[["csfile", "B"]].values, styles.DCC):
            pname = csfile.replace(".dat", "") + label
            cs_dat = self.data[csfile]
            x = cs_dat[en]/B - 1 if thresholds else cs_dat[en]
            cs = (np.sqrt(2*cs_dat[en]*QE/ME)*cs_dat[csn]
                  if vsig else cs_dat[csn])
            ax.plot(x, cs, label=pname, c=color, **plot_args)

        # Setup axes settings
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\epsilon$ (eV)")
        if thresholds:
            ax.set_xlabel(r"$\frac{\epsilon}{\epsilon_{\rm ion}} - 1$")
        ax.set_ylabel(r"$\sigma(\epsilon)$")
        if vsig:
            ax.set_ylabel(r"$\sqrt{\frac{2\epsilon}{m_{\rm e}}}\sigma(\epsilon)$")
        if legend:
            ax.legend(fontsize=12)
        if isinstance(legend, dict):
            lg_args = {"fontsize": 12}
            lg_args.update(legend)
            ax.legend(**lg_args)

        return ax

    def get_deck(self):
        """Return the string formatted cross section table portion of the
        ThunderBoltz indeck."""
        s = ""
        for row in self.table.values:
            cspath = pjoin(self.cs_dir_name, row[0])
            description = " ".join(str(a) for a in row[1:] if not pd.isna(a)).rstrip() + "\n"
            s += f"CS {cspath} {description}"
        return s

    def from_LXCat(self, fname):
        """Load cross section data from an LXCat .txt file"""
        if not os.path.exists(fname):
            raise IOError(f"LXCat file {fname} does not exist.")
        if self.data or len(self.table):
            raise RuntimeError(f"LXCat file was passed when cross section "
                                "data already exists")
        df = parsing.read_LXCat_cs(fname)
        # Loop through each process
        for process_string, pdf in df.groupby("process_string"):
            # Derive the names of each process from the process string
            name = process_string.replace(" ", "")
            # TODO: implement LXCat species interpretation for TB species indexing
            # parse_process_string(name)
            # Make sure the threshold matches
            if not len(pdf.threshold.unique()) == 1:
                raise RuntimeError("Process string does not uniquely identify process")
            cs_dat = (pdf[["Energy (eV)", "Cross section (m2)"]].copy()
                .sort_values("Energy (eV)")
                .reset_index(drop=True)
                .rename(columns={"Cross section (m2)": "Cross Section (m^2)"})
            )
            # Create thunderboltz.kinetic.Process to enfore threshold rules
            ptype = pdf.process_type.unique()[0]
            B = pdf.threshold.unique()[0]
            B = 0.0 if pd.isna(B) else B
            if ptype not in lxcat_map:
                raise NotImplementedError(f"There is no current {ptype} "
                        "process implemented in ThunderBoltz")
            self.add_process(Process(lxcat_map[ptype],
                0, 1, 0, 1, threshold=B, cs_data=cs_dat, name=name))

        return self # So it's easier to load directly from LXCat

### Automated cross section constructions.
def He_TB(n=4, egen=True, analytic_cs=True, eadf="default", ECS=None,
          nsamples=250, mix_thresh=300., fixed_background=True):
    """Generate parameterized He cross section sets in the ThunderBoltz format.
    The data is from Igor Bray and Dmitry V Fursa 2011
    J. Phys. B: At. Mol. Opt. Phys. 44 061001.

    Args:
        n (int): Include CCC excitation processes from the ground state to
            (up to and including) states with principle quantum number ``n``.
        egen (bool): Allow secondary electron generation for the ionization model.
        analytic_cs (str or bool): use either tabulated data, analytic fits,
            or a mix of both. Options are ``False``, ``True``, or ``"mixed"``.
        eadf (str) Elastic angular distribution function model. Options are
            ``"default"``, or ``"He_Park"``.
        ECS (str or None): The total elastic cross section model. Options are
            ``"ICS"`` or ``"MTCS"``, default is ``"ICS"`` if an anisotropic
            angular distribution function is used and ``"MTCS"`` if an isotropic
            angular distribution function is used.
        nsamples (int): The number of tabulated cross section values for analytic
            sampling.
        mix_thresh (float): If ``analytic_cs`` is ``"mixed"``, use numerical data
            at energies lower than this threshold value (in eV), and use analytic
            data at higher energies.
        fixed_background (bool): Flag to append "FixedParticle2" to each of the
            reaction types in the indeck.

    Returns:
        Tuple[dict,dict]: The :class:`~.thunderboltz.CrossSections` object for
            Helium and the dictionary of ThunderBoltz parameters suitable
            for the cross section model.
    """

    # Determine elastic flag
    elastic_only = False
    if n==0:
        elastic_only = True
        n = 1
    # Select elastic integrated cross section type
    if ECS is None:
        ECS = "MTCS" if eadf == "default" else "ICS"
    # Load in cross sections from CCC data
    cmod = He_CCC(n, ECS=ECS)
    if analytic_cs:
        # Use continuous fits
        cmod_aly = make_He_aly_cs(n, ECS, nsamples=nsamples)
        if analytic_cs == "mixed":
            cmod_aly = split_cs_model(cmod, cmod_aly, mix_thresh, ECS)
        # overwrite original mod
        cmod = cmod_aly

    # Make sure threshold format is correct
    cmod = fix_thresh(cmod)

    # Create CS object from cmod dict
    itype = "Ionization" if egen else "IonizationNoEGen"
    tbcs = mod_to_tbcs(cmod, ion_type=itype, elastic_only=elastic_only)
    tbcs.set_fixed_background(fixed_background)

    # Specific TB settings for He background gas
    # These cannot be overwritten
    tbp = {
        "SP": 2,
        "NP": [100000, 10000],
        "TP": [0.0, 0.0259],
        "VV": [0.0, 0.0],
        "QP": [-1.0, 0.0],
        "MP": [5.4857e-4, 4.0],
        "FV": [1000, 150000, 0],
    }

    return tbcs, tbp

# Helper functions
def mod_to_tbcs(mod, ion_type="Ionization", elastic_only=False):
    """Convert dictionary cross section format to CrossSections Object."""
    sind = {"r1": 0, "r2": 1, "p1": 0, "p2": 1}
    # Elastic CS, assume fixed particles
    cmod = mod["cross_sections"]
    elastic_type = cmod["elastic"].columns[1]
    # cs_data init with elastic data
    df = pd.DataFrame([dict(csfile=f"{elastic_type}.dat", B=0.,
        rtype="ElasticFixedParticle2", **sind)])
    if not elastic_only:
        # Excitation cross sections
        for exc, B in mod["exc_thresh"].items():
            dfrow = pd.DataFrame([dict(csfile=f"{exc}.dat", B=B,
                rtype="InelasticFixedParticle2", **sind)])
            df = pd.concat((df, dfrow), ignore_index=True)
        # Ionization cross section
        B = mod["ion_thresh"]
        ion_row = pd.DataFrame([dict(csfile="ionz.dat", B=B, rtype=ion_type, **sind)])
        df = pd.concat((df, ion_row), ignore_index=True)

    # Convert type
    df = df.astype({"r1": int, "r2": int, "B": float, "p1": int, "p2": int})
    # Reorder columns
    df = df[["csfile", "r1", "r2", "rtype", "B", "p1", "p2"]].copy()

    # Create cs data dictionary
    # Elastic data
    cs_dat = {f"{elastic_type}.dat": cmod["elastic"].rename(
        columns={elastic_type: "Cross Section (m^2)"})}

    # Excitation data
    if not elastic_only:
        for name, exc_dat in cmod.items():
            if name in ["elastic", "ionization"]:
                continue
            c2 = exc_dat.columns[1]
            cs_dat.update({f"{name}.dat": exc_dat.rename(
                columns={c2: "Cross Section (m^2)"})})

        # Ionization data
        cs_dat.update({f"ionz.dat": cmod["ionization"].rename(
            columns={"TICS (m^2)": "Cross Section (m^2)"})})

    # Create CrossSections object
    cs = CrossSections()
    cs.table = pd.concat((cs.table, df), ignore_index=True)
    cs.data = cs_dat
    return cs

def He_CCC(nmax=3, ECS="MTCS"):
    """Helium package cross section model from CCC data.
    Includes excitation cross sections up to nmax"""
    # Create tables
    # Import excitation cross sections and threshold energies.
    ccc_path = pjoin(data_path, "CCC")
    tcs = pd.read_csv(pjoin(ccc_path, "He_DCS_tcs.csv"))
    potl = pd.read_csv(pjoin(ccc_path, "He_DCS_potl.csv"))

    # Manipulate extract for bolsig.
    potl = potl[["Energy (eV)"] + [c for c in
        potl.columns if "(cm^2)" in c]]
    potl = convert(potl, "cm^2", "m^2", drop=True)
    tcs = convert(tcs, "cm^2", "m^2")

    # Threshold energy from sample DCS file.
    be = pd.read_csv(pjoin(ccc_path, "He_thresholds.csv"), header=None,
            names=["transition", "Energy (Ryd)"])
    be = convert(be, "Ryd", "eV")

    # ELASTIC
    eltable = potl[["Energy (eV)", ECS]]
    # Calculate electron-target mass ratio
    He_m = 2*(pc["electron mass"][0]) + pc["alpha particle mass"][0]
    etmr = pc["electron mass"][0] / He_m

    # IONIZATION
    iotable = tcs[["Energy (eV)", "TICS (m^2)"]]
    io_thr = -be.loc[be.transition=="s1S", "Energy (eV)"].values[0]

    # EXCITATION
    # Select the desired excitation states.
    n_states = limit_excitations(potl, nmax)
    n_thresholds = exthresh_from_BE(be, n_states)
    # Binding energies obtained from CCC data: Igor Bray and Dmitry V Fursa 2011 J. Phys. B: At. Mol. Opt. Phys. 44 061001

    # Assemble and name cross section tables.
    cs_tables = {
        "elastic": eltable,
        "ionization": iotable,
    }
    cs_tables.update({k: potl[["Energy (eV)", k]]
        for k in n_states})

    # Assemble model with all info needed for BOLSIG.
    model = {
        "format": "TABULAR", # Indicate that cross sections are stored as DataFrames.
        "name": "He_CCC_n<={}".format(nmax),
        "cross_sections": cs_tables,  # Columns of energy, CS for every interaction.
        "exc_thresh": n_thresholds, # threshold for all inelastic processes.
        "ion_thresh": io_thr, # eV
        "etmr": etmr,
    }
    return model

def elastic_mtcs_LANL(e, sigt):
    """Takes incident energy and ICS(e) function and returns the
    momentum transfer cross section for He.

    This functional form is obtained from Park et al.
    https://iopscience.iop.org/article/10.1088/1361-6595/ac781f
    """
    ps = [0.283, 0.667, 0.0307, 16.971, 6.59, 29.02, 0.0258, 0.295, 0.00328, 1.794]
    a, laf, lab = ps[:4], ps[4:7], ps[7:]
    # Derived parameters
    C = a[0] + a[1]*np.tanh(a[2]*(e - a[3]))
    etaF = (laf[0] / (e + laf[1])) + laf[2]
    etaB = (lab[0] / (e + lab[1])) + lab[2]
    lF = (etaF * (etaF + 1)) / np.pi
    lB = (etaB * (etaB + 1)) / np.pi
    mF = (2*etaF + 1)
    mB = (2*etaB + 1)
    # Evaluation
    Phi = lambda l, m, n: ((l/(n**2)) * (-1 + ((m+n)/(m-n)) + np.log((m-n)/(m+n))))
    sigm = 2*np.pi*sigt*(C*Phi(lF, mF, -1) + (1-C)*(Phi(lB, mB, 1)))
    return sigm

def elastic_ics_LANL(e):
    """Takes incident energy and computes ICS (m^2) using fitting function.

    This functional form is obtained from Park et al.
    https://iopscience.iop.org/article/10.1088/1361-6595/ac781f
    """
    b = [1.47183e-18, 6.23206e-18, -1.30586e-18, 7.52062e-1,
         1.92954e-1, 2.42524, 4.49405, 1.25050e2]
    t1 = b[0] + (b[1]/(e+1)) + (b[2]/(e**2 + 1))
    t2 = (e**b[3] + np.exp(-b[4]*e))**b[5]
    t3 = (e+1)**3 + b[6]*(e+1)**2 + b[7]*(e+1)
    return t1*(t2/t3)

def fix_thresh(mod):
    """Fix the threshold value such that cross sections are 0
    at their threshold."""

    en = "Energy (eV)"
    csn = "cross_sections"
    el = mod[csn]["elastic"]
    ECS = el.columns[1]
    # Check elastic cross section to see if it has data at 0 eV.
    if el[en].values[0] != 0.:
        mod[csn]["elastic"] = pd.concat((pd.DataFrame({en: [0.], ECS:
            el.iloc[[0], 1]}), el), ignore_index=True)
    # Fix thresh on ionization cross section.
    ionn = "TICS (m^2)"
    B = mod["ion_thresh"]
    ion = mod[csn]["ionization"]
    mod[csn]["ionization"] = pd.concat((pd.DataFrame({en: [0.0,
        B], ionn: [0.0, 0.0]}), ion[ion[en] > B]), ignore_index=True)

    # Fix thresh on excitation cross sections.
    for proc in mod[csn]:
        if proc in ["elastic", "ionization"]:
            continue
        thr = mod["exc_thresh"][proc]
        inel = mod[csn][proc]
        mod[csn][proc] = pd.concat((pd.DataFrame({en: [0.0,
            thr], proc: [0.0, 0.0]}), inel[inel[en] > thr]), ignore_index=True)
    return mod


def make_He_aly_cs(n=3, ECS="ICS", nsamples=248):
    """Construct a grid from analytic He cross section fits proposed by
    Ralchenko et al. https://doi.org/10.1016/j.adt.2007.11.003."""

    # Get thresholds
    mod = He_CCC(n, ECS=ECS)
    exth = mod["exc_thresh"]

    # Named energy column
    en = "Energy (eV)"
    ionn = "TICS (m^2)"

    # Read in cs parameters from Ralchenko
    csp = pd.read_csv(pjoin(data_path, "Ralchenko_He_CS.csv"))
    # Only use transitions from ground state
    csp = csp[csp.i == "s1S"].copy()

    def sample_exc(name, e):
        """Sample excitation cross section with 'name' at energy 'e' (thresholds)"""
        # params
        rp = csp[csp.f == name].copy()
        a = rp.loc[:,"A1":"A6"].values[0]
        # Employ Ralchenko model
        if rp.table.values[0] == "dipole-allowed":
            return (a[0]*np.log(e) + a[1] + a[2]/e + a[3]/e**2 + a[4]/e**3)*(e + 1)/(e+a[5])
        elif rp.table.values[0] == "dipole-forbidden":
            return (a[0] + a[1]/e + a[2]/e**2 + a[3]/e**3) * e**2/(e**2 +a[4])
        elif rp.table.values[0] == "spin-forbidden":
            return (a[0] + a[1]/e + a[2]/e**2 + a[3]/e**3) / (e**2 + a[4])

    def sample_exc_set(name):
        """Sample a set of energies (eV) of excitation cross section with transition 'name'"""
        # Get threshold energy
        B = exth[name]
        # Make energies in eV
        es = np.logspace(-4, np.log10(10**6 - B), nsamples) + B
        # Calculate omega with energy in thresholds
        omega = np.vectorize(lambda e: sample_exc(name, e))(es/B)
        # Add pre factor to convert to cross section in m^2
        prefac = np.pi*a_0**2*(eVRyd/es)
        # Create DataFrame for this data
        dat = pd.DataFrame({en: es, name: prefac*omega})
        # Add extrapolation near 0 energy
        return pd.concat((pd.DataFrame({en: [0, B], name: [0, 0]}),
                          dat[dat[en] > B]), ignore_index=True)

    # Update elastic model
    # Add elastic model, include 0 eV data point
    es = np.concatenate(([0], np.logspace(-4, 6, nsamples)))
    el_df = pd.DataFrame({en: es, ECS: elastic_ics_LANL(es)})
    if ECS=="MTCS":
        # Use ICS to compute MTCS from LANL model
        el_df[ECS] = elastic_mtcs_LANL(es, el_df[ECS])
    # Add one point at 0 incident energy with the same value of cross section as first NZ point
    mod["cross_sections"]["elastic"] = el_df

    # Update ionization model
    a = csp.loc[csp.table=="ionization", "A1":"A6"].values[0]
    # Make a[1] = A1 for clarity
    a = np.concatenate(([0], a))
    B = mod["ion_thresh"]
    # Do not include 0 eV data point
    es = np.logspace(-4, np.log10(10**6 - B), nsamples) + B
    TICS = (1e-17/(B*es))*(a[1]*np.log(es/B) + sum(a[i]*(1-B/es)**(i-1) for i in range(2,7)))
    ion_df = pd.DataFrame({en: es, ionn: TICS})
    # Add a point at 0,0 and
    mod["cross_sections"]["ionization"] = pd.concat((
        pd.DataFrame({en: [0, B], ionn: [0, 0]}), ion_df[ion_df[en]>B]), ignore_index=True)

    # Update excitation cross sections.
    for name in mod["cross_sections"]:
        if name in ["ionization", "elastic"]:
            continue
        mod["cross_sections"][name] = sample_exc_set(name)

    return mod

def split_cs_model(m1, m2_, thresh, ECS):
    """Split two cross section models such that `m1` provides energies lower
    than `thresh` and `m2` provides energies higher than `thresh`."""
    en = "Energy (eV)"
    # Copy right mod
    m2 = copy.deepcopy(m2_)

    def split_cs(key):
        col = key
        if key == "elastic":
            col = ECS
        elif key == "ionization":
            col = "TICS (m^2)"
        m1e = m1["cross_sections"][key]
        m1cut = m1e[m1e[en] < thresh].copy()
        m2e = m2["cross_sections"][key]
        m2cut = m2e[m2e[en] > thresh].copy()
        m2["cross_sections"][key] = pd.concat((m1cut, m2cut),
                ignore_index=False).reset_index(drop=True)
    # Swap elastic model
    split_cs("elastic")
    split_cs("ionization")
    for proc in m2["cross_sections"]:
        if proc in ["elastic", "ionization"]:
            continue
        split_cs(proc)

    return m2

def parse_process_string(s):
    """Find species stings from process string."""
    # remove unnecessary eV tags
    s = re.sub(r"\([^()]*eV\)", "", s)

def limit_excitations(df, nmax=3):
    """Ignore transitions with principle quantum numbers
    above `nmax`."""
    # Isolate relevant excitations from CCC He data
    complete_states = [c for c in df.columns
        if not df[c].isnull().values.all()]
    # singlet and triplet states.
    st_states = [s for s in complete_states if s[0] in ["s", "t"]]
    # Quantum number of allowed transitions.
    n = np.arange(2, nmax+1)
    n_states = [s for s in st_states if int(s[1]) in n]
    return n_states

def exthresh_from_BE(be, states):
    """Create mapping from state transition name to excitation threshold"""
    smap = lambda s: s in states
    gs_be = be.loc[be.transition=="s1S", "Energy (eV)"].values
    exthresh = be[be.transition.map(smap)].copy()
    exthresh.loc[:, "thresh"] = exthresh["Energy (eV)"] - gs_be
    return {k: v for k, v in exthresh[["transition", "thresh"]].values}

def convert(df, u1, u2, inv=False, drop=False, add=False):
    """Convert easily between units with a labeled DataFrame.
    Columns in the format `<name> (<unit>)` will be converted."""
    umap = {
        # Area in one m^2
        "m^2": 1,
        "m2": 1,
        "Angs^2": 1e20,
        "10^-22m^2": 10e22,
        "cm^2": 100**2,
        "10^-16cm^2": 1e16*100**2,
        "a_0^2": 1.8897259e10**2,
        "Mb": 1e22,
        "kb": 1e25,
        # Energy in one eV
        "eV": 1.,
        "Ryd": 1/13.6056980659,
    }
    for col in df.columns:
        if u1 in col:
            if drop:
                # Convert but remove the unit tag.
                ncol = col.replace(" "+col.split()[1], "")
            else:
                # Change unit tag to new unit
                ncol = col.replace(u1, u2)
            if inv:
                # Convert inverse units
                conversion = umap[u1]/umap[u2]
            else:
                # Convert units
                conversion = umap[u2]/umap[u1]
            if add:
                # Make a new column
                df.loc[:, ncol] = df[col] * conversion
            else:
                # Rename the old column
                df.loc[:, col] *= conversion
                df = df.rename(columns={col: ncol})
    return df
