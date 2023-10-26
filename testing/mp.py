"""ThunderBoltz multiprocessing tests."""
import multiprocessing as mp
import os
from os.path import join as pjoin
from os.path import expandvars as expand
from pathlib import Path
import pickle
import shutil
import time

import numpy as np
import pandas as pd

import pytb
from pytb import ThunderBoltz as TB
from pytb.parallel import DistributedPool
from pytb.parallel import SlurmManager
import testing
from testing.utils import setup_
from testing.utils import He_settings
from testing.utils import TEST_OUTPUT_DIR

def test_simple_mp():
    """Run a simple He calculation with different inputs
    across all the cores of the computer."""
    # Get number of cores
    cores = mp.cpu_count()
    # Specify a set of reduced field values
    Efields = 10 + np.arange(cores)*10
    # Create a base TB object
    tb = TB(NS=1001, **He_settings())
    # Open up a multiprocessing context with this ThunderBoltz object
    with DistributedPool(tb) as pool: # Loop through the conditions
        for field in Efields:
            path = setup_(subdir=f"{field}vpm")
            # Send any updates through submit
            pool.submit(directory=path, E=field)
    # Exit from with block will wait for all processes to finish
    p = pjoin(*Path(path).parts[:-1])
    tbs = pytb.tb.query_tree(p)
    tss = pd.concat([tb.get_timeseries() for tb in tbs], ignore_index=True)
    # Ensure all the steps of all processes get collected from reader
    assert len(tss) == 11*cores

def test_slurm_submission():
    """Run a simple slurm script through the SlurmManager."""
    # Run a mockup of this test if there is no slurm + compatible cpu count on this device
    if shutil.which("sbatch"):
        # (this assumes scratch directory environment variable is set)
        mock = False
        base_dir = expand("$SCR")
    else:
        mock = True
        base_dir = TEST_OUTPUT_DIR
    path = setup_(base_directory=base_dir)
    # Create base ThunderBoltz Object
    tb = TB(NS=201, indeck=pytb.input.He_TB)
    # Setup (5/4 * ncores) calculations (e.g. require 2 nodes, 1 with 36 cores, # 1 with 9 cores)
    ncores = mp.cpu_count()
    njobs = int(5/4 * ncores)
    Efields = 100 + np.arange(njobs)*10
    # Create base directory for calculations
    # Configure SLURM parameters
    slurm_options = {
        "account": "xd",
        "time": 1, # 1 minutes
        "job-name": "test_slurm",
        "ntasks-per-node": ncores, # Use all cores
        "qos": "debug",
        "reservation": "debug",
    }
    # Add modules required for compute nodes
    mods = ["python", "gcc"]
    with SlurmManager(tb, path, modules=mods, mock=mock, **slurm_options) as slurm:
        for field in Efields:
            subpath = setup_(subdir=f"{field}vpm", base_directory=base_dir)
            slurm.submit(directory=subpath, E=field)
    # Exiting the slurm context will call sbatch for each set of tasks
    if not mock:
        # Wait for slurm to finish
        while slurm.has_active():
            time.sleep(1)

    # Read the resulting calculations
    tbs = pytb.tb.query_tree(path)
    tss = pd.concat([tb.get_timeseries() for tb in tbs], ignore_index=True)
    fields = sorted([tb.tb_params["E"] for tb in tbs])
    # Ensure all the steps of all processes get collected from reader
    assert np.all(fields == Efields)
    assert len(tss) == 3*njobs
