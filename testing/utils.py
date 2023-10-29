"""Helper methods for unit test functions."""
import inspect
import os
from os.path import join as pjoin
from os.path import expandvars as expand
import shutil

import thunderboltz as tb
from thunderboltz import ThunderBoltz as TB
from thunderboltz.parallel import SlurmManager

TEST_OUTPUT_DIR = pjoin("testing", "output")

def setup_(i=0, subdir=None, base_directory=None, dry=False):
    """Call from a run function such that it makes a directory
    in simulation/output with the name of that test function. Increment
    `i` if this function is being wrapped in another intermediary."""
    if base_directory is None:
        base_directory = TEST_OUTPUT_DIR
    # Create the path to the test file
    fpath = pjoin(base_directory, inspect.stack()[1+i][3])
    if subdir is not None:
        fpath = pjoin(fpath, str(subdir))
    # Just return path if dry run
    if dry: return fpath
    # Overwrite directory if necessary
    if os.path.isdir(fpath):
        shutil.rmtree(fpath)
    os.makedirs(fpath)
    return fpath

def array_test(key, vals, **fixed):
    """Take a parameter name and several values, and run 100 steps
    of ThunderBoltz, ensuring the error stack is empty and data
    has no NaNs."""
    for v in vals:
        calc = TB(L=1e-7, NS=101, NP=[100, 10],
                **He_settings(1), **{key:v}, **fixed)
        calc.run()
        # Read files
        calc.get_timeseries()
        # error stack should be empty
        assert not calc.err_stack

def He_settings(i=0):
    """Use automated helium input deck with auto file writing.
    Increment `i` if this function is being wrapped in another
    intermediary."""
    return {"directory": setup_(1+i), "indeck": tb.input.He_TB}

def slurm_debug(calc, path=None, nodes=1, time=60, **kwargs):
    """Setup a default debug slurm context"""
    if path is None:
        path = setup_(1, base_directory=expand("$SCR"))
    slurm_options = {
        "time": time,
        "job-name": "slm_debug",
        "ntasks-per-node": nodes,
        "qos": "debug",
        "reservation": "debug",
    }
    mods = ["python", "gcc"]
    return SlurmManager(calc, path, mods=mods, **slurm_options, **kwargs)
