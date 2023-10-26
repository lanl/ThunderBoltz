"""Test pytb.input.CrossSections class."""
import copy
import os
from os.path import join as pjoin

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import pytb
from pytb import ThunderBoltz as TB
from pytb.input import CrossSections
from pytb.kinetic import Process
import testing
from testing.utils import setup_
from testing.utils import He_settings

##### UNIT TESTS #####

def test_add_process():
    """Make sure a sample added process gets properly
    integrated into the CrossSections object."""
    # Create process with name
    name = "A+B_v=0->A+B_v=1"
    process = Process("Inelastic", threshold=10.0,
        cs_func=lambda x: max(1-.01*x, 0), name=name)
    # Create CrossSections from CCC data.
    cs, _ = pytb.input.He_TB()
    N = len(cs.table)
    cs.add_process(process)
    assert len(cs.table) == N + 1
    assert name+".dat" in cs.data

def test_set_tracking():
    """Ensure that ion tracking works properly."""
    # Setup an N2 calculation.
    tb = TB(NS=200, indeck=pjoin("indecks", "ac_field"),
            directory=setup_())
    # Set one random reaction to be between particle 0
    rp = ["r1", "r2", "p1", "p2"]
    tb.cs.table.loc[3, rp] = 0
    tb.set_fixed_tracking()
    assert tb.cs.table[rp].max().max() == 1

def test_cs_too_many_species():
    """Ensure that a species index used in the cross section
    table which."""
    tb = TB(NS=200, **He_settings())
    tb.cs.table.loc[0,"p2"] = 2
    with pytest.raises(IndexError):
        tb.run()

def test_read_LXCat_cs():
    """Test parsing the LXCat text"""
    df = pytb.parsing.read_LXCat_cs(pjoin("lxcat", "assortment.txt"))
    keys = ["process_type", "threshold", "db", "target"]
    assert len(df.drop_duplicates(keys)[keys]) == 264

def test_from_LXCat():
    """Try running a program from LXCat input."""
    tb = TB(NS=200, directory=setup_())
    tb.cs.from_LXCat(pjoin("lxcat", "trinity_N2.txt"))
    tb.run()

def test_He_settings():
    """Make sure the He settings are loaded properly from pytb.He_TB"""
    tb = TB(NS=200, **He_settings())
    assert tb.tb_params["MP"][1] == 4

def test_plot_cs():
    """Test the plot created by CrossSections objects."""
    tb = TB(NS=200, **He_settings())
    tb.plot_cs()
    tb.plot_cs(vsig=True)
    tb.plot_cs(thresholds=True)
    tb.plot_cs(vsig=True, thresholds=True)

def _test_indexed_naming():
    """Make sure that naming is automatically handled
    when not provided alongside CS objects read in from
    ThunderBoltz simulation directories."""
    # TODO: do this for N2 gas
    pass
