"""Test thunderboltz.input.CrossSections class."""
import copy
import os
from os.path import join as pjoin

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import thunderboltz as tb
from thunderboltz import ThunderBoltz as TB
from thunderboltz.input import CrossSections
from thunderboltz.kinetic import Process
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
    cs, _ = tb.input.He_TB()
    N = len(cs.table)
    cs.add_process(process)
    assert len(cs.table) == N + 1
    assert name+".dat" in cs.data

def test_set_tracking():
    """Ensure that ion tracking works properly."""
    # Setup an N2 calculation.
    calc = TB(NS=200, indeck=pjoin("indecks", "ac_field"),
            directory=setup_())
    # Set one random reaction to be between particle 0
    rp = ["r1", "r2", "p1", "p2"]
    calc.cs.table.loc[3, rp] = 0
    calc.set_fixed_tracking()
    assert calc.cs.table[rp].max().max() == 1

def test_cs_too_many_species():
    """Ensure that a species index used in the cross section
    table which."""
    calc = TB(NS=200, **He_settings())
    calc.cs.table.loc[0,"p2"] = 2
    with pytest.raises(IndexError):
        calc.run()

def test_read_LXCat_cs():
    """Test parsing the LXCat text"""
    df = tb.parsing.read_LXCat_cs(pjoin("lxcat", "assortment.txt"))
    keys = ["process_type", "threshold", "db", "target"]
    assert len(df.drop_duplicates(keys)[keys]) == 264

def test_from_LXCat():
    """Try running a program from LXCat input."""
    calc = TB(NS=200, directory=setup_())
    calc.cs.from_LXCat(pjoin("lxcat", "trinity_N2.txt"))
    calc.run()

def test_He_settings():
    """Make sure the He settings are loaded properly from He_TB."""
    calc = TB(NS=200, **He_settings())
    assert calc.tb_params["MP"][1] == 4

def test_plot_cs():
    """Test the plot created by CrossSections objects."""
    calc = TB(NS=200, **He_settings())
    calc.plot_cs()
    calc.plot_cs(vsig=True)
    calc.plot_cs(thresholds=True)
    calc.plot_cs(vsig=True, thresholds=True)

def _test_indexed_naming():
    """Make sure that naming is automatically handled
    when not provided alongside CS objects read in from
    ThunderBoltz simulation directories."""
    # TODO: do this for N2 gas
    pass
