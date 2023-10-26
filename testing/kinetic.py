"""Unit testing for the pytb.kinetic module."""
import copy

import numpy as np
from numpy import vectorize as vec
import pandas as pd

import pytb
from pytb.kinetic import Process

##### UNIT TESTS #####

def test_auto_grid_simple():
    """Make sure simple auto grid is selected for a Heaviside operator
    cross section."""
    ics = lambda e: 1
    p = Process("Inelastic", threshold=1., cs_func=ics, name="A+B->A+C")
    assert len(p.data) == 4
    # Make sure it is chosen for a linear cs as well
    ics = lambda e: 1+0.01*e
    p = Process("Inelastic", threshold=1., cs_func=ics, name="A+B->A+C")
    assert len(p.data) == 4

def test_elastic_near_0():
    """Make sure the elastic cross section if sampled correctly."""
    ics = lambda e: 1
    p = Process("Elastic", cs_func=ics, name="A+B->A+C")
    en = "Energy (eV)"
    cs = p.data
    # There should only be unique energy points
    assert len(cs) == len(cs[en].unique()) == 2

def test_nonlinear_elastic():
    """Make sure an elastic cross section with a very small
    threshold-like shelf is properly sampled."""
    # To pass, this must be interpreted as float in Process caller,
    # even if 0 is returned
    ics = lambda e: 2e-20 if e > 1 else 0
    p = Process("Elastic", cs_func=ics)
    assert len(p.data) > p.nsamples

def test_auto_grid_dense():
    """Make sure dense grid is selected for a nonlinear cross section."""
    ics = lambda e: np.log(e)/e if e > 10 else 0
    p = Process("Inelastic", threshold=10., cs_func=ics, name="A+B->A+C")
    assert len(p.data) == 252

def test_cs_data_input():
    """Test input of tabulated cross section data."""
    cols = ["Energy (eV)", "Cross Section (m^2)"]
    vals = [
        [0.0000, 0.0],
        [1.0000, 0.0],
        [1.0001, 1.0],
        [1000000.0000, 1.0],
    ]
    df = pd.DataFrame(vals, columns=cols)
    p = Process("Inelastic", threshold=1., cs_data=df, name="A+B->A+C")
    cs1 = copy.deepcopy(p.data)
    dat = np.array(vals)
    p = Process("Inelastic", threshold=1., cs_data=dat, name="A+B->A+C")
    cs2 = copy.deepcopy(p.data)
    p = Process("Inelastic", threshold=1., cs_data=vals, name="A+B->A+C")
    cs3 = copy.deepcopy(p.data)
    # These should be the same
    assert (cs1 == cs2).all().all() and (cs2 == cs3).all().all()

def test_threshold_cutoff():
    """Cross section should be 0 at and below the threshold, even
    when tabulated data is passed in."""
    vals = [
        [0.0000, 0.0],
        [1.0000, 0.0],
        [1.0001, 1.0],
        [1000000.0000, 1.0],
    ]
    en = "Energy (eV)"
    csn = "Cross Section (m^2)"
    p = Process("Inelastic", threshold=2., cs_data=vals, name="A+B->A+C")
    assert p.data.loc[p.data[en] <= 2.0, csn].sum() == 0
