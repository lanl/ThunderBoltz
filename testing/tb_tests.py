"""Unit testing for thunderboltz package."""
import os
from os.path import join as pjoin
from os.path import expandvars as expand
import shutil
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import thunderboltz as tb
from thunderboltz import ThunderBoltz as TB
from thunderboltz import query_tree
from thunderboltz.kinetic import Process
import testing
from testing.utils import array_test
from testing.utils import setup_
from testing.utils import He_settings
from testing.utils import slurm_debug

def test_get_timeseries():
    """Simplest TB setup."""
    calc = TB(NS=300, **He_settings())
    calc.run()
    ts = calc.get_timeseries()
    assert len(ts) == 4

def test_short_run():
    """One line of time series."""
    calc = TB(NS=10, **He_settings())
    calc.run()
    ts = calc.get_timeseries()
    assert len(calc.counts) == 1

def test_dryrun():
    """Run setup operations without running simulation."""
    calc = TB(NS=200, **He_settings())
    calc.run(dryrun=True)
    assert calc._proc == None

def test_proper_thresh():
    """Make sure threshold from auto indeck is correct."""
    calc = TB(NS=200, **He_settings())
    assert round(calc.cs.table.loc[11, "B"], 3) == 24.449

def test_stdout():
    """Run simple setup with stdout."""
    calc = TB(NS=200, out_file="_stdout_", **He_settings())
    calc.run()

def test_user_src():
    """Test compilation from repo level source directory."""
    calc = TB(NS=200, **He_settings())
    calc.run(src_path=pjoin("src", "thunderboltz", "cpp"))
    ts = calc.get_timeseries()
    assert len(ts) == 3

def _test_user_bin():
    """Test simulation from user input binary file."""
    calc = TB(NS=200, **He_settings())
    calc.run(bin_path=pjoin("bin", "thunderboltz.bin"))
    ts = calc.get_timeseries()
    assert len(ts) == 3

def test_too_many_part():
    """Make sure too many particles registers into outfile."""
    calc = TB(NS=2001, OS=1000, NP=[300000, 10], DT=1e-9, egen=True,
            # live=True,
            # live_rate=True,
            MEM=0.025, **He_settings())
    calc.run()
    calc.get_timeseries()
    assert "particle overload" in calc.err_stack

def test_ss_not_enough_steps():
    """Try to get steady state statistics with very few steps.
    Ensure that this raises a RuntimeError."""
    calc = TB(NS=200, **He_settings())
    calc.run()
    # Make sure a runtime error is raised
    with pytest.raises(RuntimeError):
        sst = calc.get_ss_params()

def test_ss_enough_steps():
    """Try to get steady state statistics on a run
    with very few steps"""
    calc = TB(NS=1201, NP=[10, 2], **He_settings())
    calc.run()
    sst = calc.get_ss_params(fit=False)
    assert len(sst) == 1

def test_autostep():
    """Test the autostep feature."""
    Ered = 10
    calc = TB(NS=200, Ered=Ered, pct_ion=10, DE=0.1, autostep=True, **He_settings())
    n = calc.tb_params["NP"][1]/calc.tb_params["L"]**3
    # Enforce the conversion E[V/m] = 1e-21*(E/n)[Td]*n[m^-3]
    assert abs(calc.tb_params["E"]) == 1e-21*Ered*n

def test_N2_time_dependent_indeck():
    """Ensure loading from ThunderBoltz input files works
    correctly."""
    calc = TB(NS=200, indeck=pjoin("indecks", "ac_field"),
            directory=setup_())
    calc.run()
    calc.get_timeseries()

def test_fixed_bg():
    """Test a fixed background gas, ensure particles mean energies are
    unchanging"""
    calc = TB(NS=301, NP=[100, 10], **He_settings())
    calc.run()
    calc.get_timeseries()
    pt = calc.particle_tables
    # Ensure the mean energy of neutrals does not change
    assert pt[1]["Mi"].diff().abs().sum() == 0

def test_dynamic_bg():
    """Test a dynamic background gas, ensure particles mean energies are
    changing"""
    calc = TB(fixed_background=False, L=1e-8, NS=501,
            NP=[1000, 1], **He_settings())
    calc.run()
    calc.get_timeseries()
    pt = calc.particle_tables
    # Ensure the mean energy of neutrals does change
    assert pt[1]["Mi"].diff().abs().sum() != 0

def test_eadf_options():
    """Ensure all the elastic angular differential functions run"""
    # All the valid models should work
    array_test("eadf", ["He_Park", "default"])
    # An unimplemented model should fail
    with pytest.raises(KeyError):
        array_test("eadf", ["The_Perfect_Model"])

def test_eesd_options():
    """Ensure all the electron impact ionization electron energy sharing
    energy distribution functions run."""
    # All the valid models should work
    array_test("eesd", ["equal", "uniform", "default"],
        TP=[30, .0235])

    # An unimplemented model should fail
    with pytest.raises(KeyError):
        array_test("eesd", ["The_Perfect_Model"])

def test_build_from_Process():
    """Ensure that ThunderBoltz runs when its cross section data
    is built using the thunderboltz.kinetic module."""
    # Make single cross section
    cs_func = lambda e: 1e-20*max(10 - .001*e, 0)
    p = Process("Elastic", cs_func=cs_func, name="procA")
    calc = TB(out_file="_stdout_", NS=200, NP=[1000, 1], **He_settings())
    calc.cs.add_process(p)
    calc.run()

def test_reduced_params():
    """Ensure density-reduced field parameters are calculated with the
    proper units."""
    # First with two species
    # Default charge and mass parameters
    p = {"QE": [-1, 0], "ME": [1, 1]}
    L = 1e-6
    calc = TB(NP=[1, 2], L=L, Ered=2, Bred=[0, 2.5, 0], **p)
    n = 2/L**3
    # Enforce the conversion E[V/m] = 1e-21*(E/n)[Td]*n[m^-3]
    assert abs(calc.tb_params["E"]) == 1e-21*2*n
    # Enforce the conversion B_i[Tesla] = 1e-27*(B/n)[Hx]*n[m^-3]
    assert abs(calc.tb_params["B"][1]) == 1e-27*2.5*n
    # Now with one species
    p = {"QE": [-1], "ME": [1]}
    L = 1e-6
    calc = TB(NP=[2], L=L, Ered=2, Bred=[0, 2.5, 0], **p)
    n = 2/L**3
    # Enforce the conversion E[V/m] = 1e-21*(E/n)[Td]*n[m^-3]
    assert abs(calc.tb_params["E"]) == 1e-21*2*n
    # Enforce the conversion B_i[Tesla] = 1e-27*(B/n)[Hx]*n[m^-3]
    assert abs(calc.tb_params["B"][1]) == 1e-27*2.5*n

def test_differential_param_changes():
    """ThunderBoltz object should be able to run a calculation with one
    set of parameters, then change settings, then run again properly."""
    calc = TB(NS=200, **He_settings())
    m = {"He_Park": "Park", "equal": "Equal", "uniform": "Uniform"}
    for eadf, eesd in [["He_Park", "default"], ["default", "equal"], ["default", "uniform"]]:
        calc.set_(eadf=eadf, eesd=eesd)
        calc.run()
        deck = calc.cs.get_deck()
        if eadf != "default":
            assert m[eadf] in deck
        if eesd != "default":
            assert m[eesd] in deck

def test_indeck_changes():
    """ThunderBoltz object should be able to run a calculation with one
    set of parameters, then change the indeck, then run again properly."""
    calc = TB(NS=200, **He_settings())
    calc.run()
    assert len(calc.cs.table) == 12
    # Remake dir
    setup_()
    # Drastic change of parameters
    calc.set_(indeck=pjoin("indecks", "ac_field"))
    calc.run()
    assert len(calc.cs.table) == 29

def test_to_primitive():
    """Test that numpy type conversion works."""
    a = {"hp": {"a": 3, "b": "asdf", "c": np.arange(3),
        "d": np.int64(3), "e": {"a": np.arange(3)}}}

    a_native = {"hp": {"a": 3, "b": "asdf", "c": [0., 1., 2.],
        "d": 3., "e": {"a": [0., 1., 2.]}}}

    tb.tb.to_primitive(a)
    assert a == a_native

def test_plot_timeseries():
    """Make sure tb object can generate a time series plot with
    the correct structure."""
    # Set up simple calculation.
    calc = TB(NS=1601, **He_settings())
    calc.run()
    calc.plot_timeseries(stamp=["directory"])

@pytest.mark.skipif(not shutil.which("sbatch"), reason='slurm not available on device')
def test_particle_overload_exit():
    """Make sure there is a particle dump
    after particle overload."""
    # Bypass
    cs = tb.input.CrossSections()
    p = Process("Ionization", threshold=1, cs_func=lambda e: 5e-20,
                name="TICS")
    cs.add_process(p)
    n = 1000
    calc = TB(NS=4001,
            NP=[n, 100],
            TP=[0.2, 0],
            MP=[1e-5, 1],
            QP=[-1, 1],
            SP=2,
            MEM=2000*72 / 2**30,
            cs=cs,
    )
    with slurm_debug(calc, account="xd") as slurm:
        directory = setup_(subdir=f"task1", base_directory=expand("$SCR"))
        slurm.submit(directory=directory)

    # Print while pending
    while slurm.has_active():
        time.sleep(2)
    vdfs = calc.read_vdfs(directory=directory)
    times = vdfs.step.unique()
    assert len(times) == 1

@pytest.mark.skipif(not shutil.which("sbatch"), reason='slurm not available on device')
def _test_mem_bounds():
    """Ensure memory is allocated properly, and no leaks near the
    end of large particle runs."""
    cs = tb.input.CrossSections()
    p = Process("Ionization", threshold=.1, cs_func=lambda e: .1e-20,
                name="TICS")
    cs.add_process(p)
    n = 100000
    mem = n*72/2**30
    calc = TB(NS=40000001,
            NP=[n, 100],
            TP=[10.0, 0.0],
            E = 0,
            MP=[5.4857e-4, 1e4],
            eesd="uniform",
            QP=[-1, 1],
            SP=2,
            MEM=mem*1.1,
            cs=cs,
    )
    with slurm_debug(calc, account="xd") as slurm:
        directory = setup_(subdir=f"task1", base_directory=expand("$SCR"))
        slurm.submit(directory=directory)

    while slurm.has_active():
        try:
            ts = tb.read_timeseries(directory=directory)
            print(ts.tail(), flush=True)
        except:
            pass
        time.sleep(2)

def test_energy_conv():
    """Ensure memory is allocated properly, and no leaks near the
    end of large particle runs."""
    m = 5.48579909067e-4
    calc = TB(NS=200,
            OS=1,
            TP=[80, 0],
            L=8.00e-7, # Force many collisions
            MP=[5.4857e-4, 1e11],
            E=-0.,
            # out_file="_stdout_",
            NP=[10000, 100],
            **He_settings())
    # Change the reaction table to have no threshold values
    calc.cs.table.loc[:,"B"] = 0.0
    # Remove some processes
    # calc.cs.table = calc.cs.table.iloc[-1:,:]
    calc.run()
    ts = calc.get_timeseries()
    ts = pd.concat((ts, calc.particle_tables[0]), axis=1)
    keys = ["MEe", "KEe", "E", "Ni"]
    diffs = ts[keys].select_dtypes(
        include=np.number).diff().values[1:,:]
    assert not np.any(diffs.flatten())

def test_compile_debug():
    """Test that at zero field, frozen background, 0 thresholds,
    mean energy is conserved for all the cross sections."""
    calc = TB(NS=200000000001,
            Ered=500,
            NP=[1e9, 1000],
            MEM=90,
            **He_settings())
    calc.compile_debug()

def test_query_tree():
    """Make sure query tree is working properly."""

    # Setup a run with various kinds of data in various places.
    # Setup base object
    base_path = setup_()

    calc = TB(
        directory=base_path,
        NS=300,
        L=1e-6,
        NP=[1000, 100],
        indeck=tb.input.He_TB,
        FV=[99, 100, 0],
    )

    # Test round 1
    for Ered in [100, 200]:
        path = setup_(subdir=f"constB10/{Ered}Td")
        calc.set_(directory=path, Ered=Ered, Bred=[0,0,10])
        calc.run()

    for Ered in [100, 200]:
        path = setup_(subdir=f"constB20/{Ered}Td")
        calc.set_(directory=path, Ered=Ered, Bred=[0,0,20])
        calc.run()

    # Test round 2 (potentially run at a different time from round 1)
    i = 0
    for Ered in [100, 200]:
        for Bred in [10, 20]:
            path = setup_(subdir=f"job{i}")
            calc.set_(directory=path, Ered=Ered, Bred=[0,0,Bred])
            calc.run()
            i+=1

    # Test defaults
    calcs = query_tree(base_path)
    assert len(calcs) == 8
    assert all(isinstance(calc, TB) for calc in calcs)

    # Test the name requirements
    calcs = query_tree(base_path, name_req=lambda s: "B10" in s)
    assert len(calcs) == 2
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, name_req=lambda s: "const" in s)
    assert len(calcs) == 4
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, name_req=lambda s: "job" in s)
    assert len(calcs) == 4
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, name_req=lambda s: "query" in s)
    assert len(calcs) == 8
    assert all(isinstance(calc, TB) for calc in calcs)


    # Test the param requirements
    calcs = query_tree(base_path, param_req={"Ered": 100})
    assert len(calcs) == 4
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, param_req={"Ered": 200})
    assert len(calcs) == 4
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, param_req={"Bred": [0, 0, 20]})
    assert len(calcs) == 4
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, param_req={"Bred": [0, 0, 20], "Ered": 100})
    assert len(calcs) == 2
    assert all(isinstance(calc, TB) for calc in calcs)

    calcs = query_tree(base_path, param_req={"NS": 300})
    assert len(calcs) == 8
    assert all(isinstance(calc, TB) for calc in calcs)


    # Test the callback behavior
    df = query_tree(base_path, callback=TB.get_timeseries)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4*8
    vdfs = query_tree(
        base_path,
        callback=lambda calc: calc.get_vdfs(steps="all")
    )
    assert isinstance(vdfs, pd.DataFrame)
    assert len(vdfs) == 3*1000*8
    assert len(vdfs.step.unique()) == 3

    # Callback with concatenation
    def callback_mult_df(calc):
        return calc.get_timeseries(), calc.get_vdfs(steps="all")

    df, vdfs = query_tree(
        base_path,
        callback=callback_mult_df,
    )
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4*8
    assert isinstance(vdfs, pd.DataFrame)
    assert len(vdfs) == 3*1000*8
    assert len(vdfs.step.unique()) == 3

    # Callback without concatenation
    calc_ex = query_tree(
        base_path,
        callback=callback_mult_df,
        agg=False,
    )
    assert len(calc_ex[0]) == 8
    assert len(calc_ex[1]) == 8
    assert all(isinstance(ts, pd.DataFrame) and isinstance(vdf, pd.DataFrame)
               for ts, vdf in zip(*calc_ex))



def _test_get_distributions():
    """Test all the energy and velocity distribution parsers."""
    path = setup_(dry=True)
    ns = 5000
    fv = [ns-101, 50, 0]
    ndumps = (ns-fv[0])//fv[1]  + 1
    if False:
        calc = TB(
            NS=ns,
            Ered=200,
            # autostep=True,
            pct_ion=10,
            OS=100,
            TP=[70, 0],
            MP=[1e-5, 1e10],
            L=1.00e-7, # Force many collisions
            FV=fv,
            NP=[100000, 100],
            indeck=tb.input.He_TB,
            directory=path,
            live=True,
            live_rate=True,
        )
        setup_()
        calc.run(std_banner=True)

    else:
        calc = tb.read(path, read_cs_data=True)

    calc.plot_timeseries()
    edfs = calc.get_edfs(steps="all")
    ddist = calc.describe_dist()
    assert len(ddist) == 8*ndumps
    fig_path = pjoin(path, "figures")
    os.makedirs(fig_path, exist_ok=True)

    # calc.plot_vdfs(steps="all", save=fig_path)
    # calc.plot_timeseries(save=fig_path, stamp=["Ered", "E", "B", "directory"])
    calc.plot_edf_comps(steps="all", bins=100)
    # os.system(f"open {fig_path}/*.pdf")
    # os.system(f"open testing/output/test_get_distributions/figures/test_get_distributions.pdf")

def _test_live_plot():
    """Simplest TB setup with live plotting"""
    calc = TB(NS=3001, **He_settings())
    calc.run(std_banner=True, monitor=True, live=True, live_rate=True)
    ts = calc.get_timeseries()
    assert len(calc.ts_fig.axes) == 3
    assert len(calc.rate_fig.axes) == 1

def _test_He_transport():
    """Produce steady state electron transport parameters for Helium at various
    fields with various electron energy sharing distribution models."""
    path = setup_()
    # Create ThunderBoltz object
    calc = TB(
        indeck=tb.input.He_TB,
        out_file="thunderboltz",
        eesd="default",
        eadf="default",
        Ered=1500,
        NS=100000001,
        L=1.46e-6,
        DE=.1,
        OS=100,
        directory=path,
        DT=2.0485560548363354e-09,
        # NP=[5000000, 167],
        autostep=True,
        pct_ion=10,
        MEM=5,
        # FV=[10000, 100, 0],
        # TP=[200, 0],
        n=3,
        SE=1, # Slurm Exit
        egen=True,
        monitor=True,
        live=True,
        live_rate=True,
    )
    calc.run(std_banner=True)

def _test_e_beam():
    path = setup_()
    cs = tb.CrossSections()
    p = Process("Ionization", threshold=100.,
            cs_func=lambda e: 20e-20, name="ionz")
    cs.add_process(p)

    calc = TB(
        cs=cs,
        # indeck=tb.input.He_TB,
        out_file="thunderboltz",
        eesd="equal",
        eadf="default",
        E=100,
        NS=30,
        # NS=300000,
        VV=[150, 0],
        TP=[0, 0],
        L=1.46e-6,
        DE=.1,
        OS=1,
        # Nmin=1000,
        directory=path,
        DT=2.0485560548363354e-8,
        NP=[17, 16700],
        MEM=2,
        SP=2,
        QP=[-1.0, 0.0],
        MP=[5.4857e-4, 4.0],
        # TP=[200, 0],
        n=3,
        SE=0, # Slurm Exit
        egen=True,
        monitor=True,
        # live=True,
        # live_rate=True,
    )

    calc.run(std_banner=True)

def _test_e_beam_cp():
    """Test e beam with inelastic particle change"""
    path = setup_()
    cs = tb.CrossSections()
    p = Process("InelasticChangeParticle2", 1, 0, 1, 2,
            cs_func=lambda e: 10e-20, name="elastic")
    cs.add_process(p)

    calc = TB(
        cs=cs,
        # indeck=tb.input.He_TB,
        out_file="thunderboltz",
        eesd="equal",
        eadf="default",
        E=0,
        # NS=3,
        NS=30000,
        VV=[110, 0, 110],
        TP=[0, 0, 0],
        FV=[0, 1, 2],
        L=1.46e-6,
        DE=.1,
        OS=1,
        # Nmin=1000,
        directory=path,
        DT=2.0485560548363354e-9,
        NP=[1670, 167, 1670],
        MEM=2,
        SP=3,
        QP=[-1.0, 0.0, -1.0],
        MP=[5.4857e-4, 4.0, 5.4857e-4],
        # TP=[200, 0],
        n=3,
        SE=0, # Slurm Exit
        egen=True,
        monitor=True,
        # live=True,
        # live_rate=True,
    )

    print(calc.cs.table)
    calc.run(std_banner=True)
    calc.get_vdfs(particle_type=[0,1,2])
    print(calc.describe_dist(steps="all"))
    calc.describe_dist().to_csv("~/Desktop/dist.csv")

    calc.set_ts_plot_params(["Ni"])
    calc.plot_timeseries()

def _test_hot_start():
    """Produce steady state electron transport parameters for Helium at various
    fields with various electron energy sharing distribution models."""
    path = setup_(dry=True)
    # Create ThunderBoltz object
    calc = TB(
        indeck=tb.input.He_TB,
        out_file="thunderboltz",
        eesd="default",
        eadf="default",
        Ered=1500,
        NS=100000001,
        L=1.46e-6,
        DE=.1,
        OS=100,
        directory=path,
        DT=2.0485560548363354e-09,
        NP=[50000, 167],
        VV=[10, 0],
        TP=[200, 0],
        autostep=True,
        pct_ion=10,
        MEM=5/3,
        # FV=[10000, 100, 0],
        # TP=[200, 0],
        n=3,
        SE=1, # Slurm Exit
        egen=True,
        monitor=True,
        live=True,
        live_rate=True,
    )
    # calc.run(std_banner=True)
    calc_r = tb.read(path)

    calc_r.plot_timeseries()

def _test_vdf_loader():
    """Check that the VL ThunderBoltz behavior works, then check that
    API interfaces properly."""
    path = setup_(dry=True)
    # Create ThunderBoltz object
    if False:
        path = setup_()
        calc = TB(
            indeck=tb.input.He_TB,
            out_file="thunderboltz",
            eesd="default",
            eadf="default",
            Ered=1500,
            NS=1001,
            L=1.46e-6,
            DE=.1,
            OS=100,
            directory=path,
            DT=2.0485560548363354e-09,
            autostep=True,
            pct_ion=10,
            MEM=5,
            FV=[999, 1, 0],
            # TP=[200, 0],
            egen=True,
            monitor=True,
            live=True,
            live_rate=True,
        )
        calc.run(std_banner=True)
    else:
        calc = tb.read(path)

    print(calc._check_available_vdfs())

def test_vdf_init_setup():
    """Test that the velocity dump init file can be read from arbitrary places,
    specified programmatically, or read from a previous run."""
    calc = TB(
        NS=200,
        vdf_init=[pjoin("indecks", "100eV_beam.dat"), 0],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 10000

    # Create programmatic vdfs

    vdf_init_lists = [[0, 0, 100]]*10000
    vdf_init_np = np.array(vdf_init_lists)
    vdf_init_df = pd.DataFrame(
        vdf_init_np,
        columns=["vx", "vy", "vz"]
    )

    # Test pandas
    calc = TB(
        NS=200,
        vdf_init=[vdf_init_df, 0],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 10000

    # Test numpy
    calc = TB(
        NS=200,
        vdf_init=[vdf_init_np, 0],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 10000

    # Test lists
    calc = TB(
        NS=200,
        vdf_init=[vdf_init_lists, 0],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 10000

    # Check downsample
    calc = TB(
        NS=200,
        vdf_init=[pjoin("indecks", "100eV_beam.dat"), 0],
        downsample=True,
        NP=[2000, 10],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 2000

    # Check multi dump
    calc = TB(
        NS=200,
        vdf_init=[[pjoin("indecks", "100eV_beam.dat"), 0],
                  [pjoin("indecks", "100eV_beam.dat"), 1]],
        FV=[[190, 10, 0],
            [190, 10, 1]],
        **He_settings(),
    )
    calc.run(dryrun=True)
    pth = pjoin(calc.get_directory(), "velocity_init", "ParticleType0.dat")
    with open(pth, "r") as f:
        assert len(f.readlines()) == 10000

def test_vdf_init_run():
    """Test that the velocity dump init file can be read from arbitrary places,
    specified programmatically, or read from a previous run."""
    a = np.arange(1)
    b = np.stack((a, a, a), axis=1)
    calc = TB(
        NS=200,
        vdf_init=[[pjoin("indecks", "100eV_beam.dat"), 0],
                  [[[10, 10, 10]]*160, 1]],
        # vdf_init=[b,1],
        FV=[[100, 50, 0],
            [190, 1, 0],
            [100, 50, 1]],
        egen=False,
        **He_settings(),
    )
    calc.run()

def _test_plot_edfs():
    """Test the energy distribution plot."""

    if True:
        calc = TB(DT=1e-9, NS=40000, eesd="equal", FV=[15000, 1000, 0],
                NP=[10000, 1000], Ered=1000, **He_settings())

        calc.run(std_banner=True)

    calc = tb.read("testing/output/test_plot_edfs", read_cs_data=True)

    print(calc.describe_dist())
    calc.plot_edfs(steps="all", plot_cs=True)

def test_ExB():
    """In crossed electric and magnetic fields, the particles should drift in the
    ExB direction, regardless of the sign of charge"""
    calc = TB(DT=1e-4, NP=[1], TP=[0], NS=1000, Ered=100,
              MP=[1], QP=[-1], Bred=[0, 1000, 0], cs=tb.CrossSections(),
              directory=setup_(1))
    calc.run()
    pt = calc.get_particle_tables()
    v1 = pt[0].Vxi.mean()

    calc = TB(DT=1e-4, NP=[1], TP=[0], NS=1000, Ered=100,
              MP=[1], QP=[1], Bred=[0, 1000, 0], cs=tb.CrossSections(),
              directory=setup_(1))
    calc.run()
    pt = calc.get_particle_tables()
    v2 = pt[0].Vxi.mean()

    assert np.sign(v1) == np.sign(v2) == 1

def test_E_pulse():
    """Create a Gaussian model of an E field pulse."""
    calc = TB(NS=3000, DT=1e-10, OS=10, EP=[2e-7, 1e-8], **He_settings())
    calc.run()
    calc.plot_timeseries(series=["E", "k_ion", "k_8"])
    ts = calc.get_timeseries()
    # assert len(ts) == 4
