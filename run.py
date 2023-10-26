"""Setup and run the simulations found in the ThunderBoltz paper."""
import inspect
import os
from os.path import join as pjoin
from os.path import expandvars as expand
import shutil

import matplotlib.pyplot as plt
import pandas as pd

import pickle
import pytb
from pytb.kinetic import Process
from pytb.parallel import SlurmManager
from pytb.parallel import DistributedPool
import visualize

# Physical Constants
ME = 9.1093837e-31 # kg
QE = 1.60217663e-19 # C or J/eV
AMU_TO_KG = 1.6605e-27 # kg / amu

def onsager_relation():
    """Test a simple system with a closed-form equilibrium condition and
    compare the results."""
    # The species are indexed in the following order
    sp = ["X", "A", "B", "C"]
    # Initial density conditions for species X, A, B, C (m^-3)
    n_X = 2e24
    n_A = 1e22
    n_B = n_C = 0
    # Cell size (m)
    L = 1e-6

    # Initialize cross sections
    cross_sections = pytb.CrossSections()
    # Make the processes
    # Some processes have a threshold energy of 1eV
    activated = lambda cs: (lambda e: cs if e > 1 else 0)
    # Others are constant
    const = lambda cs: (lambda e: cs)

    # Create general process constructors

    # These have varying species and cross sections
    chem_proc = lambda i, j, cs_func, thresh: Process("InelasticChangeParticle2",
        0, i, 0, j, threshold=thresh, cs_func=cs_func, name=f"{sp[i]}+X->{sp[j]}+X")
    # These have constant cross sections and unchanging species.
    nonchem_proc_1 = lambda i, j: Process("Elastic",
        i, j, i, j, cs_func=const(1e-20), name=f"{sp[i]}+{sp[j]}->{sp[i]}+{sp[j]}")
    nonchem_proc_100 = lambda i, j: Process("Elastic",
        i, j, i, j, cs_func=const(1e-18), name=f"{sp[i]}+{sp[j]}->{sp[i]}+{sp[j]}")

    cross_sections.add_processes([
        # A <-> B
        chem_proc(1, 2, activated(1e-20), 1.0),
        chem_proc(2, 1, const(1e-20), -1.0),
        # B <-> C
        chem_proc(2, 3, activated(2e-20), 1.0),
        chem_proc(3, 2, const(2e-20), -1.0),
        # A <-> C
        chem_proc(1, 3, activated(2e-20), 1.0),
        chem_proc(3, 1, const(3e-20), -1.0),
        # Non reacting processes (every combo of species thermally interact)
        nonchem_proc_1(0, 0),
        nonchem_proc_100(0, 1),
        nonchem_proc_100(0, 2),
        nonchem_proc_100(0, 3),
        nonchem_proc_100(1, 1),
        nonchem_proc_100(1, 2),
        nonchem_proc_100(1, 3),
        nonchem_proc_100(2, 2),
        nonchem_proc_100(2, 3),
        nonchem_proc_100(3, 3),
    ])

    # Optionally print cross section info
    # print(cross_sections.table)
    # print(cross_sections.data)

    # Generate the ThunderBoltz wrapper object
    tb = pytb.ThunderBoltz(
        cs=cross_sections,
        duration=.2e-6, # run .2 millisecond sim
        DT=2e-11,       # Seconds per step
        E=0,            # Field strength (V/m)
        L=L,            # Cell length
        # Number of each type of particle
        NP=[n*L**3 for n in [n_X, n_A, n_B, n_C]],
        # Temperature of each species (eV)
        TP=[1., 1., 0, 0],
        # Masses (amu)
        MP = 4*[14.],
        # Charges (elem.)
        QP = 4*[0],     # all neutral
        VS=100,
        # Produce simulation files in simulation directory
        directory=setup_(),
    )

    # Run the simulation
    # File writing will be invoked here
    tb.run(std_banner=False, live=False)

    # Plot results once finished
    visualize.plot_onsager()

def ikuta_sugai():
    """Test electron transport in crossed electric
    and magnetic fields."""
    # electron and gas particle mass (amu)
    m_e = 5.4857e-4
    m = m_e * 100
    # Particle density (m^-3)
    n = 1e23
    # Cell size (m)
    L = 1e-6
    # Initial flow velocity
    tb = pytb.tb.ThunderBoltz(
        NS=200000, # Number of steps
        DT=5e-11, # Seconds per step
        Ered=1, # E-field strength [Td]
        L=L, # Cell length
        # Number of each particle
        NP=[6000, n*L**3],
        # Temp (eV)
        TP = [0.0, 0.0],
        # Masses (amu)
        MP = [m_e, m],
        # Charges (elem.)
        QP = [-1, 0],
    )

    # Add a constant process
    tb.cs.add_process(Process("ElasticFixedParticle2",
        cs_func = lambda e: 1e-20, name="IkutaElastic"))

    # Loop through B/n values
    with DistributedPool(tb) as pool:
        for Bred in [0, 10, 25, 50]:
            # Use new directory and assign new reduced B-field.
            pool.submit(directory=setup_(subdir=f"{Bred}Hx"), Bred=[0, Bred, 0])

    visualize.plot_ikuta_sugai()

def He_transport():
    """Produce steady state electron transport parameters for Helium at various
    fields with various electron energy sharing distribution models."""
    path = setup_(scratch=False)
    # Create ThunderBoltz object
    tb = pytb.ThunderBoltz(
        indeck=pytb.input.He_TB,
        eesd="uniform",
        eadf="default",
        Ered=1500,
        NS=15000,
        L=1.46e-6,
        DE=.1,
        MEM=15,
        Nmin=100,
        autostep=True,
        pct_ion=10,
        n=3,
        SE=1, # Slurm Exit
        egen=True,
        monitor=True,
    )
    # Run the following reduced fields
    Ereds = [1, 10, 50, 70, 100, 150, 200, 272.8, 323, 407, 500, 600, 704, 823, 1000, 1500]
    job_counter = 0
    with DistributedPool(tb) as pool:
        for eesd in ["default", "equal"]:
            for Ered in Ereds:
                # Setup path in scratch directory
                path = setup_(subdir=f"JOB{job_counter:03d}", scratch=False)
                job_counter += 1
                pool.submit(directory=path, eesd=eesd, Ered=Ered)

# Helpers
SIM_OUTPUT_DIR = pjoin("simulations")

def setup_(i=0, subdir=None, scratch=False):
    """Call from a run function such that it makes a directory
    in simulation/output with the name of that test function. Increment
    `i` if this function is being wrapped in another intermediary."""
    base_directory = SIM_OUTPUT_DIR
    if scratch:
        base_directory = expand("$SCR")
    # Create the path to the test file
    fpath = pjoin(base_directory, inspect.stack()[1+i][3])
    if subdir is not None:
        fpath = pjoin(fpath, str(subdir))
    # Overwrite directory if necessary
    if os.path.isdir(fpath):
        shutil.rmtree(fpath)
    os.makedirs(fpath)
    return fpath
