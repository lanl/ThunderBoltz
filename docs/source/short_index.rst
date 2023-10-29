.. _api:

====================================
ThunderBoltz |release| Documentation
====================================

This documentation includes instruction on installation
and compilation of ThunderBoltz source code and descriptions
of the input and output simulation parameters involved.

Installation
------------

For now, the code must be downloaded from a repository.
Use the following command to clone the code into a local repository.

.. code-block:: bash

   git clone git@gitlab.com/Mczammit/thunderboltz.git

You may need to set up SSH keys in order to access gitlab. See the
`Gitlab SSH Guide <https://docs.gitlab.com/ee/user/ssh.html>`_ to
set up access to GitLab repositories.


The basic ThunderBoltz functionality is available either
as an executable in ``bin/thunderboltz.bin`` or can be compiled from the
source in ``src/thunderboltz``. 


Compilation
-----------

ThunderBoltz requires a g++ of clang compiler and should be compiled
from source directories as

.. code-block:: bash

   g++ -std=c++17 -o thunderboltz.bin DSMC0D.cpp

Then run with

.. code-block:: bash

   ./thunderboltz.bin inputfile.in

to use a manually constructed indeck file. The code is
maintained with the standard `-Wall` and `-Werror`
compiler options.

Here is an example of how to run a simple ThunderBoltz calculation.

.. code-block:: bash

   # Make a directory for testing
   mkdir example_sim
   cd example_sim
   # Copy the source over
   cp ../src/thunderboltz/* .
   # Copy example input files over
   cp -r ../indecks/N2/* .
   # Compile
   g++ -std=c++17 -o thunderboltz.bin DSMC0D.cpp
   # Run
   ./thunderboltz.bin N2vib.in

Simulation Parameters
---------------------

Input Parameters
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: api/

   thunderboltz.parameters.TBParameters

.. _output_params:

Output Parameters
~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: api/

   thunderboltz.parameters.OutputParameters
   thunderboltz.parameters.ParticleParameters


.. _memory:

Electron Growth and Memory Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the ionization model and field strength,
ThunderBoltz may generate a large number of electrons.
In these cases, the appropriate amount of memory must be
allocated. The correct amount will be allocated automatically
in scenarios where no ionization process is used,
or when the ``IonizationNoEgen`` model is used. This amount
will be allocated based on the sum of all ``NP`` elements
times 4.

However, in scenarios where there is significant electron generation,
i.e. at high :math:`E` fields with the ``Ionization`` model on,
the default memory settings are not sufficient and the simulation
will exit with the error "Too many particles!". To prevent this
specify the ``MEM`` flag in the indeck. MEM will accept any
float representing the number of gigabytes to be made available
to the particle arrays.

.. warning::

   If the value of ``MEM`` is more than the actual number of
   available GB, then the simulation will still run, but will
   exit with a segmentation fault once too many particles are
   created.

.. warning::

   When using multiple cores on the same machine / node, ensure
   that each process has enough memory requested and that
   the sum of memory requests does not exceed the available
   pool of RAM.


