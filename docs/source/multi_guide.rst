=============================
Running Multiple Calculations
=============================

In Sequence
~~~~~~~~~~~

You can change simulation parameters in the :class:`~.thunderboltz.ThunderBoltz`
object and run the program again in a new directory. Use
the :meth:`~.thunderboltz.ThunderBoltz.set_` method to update the desired parameters.

Suppose you wanted to run several calculation at various
field values. To do this, loop through the field values,
create new directories for the new calculation and run the
object like so:

.. code-block:: python

   import os
   import thunderboltz as tb

   # Make a base directory for this ensemble of simulations
   os.makedirs("multi_sim")

   calc = ThunderBoltz(indeck=tb.input.He_TB)

   fields = [10, 100, 500]
   # Loop through the field values
   for field in fields:
       # Create a new directory for this calculation
       subdir = os.path.join("multi_sim", f"{field}Td")
       os.makedirs(subdir)
       calc.set_(Ered=field, directory=subdir)
       # Run the calculation
       calc.run()

Each call to :meth:`~..ThunderBoltz.run` will block until the
corresponding simulation is finished.

In Parallel
~~~~~~~~~~~

Now suppose you would like to take advantage of multiple cores to run
several ThunderBoltz calculations at once. Though the internal kinetic
code is not (yet) parallelized, the python interface can run several
ThunderBoltz subprocesses in parallel like so:

.. code-block:: python

   import os
   import thunderboltz as tb

   # Make a base directory for this ensemble of simulations
   base_path = "multi_sim_parallel"
   os.makedirs(base_path)

   # Create the base object for the calculation
   calc = ThunderBoltz(indeck=tb.input.He_TB)

   fields = [10, 100, 500]

   # This time use the DistributedPool context,
   # passing the ThunderBoltz object like so
   with DistributedPool(calc) as pool:
       # Loop through the field values
       for field in fields:
           # Create a new directory for this calculation
           subdir = os.path.join(base_path, f"{field}Td")
           os.makedirs(subdir)

           # Rather than running with the ``ThunderBoltz`` object,
           # submit the changes to the pool, and it will automatically
           # run each each submitted calculation in parallel.
           pool.submit(Ered=field, directory=subdir)

   # The DistributedPool context will wait for all the jobs to finish
   # before continuing execution outside the 'with' block.

.. warning::

   The forking process used to run multiple simulations has thusfar
   only been tested on UNIX/LINUX operating systems.

.. warning::

   Ensure there is enough simultaneous memory for all jobs when running
   them in parallel. See the section on
   :ref:`Electron Growth and Memory Management <memory>`.


With a Job Manager
~~~~~~~~~~~~~~~~~~

If HPC resources are available to the user, the python API
includes a job manager compatible with the
`SLURM <https://slurm.schedmd.com/documentation.html>`_ protocol.
The :class:`~.thunderboltz.parallel.SlurmManager` context allows
for many different calculations to be split up among compute nodes,
and further distributed across cores. Use it as follows:

.. code-block:: python

   import os
   import thunderboltz as tb

   # Make a base directory for this ensemble of simulations
   base_path = "multi_sim_slurm"
   os.makedirs(base_path)

   # Create the base object for the calculation
   calc = ThunderBoltz(indeck=tb.input.He_TB)

   fields = [10, 100, 500]

   # Configure SLURM parameters for your job
   slurm_options = {
       "account": "my_account",
       "time": 100, # in minutes
       "job-name": "test_slurm",
       "ntasks-per-node": 8, # Specify number of cores to use
       "qos": "debug",
       "reservation": "debug",
   }

   # Use the SlurmManager Context, just like the DistributedPool context,
   # but also give it your SLURM options.
   with SlurmManager(calc, base_path, **slurm_options) as slurm:
       # Loop through the field values
       for field in fields:
           # Create a new directory for this calculation
           subdir = os.path.join(base_path, f"{field}Td")
           os.makedirs(subdir)
           # Use the slurm manager the same way as the pool, it will
           # handle node and core allocation internally.
           slurm.submit(Ered=field, directory=subdir)

See `here <https://docs.python.org/3/reference/expressions.html#calls>`_
for an explanation of the ``**`` (unpacking) operator used
in the previous example.

.. note::

   This job manager currently only works for clusters that either
   already have the gcc and python requirements installed on each
   compute node, or clusters that use the
   `Module System <https://hpc-wiki.info/hpc/Modules>`_ to load
   functionality.

   The default behavior is to accomodate the module system as it
   is common on most HPC machines. If you wish to avoid writing
   ``module load`` commands in the SLURM script, simply specify
   ``modules=[]`` in the ``SlurmManager`` constructor.

.. warning::

   Ensure there is enough memory for all parallel jobs when running
   them in parallel.
