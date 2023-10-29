=====================
Simulation Parameters
=====================


Input Parameters
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: api/

   thunderboltz.parameters.TBParameters
   thunderboltz.parameters.WrapParameters

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
specify the ``MEM`` flag in the :class:`~.thunderboltz.ThunderBoltz` constructor:

.. code-block:: python

   import thunderboltz as tb

   calc = tb.ThunderBoltz(
        # For example, using the Helium model.
        indeck=tb.input.He_TB,
        # This will turn on electron generation for the Helium model
        # i.e. this will ensure the "Ionization" collision model is
        # used in the generated indeck.
        egen=True,
        # Now we must set the MEM flag, since we will be generating
        # a lot of electrons.
        MEM = 10, # in GB
   )

MEM will accept any float representing the number of gigabytes
to be made available to the particle arrays.

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


