
========================
Preparing Cross Sections
========================

There are a variety of ways to specify cross sections with the
ThunderBoltz interface. In the :doc:`Quick Start Guide <quickstart>`,
we used a built-in Helium cross section model. A more general approach
to preparing cross sections is with the :class:`~.thunderboltz.CrossSections`
object.

Initializing the CrossSections Object
-----------------------------------------

There are three main ways to initialize a :class:`~.thunderboltz.CrossSections`
object:

#. With cross section data from another ThunderBoltz run

   .. code-block:: python

      from thunderboltz import CrossSections
      # Just specify the path to the simulation directory of a
      # different ThunderBoltz run.
      cross_sections = CrossSections(input_path="path/to/thunderboltz_sim_dir")

   Refer to the :class:`~.thunderboltz.CrossSections` section of the
   :doc:`API Reference <ref>` to ensure the
   simulation data is set up correctly for interpretation by
   :class:`~.thunderboltz.CrossSections`.

#. By reading from an `LXCat <https://nl.lxcat.net>`_ text file extract.

   .. code-block:: python

      from thunderboltz import CrossSections
      # First initialize an empty cross sections object
      cross_sections = CrossSections()
      # Then reference a text file extract from LXCat
      cross_sections.from_LXCat("path/to/LXCat_data.txt")

   .. note::

      For now, the LXCat parser assumes two species electron-gas
      systems where all processes are between electrons and
      gas macroparticles. If you wish to use LXCat data for other
      purposes, you can alter the species indices to your liking
      via ``CrossSections.table`` after loading in LXCat data.

#. By programmatically generating cross section data in python.

   This approach involves the :class:`~.thunderboltz.Process` object.

   .. code-block:: python

      from thunderboltz import CrossSections
      from thunderboltz import Process

      # Initialize an empty cross sections object
      cross_sections = CrossSection()

      # Next make a few processes

      # You can pass arbitrary tabulated data like so
      elastic_data = [
        # [eV], [m^2]
          [0.0, 2e-20],
          [0.001, 2.1e-20],
          [.01, 3e-20],
          [10.0, 1e-19],
          [1000, 1e-18],
          [10000, 2e-19],
      ]
      elastic_process = Process(
          "Elastic", # The type of process
          r1=0, # The first reactant species index
          r2=1, # The second reactant species index
          p1=0, # The first product species index
          p2=1, # The second product species index
          cs_data=elastic_data,
          # This will determine the name of the
          # written cross section file and ideally should
          # be unique.
          name="elastic_example",
      )
      # You can also pass data frames, or ndarrays if that is
      # preferable

      # Or, use an analytic form defined with a python
      # function.
      import numpy as np # Import math functionality
      def inelastic_model(energy, parameter):
          # It's okay to have conditional statements
          if energy < 5:
              return parameter

          # And nonlinear functions
          return parameter*np.log(energy)/energy

      # You can parameterize your model
      cs_mod_1 = lambda e: inelastic_model(e, 1e-20)
      cs_mod_2 = lambda e: inelastic_model(e, 2e-20)
      cs_mod_3 = lambda e: inelastic_model(e, 3e-20)

      # And create multiple cross sections
      inelastic_1 = Process(
          "Inelastic", threshold=1., cs_func=cs_mod_1, name="inelastic1")
      inelastic_2 = Process(
          "Inelastic", threshold=1., cs_func=cs_mod_2, name="inelastic2")
      inelastic_3 = Process(
          "Inelastic", threshold=1., cs_func=cs_mod_3, name="inelastic3")

      # Finally, you can create processes with differential cross section
      # models, if they are available in your ThunderBoltz version.
      ionization = Process("Ionization", threshold=10.,
          cs_func=lambda e: 1e-19*np.log(e)/e,
          # This, for example, will add the equal energy sharing condition
          differential_process="equal",
          name="ionization")


      # You can add your process to the CrossSections object one at a time
      cross_sections.add_process(elastic_process)
      # Or all at once
      cross_sections.add_processes(
          [inelastic_1, inelastic_2, inelastic_3, ionization]
      )

   .. note::

         It is important to explicitly specify threshold values for
         inelastic and superelastic processes because their values will
         not be inferred from the cross section data.



Viewing Your Cross Sections
---------------------------
When parsing data from external sources, it is important to ensure
that the correct data is being used in the intended context for the
simulation. You can view the reaction table for the model by
printing out the ``table`` attribute.

.. code-block:: python

   print(cross_section.table)

And you can view the cross section data associated with each process
by printing out the ``data`` attribute.

.. code-block:: python

   print(cross_section.data)

To view a plot of the cross section data, use the
:meth:`~thunderboltz.CrossSections.plot_cs` method.

.. code-block:: python

    cross_section.plot_cs()

    # Remember to show the plot at the end of plotting scripts
    # Make sure to include the import statement "import matplotlib.pyplot as plt"
    plt.show()

See the API reference for plotting related quantities with the
:meth:`~thunderboltz.CrossSections.plot_cs` method.

Attaching the CrossSections bject
--------------------------------------

Finally, attach the ``CrossSections`` object to the main ThunderBoltz
object using the ``cs`` keyword to use the cross section model within it.

.. code-block:: python

   calc = ThunderBoltz(
       # ...
       cs=cross_sections,
       # ...
   )

   calc.run()
   # ...
