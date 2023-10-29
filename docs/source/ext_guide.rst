Extracting Results
==================

Reading a Single Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After a calculation is finished, you can easily read the output
data using the python API like so:

.. code-block:: python

   import thunderboltz as tb

   # Pass the location of the simulation directory to be read
   calc = tb.read("path/to/previous/simulation_directory")

A single :class:`~.thunderboltz.ThunderBoltz` object will be returned, with which
you can easily `export <#exporting-data>`_ or `plot <#plotting-results>`_
output data.

When reading calculations in this way, you may or may not want to extract
cross section data from the simulation directory as well. To save on
runtime, cross section data is not read in by default. However, if you wanted
to read in cross section data from an old calculation and reuse that
data for other purposes, you can use the ``read_cs_data`` argument:

.. code-block:: python

   import thunderboltz as tb

   calc = tb.read("path/to/previous/simulation_directory", read_cs_data=True)

   # You will now see the cross section data has been loaded
   print(calc.cs.data)

Reading Many Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~

When running many calculations in various directories, it can be convenient to
read all of the output data at once. Imagine a directory structure like this:
::
   path/to/base_path
    /———sim1
        /———indeck_file.in
        /———cross_sections
        /———thunderboltz.out
        ...
    /———sim2
        ...
    /———sim3
        ...
    ...

where several ThunderBoltz calculations are stored in one base directory
located at ``path/to/base_path``. You locate and extract all relevant
ThunderBoltz data out of a directory tree using the
:py:func:`~.thunderboltz.tb.query_tree` function:

.. code-block:: python

   import thunderboltz as tb

   calcs = tb.query_tree("path/to/base_path")

   # Now you can access each of the simulation objects
   # separately. For example:

   # View the time series data from the first read calculation.
   print(calcs[0].get_timeseries)

   # Plot the last velocity dump data from the fourth read calculation.
   calcs[3].plot_vdfs()

See the :py:func:`~.thunderboltz.tb.query_tree` API reference to learn about
options for filtering criteria and automatically merging data from several
calculations.

As with the `single calculation <#reading-a-single-calculation>`_ case,
you can request the cross section data by providing the ``read_cs_data``
argument:

.. code-block:: python

   import thunderboltz as tb

   calcs = tb.query_tree("path/to/previous/simulation_directory", read_cs_data=True)

   # Now each of the calculations will have cross section model data attached to them.
   # For example, this will print the collision table for the 3rd read in calculation.
   print(calcs[2].cs.table)

Accessing Data
~~~~~~~~~~~~~~

Either after a calculation has finished, or after reading output data as shown
above, all data can be extracted from the :class:`~.thunderboltz.ThunderBoltz`
object:

Time-dependent data for the attributes found in
:class:`~.thunderboltz.parameters.OutputParameters` can be accessed with
:meth:`~.thunderboltz.ThunderBoltz.get_timeseries`:

.. code-block:: python

   data = tb.get_timeseries()

Time-averaged data for the attributes found in
:class:`~.thunderboltz.parameters.OutputParameters` can be accessed with
:meth:`~.thunderboltz.ThunderBoltz.get_ss_params`:

.. code-block:: python

   data = tb.get_ss_params()

This method will also compute standard deviations over the
steady-state interval for each parameter in a new column
with a "_std" suffix added to the column name.

.. warning::
   Currently, the last quarter of the timeseries data is assumed to be
   in steady-state by default when calculating these steady-state parameters.
   Please verify that this is true by viewing the figures produced by
   :py:func:`plot_timeseries`. Otherwise, run the simulation for longer,
   or provide your own appropriate criteria via the ``ss_func`` option when
   calling :meth:`~.thunderboltz.ThunderBoltz.get_ss_params`.

Output parameters for the attributes found in
:class:`~.thunderboltz.parameters.ParticleParameters` can be accessed with
:meth:`~.thunderboltz.ThunderBoltz.get_particle_tables`:

.. code-block:: python

   data = tb.get_particle_tables()

   # For example, this will write the mean energy, and
   # each of the mean displacement components to a csv
   # called "R_export.csv"
   data.to_csv("R_export.csv", index=False)

Exporting Data
~~~~~~~~~~~~~~

Once data is in the form of a :class:`~.pandas.DataFrame`, it is easy
to export it to other formats. See the
`Pandas I/O Guide <https://pandas.pydata.org/docs/user_guide/io.html>`_
for extensive options for converting from the ``DataFrame`` object.
The simplest option is to convert the data to a csv:

.. code-block:: python

   # This will write the data into a new file called "my_new_file.csv"
   data.to_csv("my_new_file.csv", index=False)

.. note::

   When exporting data to the csv format from a pandas
   DataFrame, it is usually most convenient to pass ``index`` = ``False``
   to prevent :meth:`~.pandas.DataFrame.to_csv` from writing
   the index (usually just an enumeration of the rows) into the
   first column of the csv.


Plotting Results
~~~~~~~~~~~~~~~~

The :class:`~thunderboltz.ThunderBoltz` API offers functions for
automatically plotting results. See the documentation for the following
functions

.. autosummary::

   thunderboltz.ThunderBoltz.plot_timeseries
   thunderboltz.ThunderBoltz.plot_rates
   thunderboltz.ThunderBoltz.plot_edf_comps
   thunderboltz.ThunderBoltz.plot_edfs
   thunderboltz.ThunderBoltz.plot_cs

These functions will plot the data into :class:`~.matplotlib.figure.Figure` objects,
but in order to see the plots in a GUI, you must import the plotting library
and include the line ``plt.show()`` after calling plotting methods like so:

.. code-block:: python

   import thunderboltz as tb

   # This will import the plotting library
   import matplotlib.pyplot as plt

   # Either read in data, or run calculations
   calc = tb.read("path/to/simulations_to_plot", read_cs_data=True)

   # Call plotting methods
   calc.plot_cs()

   # Show the plots and load a GUI
   plt.show()

Alternatively, you may specify a directory within which to
save a pdf file of the plot when calling any ``ThunderBoltz.plot_*`` method.
For example:

.. code-block:: python

   calc.plot_cs(save="path/to/figure_directory")

