.. _api:

====================================
ThunderBoltz |release| Documentation
====================================

This documentation includes some simple tutorials for using
the ThunderBoltz plasma simulation package and a complete
public API reference.

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
source in ``src/thunderboltz``. To install the Python interface, run the
``install.sh`` script from the root directory to install the python package.

.. code-block:: bash

   ./install.sh

.. warning:: 

   This will upgrade `pip <https://pypi.org/project/pip/>`_ and install specific versions of python packages,
   so create an environment if you are concerned with python package overwrite.


Tutorials
---------

See the :doc:`Quick Start Guide <quickstart>` to go through a brief tutorial
setting up a simple calculation and interpreting the results.

See :doc:`Preparing Cross Sections <cs_guide>` for a tutorial on the various
ways to obtain and manipulate cross sections for ThunderBoltz simulations.

See :doc:`Running Multiple Calculations <multi_guide>` for a quick guide on
how to vary simulation parameters and easily run simulations in
parallel.

See :doc:`Extracting Results <ext_guide>` for details on easily parsing,
plotting, and exporting data from many ThunderBoltz simulations at once.

.. toctree::
   :hidden:
   :maxdepth: 2

   self
   quickstart
   cs_guide
   multi_guide
   ext_guide

Benchmark Testing
-----------------

See :doc:`Benchmark Testing <bm>` to run code that reproduces
the results found in the paper.

.. toctree::
   :hidden:
   :maxdepth: 2

   bm

Simulation Parameters
---------------------

Review the :doc:`Simulation Parameters <params>` for information on
the default behavior of the code, the available input options, and
details regarding output parameter definitions and interpretations.

.. toctree::
   :hidden:
   :maxdepth: 2

   params


API Reference
-------------

See the full :doc:`API Reference <ref>` for full documentation
of the ThunderBoltz programming interface.

.. toctree::
   :hidden:
   :maxdepth: 2

   ref
