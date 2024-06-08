=================
Benchmark Testing
=================

There are three benchmark tests available in the
`repository <https://gitlab.com/Mczammit/thunderboltz>`_.
These benchmark simulations are described in detail in Sect. III of the ThunderBoltz paper.
The resulting calculations and figures from these benchmark tests can be compared directly
to the figures given in the paper. Each can be imported from the ``run.py`` python module.

Onsager Relation
----------------

The Onsager relation [#]_ predicts the kinetic rates of the following
chemical reactions between arbitrary heavy particles,

.. math::

   A \rightleftarrows B \rightleftarrows C \rightleftarrows A.


Based on the equilibrium condition, :math:`n_i k_{ij}=n_jk_{ji}`,
the rate constants :math:`k_{ij}` have the analytic solution

.. math::

    k(T) = d^2\sqrt{\frac{8 \pi k_{\textrm{B}}T}{m_r}} e^{-E_a/k_{\textrm{B}}T}.

To run this system in in ThunderBoltz, run the prepared function either
in a python script:

.. code-block:: python

   # Run this from within the repository root
   import run
   run.onsager_relation()

or directly from the command line:

.. code-block:: bash

   python -c "import run; run.onsager_relation()"

The resulting calculation will automatically run in the
directory ``simulations/onsager_relation``.

Once the simulation has finished, run the following on the command line
(or in a python script) to view a time evolution of the species densities,
reaction rates, and absolute rates.

.. code-block:: bash

   python -c "import visualize; visualize.plot_onsager()"

This will automatically save a pdf of the plot in the ``simulations``
directory.

Ikuta-Sugai
-----------

The Ikuta-Sugai benchmark problem tests electron transport
in crossed electric and magnetic fields.

To run this system in ThunderBoltz and compare it to the analytic
theory presented by Ness [#]_, run the prepared function either in a python
script of directly from the command line:

.. code-block:: bash

   python -c "import run; run.ikuta_sugai()"

Once the simulation has finished, run the following command to view
the effect of the magnetic field on the average velocity moments
and mean energy of the particles in comparison to Ness:

.. code-block:: bash

   python -c "import visualize; visualize.plot_ikuta_sugai()"

This will automatically save a pdf of the plot in the ``simulations``
directory.



He Transport
------------

Here we generate comparisons of bulk and flux electron mobility, :math:`\mu N`, and
Townshend ionization coefficient, :math:`\alpha / N`, at various reduced fields.
We compare ThunderBoltz results to the two-term Boltzmann equation solver, BOLSIG,
as well as some swarm experiments.

To simulate this system in ThunderBoltz run the prepared function either in a
python script or directly from the command line:

.. code-block:: bash

   python -c "import run; run.He_transport()"

Once the simulation has finished, run the following command to view
the reduced Townshend ionization coefficient and the reduced electron mobility
as a function of reduced electric field:

.. code-block:: bash

   python -c "import visualize; visualize.plot_He_transport()"

This will automatically save a pdf of the plot in the ``simulations``
directory. To view a plot comparing the individual reaction rate coefficients of
ThunderBoltz and BOLSIG, run the following:

.. code-block:: bash

   python -c "import visualize; visualize.rate_comp()"


.. [#] Light, J. C., Ross, J., & Shuler, K. E. (1969). Rate coefficients,
       reaction cross sections and microscopic reversibility. Kinetic
       Processes in Gases and Plasmas, 314, 281.

.. [#] K F Ness 1994 J. Phys. D: Appl. Phys. 27 1848.

