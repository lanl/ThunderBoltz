"""Main class for kinetic processes and their cross sections."""
import copy
import warnings

import numpy as np
from numpy import vectorize as vec
import pandas as pd


class Process(object):
    r"""A reaction process determined by reaction and product indices,
    a process type, a potential threshold, and a corresponding cross
    section specification.

    process_type: str
        Elastic | Inelastic | Ionization
    r1,r2,p1,p2: int (0,1,0,1)
        The indices of the reactants and products.
    threshold: float
        The threshold value for the process (e.g. the binding
        energy of an ionization process).
    cs_func: callable
        Function that returns the cross section for this process in
        :math:`\mathrm{m}^2` given an incident electron energy in
        eV (center of mass frame).
    cs_data: 2-D Array-Like
        Tabular cross section data with columns of energy (eV) and cross
        section (:math:`\mathrm{m}^2`).
    """
    # sample grid default range (orders of magnitude)
    SAMPLE_MIN = -4
    SAMPLE_MAX = 6

    def __init__(self, process_type, r1=0, r2=1, p1=0, p2=1,
                 threshold=0.0, cs_func=None, cs_data=None,
                 name=None, differential_process=None,
                 nsamples=250):
        self.process_type = process_type
        self.r1, self.r2, self.p1, self.p2 = r1, r2, p1, p2
        self.threshold = threshold
        self.differential_process = differential_process
        self.differential_parameters = None
        self.name = name
        self._cs_func = cs_func
        self._nsamples = nsamples
        if callable(cs_func):
            self.auto_sample()
        if cs_data is not None:
            self.data = self.to_cs_frame(cs_data)
            # Maintain threshold requirements
            self.zero_below_thresh()

        if cs_data is not None and cs_func:
            msg = "Both tabulated data and analytic functions were provided "
            msg += f"for process{self.name}"
            msg += "Tabulated data will be used"
            raise warnings.warn(msg)

        if "Elastic" in process_type and threshold > 0:
            raise ValueError("Process with a threshold > 0 cannot be elastic.")

    def to_df(self):
        """Convert to properly formatted pandas DataFrame.

        Returns:
            (:class:`pandas.DataFrame`): The process information.
        Raises:
            RuntimeError: if no data is available to produce a DataFrame
        """
        cols = ["csfile", "r1", "r2", "rtype", "B", "p1", "p2", "model_name", "params"]
        vals = [self.name+".dat", self.r1, self.r2, self.process_type,
                self.threshold, self.p1, self.p2, self.differential_process,
                self.differential_parameters]
        # Make sure cs_data is set
        if self.data is None:
            raise RuntimeError("Cross section data not provided")
        return pd.DataFrame([vals], columns=cols)

    def add_differential_parameters(self, name, params):
        """Typically differential processes require analytic forms
        due to the difficulty of extrapolation in several dimensions.
        Add free parameters into an analytic differential model here.

        Args:
            name: The name of the model for this differential process.
            params: The free parameters required for this differential model.
        """
        self.differential_process = name
        self.differential_parameters = params

    def sample_cs(self, e_points=None, grid_type="log dense", nsamples=None):
        """Sample self.cs_func on a grid of energies.

        Args:
            e_points (ArrayLike):
                Explicit energy (eV) grid points on which to sample.
            grid_type (str):

                =============== =============================================
                ``"log dense"`` sample ``nsamples`` near threshold up to 1MeV.
                   ``"simple"`` sample 0 eV - threshold - 1 MeV
                       ``None`` Behavior will automatically be determined
                                by ``auto_sample``.
                =============== =============================================

            nsamples (int):
                Override self.nsamples for this sampling call.
        """
        # Make sure there is a function to sample from.
        if not callable(self.cs_func):
            raise RuntimeError(
                "Cross section function is not callable/specified.")

        if "Elastic" in self.process_type and self.threshold > 0:
            raise RuntimeError(
                "Positive non-zero threshold passed for elastic cross section.")

        if nsamples == None and "dense" in grid_type:
            nsamples = self.nsamples
        if grid_type == None:
            grid_type == self.grid_type

        # Convert to min energy required by reaction (i.e. 0 if the threshold
        # is negative)
        B = max(0, self.threshold)
        if e_points:
            # Convert to numpy array
            e_points = np.array(e_points)
        else:
            # Create a energy points from preset grid based on the threshold.
            if grid_type == "log dense":
                e_points = np.logspace(self.SAMPLE_MIN,
                        np.log10(10**self.SAMPLE_MAX - B), nsamples) + B
            elif grid_type == "simple":
                if B > 0:
                    e_points = np.array([B, B+10**self.SAMPLE_MIN,
                        10**self.SAMPLE_MAX])
                else:
                    e_points = np.array([0, 10**self.SAMPLE_MAX])

        # Call function over grid (make sure to interpret as float)
        cross_section = vec(self.cs_func, otypes=[np.float64])(e_points)
        self.data = pd.DataFrame(np.stack((e_points, cross_section)).T,
                columns=["Energy (eV)", "Cross Section (m^2)"])
        # Sort data by energy
        self.data = self.data.sort_values(
            "Energy (eV)").reset_index(drop=True)
        # Clean up table near threshold
        self.zero_below_thresh()

    def zero_below_thresh(self):
        """Enforce the ThunderBoltz required cross section format.
        For processes with a non-zero threshold, include zero valued
        points at 0 eV and at threshold energy.

        Raises:
            RuntimeError: if there is no cross section data to format.
        """
        if self.data is None:
            raise RuntimeError(
                "Tabulated cross section data is not available for adjustment.")
        # Aliases for column names and threshold energy
        en = "Energy (eV)"
        csn = "Cross Section (m^2)"
        B = self.threshold

        if "Elastic" in self.process_type or B <= 0:
            # Extrapolate lowest cross section point to zero, if not available already
            if self.data.values[0,0] > 0:
                csv = self.data.values[0,1]
                self.data = pd.concat(
                    (pd.DataFrame([{en:0,csn:csv}]), self.data),
                    ignore_index=True)
            return

        # Otherwise, remove sampled data at and below the threshold
        self.data = self.data.loc[self.data[en] > B].copy()
        # Then add points at 0 eV and the threshold
        self.data = pd.concat((pd.DataFrame({en:[0,B], csn:[0,0]}),
            self.data), ignore_index=True)

    def require_cs(self):
        """Ensure there is some kind of cross section data associated with this
        process.

        Raises:
            RuntimeError: if there is no cross section data available.
        """
        if self.data is None:
            if not callable(self.cs_func):
                raise RuntimeError(
                    f"No cross section data available for process {self.name}")
            # Sample from analytic formula
            self.auto_sample()

    def to_cs_frame(self, a):
        """Convert any kind of two-dimensional data to a pandas
        DataFrame with columns `Energy (eV)` and `Cross Section (m^2)`
        """
        cols = ["Energy (eV)", "Cross Section (m^2)"]
        if isinstance(a, pd.DataFrame):
            return pd.DataFrame(a.values, columns=cols)
        return pd.DataFrame(a, columns=cols)

    def auto_sample(self):
        """Check if the cross section needs a dense grid or not by comparing
        simple and dense grids."""
        # Only need to check a few samples
        self.sample_cs(grid_type="log dense", nsamples=10)
        dense = copy.deepcopy(self.data)
        self.sample_cs(grid_type="simple")
        simple = copy.deepcopy(self.data)
        # Check difference through linear interpolation.
        cols = ["Energy (eV)", "Cross Section (m^2)"]
        dx, dy = dense[cols].values.T
        sx, sy = simple[cols].values.T
        s_interp = np.interp(dx, sx, sy)
        # Divide by zero mask
        msk = dy != 0
        is_diff = np.sum(np.abs((s_interp - dy)[msk]/dy[msk])) > 1e-6
        # Data is simple here, but if there is a discrepancy, use dense data
        if is_diff:
            self.sample_cs(grid_type="log dense")

    @property
    def nsamples(self):
        return self._nsamples

    @property
    def cs_func(self):
        return self._cs_func

    @nsamples.setter
    def nsamples(self, ns):
        self._nsamples = ns
        # Re-sample cross section if there is a callable cross section function
        if callable(self.cs_func):
            self.auto_sample()

    @cs_func.setter
    def cs_func(self, f):
        self._cs_func = f
        if callable(f):
            self.auto_sample()
