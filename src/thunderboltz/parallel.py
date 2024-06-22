"""multiprocessing tools for multi core and multi node systems."""
import inspect
import copy
import multiprocessing as mp
import os
from os.path import join as pjoin
from pathlib import Path
import re
import shlex
import subprocess
import time
import warnings

import numpy as np
import pickle

#: SLURM parameters
slurm_options = {
    # The HPC user account
    "account": None, # required
    "job-name": None, # required
    "ntasks-per-node": None, # required, number of cores in each compute node
    "qos": None, # Quality of service
    "reservation": None, # Request particular resources
    "time": 60, # job minutes
    "array": None, # option for multi-task submission.
    "bb": None, "bbf": None, # burst buffer
    "begin": None, # scheduled runtime
    "chdir": None, # working directory of cluster process
    "clusters": None, # comma separated string of clusters
    "comment": None, # slurm script comment
    "constraint": None, # more constraints
    "deadline": None,
    "error": "job.err",
    "output": "job.out",
    "mail-user": None,
    "mail-type": None,
    # name / cpu / task / node allocation left to Manager
}

class MPRunner(object):
    """Interface for any kind of calculation that is run-able
    and set-able can be compatible with multiprocessing utilities like
    SlurmManager and DistributedPool."""
    def set_(self, **state_options):
        """Update the internal state of the object being run.

        Args:
            **state_options: Keywords arguments corresponding to attributes of the object
                being updated.
        """
        raise NotImplementedError("MPRunner must implement set_(self, **state_options) function.")

    def run(self, **run_options):
        """Run the calculation with the current settings.

        Args:
            **run_options: Keywords arguments that modify the nature of the
                way the program runs.
        """
        raise NotImplementedError("MPRunner must implement run(self, **run_options) function.")

    def get_directory(self):
        """Return the directory in which the calculation is occurring.

        Returns
            (str): The path of the directory in which the program is being run.
        """
        raise NotImplementedError("MPRunner must implement get_directory(self) function.")

    def to_pickleable(self):
        """Returns a pickle-able portion of the object sufficient to run
        the calculations."""
        raise NotImplementedError("MPRunner must implement to_picklable(self) function.")

class DistributedPool(object):
    """A multiprocessing Pool context for running
    calculations among cores with different settings.

    Args:
        runner (MPRunner): The calculation runner.
        processes (int): The number of cores to divide up the work.
    """

    def __init__(self, runner: MPRunner, processes=None):
        # Always use the fork context
        mpc = mp.get_context("fork")
        self.runner = runner
        self.pool = mpc.Pool(processes=processes)

    def submit(self, run_args={}, **set_args):
        """Submit a single job with updated key words to the pool.

        Args:
            run_args (dict): Keyword arguments to be passed to
                :meth:`~.thunderboltz.parallel.MPRunner.run`.
        **set_args: Keyword arguments passed to
            :meth:`~.thunderboltz.parallel.MPRunner.set_`
            before calling :meth:`~.thunderboltz.parallel.MPRunner.run`.
        """
        # Copy the state of the runner (this ensures that
        # the parameter state does not get overwritten while
        # the processors are all busy)
        self.runner_child = copy.deepcopy(self.runner)
        # Keyword args will be sent to the set function
        self.runner_child.set_(**set_args)
        # run_args are sent separately with the caller
        self.pool.apply_async(
            self.runner_child.run, (), run_args,
            error_callback=self.err_callback)

    def __enter__(self):
        """Pipe the context input back to caller."""
        self.pool.__enter__()
        return self

    def err_callback(self, err):
        """Print out errors that subprocesses encounter."""
        path = self.runner.get_directory()
        warnings.warn(f"Error in process: PID={os.getpid()} "
                      f"path={path}\n\n{str(err)}", RuntimeWarning)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Safely wait for all processes to finish and exit the pool."""
        # Stop accepting jobs
        self.pool.close()
        # Wait for all jobs to finish
        self.pool.join()
        # Terminate multiprocessing pool
        self.pool.__exit__(exc_type, exc_val, exc_tb)

class SlurmManager(object):
    """A python context interface for the common Slurm HPC job manager
    to run more several intensive calculations on large clusters.
    See https://slurm.schedmd.com/sbatch.html.

    Args:
        runner (MPRunner): The calculation runner.
        directory (str or None) The simulation directory. Default is
            the current working directory.
        modules (list[str]): A list of modules to be loaded by the
            HPC module system.
        mock (bool): Option to test scripts without calling a slurm
            manager.
        **options: Additional keyword arguments will be interpreted as
            SLURM parameters.

    Note:
        This job manager currently only works for clusters that either
        already have the gcc and python requirements installed on each
        compute node, or clusters that use the
        `Module System <https://hpc-wiki.info/hpc/Modules>`_ to load
        functionality.
 
        The default behavior is to accommodate the module system as it
        is common on most HPC machines. If you wish to avoid writing
        ``module load`` commands in the SLURM script, simply specify
        ``modules=[]`` in the constructor.
    """
    _pckl_file = "job.pckl"
    def __init__(self, runner: MPRunner, directory=None, modules=["python", "gcc"],
            mock=False, **options):
        if directory is None:
            directory = os.getcwd()
        # Make group directory absolute
        if not os.path.isabs(directory):
            directory = str(Path(directory).resolve())
        #: The simulation directory
        self.directory = directory
        #: The list of modules to be loaded by the HPC module system.
        self.modules = modules
        #: The :class:`~.thunderboltz.parallel.MPRunner` object.
        self.runner = runner
        # Option to test the scripts without calling a slurm manager
        self.mock = mock
        # These will hold all the updates to both slurm and the runner
        self.input_sets = []
        #: Store references to the slurm job numbers after jobs are submitted
        self.job_ids = []
        #: The SLURM sbatch options
        self.options = slurm_options.copy()
        self.set_(**options)
        # Make sure cores are specified
        msg = (f"{self.__class__.__name__} requires "
                "the number of cores available at each node "
                "to be passed as ntasks-per-node.")
        if not "ntasks-per-node" in self.options:
            raise RuntimeError(msg)
        if not self.options["ntasks-per-node"]:
            raise RuntimeError(msg)

        # Assume fixed number of cores, store value
        self.ncores = self.options["ntasks-per-node"]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit and run sbatch for all the inputs"""
        ntasks = len(self.input_sets)
        nodes_required = int(np.ceil(ntasks / self.ncores))
        for start in range(nodes_required):
            # Isolate the inputs for this node
            inputs = self.input_sets[
                int(start*self.ncores):int(min((start+1)*self.ncores, ntasks))
            ]
            # Take any slurm updates from first input
            self.set_(**inputs[0])
            # Create a directory for slurm files of this node
            node_dir = pjoin(self.directory, f"NODE_{start+1}")
            os.makedirs(node_dir)
            with visit(node_dir):
                # Write pickle of input at the directory level
                with open(self._pckl_file, "wb") as f:
                    pickle.dump((inputs, self.runner), f)
                # Write the slurm script
                self.write_slurm_script()
                # Execute slurm call
                self.sbatch()

    def submit(self, run_args={}, **settings):
        """Add a set of parameter updates to the job queue.
        Slurm is not invoked until the context is exited.

        Args:
            run_args (dict): Keyword arguments to be passed to
                :meth:`~.thunderboltz.parallel.MPRunner.run`.
            **settings: Keyword arguments passed to the
                :meth:`~.thunderboltz.parallel.MPRunner.set_`.
                before calling :meth:`~.thunderboltz.parallel.MPRunner.run`.
        """
        # Resolve relative paths
        if not os.path.isabs(settings["directory"]):
            settings["directory"] = str(Path(settings["directory"]).resolve())

        combined_settings = {"run_args": run_args}
        combined_settings.update(**settings)
        self.input_sets.append(combined_settings.copy())

    def sbatch(self):
        """Call slurm with current settings."""
        if self.mock:
            # Just call the script and block
            self.mock_run()
            return
        # Otherwise, submit to slurm
        cmds = ["date | tee -a jobnum", "sbatch job_script.sh | tee -a jobnum"]
        subprocess.run(cmds[0], shell=True)
        sub_string = subprocess.check_output(cmds[1], shell=True).decode()
        self.job_ids.append(int(sub_string.split()[-1]))

    def set_(self, **options):
        """Update slurm manager options.

        Args:
            **options: SLURM settings.
        """
        self.options.update({k:v
            for k, v in options.items()
            if k in slurm_options})

    def write_slurm_script(self, path=None, script_name=None):
        """Write the SLURM batch script."""
        if script_name is None:
            script_name = "job_script.sh"
        if path is None:
            path = os.getcwd()
        with open(pjoin(path, script_name), "w") as f:
            f.write("#!/bin/bash\n")
            for k, v in self.options.items():
                if not v:
                    continue
                f.write(f"#SBATCH --{k}={v}\n")

            if not self.mock:
                # Only load modules if this is a real slurm run
                f.write("\n")
                for mod in self.modules:
                    f.write(f"module load {mod}\n")
                f.write("\n\n")

            python_script = self.process_batch_script()
            # Declare pyscript var
            f.write("SCRIPT=$(cat<<END\n")
            f.write(python_script)
            f.write("END\n")
            f.write(")\n\n")
            f.write("python -c \"$SCRIPT\"\n")

    def process_batch_script(self):
        """Inspect the batch script below and process it for use
        in sbatch."""
        source = inspect.getsource(self.batch_script)
        source = "\n".join(source.split("\n")[4:])
        # Inject ncores
        source = source.replace("ncores = None", f"ncores = {self.ncores}")
        # Inject pckl name
        source = source.replace("self._pckl_file", f"\"{self._pckl_file}\"")
        # Remove indents
        source = re.sub("^[^\S\n]{8}", "", source, flags=re.MULTILINE)
        return source

    def batch_script(self):
        """The SLURM job script. This does not get called in the
        parent process, but instead the source code is invoked in the
        sbatch script/command for subprocess startup."""
        import pickle
        from thunderboltz.parallel import DistributedPool
        from thunderboltz.parallel import slurm_options
        # number of cores and directory will be injected in sbatch
        ncores = None
        # Load in batch parameters
        with open(self._pckl_file, "rb") as f:
            updates, runner = pickle.load(f)
        # Run settings through distruted processing interface
        with DistributedPool(runner, processes=ncores) as pool:
            for update in updates:
                # Separate submit args
                submit_args = {k: v for k, v in update.items()
                    if k not in slurm_options}
                pool.submit(**submit_args)

    def mock_run(self):
        """Act as a compute node and test the job scripts sequentially."""
        os.system(f"chmod +x job_script.sh")
        os.system(f"./job_script.sh")

    def has_active(self):
        """Check whether any submitted jobs are still pending or running.

        Returns:
            (bool): ``True`` if there are still jobs that are pending or
                running. ``False`` otherwise.
        """
        has_active = False
        for job_id in self.job_ids:
            is_active = False
            try:
                cmd = f"scontrol show job {job_id}"
                sub_string = subprocess.check_output(
                    shlex.split(cmd)).decode()
                is_active = "PENDING" in sub_string or "RUNNING" in sub_string
            except subprocess.CalledProcessError as e:
                is_active = False
            has_active = has_active or is_active
        return has_active

    def has_pending(self):
        """Check whether any submitted jobs are still pending.

        Returns:
            (bool): ``True`` if there are still jobs that are pending.
                ``False`` otherwise.
        """
        has_pending = False
        for job_id in self.job_ids:
            is_active = False
            try:
                cmd = f"scontrol show job {job_id}"
                sub_string = subprocess.check_output(
                    shlex.split(cmd)).decode()
                is_pending = "PENDING" in sub_string
            except subprocess.CalledProcessError as e:
                is_pending = False
            has_pending = has_pending or is_pending
        return has_pending

    def join(self):
        """Wait for all slurm jobs to finish."""
        while self.has_active():
            time.sleep(2)

class visit(object):
    """Directory environment context manager."""
    def __init__(self, path):
        self.original_path = os.getcwd()
        self.path = path

    def __enter__(self):
        os.chdir(self.path)
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.original_path)

def setup_(i=0, subdir=None, base_directory=None):
    """Call from a run function such that it makes a directory
    in simulation/output with the name of that test function. Increment
    `i` if this function is being wrapped in another intermediary."""
    if base_directory is None:
        base_directory = os.getcwd()
    # Create the path to the test file
    fpath = pjoin(base_directory, inspect.stack()[1+i][3])
    if subdir is not None:
        fpath = pjoin(fpath, str(subdir))
    # Overwrite directory if necessary
    if os.path.isdir(fpath):
        shutil.rmtree(fpath)
    os.makedirs(fpath)
    return fpath

