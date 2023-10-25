# ThunderBoltz
Direct simulation Monte Carlo ThunderBoltz C++ and Python API codes 



# ThunderBoltz

ThunderBoltz is a 0D Direct Simulation Monte-Carlo code for the calculation of transport
coefficients and 0D gas modeling. This repository contains the ThunderBoltz C++ source code (C22063)
and ThunderBoltz Python API code (O4674), which are released with a GPL v3 license. 
The overarching goal of this project is to provide
transport coefficients for electrons, ions and neutrals, and to model energy transfer
between charged particles and neutral gases. These processes are not modeled by
commonly used Boltzmann solvers such as Bolsig+ or MagBoltz since these solvers only
track electron transport (either via trajectories in the MC case, or through modeling
of the EEDF) due to collisions with a stationary neutral background.

## Installation

### Clone the repository

Clone a local repository with the code. You may need to set up
SSH keys in order to access gitlab. See [here](https://docs.gitlab.com/ee/user/ssh.html).
```
git clone git@gitlab.com/Mczammit/thunderboltz.git
```
The basic ThunderBoltz functionality will be available either
as an executable in `bin/thunderboltz.bin` or can be compiled from the
source in `src/thunderboltz`.

## Usage

For ease of use, we recommend using the Python interface to setup, run,
and process simulations ([see here](#using-the-pytb-python-wrapper)).
Alternatively, one can use the stand-alone C++
code to run calculations (below).

### Manual compilation and execution (ThunderBoltz code only)

ThunderBoltz requires a g++ of clang compiler and should be compiled
from source directories as
```
g++ -std=c++17 -o thunderboltz.bin DSMC0D.cpp
```
Then run with
```
./thunderboltz.bin inputfile.in
```
to use a manually constructed indeck file. The code is maintained with the
standard `-Wall` and `-Werror` compiler options.

Here is an example of how to run a simple ThunderBoltz calculation.

```
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
```

### Using the pytb Python wrapper

We also include a python wrapper that allows for automated compilation
and extended input/output processing. Run the `install.sh` script from the
root directory to install the python package.
```
./install.sh
```
This will upgrade `pip` and install specific versions of python packages,
so create an environment if you are concerned with python package overwrite.

See `run.py` for implementations of the tests found in the paper
using the pytb wrapper.


## Testing and code sync

Make sure a version of pytest is installed:
```
pip install pytest
```
To run unit tests on the `pytb` package, run
```
python -m pytest testing/*.py --verbose
```

Currently, pytb has its own version of the Thunderboltz source code in
`src/pytb/thunderboltz`. To update it with the main version, run
```
cp src/thunderboltz/* src/pytb/thunderboltz
```

To recompile the standalone binary version, run
```
g++ -std=c++17 src/thunderboltz/DSMC0D.cpp -o bin/thunderboltz.bin -Wall -Werror -Wsign-compare
```

## License

ThunderBoltz is distributed under a GNU GPLv3 license. See
[http://www.gnu.org/licenses/gpl-3.0.en.html](http://www.gnu.org/licenses/gpl-3.0.en.html)
for more details.

## Roadmap
- Regression testing of internal ThunderBoltz code
- Diffusion post processing
- Random deletion after ionization option
- Addition of photon transport
