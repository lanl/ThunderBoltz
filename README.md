
# ThunderBoltz

[![CI/CD](https://github.com/lanl/ThunderBoltz/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/lanl/ThunderBoltz/actions/workflows/ci-cd.yml)
[![PyPI Latest Release](https://img.shields.io/pypi/v/thunderboltz.svg)](https://pypi.org/project/thunderboltz/)
[![License - BSD 3-Clause](https://img.shields.io/pypi/l/thunderboltz.svg)](https://github.com/lanl/ThunderBoltz/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/thunderboltz/badge/?version=latest)](https://thunderboltz.readthedocs.io/en/latest/?badge=latest)


ThunderBoltz is a 0D Direct Simulation Monte-Carlo code for the calculation of transport
coefficients and 0D gas modeling. This repository contains the ThunderBoltz C++ source code (C22063)
and ThunderBoltz Python API code (O4674), which are released with a GPL v3 license. 
The overarching goal of this project is to provide
transport coefficients for electrons, ions and neutrals, and to model energy transfer
between charged particles and neutral gases. These processes are not modeled by
commonly used Boltzmann solvers such as Bolsig+ or MagBoltz since these solvers only
track electron transport (either via trajectories in the MC case, or through modeling
of the EEDF) due to collisions with a stationary neutral background.

To reference this code and work, please cite the arXiv article:
https://arxiv.org/abs/2310.07913

## Documentation 
The ThunderBoltz C++ source code manual is found in the home directory
"cpp_manual.pdf", and the Python API code manual is found in the home directory
"api_manual.pdf". The online API documentation can be found at
[here](https://thunderboltz.readthedocs.io).

## Installation

### Pip Installation
The Python API can be installed via `pip`.
```
pip install thunderboltz
```
This package ships with the C++ ThunderBoltz source and requires only a gcc compiler
to run complete calculations.

### Clone the repository

The code can also be downloaded from the repository.
```
git clone https://github.com/lanl/ThunderBoltz.git
```
You may need to set up SSH keys in order to access Github. See the
[Gitlab SSH Guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)
to set up access to Github repositories.

The basic ThunderBoltz functionality will be available either
as an executable in `bin/thunderboltz.bin` or can be compiled from the
source in `src/thunderboltz/cpp`. To install the Python interface, run the
`install.sh` script from the root directory.
```
./install.sh
```

This will upgrade [pip](https://pypi.org/project/pip)
and install specific versions of Python packages,
so create an environment if you are concerned with Python package overwrite.

## Usage

For ease of use, we recommend using the Python interface to setup, run,
and process simulations ([see here](#using-the-python-wrapper)).
Alternatively, one can use the stand-alone C++
code to run calculations (below).

### Manual compilation and execution (ThunderBoltz code only)

The C++ source files are located in `src/thunderboltz/cpp/`.
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
cp ../src/thunderboltz/cpp/* .
# Copy example input files over
cp -r ../indecks/N2/* .
# Compile
g++ -std=c++17 -o thunderboltz.bin DSMC0D.cpp
# Run
./thunderboltz.bin N2vib.in
```
Simulation parameters can be adjusted in the `N2vib.in` file for this example.

### Using the Python wrapper

We also include a Python wrapper that allows for automated compilation
and extended input/output processing. To install the Python API from
the online Python index (PYPI), run
```
pip install thunderboltz
```
To install the Python API from the
repository, run the `install.sh` script from the root directory:
```
./install.sh
```
This will upgrade `pip` and install specific versions of Python packages,
so create an environment if you are concerned with Python package overwrite.

See `run.py` for implementations of the tests found in the paper
using the API wrapper.

## Testing and code sync

Make sure a version of pytest is installed:
```
pip install pytest
```
To run unit tests on the Python `thunderboltz`,
make sure `thunderboltz` has been installed, and 
```
python -m pytest testing/*.py --verbose
```

To recompile the standalone binary version, run
```
g++ -std=c++17 src/thunderboltz/cpp/DSMC0D.cpp -o bin/thunderboltz.bin -Wall -Werror -Wsign-compare
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
- Optional particle reweigting


