import os

# Specify internal package data location
data_path = os.path.join(os.path.dirname(__file__), "data")
# Specify internal thunderboltz C++ source directory
src_path = os.path.join(os.path.dirname(__file__), "cpp")

# Allow main objects to be available on the package level 
from .tb import ThunderBoltz
from .tb import read
from .tb import query_tree
from .tb import plot_tree
from .input import CrossSections
from .kinetic import Process


