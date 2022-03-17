# Copyright 2022 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Threshold estimations for measurement-based implementation of quantum error
correcting codes using GKP qubits."""

import sys
import platform
import os

import numpy
import scipy
import networkx
import retworkx
import matplotlib
import pandas

import flamingpy.cpp.lemonpy as lp
import flamingpy.cpp.cpp_mc_loop as cmc

from ._version import __version__


__all__ = ["version", "about"]


def version():
    r"""
    Version number of FlamingPy.
    Returns:
      str: package version number
    """
    return __version__


def about():
    """Prints the installed version numbers for FlamingPy and its dependencies,
    and some system info.

    Please include this information in bug reports.
    """

    # a QuTiP-style infobox
    print(
        "\nFlamingPy is a cross-platform Python library with a variety of backends for efficient simulations of error correction in fault-tolerant quantum computers."
    )
    print("\nCopyright 2022 Xanadu Quantum Technologies Inc.\n")

    print("Platform info:               {}".format(platform.platform()))
    print("Installation path:           {}".format(os.path.dirname(__file__)))
    print("Python version:              {}.{}.{}".format(*sys.version_info[0:3]))
    print("FlamingPy version:           {}".format(__version__))
    print("Numpy version:               {}".format(numpy.__version__))
    print("Scipy version:               {}".format(scipy.__version__))
    print("NetworkX version:            {}".format(networkx.__version__))
    print("RetworkX version:            {}".format(retworkx.__version__))
    print("Matplotlib version:          {}".format(matplotlib.__version__))
    print("Pandas version:              {}".format(pandas.__version__))
    print("lemonpy shared object:       {}".format(lp))
    print("cpp_mc_loop shared object:   {}".format(cmc))
