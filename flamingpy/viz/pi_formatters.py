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
"""A series of functions and classes to help processing and visualizing
multiples of sqrt(pi)/2."""

# pylint: disable=too-many-statements,singleton-comparison,too-many-lines

import math

import numpy as np
from matplotlib.ticker import Formatter


def to_pi_string(x, tex: bool = True, d=2):
    """Convert x, a multiple of sqrt(pi)/2, to a pretty string.

    If x is not a multiple of sqrt(pi)/2, return the unmodified string
    of x with `d` integers after the decimal. If tex is True, add LaTeX
    $ signs.
    """
    remainder = math.remainder(x, np.sqrt(np.pi) / 2)
    if np.isclose(remainder, 0):
        integer = round(x / (np.sqrt(np.pi) / 2))
        pref = int(integer * ((1 - integer % 2) / 2 + integer % 2))
        x_str = (not bool(round(x))) * "0" + bool(round(x)) * (
            bool(tex) * "$"
            + (not bool(1 + pref)) * "-"
            + bool(1 - abs(pref)) * str(pref)
            + r"\sqrt{\pi}"
            + (integer % 2) * "/2"
            + bool(tex) * "$"
        )
        return x_str
    return f"{x:.{d}f}"


class PiFormatter(Formatter):
    """Formatter for axis-ticks containing multiples of sqrt(pi)/2."""

    def __init__(self, tex: bool = True, d: int = 2):
        """Initialize the formatter.

        Args:
            tex: Whether to use LaTeX formatting (i.e. adding $ around the string).
            d: Number of decimals to use.
        """
        self.tex = tex
        self.d = d

    def __call__(self, x, pos=None):
        return to_pi_string(x, tex=self.tex, d=self.d)
