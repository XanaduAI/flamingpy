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
""""Unit tests for functions in the viz module."""

import math
import numpy as np
from numpy.random import default_rng as rng

import flamingpy.utils.viz as viz


def test_to_pi_string():
    """Test for the convenience function to_pi_string."""
    # Test +- sqrt(pi) and sqrt(pi)/2.
    assert viz.to_pi_string(np.sqrt(np.pi)) == "$\\sqrt{\\pi}$"
    assert viz.to_pi_string(-np.sqrt(np.pi)) == "$-\\sqrt{\\pi}$"
    assert viz.to_pi_string(np.sqrt(np.pi) / 2) == "$\\sqrt{\\pi}/2$"
    assert viz.to_pi_string(-np.sqrt(np.pi) / 2) == "$-\\sqrt{\\pi}/2$"

    # Test random odd integer multiples of sqrt(pi)/2, eccept 1 and -1
    odd_int = (2 * rng().integers(2, 25) - 1) * (-1) ** rng().integers(2)
    assert viz.to_pi_string(odd_int * np.sqrt(np.pi) / 2) == "${}\\sqrt{{\\pi}}/2$".format(odd_int)

    #  Test random even multiples of sqrt(pi).
    even_int = odd_int + 1
    assert viz.to_pi_string(even_int * np.sqrt(np.pi)) == "${}\\sqrt{{\\pi}}$".format(even_int)

    # Check everything else converted into a str.
    rand_numb = rng().random()
    rand_d = rng().integers(2, 25)
    if not np.isclose(math.remainder(rand_numb, np.sqrt(np.pi) / 2), 0):
        assert viz.to_pi_string(rand_numb, d=rand_d) == "{:.{}f}".format(rand_numb, rand_d)

    # Test for tex=False
    assert viz.to_pi_string(-np.sqrt(np.pi) / 2, tex=False) == "-\\sqrt{\\pi}/2"
