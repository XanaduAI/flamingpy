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
"""Unit tests for passive_construct module members."""


import itertools as it

import numpy as np
from numpy.random import shuffle, default_rng as rng
import pytest

from flamingpy.codes import SurfaceCode
from flamingpy.cv.ops import invert_permutation
from flamingpy.noise import CVMacroLayer

code_params = it.product([2, 3, 4], [0.0001], [0, 0.5, 1], ["open", "periodic"], ["primal", "dual"])


@pytest.fixture(scope="module", params=code_params)
def code_and_noise(request):
    """Defines a macronode RHG lattice, the reduced lattice, and the
    delta/p-swap paramaters for use in this module."""
    d, delta, p_swap, boundaries, ec = request.param
    # The reduced lattice.
    RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)
    return delta, p_swap, RHG_code


class TestHelpers:
    """Test helper functions for macronode reduction."""

    def test_invert_permutation(self):
        """Check that permuting and then unpermuting a random array leaves it
        unchanged."""
        N = rng().integers(1, 100)
        random_array = rng().integers(0, 100, N)
        random_p = np.arange(N)
        shuffle(random_p)
        inverted = invert_permutation(random_p)
        presumed_unchanged_array = random_array[random_p][inverted]
        assert np.array_equal(random_array, presumed_unchanged_array)

    # def test_BS_network(self):
    # pass


class TestReduction:
    """Test reduction of the macronode to the canonical RHG lattice."""

    def test_reduce_macro_and_simulate(self, code_and_noise):
        """Test the reduce_macro_and_simulate function."""
        delta, p_swap, RHG_code = code_and_noise
        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        CV_macro = CVMacroLayer(RHG_code, delta=delta, p_swap=p_swap)
        CV_macro.apply_noise()
        # Check proper reduction to effective node type.
        RHG_macro = CV_macro.egraph
        RHG_reduced = CV_macro.reduced_graph
        for central_node in RHG_macro.macro_to_micro:
            micronodes = RHG_macro.macro_to_micro[central_node]
            effective_type = RHG_reduced.nodes[central_node]["state"]
            p_count = 0
            for micronode in micronodes:
                if RHG_macro.nodes[micronode]["state"] == "p":
                    p_count += 1
            expected_type = "p" if p_count == 4 else "GKP"
            assert effective_type == expected_type
