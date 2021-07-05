# Copyright 2020 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from ft_stack.passive_construct import invert_permutation, BS_network, reduce_macro_and_simulate
from ft_stack.graphstates import CVGraph
import ft_stack.RHG as RHG
import itertools as it
import pytest
import numpy as np
from numpy.random import shuffle, default_rng as rng

code_params = it.product([2, 3, 4], [0.0001], [0, 0.5, 1])


# A macronode RHG lattice, the reduced lattice, and the delta/p-swap
# paramaters for use in this module.
@pytest.fixture(scope="module", params=code_params)
def macro_RHG(request):
    d, delta, p_swap = request.param
    boundaries = "periodic"
    # The lattice with macronodes.
    RHG_macro = RHG.RHG_graph(d, boundaries=boundaries, macronodes=True, polarity=False)
    RHG_macro.index_generator()
    RHG_macro.adj_generator(sparse=True)
    # The reduced lattice.
    RHG_code = RHG.RHGCode(d, boundaries=boundaries)
    RHG_reduced = RHG_code.graph
    RHG_reduced.index_generator()
    return delta, p_swap, RHG_macro, RHG_reduced


class TestHelpers:
    """Test helper functions for macronode reduction."""

    def test_invert_permutation(self):
        N = rng().integers(1, 100)
        random_array = rng().integers(0, 100, N)
        random_p = np.arange(N)
        shuffle(random_p)
        inverted = invert_permutation(random_p)
        presumed_unchanged_array = random_array[random_p][inverted]
        # Check that permuting and then unpermuting a random array
        # leaves it unchanged.
        assert np.array_equal(random_array, presumed_unchanged_array)

    def test_BS_network(self):
        pass


class TestReduction:
    """Test reduction of the macronode to the canonical RHG lattice."""

    def test_reduce_macro_and_simulate(self, macro_RHG):
        delta, p_swap, RHG_macro, RHG_reduced = macro_RHG
        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVGraph(RHG_reduced)
        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = BS_network(4)

        reduce_macro_and_simulate(
            RHG_macro, RHG_reduced, CVRHG_reduced, bs_network, p_swap, delta
        )
        # Check proper reduction to effective node type.
        for central_node in RHG_macro.macro.nodes:
            micronodes = RHG_macro.macro.nodes[central_node]["micronodes"]
            effective_type = RHG_reduced.nodes[central_node]["state"]
            p_count = 0
            for micronode in micronodes:
                if RHG_macro.nodes[micronode]["state"] == "p":
                    p_count += 1
            expected_type = "p" if p_count == 4 else "GKP"
            assert effective_type == expected_type
