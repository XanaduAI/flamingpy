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
"""Monte Carlo simulations for estimating FT thresholds."""
from ft_stack.RHG import alternating_polarity, RHG_graph, RHGCode
from ft_stack.graphstates import CVGraph
from ft_stack.passive_construct import BS_network
import pytest
from ft_stack.simulations import ec_monte_carlo
import itertools as it

params = it.product([2, 3, 4], ["finite", "periodic"])


class TestBlueprint:
    @pytest.mark.parametrize("distance, boundaries", params)
    def test_all_GKP_high_squeezing(self, distance, boundaries):
        p_swap = 0
        delta = 0.001
        trials = 10
        RHG_code = RHGCode(distance, boundaries=boundaries, polarity=alternating_polarity)
        errors = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects=None)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors == 0


class TestPassive:
    @pytest.mark.parametrize("distance", [2, 3, 3])
    def test_all_GKP_high_squeezing(self, distance):
        p_swap = 0
        delta = 0.001
        trials = 10

        # The lattice with macronodes.
        RHG_macro = RHG_graph(distance, boundaries="periodic", macronodes=True, polarity=False)
        RHG_macro.index_generator()
        RHG_macro.adj_generator(sparse=True)
        # The reduced lattice.
        RHG_code = RHGCode(distance, boundaries="periodic")
        RHG_reduced = RHG_code.graph
        RHG_reduced.index_generator()
        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVGraph(RHG_reduced)

        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = BS_network(4)
        passive_objects = [RHG_macro, RHG_reduced, CVRHG_reduced, bs_network]
        errors = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects=passive_objects)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors == 0


# TODO: Tests that cross-check published results for blueprint and
# passive architecture.
