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
"""Unit tests for Monte Carlo simulations for estimating FT thresholds."""
import itertools as it
import pytest

from flamingpy.codes import alternating_polarity, RHG_graph, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import BS_network
from flamingpy.simulations import ec_monte_carlo


params = it.product([2, 3, 4], ["finite", "periodic"])


class TestBlueprint:
    """A class with members to test Monte Carlo simulations for FT threshold estimations for Xanadu's blueprint architecture."""

    @pytest.mark.parametrize("distance, boundaries", params)
    def test_all_GKP_high_squeezing(self, distance, boundaries):
        """Tests Monte Carlo simulations for FT threshold estimation of a system with zero swap-out probability and high squeezing."""
        p_swap = 0
        delta = 0.001
        trials = 10
        RHG_code = SurfaceCode(distance, boundaries=boundaries, polarity=alternating_polarity)
        errors_py = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects=None)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors_py == 0


class TestPassive:
    """A class with members to test Monte Carlo simulations for FT threshold estimations for Xanadu's passive architecture."""

    @pytest.mark.parametrize("distance", [2, 3, 3])
    def test_all_GKP_high_squeezing(self, distance):
        """Tests Monte Carlo simulations for FT threshold estimation of a system with zero swap-out probability and high squeezing."""
        p_swap = 0
        delta = 0.001
        trials = 10

        # The lattice with macronodes.
        RHG_macro = RHG_graph(distance, boundaries="periodic", macronodes=True, polarity=False)
        RHG_macro.index_generator()
        RHG_macro.adj_generator(sparse=True)
        # The reduced lattice.
        RHG_code = SurfaceCode(distance, boundaries="periodic")
        RHG_reduced = RHG_code.graph
        RHG_reduced.index_generator()
        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVLayer(RHG_reduced)

        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = BS_network(4)
        passive_objects = [RHG_macro, RHG_reduced, CVRHG_reduced, bs_network]
        errors_py = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects=passive_objects)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors_py == 0


# TODO: Tests that cross-check published results for blueprint and
# passive architecture.
