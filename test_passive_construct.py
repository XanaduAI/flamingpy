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
from passive_construct import invert_permutation, BS_network, reduce_macro_and_simulate
from graphstates import CVGraph, SCZ_apply
from GKP import GKP_binner, Z_err_cond
import RHG
from decoder import correct
import itertools as it
import pytest

code_params = it.product([2], [0.0001], [0])


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
    def test_invert_permutation(self):
        pass

    def test_BS_network(self):
        pass


class TestReduction:
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
        for node in RHG_reduced:
            if p_swap == 0:
                assert RHG_reduced.nodes[node]["state"] == "GKP"
            elif p_swap == 1:
                assert RHG_reduced.nodes[node]["state"] == "p"
