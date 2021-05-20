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
""""Unit tests for decoder funcions in decoder.py."""

import numpy as np
from graphstates import CVGraph
from decoder import (assign_weights,
                     CV_decoder,
                     decoding_graph,
                     matching_graph,
                     MWPM,
                     recovery,
                     check_correction,
                     correct)
from RHG import RHGCode, RHGCube
import networkx as nx

distance = 3
boundaries = "periodic"
RHG_code = RHGCode(distance=distance, boundaries=boundaries, polarity=True)
RHG_lattice = RHG_code.graph
# CV (inner) code/state
p_swap = 0
CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)

# Noise model
delta = 0.1
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Decoding options
decoder = {"inner": "basic", "outer": "MWPM"}
# Apply noise
CVRHG.apply_noise(cv_noise)
# Measure syndrome
CVRHG.measure_hom("p", RHG_code.syndrome_inds)

# RHG_lattice.draw()

weight_options = {
    "method": "unit",
    "integer": False,
    "multiplier": 1,
    "delta": delta,
    }
assign_weights(RHG_code, **weight_options)
# RHG_lattice.draw(label='weight')

CV_decoder(RHG_code)
G_dec = decoding_graph(RHG_code)
G_match = matching_graph(G_dec)
matching = MWPM(G_match, G_dec)


class TestAssignWeights:
    "Test the weight assignment in decoder.py."

    def test_unit_weights(self):
        for point in RHG_code.syndrome_coords:
            assert RHG_lattice.nodes[point]['weight'] == 1

    def test_blueprint_weights(self):
        weight_options = {
            "method": "blueprint",
            "integer": True,
            "multiplier": 100,
            "delta": delta,
            }
        assign_weights(RHG_code, **weight_options)
        for point in RHG_code.syndrome_coords:
            weight = RHG_lattice.nodes[point]['weight']
            assert isinstance(weight, int)
            assert weight >= 0


class TestDecoder:

    def test_CV_decoder(self):
        CV_decoder(RHG_code)
        assert CVRHG.bit_values

    def test_decoding_graph(self):
        assert G_dec.graph['title'] == 'Decoding Graph'
        assert "odd_cubes" in G_dec.graph
        boundary_points = G_dec.graph["boundary_points"]
        for node in G_dec.nodes - boundary_points - {'low', 'high'}:
            assert isinstance(G_dec.nodes[node]['stabilizer'], RHGCube)
        for edge in G_dec.edges:
            assert 'weight' in G_dec.edges[edge]
            if ('high' not in edge) and ('low' not in edge):
                assert 'common_vertex' in G_dec.edges[edge]

    def test_matching_graph(self):
        assert G_match.graph["title"] == "Matching Graph"
        assert "virtual_points" in G_match.graph
        virtual_points = G_match.graph["virtual_points"]
        remaining_points = G_match.nodes - virtual_points
        n_virt = len(virtual_points)
        n_stabes = len(remaining_points)
        virtual_subgraph = G_match.subgraph(virtual_points)
        stabilizer_subgraph = G_match.subgraph(remaining_points)
        assert nx.is_isomorphic(virtual_subgraph, nx.complete_graph(n_virt))
        assert nx.is_isomorphic(stabilizer_subgraph, nx.complete_graph(n_stabes))
        for edge in virtual_subgraph.edges:
            assert virtual_subgraph.edges[edge]['weight'] == 0

    def test_MWPM(self):
        # Check that the matching is perfect (the set of all the nodes
        # in the matching is the same as the set of all nodes in the
        # matching graph)
        assert not {a for b in matching for a in b} - G_match.nodes


class TestRecovery:

    def test_recovery(self):
        recovery(RHG_code, G_match, G_dec, matching)
        G_dec_new = decoding_graph(RHG_code, draw=False)
        odd_cubes = G_dec_new.graph["odd_cubes"]
        assert not odd_cubes

    def test_corection_check(self):
        pass
