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
""""Unit tests for decoding funcions in the decoder module."""
from ast import literal_eval
import io
import itertools as it
import sys

import networkx as nx
import numpy as np
import pytest

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders.decoder import (
    assign_weights,
    CV_decoder,
    recovery,
    check_correction,
)
from flamingpy.decoders.mwpm.matching import NxMatchingGraph


code_params = it.product([2, 3, 4], ["finite", "periodic"], [1, 0.1, 0.01], [0, 0.5, 1])


@pytest.fixture(scope="module", params=code_params)
def enc_state(request):
    """An RHGCode object and an encoded CVLayer for use in this module."""
    distance, boundaries, delta, p_swap = request.param
    DVRHG = SurfaceCode(distance=distance, boundaries=boundaries, polarity=alternating_polarity)
    RHG_lattice = DVRHG.graph
    # CV (inner) code/state
    CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
    # Noise model
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
    # Apply noise
    CVRHG.apply_noise(cv_noise)
    # Measure syndrome
    CVRHG.measure_hom("p", DVRHG.syndrome_inds)
    return DVRHG, CVRHG


@pytest.fixture(scope="module")
def match_data(enc_state):
    """The decoding graphs for use in this module."""
    graph = NxMatchingGraph(enc_state[0])
    return graph, graph.min_weight_perfect_matching()


class TestAssignWeights:
    "Test the weight assignment in decoder.py."

    def test_unit_weights(self, enc_state):
        """Check that all weights are 1 with the "unit" weight option."""
        assign_weights(enc_state[0], method="unit")
        for point in enc_state[0].syndrome_coords:
            assert enc_state[0].graph.nodes[point]["weight"] == 1

    def test_blueprint_weights(self, enc_state):
        """Test the weight options for the Xanadu's blueprint architecture."""
        weight_options = {
            "method": "blueprint",
            "integer": True,
            "multiplier": 100,
            "delta": enc_state[1]._delta,
        }
        assign_weights(enc_state[0], **weight_options)
        for point in enc_state[0].syndrome_coords:
            weight = enc_state[0].graph.nodes[point]["weight"]
            # Check that all the weights are positive integers given
            # the weight options above.
            assert isinstance(weight, int)
            assert weight >= 0


class TestDecoder:
    """A class that defines tests for the CV and MWPM decoder."""

    def test_CV_decoder(self, enc_state):
        """Test CV_decoder function."""
        CV_decoder(enc_state[0])
        bits = enc_state[1].bit_values()
        # Check that bit values have been computed, and that they are
        # either 0 or 1.
        assert bits
        for point in enc_state[0].syndrome_coords:
            bit_val_attribute = enc_state[0].graph.nodes[point]["bit_val"]
            assert bit_val_attribute in (0, 1)
            index = enc_state[0].graph.to_indices[point]
            # Check that there is not index mismatch between bit list
            # from CVLayer.bit_values() and bits stored inside the graph
            # attributes.
            assert bits[index] == bit_val_attribute
        # Check that all the stabilizer cubes point to nodes with bit
        # values
        for cube in enc_state[0].stabilizers:
            for point in cube.egraph:
                assert cube.egraph.nodes[point].get("bit_val") is not None

    def test_stab_graph(self, enc_state):
        """Check that edges in a stabilizer grapg contain the coordinates
        of the common vertex between neighbouring stabilizers."""
        stab_graph = enc_state[0].stab_graph
        odd_stabs = list(stab_graph.odd_parity_stabilizers())
        even_stabs = set(stab_graph.stabilizers) - set(odd_stabs)
        # Check odd and even cubes appropriately indexed.
        for stab in odd_stabs:
            assert stab.parity == 1
        for stab in even_stabs:
            assert stab.parity == 0
        for edge in stab_graph.edges():
            if {"high", "low"}.isdisjoint(edge):
                assert stab_graph.edge_data(*edge)["common_vertex"] is not None

    def test_matching_graph(self, match_data):
        """Test the structure of the matching graph."""
        graph = match_data[0]
        virtual_points = graph.virtual_points
        remaining_points = graph.graph.nodes - virtual_points
        n_virt = len(virtual_points)
        n_stabes = len(remaining_points)
        virtual_subgraph = graph.graph.subgraph(virtual_points)
        stabilizer_subgraph = graph.graph.subgraph(remaining_points)
        # Check that the graph formed between the virtual boundary excitations,
        # as well as the graph formed between the nodes corresponding to
        # stabilizers, are both complete graphs.
        assert nx.is_isomorphic(virtual_subgraph, nx.complete_graph(n_virt))
        assert nx.is_isomorphic(stabilizer_subgraph, nx.complete_graph(n_stabes))
        # Check the the edges in the graph formed between the boundary
        # excitations all have weight 0.
        for edge in virtual_subgraph.edges:
            assert virtual_subgraph.edges[edge]["weight"] == 0

    def test_MWPM(self, match_data):
        """Check that the matching is perfect (the set of all the nodes in the matching is the same as the set of all nodes in the matching graph)."""
        graph, matching = match_data[0], match_data[1]
        assert not {a for b in matching for a in b} - graph.graph.nodes


class TestRecovery:
    """A class that defines recovery and correction tests."""

    def test_recovery(self, enc_state, match_data):
        """Check that there remain no unsatisfied stabilizers after the
        recovery operation."""
        recovery(enc_state[0], match_data[0], match_data[1])
        odd_cubes = enc_state[0].stab_graph.odd_parity_stabilizers()
        assert len(list(odd_cubes)) == 0

    def test_corection_check(self, enc_state):
        """Check that error correction fails in case of at least one plane with odd parity."""
        old_stdout = sys.stdout
        new_stdout = io.StringIO()
        sys.stdout = new_stdout

        result = check_correction(enc_state[0], sanity_check=True)

        output = new_stdout.getvalue()
        sys.stdout = old_stdout

        surface_dict = literal_eval(output)
        boundaries = enc_state[0].boundaries
        if boundaries == ["periodic"] * 3:
            planes = ["x", "y", "z"]
        if tuple(boundaries) in set(it.permutations(["primal", "dual", "dual"])):
            where_primal = np.where(np.array(boundaries) == "primal")[0][0]
            planes = [["x", "y", "z"][where_primal]]
        failure_events = 0
        for plane in planes:
            plane_parities = np.array(surface_dict[plane])
            # Check that parity along a plane is conserved.
            assert np.all(plane_parities) or np.all(plane_parities ^ 1)
            if np.all(plane_parities ^ 1):
                failure_events += 1
        print(planes)
        print(surface_dict)
        print(failure_events, result)
        # Check that error correction fails in case of at least one
        # plane with odd parity.
        if failure_events:
            assert not result

    # def test_correct(self):
    # pass
