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
""""Unit tests for MatchingGraph classes in matching.py.
The networkx implementation is used as a reference.
"""
import itertools as it

import numpy as np
import pytest
import networkx as nx

from flamingpy.decoders.mwpm.matching import LemonMatchingGraph, NxMatchingGraph, RxMatchingGraph


# Test parameters
matching_graph_types = [LemonMatchingGraph, RxMatchingGraph]
num_nodes = range(4, 24, 4)


@pytest.fixture(scope="module", params=it.product(matching_graph_types, num_nodes))
def matching_graphs(request):
    """
    Return an instance of the given matching graph type with random weights
    and an edge between each of the num_nodes nodes.
    Also return the same graph as a NxMatchingGraph for comparison.
    """
    MatchingGraphType, num_nodes = request.param
    graph = MatchingGraphType()
    nx_graph = NxMatchingGraph()
    rng = np.random.default_rng()
    for edge in it.combinations(range(num_nodes), r=2):
        weight = rng.integers(0, 10)
        graph.add_edge(edge, weight)
        nx_graph.add_edge(edge, weight)
    return graph, nx_graph


def test_conversion(matching_graphs):
    """Test that different backends return the same graph as networkx."""
    graph, nx_graph = matching_graphs
    assert nx.is_isomorphic(graph.to_nx().graph, nx_graph.graph)


def test_matching_has_same_weight(matching_graphs):
    """Test that different backends return matching similar to networkx."""
    graph, nx_graph = matching_graphs
    matching = graph.min_weight_perfect_matching()
    nx_matching = graph.min_weight_perfect_matching()
    assert graph.total_weight_of(matching) == nx_graph.total_weight_of(nx_matching)
    assert len(matching) == len(nx_matching)
