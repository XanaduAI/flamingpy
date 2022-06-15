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
""""Unit tests for the graph state classes in the graphstates module."""

# pylint: disable=protected-access

import string

import networkx as nx
import numpy as np
import numpy.random as rand
import pytest

from flamingpy.codes.graphs import EGraph

# A NetworkX random graph of size N for use in this module.
N = 20


@pytest.fixture(scope="module", params=[N])
def random_graph(request):
    """A convenience function to initialize random EGraphs for use in this
    module."""
    n = request.param
    G = nx.fast_gnp_random_graph(n, 0.5)
    G_adj = nx.to_numpy_array(G)
    G_adj_sparse = nx.to_scipy_sparse_array(G)
    return EGraph(G), G_adj, G_adj_sparse


def noise_model(delta, order):
    """A convenience function for outputing a noise model dictionary."""
    return {"noise": "grn", "delta": delta, "sampling_order": order}


class TestEGraph:
    """Tests for EGraphs."""

    def test_init(self, random_graph):
        """Check that the adjacency matrix of the random graph matches the
        adjancency matrix of the EGraph."""
        E = EGraph(random_graph[0])
        E_array = nx.to_numpy_array(E)
        assert np.all(random_graph[1] == E_array)

    def test_index(self, random_graph):
        """Tests a graph with nodes from a shuffled alphabet."""
        alph = list(string.ascii_lowercase)
        rand.shuffle(alph)
        reduced_alph = alph[:N]
        label_dict = dict(zip(range(N), reduced_alph))
        H = nx.relabel_nodes(random_graph[0], label_dict)
        E = EGraph(H)
        E.index_generator()
        # assert list(E.to_indices.keys()) == sorted(reduced_alph)
        # Adjacency matrix of H with rows/columns arranged according
        # to the sorted alphabet should equal to the adjacency matrix
        # of E.
        H_adj = nx.to_numpy_array(H, nodelist=sorted(reduced_alph))
        E.adj_generator(sparse=False)
        E_adj = E.adj_mat
        assert np.array_equal(H_adj, E_adj)

    # def test_macronode(self):
    # pass

    # def test_slice_coords(self):
    # pass
