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

# pylint: disable=protected-access,no-self-use
from datetime import datetime
import logging
from multiprocessing.sharedctypes import Value
import string

import random
import networkx as nx
import numpy as np
import pytest

from flamingpy.codes.graphs import EGraph

now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)
random.seed(int_time)

# A NetworkX random graph of size N for use in this module.
N = 20
N3D = 10


@pytest.fixture(scope="module", params=[N])
def random_graph(request):
    """A convenience function to initialize random EGraphs for use in this
    module."""
    n = request.param
    G = nx.fast_gnp_random_graph(n, 0.5)
    G_adj = nx.to_numpy_array(G)
    G_adj_sparse = nx.to_scipy_sparse_array(G)
    return EGraph(G), G_adj, G_adj_sparse


@pytest.fixture(scope="module", params=[N3D])
def random_graph_3D(request):
    """A convenience function to initialize random EGraphs for use in this
    module."""
    n = request.param
    G = nx.Graph()
    random_nodes = [
        (random.randint(0, 100), random.randint(0, 100), random.randint(0, 100)) for _ in range(n)
    ]
    G.add_nodes_from(random_nodes)

    random_edges = [ed for ed in nx.non_edges(G) if random.randint(0, 10) >= 7]
    G.add_edges_from(random_edges)
    return G


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
        np.random.shuffle(alph)
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

    def test_add_qubit(self, random_graph_3D):
        """Test the add_qubit function on a random EGraph."""

        E = EGraph(random_graph_3D)
        E.adj_generator()
        max_ind_old = max(E.to_points.keys())
        n_old = E.number_of_nodes()
        E.add_qubit()
        E.add_qubit()
        assert n_old + 2 == E.number_of_nodes()
        assert max(E.to_points.keys()) == max_ind_old + 2
        assert E.adj_mat is None

        E = EGraph(random_graph_3D)
        E.index_generator()
        E.add_qubit(qubit=(100, 101, -400))
        assert (100, 101, -400) in E

        with pytest.raises(Exception) as e:
            E.add_qubit(qubit=(1, 1, 1, 1))
        assert e.type == ValueError

        with pytest.raises(Exception) as e:
            E.add_qubit(qubit=0)
        assert e.type == TypeError

        E = EGraph(random_graph_3D)
        n_edges_old = E.number_of_edges()
        E.index_generator()
        E.add_qubit(neighbors=[0, 1])
        assert n_edges_old + 2 == E.number_of_edges()

        E = EGraph(random_graph_3D)
        n_edges_old = E.number_of_edges()
        new_neighs = [neigh for neigh in E.nodes() if random.randint(0, 10) > 7]
        E.add_qubit(neighbors=new_neighs)
        assert n_edges_old + len(new_neighs) == E.number_of_edges()

    def test_remove_qubit(self, random_graph_3D):
        """Test the remove_qubit function on a random EGraph."""

        E = EGraph(random_graph_3D)
        E.index_generator()
        to_points_old = E.to_points.copy()
        n_old = E.number_of_nodes()
        E.remove_qubit(1)
        assert n_old == E.number_of_nodes() + 1
        assert E.adj_mat is None
        assert to_points_old[2] == E.to_points[2]

        E = EGraph(random_graph_3D)
        n_old = E.number_of_nodes()
        nodes = list(E.nodes())
        E.remove_qubit(nodes[0])
        E.remove_qubit(nodes[1])
        assert n_old == E.number_of_nodes() + 2
        assert E.adj_mat is None

        E = EGraph(random_graph_3D)
        with pytest.raises(Exception) as e:
            E.remove_qubit("s")
        assert e.type == TypeError

    def test_add_qubit_macro(self, random_graph_3D):
        """Test add_qubit if graph is macronized."""

    def test_remove_qubit_macro(self, random_graph_3D):
        """Test remove_qubit if graph is macronized."""

        # def test_macronode(self):
        # pass

        # def test_slice_coords(self):
        # pass
