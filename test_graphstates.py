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
""""Unit tests for the graph state classes in graphstates.py."""
import numpy as np
import scipy.sparse as sp
import numpy.random as rand
import networkx as nx
import string

from graphstates import EGraph, CVGraph, SCZ_mat, SCZ_apply

# Some abstract graphs.
n = 10
G_complete = nx.complete_graph(n)
G_c_adj = nx.to_numpy_array(G_complete)
G_c_adj_sparse = nx.to_scipy_sparse_matrix(G_complete)


class TestEGraph:
    """Tests for EGraphs."""

    def test_init(self):
        E = EGraph(G_complete)
        E_array = nx.to_numpy_array(E)
        assert np.all(G_c_adj == E_array)

    def test_index(self):
        """Test consistency of the indexing function."""
        alph = list(string.ascii_lowercase)
        rand.shuffle(alph)

        E = EGraph()
        E.add_nodes_from(alph)
        E.index_generator()
        E.adj_generator()
        ordered_nodes = [E.to_points[i] for i in range(len(E))]

        G_adj = nx.to_numpy_array(G_complete, nodelist=ordered_nodes)
        E_adj = E.adj_mat
        assert np.all(G_adj == E_adj)

    def test_slice_coords(self):
        pass

    def test_draw(self):
        pass


class TestCVHelpers:
    """Test for CVGraph helper functions."""

    def test_SCZ_mat(self):
        """Test for the SCZ_mat function."""
        SCZ = SCZ_mat(G_c_adj)
        SCZ_sparse = SCZ_mat(G_c_adj_sparse)
        assert type(SCZ) == np.ndarray
        assert type(SCZ_sparse) == sp.coo_matrix

    def test_SCZ_apply(self):
        """Test for SCZ_apply function."""
        pass


class TestCVGraph:
    """Tests for the CVGraph class."""

    def test_init(self):
        G = CVGraph(G_complete)
        G_array = nx.to_numpy_array(G.egraph)
        H = CVGraph(EGraph(G_complete))
        H_array = nx.to_numpy_array(H.egraph)
        assert np.all(G_c_adj == H_array)
        assert np.all(G_c_adj == G_array)

    def test_states(self):
        pass

    def test_apply_noise(self):
        pass

    def test_measure_hom(self):
        pass

    def test_apply_Z_probs(self):
        pass

    def test_getters(self):
        pass
