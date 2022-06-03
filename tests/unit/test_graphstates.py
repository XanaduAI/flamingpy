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

import string

import networkx as nx
import numpy as np
import numpy.random as rand
from numpy.random import default_rng as rng
import pytest
import scipy.sparse as sp

from flamingpy.codes import SurfaceCode
from flamingpy.codes.graphs import EGraph
from flamingpy.cv.ops import CVLayer, SCZ_mat, SCZ_apply
import matplotlib

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

    def test_draw(self):
        E = random_graph[0]
        f, a = E.draw()
        assert issubclass(type(f), matplotlib.figure.Figure)
        assert issubclass(type(a), matplotlib.axes.Axes)


class TestCVHelpers:
    """Tests for CVLayer helper functions."""

    @pytest.mark.parametrize(
        "sparse, expected_out_type", [(True, sp.coo_matrix), (False, np.ndarray)]
    )
    def test_SCZ_mat_sparse_param(self, random_graph, sparse, expected_out_type):
        """Tests the SCZ_mat function outputs sparse or dense arrays."""
        SCZ = SCZ_mat(random_graph[2], sparse=sparse)
        assert isinstance(SCZ, expected_out_type)

    def test_SCZ_mat(self, random_graph):
        """Tests the SCZ_mat function."""
        SCZ = SCZ_mat(random_graph[1])
        SCZ_sparse = SCZ_mat(random_graph[2])
        # Check if SCZ_mat adjusts type of output matrix based on
        # type of input.
        assert isinstance(SCZ, np.ndarray)
        assert isinstance(SCZ_sparse, sp.coo_matrix)
        # Check that structure of SCZ matrix is correct.
        for mat in (SCZ, SCZ_sparse.toarray()):
            assert np.array_equal(mat[:N, :N], np.identity(N))
            assert np.array_equal(mat[:N, N:], np.zeros((N, N)))
            assert np.array_equal(mat[N:, :N], random_graph[1])
            assert np.array_equal(mat[N:, N:], np.identity(N))

    @pytest.mark.parametrize("one_shot", [True, False])
    @pytest.mark.parametrize("n", [1, 2])
    def test_SCZ_apply(self, random_graph, one_shot, n):
        """Test SCZ matrix application."""

        adj = random_graph[1]
        SCZ = SCZ_mat(adj)
        N = adj.shape[0]
        quads_shape = [N * 2] * n
        quads = np.random.rand(*quads_shape)

        if n == 1:
            expected_quads = SCZ_mat(adj).dot(quads)
        else:
            expected_quads = SCZ.dot(SCZ.dot(quads).T).T

        new_quads = SCZ_apply(adj, quads, one_shot=one_shot)

        assert np.allclose(new_quads, expected_quads)


class TestCVLayer:
    """Tests for functions in the CVLayer class."""

    def test_empty_init(self, random_graph):
        """Test the instantiation of CVLayer from codes and EGraphs."""
        G = CVLayer(random_graph[0], states=None)
        G_array = nx.to_numpy_array(G.egraph)
        H = CVLayer(EGraph(random_graph[0]), states=None)
        H_array = nx.to_numpy_array(H.egraph)
        # Check that the _N attribute is populated with the number
        # of nodes in the random graph.
        assert G._N == len(random_graph[0])
        # Check that instantiating a CVLayer with an EGraph or with
        # a regular NetworkX graph has the same effect.
        assert np.array_equal(random_graph[1], H_array)
        assert np.array_equal(random_graph[1], G_array)
        CVRHG = CVLayer(SurfaceCode(2))
        assert CVRHG._N == len(SurfaceCode(2).graph)

    def test_all_GKP_init(self, random_graph):
        """Test the all-GKP initialization of EGraph."""
        G = CVLayer(random_graph[0])
        n = len(random_graph[0])
        for node in G.egraph:
            assert G.egraph.nodes[node]["state"] == "GKP"
        assert len(G._states["p"]) == 0
        assert len(G.p_inds) == 0
        assert np.array_equal(G._states["GKP"], np.arange(n))
        assert np.array_equal(G.GKP_inds, np.arange(n))

    @pytest.mark.parametrize("p_swap", [0, 0.99 * rng().random() + 0.01, 1])
    def test_hybridize(self, random_graph, p_swap):
        """Test whether CVLayer properly populates p-squeezed states for non-
        zero p-swap."""
        n = len(random_graph[0])
        # Test all-p case
        G = CVLayer(random_graph[0], p_swap=1)
        for node in G.egraph:
            assert G.egraph.nodes[node]["state"] == "p"
        assert len(G._states["GKP"]) == 0
        assert np.array_equal(G._states["p"], np.arange(n))
        # Test various swap-out probabilities; mean p population
        # should be close to p_swap parameter, within tolerance.
        p_list = []
        for _ in range(1000):
            G = CVLayer(random_graph[0], p_swap=p_swap)
            p_list += [len(G._states["p"]) / n]
        p_prob = sum(p_list) / 1000
        assert np.isclose(p_prob, p_swap, rtol=1e-1)

    def test_state_indices(self, random_graph):
        """Test that _states, p_inds, and GKP_inds are populated with the
        correct indices."""
        n = len(random_graph[0])
        num_ps = rng().integers(n)
        p_inds = rng().choice(n, num_ps, replace=False)
        gkp_inds = list(set(np.arange(n)) - set(p_inds))
        G = CVLayer(random_graph[0], states={"p": p_inds})
        assert np.array_equal(G._states.get("p"), p_inds)
        assert np.array_equal(G._states.get("GKP"), gkp_inds)
        assert np.array_equal(G.p_inds, p_inds)
        assert np.array_equal(G.GKP_inds, gkp_inds)

    @pytest.mark.parametrize("order", ["initial", "final", "two-step"])
    def test_apply_noise(self, random_graph, order):
        """Check _delta, _sampling_order attributes with default noise
        model."""
        G = CVLayer(random_graph[0])
        G.apply_noise()
        assert G._delta == 0.01
        assert G._sampling_order == "initial"
        # Check supplied noise model
        H = CVLayer(random_graph[0])
        delta = rng().random()
        H.apply_noise(noise_model(delta=delta, order=order))
        assert H._delta == delta
        assert H._sampling_order == order

    def test_grn_model(self, random_graph):
        """Compare expected noise objects (quadratures, covariance matrix) with
        those obtained through supplying the grn_model dictionary with
        different parameters."""
        delta = rng().random()
        model_init = noise_model(delta, "initial")
        model_fin = noise_model(delta, "final")
        model_two_step = noise_model(delta, "two-step")

        n = len(random_graph[0])
        G = CVLayer(random_graph[0])
        H = CVLayer(random_graph[0], p_swap=1)
        G.apply_noise(model_init)
        H.apply_noise(model_init)
        init_noise_all_GKP = np.full(2 * n, (delta / 2) ** 0.5, dtype=np.float32)
        init_noise_all_p = np.array(
            [1 / (2 * delta) ** 0.5] * n + [(delta / 2) ** 0.5] * n, dtype=np.float32
        )
        assert np.array_equal(G._init_noise, init_noise_all_GKP)
        assert np.array_equal(H._init_noise, init_noise_all_p)

        G.apply_noise(model_fin)
        H.apply_noise(model_fin)
        noise_cov_all_GKP = SCZ_apply(G._adj, np.diag(init_noise_all_GKP) ** 2)
        noise_cov_all_p = SCZ_apply(H._adj, np.diag(init_noise_all_p) ** 2)
        assert np.array_equal(G._noise_cov.toarray(), noise_cov_all_GKP)
        assert np.array_equal(H._noise_cov.toarray(), noise_cov_all_p)

        G.apply_noise(model_two_step)
        H.apply_noise(model_two_step)
        assert np.max(H._init_quads[:n]) <= 2 * np.sqrt(np.pi)
        assert np.min(H._init_quads[:n]) >= 0
        assert np.isclose(np.max(G._init_quads[:n]), np.sqrt(np.pi))
        assert np.isclose(np.min(G._init_quads[:n]), 0)

    @pytest.mark.parametrize("order", ["initial", "final"])
    def test_measure_hom(self, random_graph, order):
        """Test closeness of average homodyne outcomes value to 0 in the all-
        GKP high-squeezing limit."""
        n = len(random_graph[0])
        G = CVLayer(random_graph[0])
        delta = 0.0001
        G.apply_noise(noise_model(delta=delta, order=order))
        G.measure_hom("p")
        G.measure_hom("q")
        outcomes_p = []
        outcomes_q = []
        for node in G.egraph:
            outcomes_p += [G.egraph.nodes[node]["hom_val_p"]]
            outcomes_q += [G.egraph.nodes[node]["hom_val_q"]]
        assert np.isclose(sum(outcomes_p) / n, 0, atol=1e-1)
        assert np.isclose(sum(outcomes_q) / n, 0, atol=1e-1)
        # TODO: test two-step sampling

    # def test_eval_Z_probs(self, random_graph):
    #     """Test that p_phase and p_phase_cond attribute get populated when
    #     phase error probabilities are evaluated."""
    #     pass
    #     G = CVLayer(random_graph[0])
    #     G.apply_noise(noise_model(delta=rng().random(), order="final"))
    #     G.measure_hom("p")
    #     G.eval_Z_probs()
    #     G.eval_Z_probs(cond=True)
    #     for node in G.egraph:
    #         assert "p_phase" in G.egraph.nodes[node]
    #         assert "p_phase_cond" in G.egraph.nodes[node]
