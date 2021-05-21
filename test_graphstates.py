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
import pytest
from numpy.random import default_rng as rng

from graphstates import EGraph, CVGraph, SCZ_mat, SCZ_apply


@pytest.fixture(scope="module", params=[10, 20])
def G_complete(request):
    n = request.param
    G_c = nx.complete_graph(n)
    G_c_adj = nx.to_numpy_array(G_c)
    G_c_adj_sparse = nx.to_scipy_sparse_matrix(G_c)
    return G_c, G_c_adj, G_c_adj_sparse


def noise_model(delta, order):
    return {'noise': 'grn', 'delta': delta, 'sampling_order': order}


class TestEGraph:
    """Tests for EGraphs."""

    def test_init(self, G_complete):
        print(G_complete[0])
        E = EGraph(G_complete[0])
        E_array = nx.to_numpy_array(E)
        assert np.all(G_complete[1] == E_array)

    def test_index(self, G_complete):
        """Test consistency of the indexing function."""
        alph = list(string.ascii_lowercase)
        rand.shuffle(alph)

        E = EGraph()
        E.add_nodes_from(alph)
        E.index_generator()
        E.adj_generator()
        ordered_nodes = [E.to_points[i] for i in range(len(E))]

        G_adj = nx.to_numpy_array(G_complete[0], nodelist=ordered_nodes)
        E_adj = E.adj_mat
        assert np.all(G_adj == E_adj)

    def test_slice_coords(self):
        pass

    def test_draw(self):
        pass


class TestCVHelpers:
    """Test for CVGraph helper functions."""

    def test_SCZ_mat(self, G_complete):
        """Test for the SCZ_mat function."""
        SCZ = SCZ_mat(G_complete[1])
        SCZ_sparse = SCZ_mat(G_complete[2])
        assert type(SCZ) == np.ndarray
        assert type(SCZ_sparse) == sp.coo_matrix

    def test_SCZ_apply(self):
        """Test for SCZ_apply function."""
        pass


class TestCVGraph:
    """Tests for the CVGraph class."""

    def test_init(self, G_complete):
        G = CVGraph(G_complete[0], states=None)
        G_array = nx.to_numpy_array(G.egraph)
        H = CVGraph(EGraph(G_complete[0]), states=None)
        H_array = nx.to_numpy_array(H.egraph)
        assert G._N == len(G_complete[0])
        assert np.all(G_complete[1] == H_array)
        assert np.all(G_complete[1] == G_array)

    def test_states_init(self, G_complete):
        G = CVGraph(G_complete[0])
        n = len(G_complete[0])
        for node in G.egraph:
            assert G.egraph.nodes[node]['state'] == 'GKP'
        assert len(G._states['p']) == 0
        assert len(G.p_inds) == 0
        assert np.array_equal(G._states['GKP'], np.arange(n))
        assert np.array_equal(G.GKP_inds, np.arange(n))

    @pytest.mark.parametrize("p_swap", [0, rng().random(), 1])
    def test_hybridize(self, G_complete, p_swap):
        n = len(G_complete[0])
        G = CVGraph(G_complete[0], p_swap=1)
        for node in G.egraph:
            assert G.egraph.nodes[node]['state'] == 'p'
        assert len(G._states['GKP']) == 0
        assert np.array_equal(G._states['p'], np.arange(n))
        p_list = []
        for i in range(1000):
            G = CVGraph(G_complete[0], p_swap=p_swap)
            p_list += [len(G._states['p']) / n]
        p_prob = sum(p_list) / 1000
        assert np.isclose(p_prob, p_swap, rtol=1e-1)

    def test_state_indices(self, G_complete):
        n = len(G_complete[0])
        num_ps = rng().integers(n)
        p_inds = rng().choice(n, num_ps, replace=False)
        gkp_inds = list(set(np.arange(n)) - set(p_inds))
        G = CVGraph(G_complete[0], states={"p": p_inds})
        assert np.array_equal(G._states.get("p"), p_inds)
        assert np.array_equal(G._states.get("GKP"), gkp_inds)
        assert np.array_equal(G.p_inds, p_inds)
        assert np.array_equal(G.GKP_inds, gkp_inds)

    @pytest.mark.parametrize("order", ["initial", "final", "two-step"])
    def test_apply_noise(self, G_complete, order):
        # Check default noise model
        G = CVGraph(G_complete[0])
        G.apply_noise()
        assert G._delta == 0.01
        assert G._sampling_order == 'initial'
        # Check supplied noise model
        H = CVGraph(G_complete[0])
        delta = rng().random()
        H.apply_noise(noise_model(delta=delta, order=order))
        assert H._delta == delta
        assert H._sampling_order == order

    def test_grn_model(self, G_complete):
        delta = rng().random()
        model_init = noise_model(delta, "initial")
        model_fin = noise_model(delta, "final")
        model_two_step = noise_model(delta, "two-step")

        n = len(G_complete[0])
        G = CVGraph(G_complete[0])
        H = CVGraph(G_complete[0], p_swap=1)
        G.apply_noise(model_init)
        H.apply_noise(model_init)
        init_noise_all_GKP = np.full(2 * n, delta / 2, dtype=np.float32)
        init_noise_all_p = np.array([1 / (2 * delta)] * n + [delta / 2] * n, dtype=np.float32)
        assert np.array_equal(G._init_noise, init_noise_all_GKP)
        assert np.array_equal(H._init_noise, init_noise_all_p)

        G.apply_noise(model_fin)
        H.apply_noise(model_fin)
        noise_cov_all_GKP = SCZ_apply(G._adj, np.diag(init_noise_all_GKP))
        noise_cov_all_p = SCZ_apply(H._adj, np.diag(init_noise_all_p))
        assert np.array_equal(G._noise_cov.toarray(), noise_cov_all_GKP)
        assert np.array_equal(H._noise_cov.toarray(), noise_cov_all_p)
        assert G.noise_cov is not None

        G.apply_noise(model_two_step)
        H.apply_noise(model_two_step)
        assert np.max(H._init_quads[:n]) <= 2 * np.sqrt(np.pi)
        assert np.min(H._init_quads[:n]) >= 0
        assert np.isclose(np.max(G._init_quads[:n]), np.sqrt(np.pi))
        assert np.isclose(np.min(G._init_quads[:n]), 0)

    @pytest.mark.parametrize("order", ["initial", "final"])
    def test_measure_hom(self, G_complete, order):
        n = len(G_complete[0])
        G = CVGraph(G_complete[0])
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

    def test_eval_Z_probs(self, G_complete):
        G = CVGraph(G_complete[0])
        G.apply_noise(noise_model(delta=rng().random(), order="final"))
        G.measure_hom("p")
        G.eval_Z_probs()
        G.eval_Z_probs(cond=True)
        for node in G.egraph:
            assert 'p_phase' in G.egraph.nodes[node]
            assert 'p_phase_cond' in G.egraph.nodes[node]
