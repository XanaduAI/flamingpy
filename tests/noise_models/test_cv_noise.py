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

from datetime import datetime
import itertools as it
import logging

import networkx as nx
import numpy as np
from numpy.random import default_rng as rng
import pytest

from flamingpy.codes import SurfaceCode
from flamingpy.codes.graphs import EGraph
from flamingpy.noise import CVLayer

# A NetworkX random graph of size N for use in this module.
N = 20

now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)

code_params = it.product(
    [rng(int_time).integers(2, 5), rng(int_time).integers(2, 5, 3)],
    [0.0001],
    [0, 0.5, 1],
    ["open", "periodic"],
    ["primal", "dual"],
)


@pytest.fixture(scope="module", params=[N])
def random_graph(request):
    """A convenience function to initialize random EGraphs for use in this
    module."""
    n = request.param
    G = nx.fast_gnp_random_graph(n, 0.5)
    G_adj = nx.to_numpy_array(G)
    G_adj_sparse = nx.to_scipy_sparse_array(G)
    return EGraph(G), G_adj, G_adj_sparse


@pytest.fixture(scope="module", params=code_params)
def code_and_noise(request):
    """Defines a macronode RHG lattice, the reduced lattice, and the
    delta/p-swap paramaters for use in this module."""
    d, delta, p_swap, boundaries, ec = request.param
    # The reduced lattice.
    RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)
    return delta, p_swap, RHG_code


now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)


class TestCVLayer:
    """Tests for for the CVLayer class."""

    def test_empty_init(self, random_graph):
        """Test the instantiation of CVLayer from codes and EGraphs."""
        G = CVLayer(random_graph[0], delta=0.1)
        G_array = nx.to_numpy_array(G.egraph)
        H = CVLayer(EGraph(random_graph[0]), delta=0.1)
        H_array = nx.to_numpy_array(H.egraph)
        # Check that the _N attribute is populated with the number
        # of nodes in the random graph.
        assert G._N == len(random_graph[0])
        # Check that instantiating a CVLayer with an EGraph or with
        # a regular NetworkX graph has the same effect.
        assert np.array_equal(random_graph[1], H_array)
        assert np.array_equal(random_graph[1], G_array)
        CVRHG = CVLayer(SurfaceCode(2), delta=0.1)
        assert CVRHG._N == len(SurfaceCode(2).graph)

    def test_all_GKP_init(self, random_graph):
        """Test the all-GKP initialization of EGraph."""
        G = CVLayer(random_graph[0], delta=0.1)
        G.populate_states()
        n = len(random_graph[0])
        for node in G.egraph:
            assert G.egraph.nodes[node]["state"] == "GKP"
        assert len(G.states["p"]) == 0
        assert len(G.p_inds) == 0
        assert np.array_equal(G.states["GKP"], np.arange(n))
        assert np.array_equal(G.gkp_inds, np.arange(n))

    @pytest.mark.parametrize("p_swap", sorted([0, 0.99 * rng(int_time).random() + 0.01, 1]))
    def test_hybridize(self, random_graph, p_swap):
        """Test whether CVLayer properly populates p-squeezed states for non-
        zero p-swap."""
        n = len(random_graph[0])
        # Test all-p case
        G = CVLayer(random_graph[0], delta=0.1, p_swap=1)
        G.populate_states()
        for node in G.egraph:
            assert G.egraph.nodes[node]["state"] == "p"
        assert len(G.states["GKP"]) == 0
        assert np.array_equal(G.states["p"], np.arange(n))
        # Test various swap-out probabilities; mean p population
        # should be close to p_swap parameter, within tolerance.
        p_list = []
        for _ in range(1000):
            G = CVLayer(random_graph[0], delta=0.1, p_swap=p_swap)
            G.populate_states()
            p_list += [len(G.states["p"]) / n]
        p_prob = sum(p_list) / 1000
        assert np.isclose(p_prob, p_swap, rtol=1e-1)

    def test_state_indices(self, random_graph):
        """Test that _states, p_inds, and gkp_inds are populated with the
        correct indices."""
        n = len(random_graph[0])
        num_ps = rng(int_time).integers(n)
        p_inds = rng(int_time).choice(n, num_ps, replace=False)
        gkp_inds = list(set(np.arange(n)) - set(p_inds))
        G = CVLayer(random_graph[0], delta=0.1, states={"p": p_inds})
        G.populate_states()
        assert np.array_equal(G.states.get("p"), p_inds)
        assert np.array_equal(G.states.get("GKP"), gkp_inds)
        assert np.array_equal(G.p_inds, p_inds)
        assert np.array_equal(G.gkp_inds, gkp_inds)

    @pytest.mark.parametrize("order", sorted(["initial", "two-step"]))
    def test_apply_noise(self, random_graph, order):
        """Check delta, sampling_order attributes with default noise model."""
        delta = rng(int_time).random()
        G = CVLayer(random_graph[0], delta=delta, sampling_order=order)
        assert G.delta == delta
        assert G._sampling_order == order

    def test_grn_model_init(self, random_graph):
        """Compare expected noise objects (quadratures, covariance matrix) with
        those obtained through the GRN model in CVLayer with the initial
        sampling order."""
        delta = rng(int_time).random()

        n = len(random_graph[0])
        G = CVLayer(random_graph[0], delta=delta)
        H = CVLayer(random_graph[0], delta=delta, p_swap=1)
        G.populate_states(), H.populate_states()
        H._covs_sampler()
        G._covs_sampler()
        init_noise_all_GKP = np.full(2 * n, (delta / 2) ** 0.5, dtype=np.float32)
        init_noise_all_p = np.array(
            [1 / (2 * delta) ** 0.5] * n + [(delta / 2) ** 0.5] * n, dtype=np.float32
        )
        assert np.array_equal(G._init_covs, init_noise_all_GKP)
        assert np.array_equal(H._init_covs, init_noise_all_p)

        # G.apply_noise(sampling_order="final")
        # H.apply_noise(sampling_order="final")
        # noise_cov_all_GKP = SCZ_apply(G._adj, np.diag(init_noise_all_GKP) ** 2)
        # noise_cov_all_p = SCZ_apply(H._adj, np.diag(init_noise_all_p) ** 2)
        # assert np.array_equal(G.Noise_cov.toarray(), noise_cov_all_GKP)
        # assert np.array_equal(H.Noise_cov.toarray(), noise_cov_all_p)

    def test_grn_model_two_step(self, random_graph):
        """Compare expected noise objects (quadratures, covariance matrix) with
        those obtained through the GRN model in CVLayer with the final sampling
        order."""
        delta = rng(int_time).random()

        n = len(random_graph[0])
        G = CVLayer(random_graph[0], delta=delta, sampling_order="two-step")
        H = CVLayer(random_graph[0], delta=delta, p_swap=1, sampling_order="two-step")
        G.populate_states(), H.populate_states()
        G._means_sampler()
        H._means_sampler()
        assert np.max(H._init_means[:n]) <= 2 * np.sqrt(np.pi)
        assert np.min(H._init_means[:n]) >= 0
        assert np.isclose(np.max(G._init_means[:n]), np.sqrt(np.pi))
        assert np.isclose(np.min(G._init_means[:n]), 0)

    @pytest.mark.parametrize("order", sorted(["initial"]))
    def test_measure_hom(self, random_graph, order):
        """Test closeness of average homodyne outcomes value to 0 in the all-
        GKP high-squeezing limit."""
        delta = 0.0001
        n = len(random_graph[0])
        G = CVLayer(random_graph[0], delta=delta)
        G.populate_states()
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
    #     G.apply_noise()
    #     G.measure_hom("p")
    #     G.eval_Z_probs()
    #     G.eval_Z_probs(cond=True)
    #     for node in G.egraph:
    #         assert "p_phase" in G.egraph.nodes[node]
    #         assert "p_phase_cond" in G.egraph.nodes[node]


class CVMacroLayer:
    """Tests for the CVMacroLayer class."""

    def test_reduction(self, code_and_noise):
        """Test the reduction from the macronode to the canonical lattice."""
        delta, p_swap, RHG_code = code_and_noise
        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        CV_macro = CVMacroLayer(RHG_code, delta=delta, p_swap=p_swap)
        CV_macro.apply_noise()
        # Check proper reduction to effective node type.
        RHG_macro = CV_macro.egraph
        RHG_reduced = CV_macro.reduced_graph
        for central_node in RHG_macro.macro_to_micro:
            micronodes = RHG_macro.macro_to_micro[central_node]
            effective_type = RHG_reduced.nodes[central_node]["state"]
            p_count = 0
            for micronode in micronodes:
                if RHG_macro.nodes[micronode]["state"] == "p":
                    p_count += 1
            expected_type = "p" if p_count == 4 else "GKP"
            assert effective_type == expected_type
