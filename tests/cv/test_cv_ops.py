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
import logging

import networkx as nx
import numpy as np
from numpy.random import default_rng as rng
import pytest
import scipy.sparse as sp

from flamingpy.codes.graphs import EGraph
from flamingpy.cv.ops import invert_permutation, SCZ_mat, SCZ_apply

now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)


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


class TestSCZ:
    """Tests for symplectic CZ matrices."""

    @pytest.mark.parametrize(
        "sparse, expected_out_type", sorted([(True, sp.coo_array), (False, np.ndarray)])
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
        assert isinstance(SCZ_sparse, sp.coo_array)
        # Check that structure of SCZ matrix is correct.
        for mat in (SCZ, SCZ_sparse.toarray()):
            assert np.array_equal(mat[:N, :N], np.identity(N))
            assert np.array_equal(mat[:N, N:], np.zeros((N, N)))
            assert np.array_equal(mat[N:, :N], random_graph[1])
            assert np.array_equal(mat[N:, N:], np.identity(N))

    @pytest.mark.parametrize("one_shot", sorted([True, False]))
    @pytest.mark.parametrize("n", sorted([1, 2]))
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


def test_invert_permutation():
    """Check that permuting and then unpermuting a random array leaves it
    unchanged."""
    N = rng(int_time).integers(1, 100)
    random_array = rng(int_time).integers(0, 100, N)
    random_p = np.arange(N)
    rng(int_time).shuffle(random_p)
    inverted = invert_permutation(random_p)
    presumed_unchanged_array = random_array[random_p][inverted]
    assert np.array_equal(random_array, presumed_unchanged_array)
