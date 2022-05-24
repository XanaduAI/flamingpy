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
""""Unit tests for the Union Find decoder."""

# pylint: disable=too-many-function-args,no-self-use,unsubscriptable-object,consider-using-dict-items,too-few-public-methods,too-many-locals,unused-argument

import itertools as it
import random

import networkx as nx
import numpy as np
import pytest
import retworkx as rx

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders.decoder import assign_weights, CV_decoder

from flamingpy.decoders.unionfind.uf_classes import Node, Root, Support, Boundary
from flamingpy.codes.graphs.stabilizer_graph import RxStabilizerGraph
from flamingpy.codes.stabilizer import Stabilizer
from flamingpy.decoders.unionfind.algos import (
    initialize_cluster_trees,
    obtain_spanning_forest,
    union_find,
    peeling,
    union,
)


code_params = it.product(
    [2, 3, 4], ["primal", "dual"], ["open", "periodic"], [1, 0.1, 0.01], [0, 0.5, 1]
)  # distance, ec, boundaries, delta, p_swap


@pytest.fixture(scope="module", params=code_params)
def enc_state(request):
    """An RHGCode object and an encoded CVLayer for use in this module."""
    distance, ec, boundaries, delta, p_swap = request.param
    DVRHG = SurfaceCode(distance, ec, boundaries, alternating_polarity, backend="retworkx")
    RHG_lattice = DVRHG.graph
    # CV (inner) code/state
    CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
    # Noise model
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
    # Apply noise
    CVRHG.apply_noise(cv_noise)
    # Measure syndrome
    CVRHG.measure_hom("p", DVRHG.all_syndrome_inds)
    return DVRHG, CVRHG


code_params = it.product(
    [2, 3, 4], ["primal", "dual"], ["open", "periodic"], [1, 0.1, 0.01]
)  # distance, ec, boundaries, delta, p_swap


@pytest.fixture(scope="module", params=code_params)
def enc_state_swap_list(request):
    """
    In this function, lattice is generated based on a list of psqueezed states,
    which are generated randomly.
    """
    distance, ec, boundaries, delta = request.param

    if boundaries == "periodic":
        psqueezed = random.sample(range(6 * distance**3), random.randint(1, 6 * distance**3))
    else:
        # Here:
        # d(d-1)*(2d-1) [For slice with surface code] +
        # (2d-1)(2d^2-2d+1) [For the alternate lattice] = 6d^3-9d^2+5d-1
        num_qubits = (6 * distance**3) - (9 * distance**2) + (5 * distance) - 1
        psqueezed = random.sample(range(num_qubits), random.randint(1, num_qubits))

    states = {"p": np.array(psqueezed)}

    DVRHG = SurfaceCode(distance, ec, boundaries, alternating_polarity, backend="retworkx")
    RHG_lattice = DVRHG.graph
    # CV (inner) code/state
    CVRHG = CVLayer(RHG_lattice, states=states)
    # Noise model
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
    # Apply noise
    CVRHG.apply_noise(cv_noise)
    # Measure syndrome
    CVRHG.measure_hom("p", DVRHG.all_syndrome_inds)

    num_psqueezed_neighbor = {}
    for p in psqueezed:
        for key in dict(DVRHG.graph[CVRHG.to_points[p]]).keys():
            if key in set(DVRHG.all_syndrome_coords):
                num_psqueezed_neighbor[key] = (
                    num_psqueezed_neighbor[key] + 1 if key in num_psqueezed_neighbor else 1
                )
    pswaps = [
        key for key in num_psqueezed_neighbor if num_psqueezed_neighbor[key] > 1
    ]  # Modified pswaps to obtain all the qubits nodes that have at least 2 p-squeezed
    return DVRHG, CVRHG, pswaps


class TestUnionFindStructures:
    """Test the structure of UnionFind members."""

    def test_cluster_tree(self):
        """Test basic cluster tree structures."""

        # Testing by hand core functionalities
        r1, r2 = Node(), Node()
        cluster1, cluster2 = Root(r1), Root(r2)
        a, b, c, d, e, f = Node(), Node(), Node(), Node(), Node(), Node()
        cluster1.add_child(a)
        cluster1.add_child(b)
        cluster1.add_child(e)
        b.add_child(c)
        cluster2.add_child(d)

        assert cluster1.size == 5
        assert cluster1.parity == 1
        assert b.children == set([c])
        assert set(cluster1.node.children) == set([a, b, e])
        assert c.parity() == 1
        assert a.find_root() == c.find_root()
        assert e.find_root() != d.find_root()
        assert f.find_root() is None
        assert f.parity() is None

    def test_support(self, enc_state_swap_list):
        """Test the support table."""

        CV_decoder(enc_state_swap_list[0])
        assign_weights(enc_state_swap_list[0], "UF")

        for ec in enc_state_swap_list[0].ec:
            if ec == "primal":
                stabilizer_graph = enc_state_swap_list[0].primal_stab_graph
            elif ec == "dual":
                stabilizer_graph = enc_state_swap_list[0].dual_stab_graph
            support = Support(stabilizer_graph)
            half_edges = []
            for edge in stabilizer_graph.edges():
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] != -1:
                        if random.random() > 0.5:
                            half_edges += [frozenset((edge[0], edge[1]))]
                            support.grow(frozenset((edge[0], edge[1])))
            for edge in stabilizer_graph.edges():
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] == -1:
                        assert support.status[frozenset((edge[0], edge[1]))] == "grown"
                    else:
                        if frozenset((edge[0], edge[1])) in half_edges:
                            assert support.status[frozenset((edge[0], edge[1]))] == "half-grown"
                        else:
                            assert (
                                support.status[frozenset((edge[0], edge[1]))] == "empty"
                            )  # empty edges replaced by value 0 as we consider weighted edges
                    support.grow(frozenset((edge[0], edge[1])))
                    support.grow(frozenset((edge[0], edge[1])))
            for edge in stabilizer_graph.edges():
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    assert support.status[frozenset((edge[0], edge[1]))] == "grown"

    def test_support_value(self, enc_state_swap_list):
        """Test Support function assignments"""
        CV_decoder(enc_state_swap_list[0])
        assign_weights(enc_state_swap_list[0], "UF")

        for ec in enc_state_swap_list[0].ec:
            if ec == "primal":
                stabilizer_graph = enc_state_swap_list[0].primal_stab_graph
            elif ec == "dual":
                stabilizer_graph = enc_state_swap_list[0].dual_stab_graph
            support = Support(stabilizer_graph)
            with pytest.raises(
                KeyError
            ):  # Changed Value error to Key error as weights is a dictionary
                support.status[(0.32, 0.1, 0)]

    def test_boundary(self, enc_state):
        """Test the initialization of the Boundary class."""
        CV_decoder(enc_state[0])
        assign_weights(enc_state[0], "UF")
        for ec in enc_state[0].ec:
            if ec == "primal":
                stabilizer_graph = enc_state[0].primal_stab_graph
            elif ec == "dual":
                stabilizer_graph = enc_state[0].dual_stab_graph
            # index of (0,0,0)

            if ec == "primal":
                root_tuple = (0, 0, 0)
            else:
                root_tuple = (1, 1, 1)
            root_object = [
                node
                for node in stabilizer_graph.nodes()
                if (
                    isinstance(node, Stabilizer)
                    and (tuple((int(node.midpoint()[i] - 1) for i in range(3))) == root_tuple)
                )
            ][0]
            for edge in stabilizer_graph.out_edges(root_object):
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    stabilizer_graph.edge_data(edge[0], edge[1])["weight"] = -1
            support_list = Support(stabilizer_graph)
            root = Node(root_object)
            cluster = Root(root)
            bound_set = set()
            for edge in stabilizer_graph.out_edges(root_object):
                if isinstance(edge[0], Stabilizer) and isinstance(edge[1], Stabilizer):
                    for out_edge in stabilizer_graph.out_edges(edge[1]):
                        if support_list.status[frozenset((out_edge[0], out_edge[1]))] != "grown":
                            node = Node(edge[1])
                            cluster.add_child(node)
                            bound_set.add(node)
                            break
        boundary = Boundary(cluster, support_list, stabilizer_graph)
        assert set(boundary.nodes) == bound_set


class TestUnionFindFunctions:
    """Test function members of UnionFind decoder"""

    def test_union(self):
        """Test the union of trees."""

        root1, root2, root3 = Node(), Node(), Node()
        cluster1, cluster2, cluster3 = Root(root1), Root(root2), Root(root3)
        tree_list = [cluster1, cluster2, cluster3]
        a, b, c, d = Node(), Node(), Node(), Node()
        cluster1.add_child(a)
        cluster1.add_child(b)
        cluster2.add_child(c)
        cluster3.add_child(d)
        root_to_remove = root3.find_root()
        union(root3.find_root(), a.find_root())
        tree_list.remove(root_to_remove)

        assert root3.parent == root1
        assert cluster1 == root1.parent
        assert d.parent == root3

        root_to_remove = root2.find_root()
        union(root3.find_root(), root2.find_root())
        tree_list.remove(root_to_remove)

        assert tree_list[0].size == 7
        assert tree_list[0].parity == 1
        assert tree_list[0].node == root1
        assert union(a.find_root(), d.find_root()) is None

        with pytest.raises(ValueError):
            Root(5)

    def test_compression(self):
        """Test path compression."""

        root1, root2 = Node(), Node()
        cluster1, cluster2 = Root(root1), Root(root2)
        a, b, c = Node(), Node(), Node()
        cluster1.add_child(a)
        a.add_child(b)
        cluster2.add_child(c)
        union(b.find_root(), c.find_root())
        del cluster2

        assert cluster1.size == 5
        assert root1.children == set([a, b, root2])
        assert a.children == set()
        assert b.children == set()
        assert root2.children == set([c])
        assert c.children == set()
        assert c.find_root() == a.find_root()
        assert cluster1.children == set([a, b, c, root2])
        assert c.find_root() == cluster1

    def test_initialize_cluster_trees(self, enc_state):
        """Testing cluster initialization."""
        CV_decoder(enc_state[0])
        assign_weights(enc_state[0], "UF")
        for ec in enc_state[0].ec:
            if ec == "primal":
                stabilizer_graph = enc_state[0].primal_stab_graph
            elif ec == "dual":
                stabilizer_graph = enc_state[0].dual_stab_graph
            # Initializing the clusters based on erased edges and non-trivial syndrome
            _, cluster_trees, _ = initialize_cluster_trees(stabilizer_graph)

            # Generate the erasure graph
            erasure_graph = rx.PyGraph()
            stab_to_index = {}

            for edge in stabilizer_graph.edges():
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    vertices = []
                    for i in range(2):
                        if edge[i] in stab_to_index:
                            vertices.append(stab_to_index[edge[i]])
                        else:
                            vertices.append(erasure_graph.add_node(edge[i]))
                            stab_to_index[edge[i]] = vertices[i]
                    if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] == -1:
                        # edge_with_indices[2] is a dictionary containing the qubit
                        # coordinate corresponding to the edge
                        erasure_graph.add_edge(vertices[0], vertices[1], None)

            nb_clusters = len(rx.connected_components(erasure_graph))

            assert len(cluster_trees) == nb_clusters
            assert len(set(cluster_trees)) == len(cluster_trees)

            for cluster in cluster_trees:
                if cluster.node.children != set():
                    cluster_node_ids = set(
                        [cluster.node.id] + [child.id for child in cluster.node.children]
                    )
                    for node_id in cluster_node_ids:
                        # Check if at least one incident edge from another node in the
                        # cluster is erased
                        check_erasure = False
                        for edge in stabilizer_graph.out_edges(node_id):
                            if edge[1] in cluster_node_ids:
                                if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] == -1:
                                    check_erasure = True
                                    break
                        assert check_erasure is True


class PeelingDecoder:
    """A class to host methods for testing peeling decoder assignments"""

    def test_peeling(self):
        """Test the peeling decoder for the surface code with periodic boundary conditions."""
        # Test 1: Check if the peeling decoder works fine for one even cluster
        spanning_forest = rx.PyGraph()

        # Add edges from the graph used in the slides of UF decoder
        # demo in QPU Decoder Session at FTQC HW Wkshp
        spanning_forest.add_nodes_from(range(19))
        spanning_forest.add_edges_from(
            [
                (1, 3, {"common_vertex": 1}),
                (2, 3, {"common_vertex": 2}),
                (3, 5, {"common_vertex": 3}),
                (3, 4, {"common_vertex": 4}),
                (4, 6, {"common_vertex": 5}),
                (6, 7, {"common_vertex": 6}),
                (7, 9, {"common_vertex": 7}),
                (9, 10, {"common_vertex": 8}),
                (8, 10, {"common_vertex": 9}),
                (10, 11, {"common_vertex": 10}),
                (11, 12, {"common_vertex": 11}),
                (11, 15, {"common_vertex": 12}),
                (15, 16, {"common_vertex": 13}),
                (14, 15, {"common_vertex": 14}),
                (13, 14, {"common_vertex": 15}),
                (14, 17, {"common_vertex": 16}),
                (15, 18, {"common_vertex": 17}),
            ]
        )

        nontrivial_syndrome_set = set([3, 6, 7, 11])
        parity_dict = {}
        for i in range(19):
            if i in nontrivial_syndrome_set:
                parity_dict[i] = 1
            else:
                parity_dict[i] = 0

        recovery = peeling(spanning_forest, parity_dict)

        assert recovery == set([4, 5, 7, 8, 10])

        # Test 2: Check if the peeling decoder works fine for two even clusters
        # with one cluster with no non-trivial syndrome

        spanning_forest = rx.PyGraph()
        spanning_forest.add_nodes_from(range(23))
        node_to_index_dict = {
            0: (0, 0, 0),
            1: (2, 5, 0),
            2: (3, 4, 0),
            3: (3, 5, 0),
            4: (3, 6, 0),
            5: (4, 5, 0),
            6: (4, 6, 0),
            7: (5, 5, 0),
            8: (5, 6, 0),
            9: (5, 7, 0),
            10: (6, 6, 0),
            11: (7, 8, 0),
            12: (7, 9, 0),
            13: (7, 10, 0),
            14: (7, 11, 0),
            15: (8, 4, 0),
            16: (8, 5, 0),
            17: (8, 8, 0),
            18: (8, 9, 0),
            19: (8, 10, 0),
            20: (8, 11, 0),
            21: (9, 4, 0),
            22: (9, 5, 0),
        }
        # Add edges from the graph used in the slides of UF decoder demo
        # in QPU Decoder Session at FTQC HW Wkshp
        spanning_forest.add_edges_from(
            [
                (1, 3, {"common_vertex": 1}),
                (2, 3, {"common_vertex": 2}),
                (3, 4, {"common_vertex": 3}),
                (3, 5, {"common_vertex": 4}),
                (4, 6, {"common_vertex": 5}),
                (6, 8, {"common_vertex": 6}),
                (7, 8, {"common_vertex": 7}),
                (8, 9, {"common_vertex": 8}),
                (8, 10, {"common_vertex": 9}),
                (11, 12, {"common_vertex": 10}),
                (11, 17, {"common_vertex": 11}),
                (12, 13, {"common_vertex": 12}),
                (12, 18, {"common_vertex": 13}),
                (13, 14, {"common_vertex": 14}),
                (13, 19, {"common_vertex": 15}),
                (14, 20, {"common_vertex": 16}),
                (15, 16, {"common_vertex": 17}),
                (15, 21, {"common_vertex": 18}),
                (16, 22, {"common_vertex": 19}),
            ]
        )

        nontrivial_syndrome_set = set([3, 8, 12, 13])
        parity_dict = {}
        for i in range(23):
            if i in nontrivial_syndrome_set:
                parity_dict[i] = 1
            else:
                parity_dict[i] = 0

        recovery = peeling(spanning_forest, parity_dict)

        assert recovery == set([3, 5, 6, 12])


class UnionFindDecoder:
    """A class to host methods testing UnionFind decoder assignments."""

    def test_uf_erasure(self, enc_state_swap_list):
        """Test if weights for the erased edges are set to -1."""

        CV_decoder(enc_state_swap_list[0])
        assign_weights(enc_state_swap_list[0], "UF")
        erased_qubits = set()
        for node in enc_state_swap_list[0].graph.nodes():
            if (node in enc_state_swap_list[0].all_syndrome_coords) and (
                enc_state_swap_list[0].graph.nodes[node]["weight"] == -1
            ):
                erased_qubits.add(node)
        assert set(enc_state_swap_list[2]) == set(erased_qubits)

    def test_assign_weight_UF(self):
        """
        Testing the different types of weights and number of weight levels in
        different weighted UF decoders in UFdecoder.py

        The graph considered is
        (0,0,0) A'p' -
        -             -
        -              -
        -                 B (3,3,1) - - - D'p' (6,3,2)
        -              -
        (0,5,4) C   -
        """

        class MinimalCode:
            """
            A class similar to SurfaceCode class but only with
            attributes that are required for the testing.
            """

            def __init__(self, G, synd_coords):
                self.graph = G
                self.all_syndrome_coords = synd_coords

        # Creating a graph for the MinimalCode class object
        G = nx.Graph()

        G.add_node((0, 0, 0))
        G.add_node((3, 3, 1))
        G.add_node((0, 5, 4))
        G.add_node((6, 3, 2))

        G.nodes[(0, 0, 0)]["state"] = "p"
        G.nodes[(3, 3, 1)]["state"] = "gkp"
        G.nodes[(0, 5, 4)]["state"] = "gkp"
        G.nodes[(6, 3, 2)]["state"] = "p"

        G.add_edge((0, 0, 0), (0, 5, 4))
        G.add_edge((0, 0, 0), (3, 3, 1))
        G.add_edge((3, 3, 1), (6, 3, 2))
        G.add_edge((3, 3, 1), (0, 5, 4))

        # error probability p = value of "p_phase_cond" key
        G.nodes[(0, 0, 0)]["p_phase_cond"] = 0.001
        G.nodes[(3, 3, 1)]["p_phase_cond"] = 0.5
        G.nodes[(0, 5, 4)]["p_phase_cond"] = 0.1
        G.nodes[(6, 3, 2)]["p_phase_cond"] = 0.23

        synd_coords = [(0, 0, 0), (3, 3, 1), (0, 5, 4), (6, 3, 2)]

        code = MinimalCode(G, synd_coords)

        # Checking if the default UF decoder is the non-weighted UF decoder
        assign_weights(code, "UF")

        G = code.graph

        assert G.nodes[(0, 0, 0)]["weight"] == 2
        assert G.nodes[(3, 3, 1)]["weight"] == -1
        assert G.nodes[(0, 5, 4)]["weight"] == 2
        assert G.nodes[(6, 3, 2)]["weight"] == 2

    def test_UF(self):
        """
        Testing the UFdecoder
        """
        syndrome_dict = {1: 1, 2: 0, 3: 1, 4: 0, 5: 0, 6: 1, 7: 1, 8: 0, 9: 0}
        syndrome_nodes = []
        for i in range(len(syndrome_dict)):
            G = nx.Graph()
            G.add_node(i + 1, bit_val=syndrome_dict[i + 1])
            syndrome_nodes.append(Stabilizer(G))

        rx_stabilizer_graph = RxStabilizerGraph("primal")

        # The RxStabilizerGraph has two nodes 'low' and 'high' initially.
        # Thus, the node index for the rest starts from 2.
        rx_stabilizer_graph.graph.add_nodes_from(syndrome_nodes)
        rx_stabilizer_graph.graph.add_edges_from(
            [
                (2, 3, {"common_vertex": 1}),
                (3, 4, {"common_vertex": 2}),
                (2, 5, {"common_vertex": 3}),
                (3, 5, {"common_vertex": 4}),
                (4, 7, {"common_vertex": 5}),
                (5, 6, {"common_vertex": 6}),
                (6, 7, {"common_vertex": 7}),
                (5, 8, {"common_vertex": 8}),
                (6, 9, {"common_vertex": 9}),
                (7, 10, {"common_vertex": 10}),
                (8, 9, {"common_vertex": 11}),
                (9, 10, {"common_vertex": 12}),
            ]
        )

        def fn(code):
            return 2

        weight_fn = fn

        node_dict, cluster_trees, odd_clusters = initialize_cluster_trees(
            rx_stabilizer_graph, weight_fn
        )
        support = Support(rx_stabilizer_graph, weight_fn)
        boundary = {}
        for cluster in cluster_trees:
            boundary[cluster.node] = Boundary(cluster, support, rx_stabilizer_graph)
        union_find(odd_clusters, boundary, rx_stabilizer_graph, support, node_dict)

        # Constructing the spanning forest
        spanning_forest, parity_dict = obtain_spanning_forest(rx_stabilizer_graph, support)

        assert peeling(spanning_forest, parity_dict) == set([3, 5, 8])
