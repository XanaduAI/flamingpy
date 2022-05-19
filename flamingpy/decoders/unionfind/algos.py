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
"""Implementation of the Union-Find decoder, adapted from arXiv:1709.06218 and
arXiv:1703.01517 ."""

# pylint: disable=import-outside-toplevel,too-many-locals

import retworkx as rx

from flamingpy.decoders.unionfind.uf_classes import Node, Root, Support, Boundary
from flamingpy.codes.stabilizer import Stabilizer


def union(root1, root2):
    """Perform a weighted union between root1 and root2.

    Args:
        root1 (Root): root of the first cluster in the union
        root2 (Root): root of the second cluster in the union

    Returns:
        NoneType or (Root, Root): If root1 and root2 are same, returns None;
            else, the big and small root node after the union.
    """
    if root1 != root2:
        # The equal case is important here, given the use in initialize_cluster_trees
        if root1.size >= root2.size:
            big_root = root1
            small_root = root2
        else:
            big_root = root2
            small_root = root1

        big_root.size += small_root.size

        if small_root.boundary:
            big_root.boundary = small_root.boundary

        big_root.parity = (big_root.parity + small_root.parity) % 2
        big_root.node.children.add(small_root.node)
        small_root.node.parent = big_root.node

        return big_root.node, small_root.node

    return None


def initialize_cluster_trees(stabilizer_graph):
    """Initialize the cluster trees (Algo 2, step 1 in arXiv:1709.06218).

    Args:
        stabilizer_graph (StabilizerGraph): stabilizer graph that contains the
            syndrome data from measurement outcomes

    Returns:
        dict, list[Root], list[Root]: a dictionary of the nodes, a list of
            roots of the various cluster trees, and a list of roots with
            odd parity.
    """

    # Generate the erasure graph
    erasure_graph = rx.PyGraph()
    stab_to_index = {}

    stabilizer_graph_nodes = stabilizer_graph.nodes()

    for edge in stabilizer_graph.edges():
        if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
            vertices = []
            for i in range(2):
                if edge[i] in stab_to_index:
                    vertices.append(stab_to_index[edge[i]])
                else:
                    # Adding all nodes (not just erasure nodes) is important to
                    # initialize the single node clusters along with erasures
                    vertices.append(erasure_graph.add_node(edge[i]))
                    stab_to_index[edge[i]] = vertices[i]
            if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] == -1:
                # edge_with_indices[2] is a dictionary containing the qubit
                # coordinate corresponding to the edge
                erasure_graph.add_edge(vertices[0], vertices[1], None)

    # Create a dictionary of nodes for stabilizers
    node_dict = {}
    stabilizer_graph_nodes = stabilizer_graph.nodes()
    for stabilizer in stabilizer_graph_nodes:
        if stabilizer != "low" and stabilizer != "high":
            node_dict[stabilizer] = Node(stabilizer)

    # Create clusters corresponding to the connected components of the erasures
    cluster_trees = []

    erasure_graph_nodes = erasure_graph.nodes()
    odd_clusters = []
    for component in rx.connected_components(erasure_graph):
        # Assign a random vertex in the erasure component to be the root
        root_stabilizer = erasure_graph_nodes[component.pop()]
        cluster_root = Root(
            node_dict[root_stabilizer],
            parity=root_stabilizer.parity
            if isinstance(root_stabilizer, Stabilizer)
            else "boundary",
        )  # boundary nodes are represented by tuples
        for vertex in component:
            vertex_stabilizer = erasure_graph_nodes[vertex]
            union(
                cluster_root,
                Root(
                    node_dict[vertex_stabilizer],
                    parity=vertex_stabilizer.parity
                    if isinstance(vertex_stabilizer, Stabilizer)
                    else "boundary",
                ),
            )
        if cluster_root.parity:
            odd_clusters += [cluster_root]
        cluster_trees += [cluster_root]

    return node_dict, cluster_trees, odd_clusters


def union_find(odd_clusters, boundary, stabilizer_graph, support, node_dict):
    """Perform the 'find' and 'union' operations.

    Each odd cluster is grown by a half edge in all the directions; odd clusters
    that have common nodes after the growth are found; and the union operation
    is performed to merge the clusters. This operation is repeated until all
    clusters become even.

    Args:
        odd_clusters (list): list of clusters with odd parity
        boundary (Boundary): dictionary of the Boundary objects of all clusters
        stabilizer_graph (StabilizerGraph): the stabilizer graph
        support (Support): the support table
        node_dict (dict): a dictionary of nodes.

    Returns:
        NoneType
    """

    # Growing the clusters until they all become even
    while odd_clusters:
        # Growing each cluster by half an edge or by a part of an edge based
        # on the type of weight used
        fusion_list = []
        for cluster in odd_clusters:
            for node in boundary[cluster.node].nodes:
                for edge in stabilizer_graph.out_edges(node.id):
                    if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                        grown = support.grow((edge[0], edge[1]))
                        if grown:
                            fusion_list += [(edge[0], edge[1])]
        for edge in fusion_list:
            # If the new edge connects different clusters, fuse them
            root0 = node_dict[edge[0]].find_root()
            root1 = node_dict[edge[1]].find_root()
            if root0 != root1:
                big_root, small_root = union(root0, root1)
                # Updating the boundary of the new cluster
                boundary[big_root].nodes = boundary[big_root].nodes.union(
                    boundary[small_root].nodes
                )
                boundary[big_root].prune(support, stabilizer_graph)

        # Updating the roots of the odd clusters
        temp_odd_clusters = set()
        for cluster in odd_clusters:
            temp_odd_clusters.add(cluster.node.find_root())
        odd_clusters = temp_odd_clusters

        # Removing even clusters
        temp_odd_clusters = set()
        for cluster in odd_clusters:
            if cluster.parity:
                if not cluster.boundary:
                    temp_odd_clusters.add(cluster)
        odd_clusters = temp_odd_clusters


def obtain_spanning_forest(stabilizer_graph, support):
    """Obtain the spanning forest at the end of the union-find step.

    Args:
        stabilizer_graph (StabilizerGraph): the stabilizer graph
        support (Support): the support table containing the grown edges

    Returns:
        rx.PyGraph, parity_dict: a graph containing the edges grown in
            the support table and a dictionary containing the parity of all
            nodes in the spanning_forest.
    """

    spanning_forest = support.span_forest(stabilizer_graph)

    spanning_forest_nodes = spanning_forest.nodes()
    parity_dict = {}
    for tree in rx.connected_components(spanning_forest):
        nb_odd_nodes = 0
        boundary_node = None
        for node_index in tree:
            node = spanning_forest_nodes[node_index]
            if isinstance(node, tuple):
                boundary_node = node_index
                parity_dict[node_index] = 0
            elif isinstance(node, Stabilizer):
                parity = node.parity
                parity_dict[node_index] = parity
                if parity == 1:
                    nb_odd_nodes += 1
        if nb_odd_nodes % 2:
            parity_dict[boundary_node] = 1
    return spanning_forest, parity_dict


def uf_decode(code, ec):
    """Run the union-find decoding algorithm.

    Args:
        code (code): the code object to decode
        ec (str): the error complex ("primal" or "dual")

    Returns:
        rx.PyGraph, parity_dict: a graph containing the edges grown in
            the support table and a dictionary containing the parity of all
            nodes in the spanning_forest.
    """

    # Obtain the stabilizer graph
    stabilizer_graph = getattr(code, ec + "_stab_graph")

    # Initializing the clusters based on erased edges and non-trivial syndrome
    node_dict, cluster_trees, odd_clusters = initialize_cluster_trees(stabilizer_graph)
    support = Support(stabilizer_graph)

    boundary = {}
    for cluster in cluster_trees:
        boundary[cluster.node] = Boundary(cluster, support, stabilizer_graph)
    union_find(odd_clusters, boundary, stabilizer_graph, support, node_dict)

    # Constructing the spanning forest
    spanning_forest, parity_dict = obtain_spanning_forest(stabilizer_graph, support)

    return spanning_forest, parity_dict


def trim_forest(spanning_forest, leaf, parity_dict, recovery):
    """Trim leaves in spanning_forest.

    Args:
        spanning_forest (rx.PyGraph): a graph containing the cluster edges
        leaf (int): index of a leaf node in spanning_forest
        parity_dict (dict): dictionary of parity of the nodes in the spanning
            forest
        recovery (set): set of recovery edges that need to be updated.
    Returns:
        NoneType
    """
    edges = list(spanning_forest.out_edges(leaf))
    if edges:
        edge = edges[0]
    else:
        return
    edge_qubit_index = spanning_forest.get_edge_data(edge[0], edge[1])["common_vertex"]
    spanning_forest.remove_edge(edge[0], edge[1])
    if edge[0] == leaf:
        new_leaf = edge[1]
    else:
        new_leaf = edge[0]
    if parity_dict[leaf] == 1:
        recovery.add(edge_qubit_index)
        parity_dict[leaf] = 0
        parity_dict[new_leaf] ^= 1
    if spanning_forest.degree(new_leaf) == 1:
        trim_forest(spanning_forest, new_leaf, parity_dict, recovery)


def peeling(spanning_forest, parity_dict):
    """Runs the peeling decoding algorithm.

    Args:
        spanning_forest (rx.PyGraph): graph containing the spanning forest,
        parity_dict (dict): dictionary of parity of the nodes in the spanning
            forest

    Returns:
        set[tuples]: the nodes (representing qubits) be fed into the recovery
            (i.e. whose bit values must be flipped).
    """

    recovery_set = set()
    leaves = [
        node for node in range(len(spanning_forest.nodes())) if (spanning_forest.degree(node) == 1)
    ]

    for leaf in leaves:
        trim_forest(spanning_forest, leaf, parity_dict, recovery_set)

    return recovery_set


def uf_decoder(code, ec, **kwargs):
    """Run the full Union-Find and peeling decoder on code.

    Args:
        code (SurfaceCode): the code class to decode and correct
        ec (string): the error complex ("primal" or "dual")

    Returns:
        set[tuples]: the nodes (representing qubits) be fed into the recovery
            (i.e. whose bit values must be flipped).
    """
    spanning_forest, parity_dict = uf_decode(code, ec)
    recovery_set = peeling(spanning_forest, parity_dict)

    if kwargs.get("draw"):
        from flamingpy.utils.viz import draw_decoding

        dec_objects = {"recovery_set": recovery_set}
        draw_decoding(code, ec, dec_objects, kwargs.get("drawing_opts"))

    return recovery_set
