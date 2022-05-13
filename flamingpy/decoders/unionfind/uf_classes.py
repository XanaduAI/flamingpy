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
"""A collection of classes used for UnionFind decoder"""

# pylint: disable=too-few-public-methods

import retworkx as rx


class Node:
    """A class for nodes of cluster trees.

    Has the functionality to add children to the node and find the root and
    parity of cluster tree to which the node belongs.

    Args:
        node_id (hashable object): the node id.

    Attributes:
        parent (Node): the parent of the node
        children (set[Node]): the children of the node
        id (hashable object): the node id.
    """

    def __init__(self, node_id=None):
        self.parent = self
        self.children = set()
        self.id = node_id

    def add_child(self, node):
        """Add node as a child of self."""
        self.children.add(node)
        node.parent = self
        root = self.find_root()
        root.size += 1

    def find_root(self):
        """Find the root of the node in the cluster."""
        # From the current node, go upstream in the tree until the root is
        # found.
        node = self
        traversed_nodes = set()
        while not isinstance(node.parent, Root):
            if node.parent == self:
                return None
            temp_node = node
            traversed_nodes.add(node)
            node = node.parent
            node.children.remove(temp_node)
        for traversed_node in traversed_nodes:
            traversed_node.parent = node
        node.children = node.children.union(traversed_nodes)
        return node.parent

    def parity(self):
        """Find the parity of the cluster."""
        if self.find_root():
            return self.find_root().parity

        return None


class Root:
    """A class for roots of cluster trees.

    Has the functionality to add children, find the root, and find children.

    Args:
        node (Node): a node to be associated with the Root object
        parity (int): 0 or 1 -- the parity of the cluster tree.

    Attributes:
        node (Node): a node object associated with the Root object
        size (int): size of the cluster tree
        parity (int): 0 or 1 -- parity of the cluster tree
        boundary (Node): a boundary node of the cluster tree
    """

    def __init__(self, node, parity=1):
        if not isinstance(node, Node):
            raise ValueError("A Root must be initiated with a Node.")
        self.node = node
        self.node.parent = self
        self.size = 1
        if parity == "boundary":
            self.parity = 0
            self.boundary = node
        else:
            self.parity = parity
            self.boundary = None

    def add_child(self, node):
        """Add node as a child to self."""
        self.node.add_child(node)

    @property
    def children(self):
        """Find the children of the Root object."""
        return self.node.children


class Support:
    """Support data structure for the Union-Find decoder.

    Has the functionality to grow cluster edges and obtain the spanning forest
    from the clusters.

    Args:
        stabilizer_graph (StabilizerGraph): a graph of code Stabilizer objects.

    Attributes:
        status (dict): a support table implemented through a dictionary of the
            form {edge_index: status}.
    """

    def __init__(self, stabilizer_graph):
        self.status = {}

        for edge in stabilizer_graph.edges():
            if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                if stabilizer_graph.edge_data(edge[0], edge[1])["weight"] == -1:
                    self.status[frozenset((edge[0], edge[1]))] = "grown"
                else:
                    self.status[frozenset((edge[0], edge[1]))] = "empty"

    def grow(self, edge):
        """Grow edge by half an edge if not fully grown."""
        edge_status = self.status[frozenset(edge)]
        if edge_status == "empty":
            self.status[frozenset(edge)] = "half-grown"
        elif edge_status == "half-grown":
            self.status[frozenset(edge)] = "grown"
            return edge

    def span_forest(self, stabilizer_graph=None):
        """Obtain the spanning forest from the grown clusters."""
        spanning_forest = rx.PyGraph()
        stab_index_to_node = {}
        while self.status != {}:
            edge, status = self.status.popitem()
            edge = list(edge)
            if status == "grown":
                vertices = []
                for i in range(2):
                    if edge[i] in stab_index_to_node:
                        vertices.append(stab_index_to_node[edge[i]])
                    else:
                        vertices.append(spanning_forest.add_node(edge[i]))
                        stab_index_to_node[edge[i]] = vertices[i]

                if stabilizer_graph:
                    common_vertex = stabilizer_graph.edge_data(edge[0], edge[1])["common_vertex"]
                    spanning_forest.add_edge(
                        vertices[0], vertices[1], {"common_vertex": common_vertex}
                    )
                else:
                    spanning_forest.add_edge(vertices[0], vertices[1], None)

        return rx.minimum_spanning_tree(spanning_forest)


# Small helper function for the Boundary class.
def find_child(node, visited):
    """Find the next node to visit."""
    for child in node.children:
        if not child in visited:
            return find_child(child, visited)
    return node


class Boundary:
    """Boundary data structure for the Union-Find decoder.

    Has the functionality to maintain a set of boundary nodes of a cluster.

    Args:
        cluster (Root): Root of the cluster tree associated with the boundary object
        support (Support): Support object of the Union-Find decoder
        stabilizer_graph (StabilizerGraph): a graph of code Stabilizer objects.

    Attributes:
        nodes (set[Node]): the nodes on the boundary,
    """

    def __init__(self, cluster, support, stabilizer_graph):
        # Find the boundary of the initial cluster.
        self.nodes = set()

        # loop over all nodes in the cluster, and add them if they are in the boundary
        visited_nodes = []
        node = cluster.node
        while True:
            node = find_child(node, visited_nodes)
            visited_nodes += [node]

            for edge in stabilizer_graph.out_edges(node.id):
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    if support.status[frozenset((edge[0], edge[1]))] != "grown":
                        self.nodes.add(node)
                        break
            if node == cluster.node:
                break
            node = node.parent

    def prune(self, support, stabilizer_graph):
        """Remove nodes from the boundary that are not in a boundary anymore.

        Args:
            support (Support): the support object of the Union-Find decoder
            stabilizer_graph (StabilizerGraph): a graph of code Stabilizer
                objects.

        Returns:
            None
        """
        nodes_to_keep = set()
        for node in self.nodes:
            for edge in stabilizer_graph.out_edges(node.id):
                if (edge[0] not in {"low", "high"}) and (edge[1] not in {"low", "high"}):
                    if support.status[frozenset((edge[0], edge[1]))] != "grown":
                        nodes_to_keep.add(node)
        self.nodes = nodes_to_keep
