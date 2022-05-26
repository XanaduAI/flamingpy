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
"""The stabilizer graph interface and some implementations."""

# pylint: disable=import-outside-toplevel

from abc import ABC
import itertools as it

import networkx as nx
import networkx.algorithms.shortest_paths as sp

import retworkx as rx


class StabilizerGraph(ABC):
    """An abstraction for a stabilizer graph.

    Stablizer graphs are a intermediate representations of qubit codes used by
    decoders such as minimum weight perfect matching (MWPM).

    The nodes of a stabilizer graph are every stabilizer element in a code and
    every boundary point (coming from the type of boundary determined
    by the error complex -- primal or dual). Two stabilizers (or a stabilizer
    and a boundary point) sharing a vertex are connected by an edge whose
    weight is equal to the weight assigned to that vertex in the code. The
    weight is computed on the fly while computing the shortest path.

    The main function of a stabilizer graph is to be fed into the decoder.
    In the case of MWPM, it is used to compute the shortest paths between
    many nodes in order to construct the matching graph.

    Parameters:
        ec (str, optional): the error complex ("primal" or "dual"). Determines whether
            the graph is generated from primal or dual stabilizers in the code.
        code (SurfaceCode, optional): the code from which to initialize the graph.

    Note:
        Both ec and code must be provided to initialize the graph with the proper edges.
        If one of them is not provided, the graph is left empty.

    Attributes:
        stabilizers (List[Stabilizer]): All the stabilizer nodes of the graph.

        {ec}_low_bound_points (List[Tuple[Int, Int, Int]]): All the points on
            the code boundary with smaller coordinates. 'ec' can be 'primal'
            or 'dual'.

        {ec}_high_bound_points (List[Tuple[Int, Int, Int]]): All the points on
            the code boundary with larger coordinates. 'ec' can be 'primal'
            or 'dual'.

    Note:
        One has first to assign weights and bit values to the nodes of the
        lattice. This can be achieved, e.g,. by applying CV noise, conducting
        homodyne measurements, computing the phase error probabilities, and
        translating the outcomes.
    """

    # pylint: disable=too-many-public-methods
    def __init__(self, ec=None, code=None):
        self.add_node("low")
        self.add_node("high")
        self.stabilizers = []
        self.low_bound_points = []
        self.high_bound_points = []
        if code is not None and ec is not None:
            self.add_stabilizers(getattr(code, ec + "_stabilizers"))
            bound_points = getattr(code, ec + "_bound_points")
            mid = int(len(bound_points) / 2)
            # All points connected to the boundary slice (as determined by ec)
            # with smaller coordinates.
            self.add_low_bound_points(bound_points[:mid])
            # All points connected to the boundary slice (as determined by ec)
            # with larger coordinates.
            self.add_high_bound_points(bound_points[mid:])
            self.connect_nodes()

    def add_node(self, node):
        """Insert a node into the stabilizer graph.

        This should not distinguish between stabilizer and boundary points.

        Returns:
            The updated stabilizer graph.
        """
        raise NotImplementedError

    def nodes(self):
        """Return an iterable of all nodes in the graph."""
        raise NotImplementedError

    def add_edge(
        self,
        node1,
        node2,
        common_vertex=None,
    ):
        """Insert a node into the graph and return the updated stabilizer
        graph.

        This should not distinguish between stabilizer and boundary points.

        Parameters:
            node1 (Stabilizer or Tuple[Int, Int, Int]): The first node of the edge.
            node2 (Stabilizer or Tuple[Int, Int, Int]): The second node of the edge.
            common_vertex (optional): The vertex shared by the two nodes in the
                corresponding code.

        Returns:
            The updated stabilizer graph.
        """
        raise NotImplementedError

    def edge_data(self, node1, node2):
        """Return a view of the edge data as a dict.

        Parameters:
            node1 (Stabilizer or Tuple[Int, Int, Int]): The first node of the edge.
            node2 (Stabilizer or Tuple[Int, Int, Int]): The second node of the edge.

        Raises:
            KeyError if there is no edge between the given nodes.
        """
        raise NotImplementedError

    def edges(self):
        """Return an iterable of node pairs corresponding to the edges of the
        graph."""
        raise NotImplementedError

    def neighbors(self, node):
        """Return the neighbours of node."""
        raise NotImplementedError

    def out_edges(self, node):
        """Return the edges incident to node."""
        raise NotImplementedError

    def shortest_paths_without_high_low(self, source):
        """Compute the shortest path from source to every other node in the
        graph, except the 'high' and 'low' connector.

        This assumes that the edge weights are asssigned.

        Note: a path can't use the 'high' and 'low' node.


        Arguments:
            source: The source node for each path.
        Returns:
            (dict, dict): The first dictionary maps a target node to the weight
                of the corresponding path. The second one maps a target node to
                the list of nodes along the corresponding path.
        """

        raise NotImplementedError

    def shortest_paths_from_high(self):
        """Compute the shortest path from the 'high' node to every other node
        in the graph.

        This assumes that the edge weights are asssigned.

        Returns:
            (dict, dict): The first dictionary maps a target node to the weight
                of the corresponding path. The second one maps a target node to
                the list of nodes along the corresponding path.
        """
        raise NotImplementedError

    def shortest_paths_from_low(self):
        """Compute the shortest path from the 'low' node to every other node in
        the graph.

        This assumes that the edge weights are asssigned.

        Returns:
            (dict, dict): The first dictionary maps a target node to the weight
                of the corresponding path. The second one maps a target node to
                the list of nodes along the corresponding path.
        """
        raise NotImplementedError

    def add_stabilizer(self, stabilizer):
        """Add a stabilizer node to the stabilizer graph.

        Parameters:
            stabilizer (Stabilizer):
                The stabilizer to add.

        Returns:
            The updated stabilizer graph.
        """
        self.stabilizers.append(stabilizer)
        self.add_node(stabilizer)
        return self

    def add_stabilizers(self, stabilizers):
        """Add many stabilizer nodes to the stabilizer graph.

        Parameters:
            stabilizers (iterable of Stabilizer):
                The stabilizers to add.

        Returns:
            The updated stabilizer graph.
        """
        for stabilizer in stabilizers:
            self = self.add_stabilizer(stabilizer)
        return self

    def add_low_bound_point(self, point):
        """Add a boundary point with an edge to the 'low' point to the graph.

        Parameters:
            point (Tuple[int, int, int]):
                The boundary point to add.

        Returns:
            The updated stabilizer graph.
        """
        self.low_bound_points.append(point)
        self.add_node(point).add_edge("low", point)
        return self

    def add_low_bound_points(self, points):
        """Add many boundary points with edges to the 'low' point to the graph.

        Parameters:
            points (iterable of Tuple[int, int, int]):
                The boundary points to add.

        Returns:
            The updated stabilizer graph.
        """
        for point in points:
            self = self.add_low_bound_point(point)
        return self

    def add_high_bound_point(self, point):
        """Add a boundary point with an edge to the 'high' point to the graph.

        Parameters:
            point (Tuple[int, int, int]):
                The boundary point to add.

        Returns:
            The updated stabilizer graph.
        """
        self.high_bound_points.append(point)
        self = self.add_node(point).add_edge("high", point)
        return self

    def add_high_bound_points(self, points):
        """Add many boundary points with edges to the 'high' point to the
        graph.

        Parameters:
            points (iterable of Tuple[int, int, int]):
                The boundary points to add.

        Returns:
            The updated stabilizer graph.
        """
        for point in points:
            self = self.add_high_bound_point(point)
        return self

    def bound_points(self):
        """Return an iterable for all boundary points of the graph."""
        return it.chain(self.low_bound_points, self.high_bound_points)

    def has_bound_points(self):
        """Check if the graph has any boundary points."""
        return self.low_bound_points or self.high_bound_points

    def connect_nodes(self):
        """Add an edge between each pair of nodes sharing a common vertex.

        Returns:
            The updated stabilizer graph.
        """
        for (stab1, stab2) in it.combinations(self.stabilizers, 2):
            common_vertex = set(stab1.coords()) & set(stab2.coords())
            if common_vertex:
                self.add_edge(stab1, stab2, common_vertex=common_vertex.pop())
        for (stab, point) in it.product(self.stabilizers, self.bound_points()):
            if point in stab.coords():
                self.add_edge(stab, point, common_vertex=point)
        return self

    def odd_parity_stabilizers(self):
        """Return an iterable of all stabilizer nodes with an odd parity."""
        return filter(lambda stab: stab.parity % 2 == 1, self.stabilizers)

    def real_nodes(self):
        """Return an iterable of all nodes excluding the 'low' and 'high'
        points."""
        return it.chain(self.stabilizers, self.bound_points())

    def real_edges(self):
        """Returns an iterable of all edges excluding the ones connected to the
        'low' or 'high' points."""
        return filter(
            lambda edge: edge[0] not in ("low", "high") and edge[1] not in ("low", "high"),
            self.edges(),
        )

    def assign_weights(self, code):
        """Assign the weights to the graph based on the weight of the common
        vertex of each stabilizer pair of the code."""
        for edge in self.edges():
            data = self.edge_data(*edge)
            if data["common_vertex"] is not None:
                data["weight"] = code.graph.nodes[data["common_vertex"]].get("weight")
            elif "high" in edge or "low" in edge:
                data["weight"] = 0

    def to_nx(self):
        """Convert the same graph into a NxStabilizerGraph.

        This involves converting the graph representation to a networkx
        graph representation.
        """
        if isinstance(self, NxStabilizerGraph):
            return self
        nx_graph = NxStabilizerGraph()
        for edge in self.edges():
            nx_graph.add_edge(*edge, self.edge_data(*edge)["common_vertex"])
            if "weight" in self.edge_data(*edge):
                nx_graph.edge_data(*edge)["weight"] = self.edge_data(*edge)["weight"]
        return nx_graph

    def draw(self, **kwargs):
        """Draw the stabilizer graph with matplotlib.

        See flamingpy.utils.viz.draw_dec_graph for more details.
        """
        from flamingpy.utils.viz import draw_dec_graph

        return draw_dec_graph(self.to_nx(), **kwargs)


class NxStabilizerGraph(StabilizerGraph):
    """An implementation of StabilizerGraph backed by a NetworkX graph.

    See StabilizerGraph for more details.

    Attributes:
        graph (networkx.Graph): The actual graph backend.
    """

    def __init__(self, ec=None, code=None):
        self.graph = nx.Graph()
        StabilizerGraph.__init__(self, ec, code)

    def add_node(self, node):
        self.graph.add_node(node)
        return self

    def nodes(self):
        return self.graph.nodes

    def add_edge(self, node1, node2, common_vertex=None):
        self.graph.add_edge(node1, node2, common_vertex=common_vertex)
        return self

    def edge_data(self, node1, node2):
        return self.graph.edges[(node1, node2)]

    def edges(self):
        return self.graph.edges()

    def neighbors(self, node):
        return self.graph.neighbors(node)

    def out_edges(self, node):
        return self.graph.edges(node)

    def shortest_paths_without_high_low(self, source):
        subgraph = self.graph.subgraph(
            self.stabilizers + self.low_bound_points + self.high_bound_points
        )
        return nx_shortest_paths_from(subgraph, source)

    def shortest_paths_from_high(self):
        return nx_shortest_paths_from(self.graph, "high")

    def shortest_paths_from_low(self):
        return nx_shortest_paths_from(self.graph, "low")


def nx_shortest_paths_from(graph, source):
    """The NetworkX shortest path implementation."""
    (weights, paths) = sp.single_source_dijkstra(graph, source)
    del weights[source]
    del paths[source]
    return (weights, paths)


class RxStabilizerGraph(StabilizerGraph):
    """An implementation of StabilizerGraph backed by a retworkx graph.

    See StabilizerGraph for more details.

    Attributes:
        graph (retworkx.PyGraph): The actual graph backend. This graph stores
            integer indices to represent nodes.
        node_to_index (dict): The map from nodes to the corresponding indices
            in the graph backend.
        index_to_node (dict): The map from indices in the graph backend to the
            corresponding nodes.
    """

    def __init__(self, ec=None, code=None):
        self.graph = rx.PyGraph()
        self.node_to_index = {}
        self.index_to_node = {}
        StabilizerGraph.__init__(self, ec, code)

    def add_node(self, node):
        index = self.graph.add_node(node)
        self.node_to_index[node] = index
        self.index_to_node[index] = node
        return self

    def nodes(self):
        return self.graph.nodes()

    def add_edge(self, node1, node2, common_vertex=None):
        index1 = self.node_to_index[node1]
        index2 = self.node_to_index[node2]
        self.graph.add_edge(index1, index2, {"common_vertex": common_vertex})
        return self

    def edge_data(self, node1, node2):
        index1 = self.node_to_index[node1]
        index2 = self.node_to_index[node2]
        return self.graph.get_edge_data(index1, index2)

    def edges(self):
        return (
            (self.graph.get_node_data(edge[0]), self.graph.get_node_data(edge[1]))
            for edge in self.graph.edge_list()
        )

    def neighbors(self, node):
        return (
            self.index_to_node[neighbor_node]
            for neighbor_node in self.graph.neighbors(self.node_to_index[node])
        )

    def out_edges(self, node):
        return (
            (self.graph.get_node_data(edge[0]), self.graph.get_node_data(edge[1]))
            for edge in self.graph.out_edges(self.node_to_index[node])
        )

    def shortest_paths_without_high_low(self, source):
        subgraph = self.graph.copy()  # This is a shallow copy.
        # We know that nodes 0 and 1 are the 'high' and 'low' nodes.
        subgraph.remove_nodes_from([0, 1])
        return self._shortest_paths_from(subgraph, source)

    def shortest_paths_from_high(self):
        return self._shortest_paths_from(self.graph, "high")

    def shortest_paths_from_low(self):
        return self._shortest_paths_from(self.graph, "low")

    # The following methods are helpers for the shortest paths methods.

    def _shortest_paths_from(self, graph, source):
        paths = rx.graph_dijkstra_shortest_paths(
            graph, self.node_to_index[source], weight_fn=rx_weight_fn
        )
        return self._all_path_weights(paths), self._all_path_nodes(paths)

    def _all_path_weights(self, paths):
        return {
            self.index_to_node[target]: self._path_weight(path) for (target, path) in paths.items()
        }

    def _path_weight(self, path):
        weight = 0
        for e in range(len(path) - 1):
            weight += int(rx_weight_fn(self.graph.get_edge_data(path[e], path[e + 1])))
        return weight

    def _all_path_nodes(self, paths):
        return {
            self.index_to_node[target]: self._path_nodes(path) for (target, path) in paths.items()
        }

    def _path_nodes(self, path):
        nodes = []
        for index in path:
            nodes.append(self.index_to_node[index])
        return nodes


def rx_weight_fn(edge):
    """A function for returning the weight from the common vertex."""
    return float(edge["weight"])
