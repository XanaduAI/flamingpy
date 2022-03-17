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
"""The matching graph interface and some implementation."""

# pylint: disable=import-outside-toplevel

from abc import ABC
from dataclasses import dataclass
import itertools as it
from typing import Iterable, List, Tuple, TypeVar, Union

import networkx as nx
import numpy as np
import retworkx as rx

from flamingpy.decoders.mwpm import lemon


# These abstract the concept of nodes and edges of a graph
# as placeholder types.
Node = TypeVar("Node")
Edge = Tuple[Node, Node]


class MatchingGraph(ABC):
    """The matching graph base class.

    This represents a graph for which we can compute a minimum weight perfect
    matching. See flammingpy.matching.NxMatchingGraph for a reference
    implementation using the networkx library as a backend.

    Arguments:
        code (SurfaceCode, optional): If provided, the matching graph
        includes an edge for each pair of unsatisfied stabilizer nodes
        in the code together with the weight of the corresponding correction.
    """

    # The type representing the weight of an edge.
    Weight = TypeVar("Weight")

    def __init__(self, ec, code=None):
        self.virtual_points = list()
        self.ec = ec
        if code is not None:
            self.with_edges_from(code, ec)

    def add_edge(self, edge: Edge, weight: Weight, path: List[Node] = []):
        """Add an edge into the graph.

        Args:
            edge: the extremal nodes of the edges.
            weight: a positive weight for the edge.
            path: an optional list of nodes connecting the two extremal nodes.
        """
        raise NotImplementedError

    def edge_weight(self, edge: Edge) -> Weight:
        """Return the weight of the given edge or raise an exception if the
        edge is not part of the graph."""
        raise NotImplementedError

    def edge_path(self, edge: Edge) -> List[Node]:
        """Return the path for the given edge or raise an exception f the edge
        is not part of the graph."""
        raise NotImplementedError

    def min_weight_perfect_matching(self) -> List[Edge]:
        """Compute a minimum weight perfect matching for the graph.

        Returns:
            The list of edges in the matching.
        """
        raise NotImplementedError

    def total_weight_of(self, matching: Iterable[Edge]) -> Weight:
        """Return the sum of the weight of each edge in a matching.

        Args:
            The pairs of nodes in the matching.

        Returns:
            The sum of the weights.
        """
        return sum(self.edge_weight(edge) for edge in matching)

    def with_edges_from(self, code, ec):
        """Update the matching graph from the given code.

        The matching graph has as half of its nodes the odd-parity stabilizers.
        The edge connecting two nodes corresponds to the weight of the
        minimum-weight-path between the nodes in the stabilizer graph of the
        code. Additionally, each unsatisfied stabilizer is connected to a unique
        boundary point (for now from a primal bundary) located at the shortest
        weighted distance from the stabilizer. Between each other, the boundary
        points are connected by an edge of weight 0. The output graph stores
        the indices of the used boundary points under the 'used_boundary_point'
        attribute. Paths are stored under the 'paths' attribute of edges, and
        'inverse_weights' are also stored, for the benefit of maximum-weight
        matching algorithms

        Args:
            code (SurfaceCode): The code from which to build the edges.
        """
        stab_graph = getattr(code, ec + "_stab_graph")
        self = self._with_edges_between_real_odd_nodes(code, ec)
        if stab_graph.has_bound_points():
            return self._with_edges_from_low_or_high_connector(code, ec)
        else:
            return self

    def _with_edges_between_real_odd_nodes(self, code, ec):
        # Get the indices of the odd parity cubes from the stabilizer graph.
        stab_graph = getattr(code, ec + "_stab_graph")
        odd_parity_stabs = list(stab_graph.odd_parity_stabilizers())
        # Combinations of odd-parity cubes.
        odd_adjacency = {i: [] for i in odd_parity_stabs[:-1]}
        for pair in it.combinations(odd_parity_stabs, 2):
            odd_adjacency[pair[0]] += [pair[1]]
        # Find the shortest paths between odd-parity stabs.
        for stab1 in odd_parity_stabs[:-1]:
            lengths, paths = stab_graph.shortest_paths_without_high_low(stab1, code)
            for stab2 in odd_adjacency[stab1]:
                length = lengths[stab2]
                path = paths[stab2]
                # Add edge to the matching graph between the stabs, with weight
                # equal to the length of the shortest path.
                self.add_edge((stab1, stab2), length, path)
        return self

    def _with_edges_from_low_or_high_connector(self, code, ec):
        stab_graph = getattr(code, ec + "_stab_graph")
        low_lengths, low_paths = stab_graph.shortest_paths_from_low(code)
        high_lengths, high_paths = stab_graph.shortest_paths_from_high(code)
        for i, cube in enumerate(stab_graph.odd_parity_stabilizers()):
            distances = (low_lengths[cube], high_lengths[cube])
            where_shortest = np.argmin(distances)
            if where_shortest == 0:
                length = low_lengths[cube]
                full_path = low_paths[cube]
            if where_shortest == 1:
                length = high_lengths[cube]
                full_path = high_paths[cube]
            point = full_path[1]
            virtual_point = (point, i)
            path = full_path[1:]
            # Add edge to the matching graph between the cube and
            # the virtual excitation corresponding to the boundary
            # vertex, with weight equal to the length of the shortest
            # path.
            self.add_edge(
                (cube, virtual_point),
                length,
                path,
            )
            self.virtual_points.append(virtual_point)
        # Add edge with weight 0 between any two virtual excitations.
        for edge in it.combinations(self.virtual_points, 2):
            self.add_edge(edge, 0)
        return self

    def draw(self, **kwargs):
        """Draw the matching graph with matplotlib.

        See flamingpy.utils.viz.draw_dec_graph for more details.
        """
        from flamingpy.utils.viz import draw_dec_graph

        draw_dec_graph(self, **kwargs)


class NxMatchingGraph(MatchingGraph):
    """A matching graph backed by networkx.

    The edge weights can be either of type int or float.

    See the MatchingGraph class for more details.
    """

    Weight = Union[int, float]

    def __init__(self, ec, code=None):
        self.graph = nx.Graph()
        MatchingGraph.__init__(self, ec, code)

    def add_edge(self, edge: Edge, weight: Weight, path: List[Node] = []):
        self.graph.add_edge(edge[0], edge[1], weight=weight, inverse_weight=-weight, path=path)

    def edge_weight(self, edge: Edge):
        return self.graph.edges[edge]["weight"]

    def edge_path(self, edge: Edge):
        return self.graph.edges[edge]["path"]

    def min_weight_perfect_matching(self) -> List[Edge]:
        return nx.max_weight_matching(self.graph, maxcardinality=True, weight="inverse_weight")


class LemonMatchingGraph(NxMatchingGraph):
    """A matching graph class backed by the Lemon package."""

    def min_weight_perfect_matching(self) -> List[Edge]:
        return lemon.max_weight_matching(self.graph, weight="inverse_weight")

    def to_nx(self):
        """Return the same graph wrapped into a NxMatchingGraph.

        This is basically free since this class already uses networkx to
        represent a graph.
        """
        nx = NxMatchingGraph(self.ec)
        nx.graph = self.graph
        return nx


@dataclass
class RxEdgePayload:
    """The edge payload for the RxMatchingGraph.

    Args:
        weight: the minimum-weight-path in the stabilizer graph between the
            extremal nodes.
        path: the nodes in the stabilizer graph along the minimum-weight-path.
    """

    weight: int
    path: List[Node]


class RxMatchingGraph(MatchingGraph):
    """A matching graph backed by retworkx.

    The edge weights must be of type int.

    See the MatchingGraph class for more details.
    """

    Weight = int

    def __init__(self, ec, code=None):
        self.graph = rx.PyGraph(multigraph=False)
        self.node_to_index = dict()
        self.index_to_node = dict()
        MatchingGraph.__init__(self, ec, code)

    def node_index(self, node):
        """Returns the index of the corresponding node.

        If the node doesn't exist in the graph, a new index is
        generated.
        """
        if node in self.node_to_index:
            return self.node_to_index[node]
        else:
            index = self.graph.add_node(node)
            self.node_to_index[node] = index
            self.index_to_node[index] = node
            return index

    def add_edge(self, edge: Edge, weight: Weight, path: List[Node] = []):
        self.graph.add_edge(
            self.node_index(edge[0]),
            self.node_index(edge[1]),
            RxEdgePayload(weight, path),
        )

    def edge_weight(self, edge: Edge):
        return self.graph.get_edge_data(self.node_index(edge[0]), self.node_index(edge[1])).weight

    def edge_path(self, edge: Edge):
        return self.graph.get_edge_data(self.node_index(edge[0]), self.node_index(edge[1])).path

    def min_weight_perfect_matching(self) -> List[Edge]:
        matches = rx.max_weight_matching(
            self.graph,
            max_cardinality=True,
            weight_fn=lambda edge: -1 * edge.weight,
        )
        return [(self.index_to_node[pair[0]], self.index_to_node[pair[1]]) for pair in matches]

    def to_nx(self):
        """Convert the same graph into a NxMatchingGraph.

        This involves converting the retworkx graph representation to a
        networkx graph representation.
        """
        nx_graph = NxMatchingGraph(self.ec)
        for (idx0, idx1) in iter(self.graph.edge_list()):
            edge = (self.graph.nodes()[idx0], self.graph.nodes()[idx1])
            nx_graph.add_edge(edge, self.edge_weight(edge), self.edge_path(edge))
        return nx_graph
