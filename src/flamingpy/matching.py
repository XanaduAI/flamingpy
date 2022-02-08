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
""" The matching graph interface and some implementation.
"""
from abc import ABC
from dataclasses import dataclass
from typing import Iterable, List, Tuple, TypeVar, Union
import itertools as it
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import numpy as np
import retworkx as rx
from flamingpy import lemon


# These abstract the concept of nodes and edges of a graph
# as placeholder types.
Node = TypeVar("Node")
Edge = Tuple[Node, Node]


class MatchingGraph(ABC):
    """The matching graph base class.

    This represents a graph for which we can compute a minimum weight perfect matching.
    See ft_stack.matching.NxMatchingGraph for a reference implementation using the networkx
    library as a backend.
    """

    # The type representing the weight of an edge.
    Weight = TypeVar("Weight")

    def __init__(self):
        self.virtual_points = list()

    def add_edge(self, edge: Edge, weight: Weight, path: List[Node] = []):
        """Add an edge into the graph.

        Args:
            edge: the extremal nodes of the edges.
            weight: a positive weight for the edge.
            path: an optional list of nodes connecting the two extremal nodes.
        """
        raise NotImplementedError

    def edge_weight(self, edge: Edge) -> Weight:
        """Return the weight of the given edge or raise an exception
        if the edge is not part of the graph.
        """
        raise NotImplementedError

    def edge_path(self, edge: Edge) -> List[Node]:
        """Return the path for the given edge or raise an exception
        if the edge is not part of the graph.
        """
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

    def with_edges_from_dec_graph(self, dec_graph: nx.Graph):
        """Update the matching graph from the decoding graph G.

        The matching graph has as half of its nodes the odd-parity stabilizers.
        The edge connecting two nodes corresponds to the weight of the
        minimum-weight-path between the nodes in the decoding graph.
        Additionally, each unsatisfied stabilizer is connected to a unique boundary point
        (for now from a primal bundary) located at the shortest weighted distance from
        the stabilizer. Between each other, the boundary points are
        connected by an edge of weight 0. The output graph stores the
        indices of the used boundary points under the 'used_boundary_point'
        attribute. Paths are stored under the 'paths' attribute of edges,
        and 'inverse_weights' are also stored, for the benefit of maximum-
        weight-matching algorithms

        Args:
            dec_graph (networkx.Graph): the decoding graph, storing information
                about indices of odd-parity-cubes (under 'odd_cubes' graph
                attribute) and boundary points (under 'boundary_points').
        Returns:
            networkx.Graph: the updated matching graph.
        """
        # Get the indices of the odd parity cubes from the decoding graph.
        odd_parity_inds = dec_graph.graph["odd_cubes"]

        # Run the matching algorithm first without the 'high' and 'low points
        real_points = dec_graph.graph["real_points"]

        alg = sp.single_source_dijkstra
        # Combinations of odd-parity cubes.
        odd_ind_dict = {i: [] for i in odd_parity_inds[:-1]}
        odd_combs = it.combinations(odd_parity_inds, 2)
        for pair in odd_combs:
            odd_ind_dict[pair[0]] += [pair[1]]
        # Find the shortest paths between odd-parity cubes.
        for cube1 in odd_parity_inds[:-1]:
            lengths, paths = alg(dec_graph.subgraph(real_points), cube1)
            for cube2 in odd_ind_dict[cube1]:
                length = lengths[cube2]
                path = paths[cube2]
                # Add edge to the matching graph between the cubes, with weight
                # equal to the length of the shortest path.
                # TODO: Is the behavior correct for negative weights, or do I
                # want 1/weight or max_num - weight?
                self.add_edge((cube1, cube2), length, path)

        if dec_graph.graph["boundary_points"]:
            i = 0
            low_lengths, low_paths = alg(dec_graph, "low")
            high_lengths, high_paths = alg(dec_graph, "high")
            for cube in odd_parity_inds:
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
                i += 1
                self.virtual_points.append(virtual_point)
            # Add edge with weight 0 between any two virtual excitations.
            for edge in it.combinations(self.virtual_points, 2):
                self.add_edge(edge, 0)
        return self


class NxMatchingGraph(MatchingGraph):
    """A matching graph backed by networkx.

    The edge weights can be either of type int or float.
    See the MatchingGraph class for more details.
    """

    Weight = Union[int, float]

    def __init__(self):
        self.graph = nx.Graph()
        MatchingGraph.__init__(self)

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

        This is basically free since this class already use networkx
        to represent a graph.
        """
        nx = NxMatchingGraph()
        nx.graph = self.graph
        return nx


@dataclass
class RxEdgePayload:
    """The edge payload for the RxMatchingGraph.

    Args:
        weight:  the minimum-weight-path in the decoding graph between the extremal nodes.
        path : the nodes in the decoding graph along the minimum-weight-path.
    """

    weight: int
    path: List[Node]


class RxMatchingGraph(MatchingGraph):
    """A matching graph backed by retworkx.

    The edge weights must be of type int.
    See the MatchingGraph class for more details.
    """

    Weight = int

    def __init__(self):
        self.graph = rx.PyGraph(multigraph=False)
        self.node_map = dict()
        MatchingGraph.__init__(self)

    def node_index(self, node):
        """Returns the index of the corresponding node.

        If the node doesn't exist in the graph,
        a new index is generated.
        """
        if node in self.node_map:
            return self.node_map[node]
        else:
            index = self.graph.add_node(node)
            self.node_map[node] = index
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
        return self.graph.get_edge_data(self.node_index(edge[0]), self.node_index(edge[1])).weight

    def min_weight_perfect_matching(self) -> List[Edge]:
        return list(
            rx.max_weight_matching(
                self.graph,
                max_cardinality=True,
                weight_fn=lambda edge: -1 * edge.weight,
            )
        )

    def to_nx(self):
        """Convert the same graph into a NxMatchingGraph.

        This involves converting the retworkx graph representation
        to a networkx graph representation.
        """
        nx_graph = NxMatchingGraph()
        for (idx0, idx1) in iter(self.graph.edge_list()):
            edge = (self.graph.nodes()[idx0], self.graph.nodes()[idx1])
            nx_graph.add_edge(edge, self.edge_weight(edge), self.edge_path(edge))
        return nx_graph
