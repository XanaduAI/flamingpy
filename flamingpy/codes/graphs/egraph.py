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
"""A class for representing quantum graph states."""

# pylint: disable=import-outside-toplevel
from typing import Union

import warnings
import networkx as nx
import numpy as np


def macronize(can_graph, pad_boundary=False, disp=0.1):
    """Create a macronode graph out of canonical graph can_graph.

    Assume can_graph represents a 'canonical' or 'reduced' graph state.
    Replace each vertex v in can_graph with a macronode (a collection of n
    vertices, where n is the size of the neighborhood of v). This is
    achieved by replacing each pair of vertices in can_graph connected by an
    edge with a 'dumbbell'.

    Assumes that: (1) the nodes of G are three-tuples; (2) edges associated
    with periodic boundary conditions, if they exist, have the attribute
    'periodic' set to True.

    If there is a list of perfect qubits stored as a graph attribute
    (specifically, can_graph.graph["perfect_points"]), this list is updated to
    reflect the new perfect qubits in the macronode graph.

    Args:
        can_graph (nx.Graph): the graph to macronize. All graph, node, and
            edge attributes of can_graph are preserved. Will access and process
            can_graph.graph["perfect_points"] if it exists.
        disp (float, optional): how much to displace the nodes
            within each macronode from the central vertex. This number
            should be small and positive, and no larger than 0.5.
        pad_boundary (bool, optional): if True, pad each boundary
            macronode to have the same number of nodes as a bulk
            macronode. For now, this only works for the RHG lattice
            connectivity (i.e. pad to 4 nodes).

    Returns:
        nx.Graph: the macronized graph, with a macronode-to-micronode
            dictionary stored as a graph attribute.
    """
    if disp >= 0.5 or disp < 0:
        raise ValueError("Please set disp to a positive value strictly less than 0.5.")
    macro_graph = nx.Graph(**can_graph.graph)
    # The macronode-to-micronode dictionary
    macro_dict = {}
    for edge in can_graph.edges:
        old_point_1, old_point_2 = np.array(edge[0]), np.array(edge[1])
        direction_vec = old_point_2 - old_point_1
        distance = np.linalg.norm(direction_vec)
        periodic_flip = -1 if can_graph.edges[edge].get("periodic") else 1
        shortened_vec = periodic_flip * disp * direction_vec / distance
        shortened_vec = np.round(shortened_vec, -int(np.log10(disp)) + 2)
        new_point_1 = tuple(old_point_1 + shortened_vec)
        new_point_2 = tuple(old_point_2 - shortened_vec)
        macro_graph.add_node(new_point_1, **can_graph.nodes[edge[0]])
        macro_graph.add_node(new_point_2, **can_graph.nodes[edge[1]])
        macro_graph.add_edge(new_point_1, new_point_2, **can_graph.edges[edge])

        # Add to the macronode-to-micronode dictionary
        if not macro_dict.get(edge[0]):
            macro_dict[edge[0]] = []
        if not macro_dict.get(edge[1]):
            macro_dict[edge[1]] = []
        macro_dict[edge[0]].append(new_point_1)
        macro_dict[edge[1]].append(new_point_2)

    if pad_boundary:
        for node in can_graph:
            macro_size = len(macro_dict[node])
            if macro_size < 4:
                n_new = 4 - macro_size
                for i in range(n_new):
                    new_node = list(node)
                    new_node[i] = new_node[i] + 0.05
                    new_node_tup = tuple(new_node)
                    macro_graph.add_node(new_node_tup, **can_graph.nodes[node])
                    macro_dict[node].append(new_node_tup)

    macro_graph.graph["macro_dict"] = macro_dict
    old_perfect_points = can_graph.graph.get("perfect_points")
    if old_perfect_points:
        new_perfect_qubits = [macro_dict[point] for point in old_perfect_points]
        macro_graph.graph["perfect_points"] = [a for b in new_perfect_qubits for a in b]
    return macro_graph


class EGraph(nx.Graph):
    """An enhanced graph for representing quantum graph states.

    A class that builds on a NetworkX graph to better represent graph states.
    Includes indexing, drawing, and convenience methods.

    Attributes:
        macro_to_micro (dict): if macronodes is set to True, the macro_dict
            object from the underlying graph (None or a dictionary of the form
            {central coordinate of macronode: [all micronode coordinates]})
        to_indices (dict): if self.index_generator() has been run,
            a dictionary of the form {points: indices}
        to_points (dict): if self.index_generator() has been run,
            a dictionary of the form {indices: points}
        adj_mat (np.array): if self.adj_generator() has been run,
            the adjacency matrix of the graph.
    """

    def __init__(self, *args, indexer="default", macronodes=False, **kwargs):
        super().__init__(*args, **kwargs)
        self._macronodes = macronodes
        self.macro_to_micro = None
        if macronodes:
            self.macro_to_micro = self.graph.get("macro_dict")
        self.to_indices = None
        self.to_points = None
        self.adj_mat = None

    def index_generator(self):
        """Generate indices for the nodes of self.

        Set the to_indices and to_points attribute of self with points-
        to-indices and indices-to-points dictionaries, respectively.
        Indices are generated using the built-in sorted function (in
        case of macronodes, the central/integer coordinates are sorted).
        """
        # TODO: User-specified index mapping?
        N = self.order()
        if self.to_indices is not None:
            return self.to_indices
        if not self._macronodes:
            ind_dict = dict(zip(sorted(self.nodes()), range(N)))
        elif self._macronodes:
            sorted_macro = sorted(self.macro_to_micro)
            points = []
            for vertex in sorted_macro:
                points += self.macro_to_micro[vertex]
            ind_dict = {points[i]: i for i in range(N)}
        self.to_indices = ind_dict
        self.to_points = {index: point for point, index in ind_dict.items()}
        return ind_dict

    def adj_generator(self, sparse=True):
        """Return the correctly indexed adjacency matrix of the graph and set
        the self.adj_mat attribute.

        Calling the NetworkX adjacency matrix methods with default
        options may create a mismatch between the indices of the
        rows/columns of the matrix and the indices generated by
        self.index_generator(). This function demands that the indices
        match.
        """
        if self.adj_mat is not None:
            return self.adj_mat
        if self.to_points is None:
            self.index_generator()
        sorted_nodes = [self.to_points[i] for i in range(self.order())]
        # TODO: Reconsider data type for more intricate weights.
        if sparse:
            adj = nx.to_scipy_sparse_array(self, nodelist=sorted_nodes, dtype=np.int8)
        else:
            adj = nx.to_numpy_array(self, nodelist=sorted_nodes, dtype=np.int8)
        self.adj_mat = adj
        return adj

    def slice_coords(self, plane, number):
        """Obtain all the coordinates in an x, y, or z slice of self.

        Args:
            plane (str): 'x', 'y', or 'z', denoting the slice direction
            number (int): the index of the slice. The allowable range is from 0
                to the total number of slices in the given direction.

        Returns:
            list of tuples: the coordinates of the slice.
        """
        plane_dict = {"x": 0, "y": 1, "z": 2}
        plane_ind = plane_dict[plane]
        coords = [point for point in self.nodes if point[plane_ind] == number]
        return coords

    def macronize(self, pad_boundary=False, disp=0.1):
        """Return a new, macronized version of self.

        See egraph.macronize for more details.
        """
        return EGraph(macronize(self, pad_boundary, disp), macronodes=True)

    def draw(self, **kwargs):
        """Draw the graph state with Matplotlib.

        See flamingpy.utils.viz.draw_EGraph for more details.
        """
        from flamingpy.utils.viz import draw_EGraph

        fig, ax = draw_EGraph(self, **kwargs)
        return fig, ax

    def draw_adj(self, **kwargs):
        """Draw the adjacency matrix with matplotlib.

        See flamingpy.utils.viz.plot_mat_heat_map for more details.
        """
        from flamingpy.utils.viz import plot_mat_heat_map

        adj = self.adj_generator(sparse=False)
        return plot_mat_heat_map(adj, **kwargs)

    def add_node(self, node: tuple, **kwargs):
        """Overwrite add_node function from nx.Graph to raise a warning
        whenever add_qubit should be used."""
        if not (self.to_indices is None and self.adj_mat is None):
            warnings.warn(
                "EGraph attributes are already populated. To add a"
                " node, please use the add_qubit function."
            )
        super().add_node(node, **kwargs)

    def remove_node(self, node: tuple):
        """Overwrite remove_node function from nx.Graph to raise a warning
        whenever remove_node should be used."""
        if not (self.to_indices is None and self.adj_mat is None):
            warnings.warn(
                "EGraph attributes are already populated. To remove a"
                " node, please use the remove_qubit function."
            )
        super().remove_node(node)

    def add_qubit(
        self,
        qubit: Union[None, tuple] = None,
        neighbors: Union[None, list] = None,
        macro: Union[None, tuple] = None,
    ) -> None:
        """Add qubit to EGraph at position with neigbours in
        existing_neighbours.

        Args:
            qubit (3D tuple, or None): qubit to add. If it is a 3D tuple,
            the qubit is added in that position. If qubit is None, the qubit is
            positioned one unit further than the maximum position in the z
            direction in position (0,0,z_max + 1).

            neighbours (list[int], list[tuple], or None): neighbors of qubit specified
            with indices or positions, respectively.
        """

        # Makes sure that input qubit value and type is supported
        if isinstance(qubit, tuple):
            if len(qubit) != 3:
                raise ValueError(
                    "The position should be a 3D tuple, but it was given a "
                    + f"{len(qubit)}D tuple"
                )
            if qubit in self:
                raise ValueError(
                    f"Could not add qubit because there is already a node in position {qubit}."
                )
        elif qubit is None:
            z_max = max(map(lambda tup: tup[2], self.nodes()))
            qubit = (0, 0, z_max + 1)
        else:
            raise TypeError(
                "Qubit type not supported. Excepted 3D tuple or None, but was"
                + f" given {type(qubit)}"
            )

        # Check that if macro is not None, then self should be macronized EGraph
        if isinstance(macro, tuple) and self.macro_to_micro is None:
            raise ValueError(
                "Cannot add a qubit to macronode because current EGraph is not macronized."
            )

        # Add node
        self.add_node(qubit)

        # Update dictionaries when adding qubit
        self._update_attributes_add_qubit(qubit, macro)

        # Update neighbors
        if neighbors is not None:
            if len(neighbors) != 0:
                if isinstance(neighbors[0], int):
                    neighborhood = [(self.to_points[ind], qubit) for ind in neighbors]
                elif isinstance(neighbors[0], tuple):
                    neighborhood = [(tup, qubit) for tup in neighbors]
                else:
                    raise TypeError(
                        "Unsupported type of neighbors. Expected tuples but"
                        f" {type(neighbors[0])} was given. Type ints can only be used when "
                        "macro is False."
                    )

                self.add_edges_from(neighborhood)

        self.adj_mat = None

    def _update_attributes_add_qubit(self, qubit, macro):
        """Update self.macro_to_micro, self.to_points, and self.to_indices."""
        # NOTE: this method was created to reduce complexity of self.add_qubit()

        # Add qubit to macro_to_micro dictionary
        if isinstance(macro, tuple):
            if macro in self.macro_to_micro:
                self.macro_to_micro[macro].append(qubit)
            else:
                self.macro_to_micro[macro] = [qubit]

        # Update dictionaries
        if self.to_indices is not None:
            new_index = max(self.to_points.keys()) + 1
            self.to_points[new_index] = qubit
            self.to_indices[qubit] = new_index

    def remove_qubit(self, qubit: Union[tuple, int]) -> None:
        """Remove qubit from EGraph.

        Args:
            qubit (3D tuple or int) : If 3D tuple, remove qubit at that
            position. For not None self.to_points and not None self.to_indices:
            If int, remove the qubit at that index. If None, remove
            the last qubit.
        """
        # Remove qubit if dictionaries are not initialized
        if self.to_indices is None:
            if isinstance(qubit, tuple):
                self.remove_node(qubit)
            else:
                self.index_generator()

        # Remove qubit if dictionaries are initialized
        if self.to_indices is not None:
            if not isinstance(qubit, (int, tuple)):
                raise TypeError(
                    "Qubit type not supported. Excepted 3D tuple, int, or None, "
                    + f"but was given {type(qubit)}"
                )
            qubit = self.to_points[qubit] if isinstance(qubit, int) else qubit
            self.remove_node(qubit)
            self.to_indices.pop(qubit)
            self.to_points = {v: k for k, v in self.to_indices.items()}

        # Remove qubit from macro_to_micro dict if EGraph is macronized
        if self.macro_to_micro is not None:
            for k in set(self.macro_to_micro.keys()):
                if qubit in self.macro_to_micro[k]:
                    self.macro_to_micro[k].remove(qubit)

                    # If macronode is empty, then delete it from macro_to_micro
                    if not self.macro_to_micro[k]:
                        del self.macro_to_micro[k]

        self.adj_mat = None
