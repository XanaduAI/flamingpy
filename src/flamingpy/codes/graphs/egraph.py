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
"""A class for representing qubit code graph states."""
import networkx as nx
import numpy as np
import scipy.sparse as sp


class EGraph(nx.Graph):
    """An enhanced graph based on a NetworkX Graph.

    A class for adding some functionality to a NetworkX graph
    with some short-hand/convenience methods.

    Attributes:
        indexer (str): method for indexing the nodes; 'default' for
            Python's sorted function; 'macronodes' for rounding
            micronodes to integers, sorting those, and
            furthermore sorting the micronodes within each macronodes,
            all using Python's 'sorted'.
        to_indices (dict): if self.index_generator() has been run,
            a dictionary of the form {points: indices}
        to_points (dict): if self.index_generator() has been run,
            a dictionary of the form {indices: points}
        adj_mat (np.array): if self.adj_generator() has been run,
            the adjacency mtrix of the graph.
    """

    def __init__(self, *args, indexer="default", macronodes=False, **kwargs):
        """Initialize an EGraph (itself an NetworkX graph)."""
        super().__init__(*args, **kwargs)
        self.indexer = indexer
        self._macronodes = macronodes
        if macronodes:
            self.macro = nx.Graph()
        self.to_indices = None
        self.to_points = None
        self.adj_mat = None

    def index_generator(self):
        """Return a relabelled graph with indices as labels.

        Point tuples are stored in the 'pos' attribute of the new graph.
        Use the default sort as the index mapping.
        """
        # TODO: Let user specify index mapping.
        # TODO: SortedDict implementation.
        N = self.order()
        if self.to_indices is not None:
            return self.to_indices
        if self.indexer == "default" and not self._macronodes:
            ind_dict = dict(zip(sorted(self.nodes()), range(N)))
        elif self._macronodes:
            macro_graph = self.macro
            for node in self.nodes():
                rounded = tuple(np.round(node).astype(int))
                macro_graph.nodes[rounded]["micronodes"].append(node)
            sorted_macro = sorted(macro_graph)
            points = []
            for vertex in sorted_macro:
                points += self.macro.nodes[vertex]["micronodes"]
            ind_dict = {points[i]: i for i in range(N)}
        self.to_indices = ind_dict
        self.to_points = {index: point for point, index in ind_dict.items()}
        return ind_dict

    def adj_generator(self, sparse=True):
        """Return the adjacency matrix of the graph.

        Indices correspond to sorted nodes.
        """
        if self.adj_mat is not None:
            return self.adj_mat
        if self.to_points is None:
            self.index_generator()
        # TODO: SortedDict implementation.
        sorted_nodes = [self.to_points[i] for i in range(self.order())]
        # TODO: New data type in case of fancier weights.
        if sparse:
            adj = nx.to_scipy_sparse_matrix(self, nodelist=sorted_nodes, dtype=np.int8)
        else:
            adj = nx.to_numpy_array(self, nodelist=sorted_nodes, dtype=np.int8)
        self.adj_mat = adj
        # TODO: Heat map?
        return adj

    def slice_coords(self, plane, number):
        """Obtain all the coordinates in an x, y, or z slice.

        Args:
            plane (str): 'x', 'y', or 'z', denoting the slice direction
            number (int): the index of the slice

        Returns:
            list of tuples: the coordinates of the slice.
        """
        plane_dict = {"x": 0, "y": 1, "z": 2}
        plane_ind = plane_dict[plane]
        coords = [point for point in self.nodes if point[plane_ind] == number]
        return coords
