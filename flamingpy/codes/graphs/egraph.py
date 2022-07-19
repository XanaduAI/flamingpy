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

    def is_lc_equivalent(self, graph2, clifford_form="tensor"):
        """Checks if two EGraph objects are LC equivalent by finding a local clifford operation on
        n qubits, where n is the number of nodes of the graphs.
        Args:
            graph2: An EGraph object to test equivalence with
            clifford_form: A string describing the output form of Local Clifford operation when it
            exists. Default string is 'tensor' which returns a list of length n of 2x2 numpy arrays
            corresponding to single qubit tensor factors. 'global' returns a single 2nx2n numpy
            array corresponding to the global operator acting on all n qubits.
        Returns:
            (equivalent, clifford): A tuple where 'equivalent' is a boolean.
            If equivalent is True, 'clifford' is the local clifford output according to
            'clifford_form' specification.
        """
        # get adjacency matrices of input graphs
        self.adj_generator(sparse=False)
        G = self.adj_mat
        graph2.adj_generator(sparse=False)
        H = graph2.adj_mat

        # handle input and edge cases
        #
        # adjacency matrices of input graphs must have type 'np.ndarray'
        if not isinstance(G, np.ndarray) or not isinstance(H, np.ndarray):
            raise ValueError(
                "Input EGraphs must have their 'adj_mat' property as type 'numpy.ndarray'."
            )
        # demand input graphs must have nonzero number of nodes.
        # check adjacency matrices for same non-empty square shape
        if np.shape(G) == (0, 0) or np.shape(H) == (0, 0):
            # raise ValueError('Input Graphs must be non-empty')
            # Mark this case as not equivalent
            return False, None
        if np.shape(G) != np.shape(H):
            # raise ValueError('Input Graphs must have same number of nodes.')
            # Mark this case as not equivalent
            return False, None
        # adjacency matrices must be square
        if np.shape(G)[0] != np.shape(G)[1]:
            raise ValueError("Input matrices must be square.")

        # construct system of constraints for two adjacency matrices G and H
        system_constraints = self.__system_constraints(G, H)

        # construct nullspace basis of system of constraints
        nullspace_basis = self.__nullspace_basis(system_constraints)

        # search nullspace for solution vectors
        solution_vector = self.__search_nullspace(nullspace_basis, np.shape(G)[0])

        # if graphs are not equivalent set clifford_output to None
        if solution_vector is None:
            equivalent = False
            clifford_output = None
        # if graphs are equivalent and clifford_form == "tensor", convert clifford_output to tensor factors
        elif clifford_form == "tensor":
            equivalent = True
            clifford_output = self.__clifford_vec_to_tensors(solution_vector)
        # if graphs are equivalent and clifford_form == "global", convert clifford_output to global form
        elif clifford_form == "global":
            equivalent = True    
            clifford_output = self.__clifford_vec_to_global(solution_vector)
        # return equivalence and clifford_output
        return equivalent, clifford_output

    def __system_constraints(self, G, H):
        """Constructs a binary system of equations that two graph adjacency matrices G and H must
        satisfy for lc equivalence.
        Args:
            G and H are n x n adjacency matrices defined for n qubit graph states with n nodes.
        Returns:
            (numpy.array): A system given as a numpy array of shape n^2 x 4n with n^2 equations in
            4n unknowns
        """
        # set the number of qubits
        n = np.shape(G)[0]
        # define empty block matrix
        M = []
        for j in range(n):
            # A block
            A = np.diag(G[j])
            # B block
            B = np.zeros((n, n))
            B[j, j] = 1
            # C block
            C = np.array([[G[i, j] * H[i, k] for i in range(n)] for k in range(n)], dtype=int)
            # D block
            D = np.zeros((n, n), dtype=int)
            for k in range(n):
                D[k, j] = H[j, k]
            # add row block to M
            M.append([A, B, C, D])
        # define numpy block matrix for system of constraints
        return np.block(M).astype(int)

    def __RREform_mod2(self, M, max_cols=None):
        """Puts a binary matrix into Row Reduced Echelon form modulo 2 up to a maximum number of
        columns given by max_cols.
        Args: 
            M (numpy.array): A numpy array
            max_cols (int): Specifies the maximum number of columns of the input array to reduce
        Returns:
            (numpy.array): A numpy array of the Row Reduced Echelon form modulo 2
        """
        # number of columns to apply row reduction
        max_cols = M.shape[1] if max_cols is None else max_cols
        # row reduced matrix
        R = M.copy()
        # pivot index
        p = 0
        # perform row reduction
        for j in range(max_cols):
            # look for pivot column `j` at or below row `p`
            index = np.nonzero(R[p:, j])[0]
            if index.size == 0:
                continue
            # pivot row
            i = p + index[0]
            # interchange `p` and `i`
            R[[p, i], :] = R[[i, p], :]
            # set pivot to 1
            R[p, :] = R[p, :] / R[p, j]
            # add zeros above and below the pivot
            index = np.nonzero(R[:, j])[0].tolist()
            index.remove(p)
            R[index, :] = (R[index, :] - np.multiply.outer(R[index, j], R[p, :])) % 2
            # increment pivot index until final row is reached
            p += 1
            if p == R.shape[0]:
                break
        return R, p

    def __nullspace_basis(self, M):
        """Constructs an array whose rows are binary basis vectors of the right null space of
        input matrix array.
        Args:
            M (numpy.array): A binary matrix
        Returns:
            (numpy.array): A numpy array whose rows are basis vectors of the right null space
        """
        # find left Null space of transposed system, which is equal to right null space
        M_transposed = M.T
        m_rows, n_cols = M_transposed.shape
        # construct augmented block matrix A = [M | I]
        A = np.concatenate((M_transposed, np.eye(m_rows, dtype=int)), axis=-1)
        # row reduce left M block of augmented matrix
        R, p = self.__RREform_mod2(A, max_cols=n_cols)
        N = R[p:, n_cols:]
        nullspace_basis, _ = self.__RREform_mod2(N)
        return nullspace_basis

    def __search_nullspace(self, basis, n):
        """Search through sums of pairs of basis vectors of the null space, and check if any
        satisfy the determinant constraints.
        Args:
            basis (numpy.array): A numpy array whose rows are basis vectors of the nullspace
        Returns: 
            clifford_output: If solution is found, clifford_output is a numpy vector specifying
            the local clifford. If no solution is found, clifford_output is None.
        """
        # number of basis vectors of null space
        d = basis.shape[0]
        # boolean for LC equivalence
        equivalent = False
        # array to store potential solutions
        sols = []
        # search through distinct pairs of basis vectors
        for i in range(d):
            for j in range(i):
                # boolean for checking solution constraints
                sat = True
                # potential solution as sum (modulo 2) of basis pair
                sol = np.vectorize(lambda x: x % 2)(np.add(basis[i], basis[j]))
                # store potential solution
                sols.append(sol)
                for k in range(n):
                    # check the determinants (modulo 2)
                    det = (sol[k] * sol[k + 3 * n] + sol[k + n] * sol[k + 2 * n]) % 2
                    if det != 1:
                        # constraints not satisfied
                        sat = False
                        # remove the false solution
                        sols.pop()
                        break
                # if solution found, set equivalent to True
                if sat:
                    # equivalence is satisfied
                    equivalent = sat
                    break
            # end search loop
            if equivalent:
                break
        # if no solution found, return (False, None)
        if not equivalent:
            return None
        # if solution found, return (True, clifford_output)
        # get clifford solution in vector form
        solution_vector = np.array(sols[0]).astype(int)
        return solution_vector
        
    def __clifford_vec_to_tensors(self, vec):
        """Converts a local clifford on n qubits in vector form to a list of n single qubit
        cliffords given by 2x2 numpy arrays.
        Args:
            vec (numpy.array): A vector of size 4n specifying a local clifford on n qubits
        Returns:
            clifford_output (list): A list of length n of 2x2 numpy arrays representing single
            qubit local cliffords
        """
        n = int(len(vec)/4)
        # initialize empty list
        local_clifford_list = []
        for i in range(n):
            single_qubit_clifford = np.array(
                [[vec[i], vec[i + n]], [vec[i + 2 * n], vec[i + 3 * n]]]
            )
            local_clifford_list.append(single_qubit_clifford)
        return local_clifford_list

    # convert solution to single global clifford given by 2nx2n numpy array
    def __clifford_vec_to_global(self, vec):
        """Converts local clifford on n qubits in vector form to a local clifford on all n qubits.
        Args:
            vec (numpy.array): A vector specifying a local clifford
        Returns:
            clifford_output (numpy.array): A single 2n x 2n numpy array representing a local
            clifford
        """
        n = int(len(vec)/4)
        blocks = [np.diag([vec[i + k * n] for i in range(n)]) for k in range(4)]
        return np.block([[blocks[0], blocks[1]], [blocks[2], blocks[3]]])        

