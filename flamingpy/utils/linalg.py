"""Helper functions for linear algebra used for testing LC equivalence."""

import numpy as np


def reduce_RREform_mod2(M, max_cols=None):
    """Put a binary matrix into row reduced echelon form modulo 2, up to a
    maximum number of columns given by max_cols.

    Args:
        M (numpy.array): the matrix to reduce.
        max_cols (int): the maximum number of columns of the input array to reduce.
    Returns:
        (numpy.array, int): the reduced matrix and the the index of pivot column.
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


def lc_constraint_system(G, H):
    """Build the constraint for LC-equivalence of two adjacency matrices.

    Construct a binary system of equations that two adjacency matrices
    G and H must satisfy for equivalence through local complementations.

    Args:
        G, H (numpy.array): n x n adjacency matrices.
    Returns:
        numpy.array: an array of shape n^2 x 4n representing a system of n^2
            binary linear equations with 4n unknowns.
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


def nullspace_basis(M):
    """Return the nullspace basis of matrix M.

    Construct an array whose rows are binary basis vectors of the right
    nullspace of input matrix array.

    Args:
        M (numpy.array): a binary matrix.
    Returns:
        numpy.array: an array whose rows are basis vectors of the right
            nullspace of M.
    """
    # find left Null space of transposed system, which is equal to right null space
    M_transposed = M.T
    m_rows, n_cols = M_transposed.shape
    # construct augmented block matrix A = [M | I]
    A = np.concatenate((M_transposed, np.eye(m_rows, dtype=int)), axis=-1)
    # row reduce left M block of augmented matrix
    R, p = reduce_RREform_mod2(A, max_cols=n_cols)
    N = R[p:, n_cols:]
    basis, _ = reduce_RREform_mod2(N)
    return basis


def search_nullspace(basis):
    """Check if the nullspace satisfies the determinant constraints.

     Search through sums of pairs of basis vectors of the nullspace and
     see if any satisfy the determinant constraints. If no solution is found,
     output None. If a solution is found, return a vector specifying the local
     Clifford in the form:

        (a_1, a_2, ..., a_n, b_1, b_2, ..., b_n, c_1, c_2, ..., c_n, d_1, d_2, ..., d_n),

     where n is the number of nodes of the graph.

    Args:
        basis (numpy.array): an array whose rows are basis vectors of the nullspace
    Returns:
        NoneType or numpy.array: None if no solution, or the solution array.
    """
    # number of basis vectors of null space
    d = basis.shape[0]
    # number of nodes/qubits
    n = int(basis.shape[1] / 4)
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
    # if no solution found, return None
    if not equivalent:
        return None
    # if solution found, return clifford in vector form
    solution_vector = np.array(sols[0]).astype(int)
    return solution_vector


def clifford_vec_to_tensors(vec):
    """Convert a local Clifford gate on n qubits to a list of n single-qubit
    Cliffords.

    Take a local Clifford operation on n qubits in the form

    (a_1, a_2, ..., a_n, b_1, b_2, ..., b_n, c_1, c_2, ..., c_n, d_1, d_2, ..., d_n),

    and return a list of single-qubit Clifford gates as a list of 2x2 arrays, where the array
    at index k corresponds to the Clifford acting on qubit k given as:

        [a_k, b_k]
        [c_k, d_k].

    Args:
        vec (numpy.array): a vector of size 4n specifying a local Clifford gate on n qubits.

    Returns:
        list[numpy.array]: n 2x2 arrays representing single-qubit local Clifford gates.
    """
    # number of nodes/qubits
    n = int(len(vec) / 4)
    # initialize empty list
    local_clifford_list = []
    for i in range(n):
        single_qubit_clifford = np.array([[vec[i], vec[i + n]], [vec[i + 2 * n], vec[i + 3 * n]]])
        local_clifford_list.append(single_qubit_clifford)
    return local_clifford_list


def clifford_vec_to_global(vec):
    """Convert a local Clifford gate on n qubits to a local Clifford gate on
    all n qubits.

    Take a vector corresponding to a local Clifford gate on n qubits,

    (a_1, a_2, ..., a_n, b_1, b_2, ..., b_n, c_1, c_2, ..., c_n, d_1, d_2, ..., d_n),

    and return a single 2n x 2n array representing a local Clifford acting on
    all n qubits, given as the block matrix:

            [A, B]
            [C, D]

        where each block is a n x n diagonal matrix:

            A = diag(a_1, a_2, ..., a_n)
            B = diag(b_1, b_2, ..., b_n)
            C = diag(c_1, c_2, ..., c_n)
            D = diag(d_1, d_2, ..., d_n).

    Args:
        vec (numpy.array): a vector of size 4n specifying a local Clifford on n qubits.
    Returns:
        numpy.array:  a 2n x 2n numpy array representing a local
            Clifford acting on all n qubits.
    """
    # number of nodes/qubits
    n = int(len(vec) / 4)
    blocks = [np.diag([vec[i + k * n] for i in range(n)]) for k in range(4)]
    return np.block([[blocks[0], blocks[1]], [blocks[2], blocks[3]]])


def are_lc_equivalent(graph1, graph2, clifford_form="tensor"):
    """Check if two EGraphs are LC equivalent, and return the Clifford
    operation if so. Implemented as in arXiv:quant-ph/0405023.

    Args:
        graph1 (EGraph): the initial graph to check Clifford equivalence with.
        graph2 (EGraph): the graph to check Clifford equivalence against.
        clifford_form: a string describing the output form of local Clifford operation, if
            it exists.

            If 'tensor' (default), produce a list of length n of 2x2 numpy arrays corresponding to
                single-qubit tensor factors.
            If 'global', return a single 2nx2n numpy array corresponding to the global operator
                acting on all n qubits.

    Returns:
        (bool, numpy.array): whether the states are LC equivalent, and if they are, the local
            Clifford output according to 'clifford_form' specification.
    """
    # get adjacency matrices of input graphs
    graph1.adj_generator(sparse=False)
    G = graph1.adj_mat
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

    # perform algorithm to search for solution
    #
    # construct system of constraints for two adjacency matrices G and H
    system = lc_constraint_system(G, H)
    # construct nullspace basis of system of constraints
    nullspace = nullspace_basis(system)
    # search nullspace for solution vector
    solution_vector = search_nullspace(nullspace)

    # if graphs are not equivalent set clifford_output to None
    if solution_vector is None:
        equivalent = False
        clifford_output = None
    # if graphs are equivalent and clifford_form == "tensor", convert clifford_output to
    # tensor factors
    elif clifford_form == "tensor":
        equivalent = True
        clifford_output = clifford_vec_to_tensors(solution_vector)
    # if graphs are equivalent and clifford_form == "global", convert clifford_output to
    # global form
    elif clifford_form == "global":
        equivalent = True
        clifford_output = clifford_vec_to_global(solution_vector)
    # return equivalence and clifford_output
    return equivalent, clifford_output
