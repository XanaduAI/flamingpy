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
"""Continuous-variable operations, states, and noise models."""

# pylint: disable=import-outside-toplevel,too-many-instance-attributes

import itertools as it

import numpy as np
from scipy.linalg import block_diag
import scipy.sparse as sp
from thewalrus.symplectic import beam_splitter, expand, rotation


def dumbbell_mat(adj, sparse=True):
    """Return a symplectic matrix for entangling sensor states.

    Give the 2N by 2N symplectic matrix for phase shifts and beamsplitter
    application that would entangle GKP sensor states. Assume that such
    an entangling link is applied between every two modes connected by an edgem
    according to the adjacency matrix adj. Assumes quadrature-like convention:

        (q_1, ..., q_N, p_1, ..., p_N).

    Args:
        adj (array): N by N binary symmetric matrix. If modes i and j are
            linked by an edge, then entry ij and ji is equal to the weight of the
            edge (1 by default); otherwise 0.
        sparse (bool): whether to return a sparse or dense array when adj
            input is a sparse array.
    Returns:
        np.array or sp.sparse.csr_matrix: 2N by 2N symplectic matrix.
            sparse if the adjacency matrix is sparse.
    """
    # Number of modes
    N = adj.shape[0]
    # 50/50 beamsplitter in the 'all q's first' convention.
    first_rot, second_rot = expand(rotation(np.pi / 4), 0, 2), expand(rotation(-np.pi / 4), 0, 2)
    entangling_circuit = first_rot @ beam_splitter(-np.pi / 4, 0) @ second_rot
    expanded_circuit = sp.identity(2 * N, format="dia")
    for (i, j) in it.combinations(range(N), 2):
        if adj[i, j]:
            expanded_circuit = sp.dia_array(expand(entangling_circuit, [i, j], N)).dot(
                expanded_circuit
            )

    if not sparse and isinstance(expanded_circuit, sp.spmatrix):
        return expanded_circuit.toarray()

    return expanded_circuit


def invert_permutation(p):
    """Invert the permutation associated with p."""
    p_inverted = np.empty(p.size, p.dtype)
    p_inverted[p] = np.arange(p.size)
    return p_inverted


def SCZ_mat(adj, sparse=True):
    """Return a symplectic matrix corresponding to CZ gate application.

    Give the 2N by 2N symplectic matrix for CZ gate application based on the
    adjacency matrix adj. Assumes quadrature-like convention:

        (q_1, ..., q_N, p_1, ..., p_N).

    Args:
        adj (array): N by N binary symmetric matrix. If modes i and j are
            linked by a CZ, then entry ij and ji is equal to the weight of the
            edge (1 by default); otherwise 0.
        sparse (bool): whether to return a sparse or dense array when adj
            input is a sparse array.
    Returns:
        np.array or sp.sparse.csr_matrix: 2N by 2N symplectic matrix.
            sparse if the adjacency matrix is sparse.
    """
    # Number of modes
    N = adj.shape[0]
    if isinstance(adj, np.ndarray):
        identity = np.eye(N, dtype=np.int8)
        zeros = np.zeros((N, N), dtype=np.int8)
        block_func = np.block
    else:
        identity = sp.identity(N, dtype=np.int8)
        zeros = sp.csr_matrix((N, N), dtype=np.int8)
        block_func = sp.bmat
    # Construct symplectic
    symplectic = block_func([[identity, zeros], [adj, identity]])

    if not sparse and isinstance(symplectic, sp.spmatrix):
        return symplectic.toarray()

    return symplectic


def SCZ_apply(adj, quads, one_shot=True):
    """Apply SCZ matrix to one- or two-dimensional array quads.

    If one-shot is True, use SCZ_mat to apply a symplectic CZ matrix to
    a matrix or vector of quadratures. Otherwise, take advantage of the
    block structure of a symplectic SCZ matrix for a more memory-
    efficient matrix multiplication.
    """
    N = quads.shape[0] // 2
    if len(quads.shape) == 1:
        if one_shot:
            new_quads = SCZ_mat(adj).dot(quads)
        else:
            old_qs = quads[:N]
            old_ps = quads[N:]
            new_quads = np.empty(2 * N, quads.dtype)
            new_quads[:N] = old_qs
            new_quads[N:] = adj.dot(old_qs) + old_ps
    if len(quads.shape) == 2:
        if one_shot:
            SCZ = SCZ_mat(adj)
            new_quads = SCZ.dot(SCZ.dot(quads).T).T
        else:
            c1, c2, c3, c4 = quads[:N, :N], quads[:N, N:], quads[N:, :N], quads[N:, N:]
            block2 = (adj.dot(c1.T)).T + c2
            block3 = adj.dot(c1) + c3
            block4 = c4 + adj.dot(c2) + (adj.dot(c3.T)).T + adj.dot(adj.dot(c1).T).T
            new_quads = np.block([[c1, block2], [block3, block4]])
    return new_quads


def splitter_symp(n=4):
    """Return the symplectic matrix of a four-splitter.

    Return the symplectic matrix of the beamsplitters connecting the four
    micronodes in each macronode. `n` refers to the total number of modes
    (so n >= 4). If n = 4, return the matrix in the 'all q's first' convention;
    otherwise, return a large block-diagonal matrix in the 'q1p1, ..., qnpn'
    convention.

    Args:
        n (int, optional): the total number of modes on which the beamsplitters
            apply (n must be >= 4).

    Returns:
        numpy.array: the sympletic matrix of the four-splitter.
    """
    # 50/50 beamsplitter in the 'all q's first' convention.
    bs5050 = beam_splitter(np.pi / 4, 0)
    bs1 = expand(bs5050, [1, 0], 4)
    bs2 = expand(bs5050, [3, 2], 4)
    bs3 = expand(bs5050, [2, 0], 4)
    bs4 = expand(bs5050, [3, 1], 4)
    bs_network = (bs4 @ bs3 @ bs2 @ bs1).astype(np.single)
    if n == 4:
        return bs_network
    if n > 4:
        # Permutation away from 'all q's first' convention for matrices of
        # with dimension 4 and the network spanning all the macronoes.
        perm_out_4 = [0, 4, 1, 5, 2, 6, 3, 7]
        bs_perm = bs_network[:, perm_out_4][perm_out_4, :]
        # Symplectic corresponding to the beasmplitter network spanning
        # the whole lattice.
        bs_full = block_diag(*[bs_perm] * (n // 4))
        return bs_full
    else:
        print("Total number of modes cannot be less than 4.")
        raise Exception
