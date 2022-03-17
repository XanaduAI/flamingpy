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
"""Functions for reducing a macronode lattice to a canonical lattice."""
import numpy as np
from scipy.linalg import block_diag

from flamingpy.cv.ops import CVLayer, SCZ_apply
from flamingpy.cv.gkp import GKP_binner, Z_err_cond


def expand(S, modes, N):
    r"""Expands a Symplectic matrix S to act on the entire subsystem.
    Args:
        S (array): a :math:`2M\times 2M` Symplectic matrix
        modes (Sequence[int]): the list of modes S acts on
        N (int): full size of the subsystem
    Returns:
        array: the resulting :math:`2N\times 2N` Symplectic matrix
    Note:
        This is from thewalrus library.
        https://github.com/XanaduAI/thewalrus/blob/master/thewalrus/symplectic.py
    """
    M = len(S) // 2
    S2 = np.identity(2 * N, dtype=S.dtype)
    w = np.array(modes)

    S2[w.reshape(-1, 1), w.reshape(1, -1)] = S[:M, :M].copy()  # X
    S2[(w + N).reshape(-1, 1), (w + N).reshape(1, -1)] = S[M:, M:].copy()  # P
    S2[w.reshape(-1, 1), (w + N).reshape(1, -1)] = S[:M, M:].copy()  # XP
    S2[(w + N).reshape(-1, 1), w.reshape(1, -1)] = S[M:, :M].copy()  # PX

    return S2


def beam_splitter(theta, phi, dtype=np.float64):
    """Beam-splitter.

    Args:
        theta (float): transmissivity parameter
        phi (float): phase parameter
        dtype (numpy.dtype): datatype to represent the Symplectic matrix
    Returns:
        array: symplectic-orthogonal transformation matrix of an interferometer with angles theta and phi
    Note:
        This is from thewalrus library.
        https://github.com/XanaduAI/thewalrus/blob/master/thewalrus/symplectic.py
    """
    ct = np.cos(theta, dtype=dtype)
    st = np.sin(theta, dtype=dtype)
    eip = np.cos(phi, dtype=dtype) + 1j * np.sin(phi, dtype=dtype)
    U = np.array(
        [
            [ct, -eip.conj() * st],
            [eip * st, ct],
        ]
    )
    return interferometer(U)


def interferometer(U):
    """Interferometer.

    Args:
        U (array): unitary matrix
    Returns:
        array: symplectic transformation matrix
    Note:
        This is from thewalrus library.
        https://github.com/XanaduAI/thewalrus/blob/master/thewalrus/symplectic.py
    """
    X = U.real
    Y = U.imag
    S = np.block([[X, -Y], [Y, X]])

    return S


def invert_permutation(p):
    """Invert the permutation associated with p."""
    p_inverted = np.empty(p.size, p.dtype)
    p_inverted[p] = np.arange(p.size)
    return p_inverted


def BS_network(n):
    """Return the symlectic matrix of the beamsplitter network.

    Return the symplectic matrix of the beamsplitters connecting four
    micronodes in each macronode out of n total micronodes. If n = 4, return
    the matrix in the 'all q's first' convention; otherwise, return a large
    block-diagonal matrix in the 'q1p1, ... qnpn' convention.
    """
    # 50/50 beamsplitter in the 'all q's first' convention.
    bs5050 = beam_splitter(np.pi / 4, 0)
    bs1 = expand(bs5050, [1, 0], 4)
    bs2 = expand(bs5050, [3, 2], 4)
    bs3 = expand(bs5050, [2, 0], 4)
    bs4 = expand(bs5050, [3, 1], 4)
    # TODO: Data type set to 'single' because there are only 0, +-0.5 entries
    # but this is really 1/2 of an array that can have dtype=np.int8,
    # so revisit this.
    bs_network = (bs4 @ bs3 @ bs2 @ bs1).astype(np.single)
    if n < 4:
        print("Too small!")
        raise Exception
    if n > 4:
        # Permutation away from 'all q's first' convention for matrices of
        # with dimension 4 and the network spanning all the macronoes.
        perm_out_4 = [0, 4, 1, 5, 2, 6, 3, 7]
        bs_perm = bs_network[:, perm_out_4][perm_out_4, :]
        # Symplectic corresponding to the beasmplitter network spanning
        # the whole lattice.
        bs_full = block_diag(*[bs_perm] * (n // 4))
        return bs_full
    return bs_network


def reduce_macro_and_simulate(RHG_macro, RHG_reduced, CVRHG_reduced, bs_network, swap_prob, delta):
    """Reduce the macronode RHG lattice to the canonical lattice.

    Take the macronode lattice EGraph, RHG_macro, and generate a
    macronode CV lattice with swap-out probaiblity swap_prob and delta
    value delta. Then, label micronodes as planets and stars, conduct
    homodyne measurements, process these measurements, and compute
    conditional phase error probabilities. Generate an canonical RHG
    lattice with effective measurement outcomes, phase error
    probabilities, and state types ('p' or 'GKP') stored as node
    attributes.
    """
    to_points = RHG_macro.to_points
    N = len(RHG_macro)

    # The hybridized CVRHG macronode lattice.
    CVRHG = CVLayer(RHG_macro, p_swap=swap_prob)
    # Noise-model
    perfect_points = RHG_macro.graph.get("perfect_points")
    if perfect_points:
        perfect_inds = [RHG_macro.to_indices[point] for point in perfect_points]
    else:
        perfect_inds = None
    noise_model = {
        "noise": "grn",
        "delta": delta,
        "sampling_order": "two-step",
        "perfect_inds": perfect_inds,
    }
    CVRHG.apply_noise(noise_model)
    # A list of permuted indices where each block of four
    # corresponds to [star, planet, planet, planet].
    permuted_inds = np.empty(N, dtype=np.int32)
    for i in range(0, N - 3, 4):
        # Indices of GKP micronodes in macronode i.
        gkps = []
        for j in range(4):
            micronode = to_points[i + j]
            if CVRHG.egraph.nodes[micronode]["state"] == "GKP":
                gkps.append(j)
        centre_point = tuple([round(i) for i in micronode])
        if gkps:
            star_ind, reduced_state = i + gkps[0], "GKP"
        else:
            star_ind, reduced_state = i, "p"
            CVRHG_reduced._states["p"] += [RHG_reduced.to_indices[centre_point]]
        # Set type of node in the reduced lattice as a p-squeezed
        # state if all micronodes are p, else GKP.
        RHG_reduced.nodes[centre_point]["state"] = reduced_state
        # Old and permuted indices of all micronodes in macronode i.
        old_inds = [i, i + 1, i + 2, i + 3]
        old_inds.pop(star_ind - i)
        new_inds = [star_ind] + old_inds
        # Associate a 'body index' (1 to 4) to each micronode,
        # with 1 being the star index and the rest being planets.
        k = 1
        for ind in new_inds:
            CVRHG.egraph.nodes[to_points[ind]]["body_index"] = k
            k += 1
        permuted_inds[[i, i + 1, i + 2, i + 3]] = new_inds

    # Indices of stars and planets.
    stars = permuted_inds[::4]
    planets = np.delete(permuted_inds, np.arange(0, N, 4))

    # Update quadrature values after CZ gate application.
    quads = CVRHG._init_quads
    quads = SCZ_apply(RHG_macro.adj_mat, quads)

    # Permute the quadrature values to align with the permuted
    # indices in order to apply the beamsplitter network.
    quad_permutation = np.concatenate([permuted_inds, N + permuted_inds])
    permuted_quads = quads[quad_permutation]
    for i in range(0, N - 3, 4):
        q_inds = np.array([i, i + 1, i + 2, i + 3])
        p_inds = q_inds + N
        updated_qs = bs_network[:4, :4] @ permuted_quads[q_inds]
        updated_ps = bs_network[4:, 4:] @ permuted_quads[p_inds]
        permuted_quads[q_inds] = updated_qs
        permuted_quads[p_inds] = updated_ps

    unpermuted_quads = permuted_quads[invert_permutation(quad_permutation)]
    # Measure stars in p, planets in q.
    CVRHG.measure_hom(quad="p", inds=stars, updated_quads=unpermuted_quads)
    CVRHG.measure_hom(quad="q", inds=planets, updated_quads=unpermuted_quads)

    def neighbor_of_i(i, j):
        """Return the neighbor of the ith micronode, jth macronode.

        Micronode i is adjacent to a neighbor with a body index (1 for
        star, 2, 3, 4 for planets). Return the vertex and the body index
        of the neighbor to help the processing rules. If there is no
        such neighbor, return None.
        """
        # Index of ith micronode in the jth macronode.
        ith_index = permuted_inds[j + i - 1]
        ith_vertex = to_points[ith_index]
        # Vertex of the neighbor of the ith micronode.
        ith_adjacency = list(CVRHG.egraph[ith_vertex])
        if ith_adjacency:
            ith_neighbor = list(CVRHG.egraph[ith_vertex])[0]
            ith_body_index = CVRHG.egraph.nodes[ith_neighbor]["body_index"]
            return ith_neighbor, ith_body_index
        else:
            return

    def m(vertex):
        """Measurement outcomes in the macronode containing vertex.

        Return the values of the homodyne measurements of the macronode
        containing vertex. Note we are only interested in q-homodyne
        outcomes; the returned list is of the form [0, 0, q2, q3, q4].
        If vertex is None, return a list of 0s, so that the processing
        is unaltered by the outcomes.
        """
        if vertex is None:
            return [0, 0, 0, 0, 0]
        else:
            meas = np.zeros(5)
            # The central node corresponding to the neighboring
            # macronode.
            central_node = tuple([round(i) for i in vertex])
            for micro in CVRHG.egraph.macro_to_micro[central_node]:
                index = CVRHG.egraph.nodes[micro]["body_index"]
                # Populate meas with the q-homodyne outcomes for
                # the planet modes.
                if index != 1:
                    meas[index] = CVRHG.egraph.nodes[micro]["hom_val_q"]
            return meas

    def Z(M, neighbor_body_index):
        """Process the homodyne outcomes for neighboring macronode i.

        Macronode j is connected to a macronode whose with an array of
        measurement outcomes M and whose micronodes have body indices
        neighbor_body_index. Use this information to process the
        measurement outcomes M.
        """
        if neighbor_body_index == 1:
            return 0
        if neighbor_body_index == 2:
            return M[2] - M[4]
        if neighbor_body_index == 3:
            return M[3] - M[4]
        if neighbor_body_index == 4:
            return M[2] + M[3]

    # sorted_homodynes = np.empty(N // 4, dtype=np.float32)
    sorted_bits = np.empty(N // 4, dtype=np.float32)
    reduced_indices = RHG_reduced.to_indices
    # Processing of homodyne outcomes and calculations of phase
    # error probabilities
    for j in range(0, N - 3, 4):
        star_index = permuted_inds[j]
        vertex = to_points[star_index]

        # Here, j corresponds to the macronode and i to to micronode.
        # i ranges from 1 to 4, to align with manuscript.

        verts_and_inds = [neighbor_of_i(i, j) for i in (1, 2, 3, 4)]
        neighbors = [tup[0] if tup else None for tup in verts_and_inds]
        body_indices = [tup[1] if tup else None for tup in verts_and_inds]

        # Array of arrays of measurement outcomes in all the
        # macronodes adjacent to j.
        m_arr = np.array([m(neighbors[i - 1]) for i in (1, 2, 3, 4)])
        # Array of processed q-homodyne outcomes from neighboring
        # macronodes of the form [0, Z(1), Z(2), Z(3), Z(4)].
        Z_arr = np.array([0] + [Z(m_arr[i - 1], body_indices[i - 1]) for i in (1, 2, 3, 4)])
        # p-homodyne outcome of the star node.
        star_p_val = CVRHG.egraph.nodes[vertex]["hom_val_p"]

        # Types of state for the four micronodes directly neighboring
        # macronode j.
        types = [RHG_macro.nodes[neighbor]["state"] if neighbor else None for neighbor in neighbors]
        # Phase error probability and number of p-squeezed states
        # among the four micronodes in the vicinity of macronode j
        p_err = 0
        num_p = 0
        outcome = 2 * star_p_val
        gkp_inds = []
        for i in (1, 2, 3, 4):
            if types[i - 1] == "p":
                num_p += 1
                outcome -= Z_arr[i]
            if types[i - 1] == "GKP":
                if delta > 0:
                    p_err += Z_err_cond(2 * delta, Z_arr[i], use_hom_val=True)
                gkp_inds += [i]
        if delta > 0:
            p_err += Z_err_cond(2 * (2 + num_p) * delta, outcome, use_hom_val=True)
        if p_err > 0.5:
            p_err = 0.5
        if p_err < 0:
            p_err = 0

        bitp = GKP_binner([outcome])[0]
        bitq = GKP_binner(Z_arr[gkp_inds].astype(np.float64)) if gkp_inds else 0

        processed_bit_val = (bitp + np.sum(bitq)) % 2

        # Update the reduced CVRHG lattice with the effective
        # homodyne value and the phase error probability.
        central_vert = tuple([round(i) for i in vertex])
        RHG_reduced.nodes[central_vert]["bit_val"] = processed_bit_val
        sorted_bits[reduced_indices[central_vert]] = processed_bit_val
        RHG_reduced.nodes[central_vert]["p_phase_cond"] = p_err
    CVRHG_reduced.bits = sorted_bits
    return
