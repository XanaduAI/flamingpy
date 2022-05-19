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

# pylint: disable=protected-access

import numpy as np
from scipy.linalg import block_diag

from flamingpy.cv.ops import CVLayer, SCZ_apply
from flamingpy.cv.gkp import GKP_binner, Z_err_cond
from thewalrus.symplectic import expand, beam_splitter


def invert_permutation(p):
    """Invert the permutation associated with p."""
    p_inverted = np.empty(p.size, p.dtype)
    p_inverted[p] = np.arange(p.size)
    return p_inverted


def symp_splitter(n):
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


def _apply_initial_noise(macro_graph, noise_layer, noise_model):
    """Apply initial noise to RHG_macro."""
    to_points = macro_graph.to_points
    N = len(macro_graph)
    # Apply a CV noise layer to the macronode lattice.
    noisy_macro_state = noise_layer(macro_graph, p_swap=swap_prob)
    # Identify perfect nodes, if they exist
    perfect_points = macro_graph.graph.get("perfect_points")
    if perfect_points:
        perfect_inds = [macro_graph.to_indices[point] for point in perfect_points]
    else:
        perfect_inds = None
    noise_model["perfect_inds"] = perfect_inds
    apply_noise(noise_model)
    return noisy_macro_state


def _permute_indices(code, macro_graph):
    """Obtain permuted indices and set type of reduced node."""
    # A list of permuted indices where each block of four
    # corresponds to [star, planet, planet, planet].
    reduced_graph = code.graph
    reduced_graph._states = []
    permuted_inds = np.empty(N, dtype=np.int32)
    for i in range(0, N - 3, 4):
        # Indices of GKP micronodes in macronode i.
        gkps = []
        for j in range(4):
            micronode = to_points[i + j]
            if macro_graph.nodes[micronode]["state"] == "GKP":
                gkps.append(j)
        centre_point = tuple(round(i) for i in micronode)
        if gkps:
            star_ind, reduced_state = i + gkps[0], "GKP"
        else:
            star_ind, reduced_state = i, "p"
            # RHG_reduced._states += [RHG_reduced.to_indices[centre_point]]
        # Set type of node in the reduced lattice as a p-squeezed
        # state if all micronodes are p, else GKP.
        reduced_graph.nodes[centre_point]["state"] = reduced_state
        # Old and permuted indices of all micronodes in macronode i.
        old_inds = [i, i + 1, i + 2, i + 3]
        old_inds.pop(star_ind - i)
        new_inds = [star_ind] + old_inds
        # Associate a 'body index' (1 to 4) to each micronode,
        # with 1 being the star index and the rest being planets.
        k = 1
        for ind in new_inds:
            reduced_graph.nodes[to_points[ind]]["body_index"] = k
            k += 1
        permuted_inds[[i, i + 1, i + 2, i + 3]] = new_inds

    return permuted_inds


def _entangle_states(macro_graph, initial_quads, symp_bs):
    """Entangle the states in the lattice."""
    quads = SCZ_apply(macro_graph.adj_mat, initial_quads)
    # Permute the quadrature values to align with the permuted
    # indices in order to apply the beamsplitter network.
    quad_permutation = np.concatenate([permuted_inds, N + permuted_inds])
    permuted_quads = quads[quad_permutation]
    for i in range(0, N - 3, 4):
        q_inds = np.array([i, i + 1, i + 2, i + 3])
        p_inds = q_inds + N
        updated_qs = symp_bs[:4, :4] @ permuted_quads[q_inds]
        updated_ps = symp_bs[4:, 4:] @ permuted_quads[p_inds]
        permuted_quads[q_inds] = updated_qs
        permuted_quads[p_inds] = updated_ps
    return permuted_quads, quad_permutation


def _measure_syndrome(noisy_macro_state, permuted_inds, entangled_quads, quad_permutation):
    """Measure the syndrome of noisy_macro_state."""
    unpermuted_quads = entangled_quads[invert_permutation(quad_permutation)]    
    # Indices of stars and planets.
    stars = permuted_inds[::4]
    planets = np.delete(permuted_inds, np.arange(0, N, 4))
    # Update quadrature values after CZ gate application.
    # Measure stars in p, planets in q.
    noisy_macro_state.measure_hom(quad="p", inds=stars, updated_quads=unpermuted_quads)
    noisy_macro_state.measure_hom(quad="q", inds=planets, updated_quads=unpermuted_quads)
    

def _neighbor_of_micro_i_macro_j(i, j, macro_graph, permuted_inds):
    """Return the neighbor of the ith micronode, jth macronode in macro_lattice.

    Micronode i is adjacent to a neighbor with a body index (1 for
    star, 2, 3, 4 for planets). Return the vertex and the body index
    of the neighbor to help the processing rules. If there is no
    such neighbor, return None.
    """
    to_points = macro_graph.to_points
    # Index of ith micronode in the jth macronode.
    ith_index = permuted_inds[j + i - 1]
    ith_vertex = to_points[ith_index]
    # Vertex of the neighbor of the ith micronode.
    ith_adjacency = list(macro_graph[ith_vertex])
    if ith_adjacency:
        ith_neighbor = list(macro_graph[ith_vertex])[0]
        ith_body_index = macro_graph[ith_neighbor]["body_index"]
        return ith_neighbor, ith_body_index
    return None


 def _hom_outcomes(vertex, macro_graph):
        """Measurement outcomes in the macronode containing vertex.

        Return the values of the homodyne measurements of the macronode
        containing vertex. Note we are only interested in q-homodyne
        outcomes; the returned list is of the form [0, 0, q2, q3, q4].
        If vertex is None, return a list of 0s, so that the processing
        is unaltered by the outcomes.
        """
        if vertex is None:
            return [0, 0, 0, 0, 0]
        meas = np.zeros(5)
        # The central node corresponding to the neighboring
        # macronode.
        central_node = tuple(round(i) for i in vertex)
        for micro in macro_graph.macro_to_micro[central_node]:
            index = macro_graph.nodes[micro]["body_index"]
            # Populate meas with the q-homodyne outcomes for
            # the planet modes.
            if index != 1:
                meas[index] = macro_graph.nodes[micro]["hom_val_q"]
        return meas


def _process_outcomes(hom_outcomes, neighbor_body_index):
    """Process the homodyne outcomes for neighboring macronode i.
    
    Macronode j is connected to a macronode whose with an array of
    measurement outcomes M and whose micronodes have body indices
    neighbor_body_index. Use this information to process the
    measurement outcomes M.
    """
    if neighbor_body_index == 1:
        return 0
    if neighbor_body_index == 2:
        return hom_outcomes[2] - hom_outcomes[4]
    if neighbor_body_index == 3:
        return hom_outcomes[3] - hom_outcomes[4]
    if neighbor_body_index == 4:
        return hom_outcomes[2] + hom_outcomes[3]
    return None


def _reduce_jth_macronode(j, macro_graph, permuted_inds):
    """Reduce macronode with index j."""
    to_points = macro_graph.to_points
    star_index = permuted_inds[j]
    vertex = to_points[star_index]

    # Here, j corresponds to the macronode and i to to micronode.
    # i ranges from 1 to 4, to align with manuscript.
    verts_and_inds = [_neighbor_of_micro_i_macro_j(i, j, macro_graph, permuted_inds) for i in (1, 2, 3, 4)]
    neighbors = [tup[0] if tup else None for tup in verts_and_inds]
    body_indices = [tup[1] if tup else None for tup in verts_and_inds]

    # Array of arrays of measurement outcomes in all the
    # macronodes adjacent to j.
    m_arr = np.array([_hom_outcomes(neighbors[i - 1], macro_graph) for i in (1, 2, 3, 4)])
    # Array of processed q-homodyne outcomes from neighboring
    # macronodes of the form [0, Z(1), Z(2), Z(3), Z(4)].
    Z_arr = np.array([0] + [_process_outcomes(m_arr[i - 1], body_indices[i - 1]) for i in (1, 2, 3, 4)])
    # p-homodyne outcome of the star node.
    star_p_val = macro_graph.nodes[vertex]["hom_val_p"]

    # Types of state for the four micronodes directly neighboring
    # macronode j.
    types = [macro_graph.nodes[neighbor]["state"] if neighbor else None for neighbor in neighbors]
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
    p_err = min(p_err, 0.5)
    p_err = max(p_err, 0)

    bitp = GKP_binner([outcome])[0]
    bitq = GKP_binner(Z_arr[gkp_inds].astype(np.float64)) if gkp_inds else 0

    processed_bit_val = (bitp + np.sum(bitq)) % 2

    # Update the reduced RHG lattice with the effective
    # homodyne value and the phase error probability.
    central_vert = tuple(round(i) for i in vertex)
    RHG_reduced.nodes[central_vert]["bit_val"] = processed_bit_val
    sorted_bits[reduced_indices[central_vert]] = processed_bit_val
    RHG_reduced.nodes[central_vert]["p_phase_cond"] = p_err
    

def reduce_macronode_graph(code, macro_graph, bs_network, noise_layer=None, noise_model=None):
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
    if noise_layer is None:
        noise_layer = CVLayer
    if noise_model is None:
        noise_model = {
        "noise": "grn",
        "delta": delta,
        "sampling_order": "two-step",
    }
    # Sample for the initial state
    noisy_macro_state = _apply_initial_noise(macro_graph, noise_layer, noise_model)
    initial_quads = noisy_macro_state._init_quads
    # Entangle macronodes
    permuted_inds = _permute_indices(code, macro_graph)
    entangled_quads, quad_permutation = _entangle_states(macro_graph, initial_quads, bs_network)
    # Measure syndrome
    _measure_syndrome(noisy_macro_state, permuted_inds, entangled_quads, quad_permutation)
    # Process homodyne outcomes and calculate of phase error probabilities
    for j in range(0, N - 3, 4):
        _reduce_jth_macronode(j, macro_graph, permuted_inds)
