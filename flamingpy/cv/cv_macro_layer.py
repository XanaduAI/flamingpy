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


class CVMacroLayer(CVLayer):
    """The CV macronode noise layer."""
    def __init__(self, *args, reduced_graph, **kwargs):
        super().__init__(*args, **kwargs)
        self.reduced_graph = reduced_graph

    def _apply_initial_noise(self, noise_model):
        """Set up the two-step noise model for macro_graph.
    
        Based on noise_layer and noise_model, populate macro_graph with states and
        sample for the initial (ideal) measurement outcomes.
    
        This function modifies macro_graph.
        """
        macro_graph = self.egraph
        perfect_points = macro_graph.graph.get("perfect_points")
        if perfect_points:
            perfect_inds = [macro_graph.to_indices[point] for point in perfect_points]
        else:
            perfect_inds = None
        noise_model["perfect_inds"] = perfect_inds
        noise_model["sampling_order"] = "two-step"
        self.apply_noise(noise_model)
    
    def _permute_indices_and_label(self):
        """Obtain permuted indices and set type of reduced node.
    
        For each macronode, permute indices so that the first encountered GKP state
        comes first, designating it as the 'star' ('central') mode. The first
        index in the resulting list and every four indices thereafter correspond
        to star modes. The rest are 'planets' ('satellite' modes).
    
        For each star, associate a 'body_index' of 1, and 2, 3, and 4 for the
        subsequent planets. Additionally, if a macronode contains at least one GKP
        state, label the reduced state as 'GKP' (otherwise 'p').
    
        This function returns a list and modified reduced_graph.
        """
        N = self._N
        macro_graph = self.egraph
        to_points = macro_graph.to_points
        # A list of permuted indices where each block of four
        # corresponds to [star, planet, planet, planet].
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
            # Set type of node in the reduced lattice as a p-squeezed
            # state if all micronodes are p, else GKP.
            self.reduced_graph.nodes[centre_point]["state"] = reduced_state
            # Old and permuted indices of all micronodes in macronode i.
            old_inds = [i, i + 1, i + 2, i + 3]
            old_inds.pop(star_ind - i)
            new_inds = [star_ind] + old_inds
            # Associate a 'body index' (1 to 4) to each micronode,
            # with 1 being the star index and the rest being planets.
            k = 1
            for ind in new_inds:
                macro_graph.nodes[to_points[ind]]["body_index"] = k
                k += 1
            permuted_inds[[i, i + 1, i + 2, i + 3]] = new_inds
        self.permuted_inds = permuted_inds
    
    
    def _entangle_states(self, symp_bs):
        """Entangle the states in the macro_graph.
    
        Apply CZ gates to (i.e. a symplectic CZ matrix to the quadratures of)
        macro_graph, based on where the edges are in the graph. Then, apply the
        four-splitter to each macronode.
    
        Return the permuted quadratures and the corresponding permutation vector.
        """
        macro_graph = self.egraph
        quads = SCZ_apply(macro_graph.adj_mat, self._init_quads)
        # Permute the quadrature values to align with the permuted
        # indices in order to apply the beamsplitter network.
        N = self._N
        quad_permutation = np.concatenate([self.permuted_inds, N + self.permuted_inds])
        permuted_quads = quads[quad_permutation]
        for i in range(0, N - 3, 4):
            q_inds = np.array([i, i + 1, i + 2, i + 3])
            p_inds = q_inds + N
            updated_qs = symp_bs[:4, :4] @ permuted_quads[q_inds]
            updated_ps = symp_bs[4:, 4:] @ permuted_quads[p_inds]
            permuted_quads[q_inds] = updated_qs
            permuted_quads[p_inds] = updated_ps
        self.permuted_quads = permuted_quads
        self.quad_permutation = quad_permutation

    
    def _measure_syndrome(self):
        """Measure the syndrome of noisy_macro_state.
    
        Conduct p-homodyne measurements on the stars (central modes) of
        noisy_macro_state and q-homodyne measurements to the planets
        (satellite modes). This effectively conducts an X-basis measurement
        on the modes of the reduced lattice.
        """
        N = self._N
        unpermuted_quads = self.permuted_quads[invert_permutation(self.quad_permutation)]
        # Indices of stars and planets.
        stars = self.permuted_inds[::4]
        planets = np.delete(self.permuted_inds, np.arange(0, N, 4))
        # Update quadrature values after CZ gate application.
        # Measure stars in p, planets in q.
        self.measure_hom(quad="p", inds=stars, updated_quads=unpermuted_quads)
        self.measure_hom(quad="q", inds=planets, updated_quads=unpermuted_quads)
    
    
    def _neighbor_of_micro_i_macro_j(self, i, j):
        """Return the neighbor of the ith micronode, jth macronode in macro_graph.
    
        Suppose micronode i (in macronode j) is adjacent to a neighbor.
        Return the vertex (tuple) and the body index of the neighbor to help
        the subsequent processing rules. If there is no such neighbor,
        return None.
        """
        macro_graph = self.egraph
        # Index of ith micronode in the jth macronode.
        ith_index = self.permuted_inds[j + i - 1]
        ith_vertex = self.to_points[ith_index]
        # Vertex of the neighbor of the ith micronode.
        ith_adjacency = list(macro_graph[ith_vertex])
        if ith_adjacency:
            ith_neighbor = list(macro_graph[ith_vertex])[0]
            ith_body_index = macro_graph.nodes[ith_neighbor]["body_index"]
            return ith_neighbor, ith_body_index
        return None
    
    
    def _hom_outcomes(self, vertex):
        """Measurement outcomes in the macronode containing vertex.
    
        Return the values of the homodyne measurements of the macronode
        containing vertex. Note we are only interested in q-homodyne
        outcomes; the returned list is of the form [0, 0, q2, q3, q4]. If
        vertex is None, return a list of 0s, so that the processing is
        unaltered by the outcomes.
        """
        macro_graph = self.egraph
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
    
    
    def _process_neighboring_outcomes(self, neighbor_hom_vals, neighbor_body_index):
        """Process homodyne outcomes for a neighboring macronode.
    
        Suppose some micronode is connected to a neighboring micronode, i. i
        has a certain body index and belongs to a macronode with measurement
        results neighbor_hom_vals. Use this information to process the
        results (i.e. find appropriate linear combinations of
        neighbor_hom_vals).
        """
        if neighbor_body_index == 1:
            return 0
        if neighbor_body_index == 2:
            return neighbor_hom_vals[2] - neighbor_hom_vals[4]
        if neighbor_body_index == 3:
            return neighbor_hom_vals[3] - neighbor_hom_vals[4]
        if neighbor_body_index == 4:
            return neighbor_hom_vals[2] + neighbor_hom_vals[3]
        return None
    
    
    def _reduce_jth_macronode(self, j):
        """Obtain the bit value and error probability for the jth macronode.
    
        Process the measurement outcomes of the jth macronode in macro_graph to
        populate the corresponding reduced node in reduced_graph with a bit value
        and conditional phase error probability.
    
        This function modifies reduced_graph.
        """
        macro_graph = self.egraph
        reuced_graph = self.reduced_graph
        delta = self._delta
        to_points = macro_graph.to_points
        star_index = self.permuted_inds[j]
        vertex = self.to_points[star_index]
    
        # Here, j corresponds to the macronode and i to to micronode.
        # i ranges from 1 to 4, to align with manuscript.
        verts_and_inds = [
            self._neighbor_of_micro_i_macro_j(i, j) for i in (1, 2, 3, 4)
        ]
        neighbors = [tup[0] if tup else None for tup in verts_and_inds]
        body_indices = [tup[1] if tup else None for tup in verts_and_inds]
    
        # Array of arrays of measurement outcomes in all the
        # macronodes adjacent to j.
        m_arr = np.array([self._hom_outcomes(neighbors[i - 1]) for i in (1, 2, 3, 4)])
        # Array of processed q-homodyne outcomes from neighboring
        # macronodes of the form [0, Z(1), Z(2), Z(3), Z(4)].
        Z_arr = np.array(
            [0]
            + [self._process_neighboring_outcomes(m_arr[i - 1], body_indices[i - 1]) for i in (1, 2, 3, 4)]
        )
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
        self.reduced_graph.nodes[central_vert]["bit_val"] = processed_bit_val
        self.reduced_graph.nodes[central_vert]["p_phase_cond"] = p_err
    
    
    def reduce(self, noise_model, symp_bs):
        """Reduce the macronode lattice macro_graph to the canonical reduced_graph.
    
        Follow the procedure in arXiv:2104.03241. Take the macronode lattice
        macro_graph, apply noise based on noise_layer and noise_model, designate
        micronodes as planets and stars, conduct homodyne measurements, process
        these measurements, and compute conditional phase error probabilities.
    
        Modify reduced_graph into a canonical lattice with effective measurement
        outcomes and phase error probabilities stored as node attributes.
    
        Args:
            macro_graph (EGraph): the macronode lattice
            reduced_graph (EGraph): the reduced lattice
            noise_layer (CVLayer): the noise layer
            noise_model (dict): the dictionary of noise parameters to be fed into
                noise_layer
            bs_network (np.array): the beamsplitter network used to entangle the
                macronodes.
    
        Returns:
            None
        """
        macro_graph = self.egraph
        # Sample for the initial state
        self._apply_initial_noise(noise_model)
        self._permute_indices_and_label()
        self._entangle_states(symp_bs)
        self._measure_syndrome()
        # Process homodyne outcomes and calculate phase error probabilities
        for j in range(0, len(macro_graph) - 3, 4):
            self._reduce_jth_macronode(j)
