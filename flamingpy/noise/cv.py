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
"""Continuous-variable noise model classes."""

# pylint: disable=import-outside-toplevel

import numpy as np
from numpy.random import default_rng

from flamingpy.codes import EGraph
from flamingpy.cv.ops import invert_permutation, SCZ_mat, SCZ_apply, splitter_symp
from flamingpy.cv.gkp import GKP_binner, Z_err_cond


class CVLayer:
    def __init__(self, code, state_sampler, homodyme_sampler, translator):
        # Set some egraph properties to attributes of self.
        # (Note: self.adjacency also generates indices in the egraph)
        self.egraph = code.egraph
        # self.adjacency = self.egraph.adj_generator(sparse=True)
        # self.num_nodes = len(self.egraph)
        self.to_points = self.egraph.to_points

        self.syndrome_inds = code.all_syndrome_inds

        # Get indices of nodes we wish to designate as ideal.
        perfect_points = self.egraph.graph.get("perfect_points")
        if perfect_points is None:
            self.perfect_inds = None
        else:
            self.perfect_inds = [self.egraph.to_indices[point] for point in perfect_points]

        # Noise related data
        self.state_sampler = state_sampler
        self.homodyme_sampler = homodyme_sampler
        self.translator = translator

    def apply_noise(self, rng=default_rng()):
        """Apply cv-level noise to the graph state.

        Identify the nodes of the EGraph with CV states, measure the syndrome
        assuming the specified noise model, and convert the simulated homodyne
        outcomes to bit values.

        This method modifies the attributes of self.egraph.
        """
        states = self.state_samplers.sample(rng)
        outcomes = self.homodyme_sampler.measure(states, self.syndrome_inds, rng)
        bit_values = self.translator(outcomes)

        self._apply_state_labels(states)
        self._apply_homodyme_outcomes(outcomes)
        self._apply_bit_values(bit_values)

    def _apply_state_labels(self, states):
        """Add state labels as node attributes to the graph state."""
        for psi in states:
            for ind in states[psi]:
                self.egraph.nodes[self.to_points[ind]]["state"] = psi

    def _apply_homodyme_outcomes(self, outcomes):
        """Add homodyme outcome as node attributes to the graph state."""
        for i, ind in enumerate(self.syndrome_inds):
            self.egraph.nodes[self.to_points[ind]]["hom_val_p"] = outcomes[i]

    def _apply_bit_values(self, bit_values):
        """Add bit values as node attributes to the graph state."""
        for i, ind in enumerate(self.syndrome_inds):
            self.egraph.nodes[self.to_points[ind]]["bit_val"] = bit_values[i]


class SqueezedStateSampler:
    """This sampler generates list of indices for p-squeezes and GKP states.

    Args:
        p_swap (float): the probability of a node being a p-squeezed state.
        size (int): the number of states to generate.
        rng (numpy.random.Generator, optional): a random number generator
            following the NumPy API. It can be seeded for reproducibility.
            By default, numpy.random.default_rng is used without a fixed seed.
    """

    def __init__(self, p_swap, size):
        self.p_swap = p_swap
        self.size = size

    def sample(self, rng=default_rng()):
        """Sample the p-squeezes and GKP states.

        Args:
            rng (numpy.random.Generator, optional): a random number generator
                following the NumPy API. It can be seeded for reproducibility.
                By default, numpy.random.default_rng is used without a fixed seed.

        Returns:
            A dict where keys are "p" or "gkp" and represent the states and the values
            are lists of the corresponding indices.
        """
        states = dict()
        states["p"] = self._generate_squeezed_indices(rng)
        states["gkp"] = self._generate_gkp_indices(states["p"])
        return states

    def _generate_squeezed_indices(self, rng=default_rng()):
        """Generate the indices of squeezed states."""
        if self.p_swap == 0:
            return np.array([])
        elif self.p_swap == 1:
            return np.arange(self.size)
        else:
            num_p = rng.binomial(self.size, self.p_swap)
            inds = rng.choice(range(self.size), size=int(np.floor(num_p)), replace=False)
            return inds

    def _generate_gkp_indices(self, used_inds):
        """Associate all indices not in self.states with GKP states."""
        remaining_inds = list(set(range(self.size)) - set(used_inds))
        return np.array(remaining_inds, dtype=int)


class InitialOrderHomodymeSampler:
    def __init__(self, code, delta, quad):
        self.delta = delta
        self.quad = quad

        self.adjacency = self.egraph.adj_generator(sparse=True)
        self.size = len(code.egraph)

        # Get indices of nodes we wish to designate as ideal.
        perfect_points = self.egraph.graph.get("perfect_points")
        if perfect_points is None:
            self.perfect_inds = None
        else:
            self.perfect_inds = [code.egraph.to_indices[point] for point in perfect_points]

    def measure(self, states, inds=None, rng=default_rng()):
        if inds is None:
            inds = range(self.size)
        means = self.means()
        covs = self.covariances(states)
        outcomes = rng.normal(means, covs)
        outcomes = SCZ_apply(self.adjacency, outcomes)
        if self.quad == "q":
            return outcomes[: self.size][inds]
        elif self.quad == "p":
            return outcomes[self.size :][inds]

    def means(self):
        """Return the means for the homodyne measurement sample."""
        return np.zeros(2 * self.size, dtype=np.float32)

    def covariances(self, states):
        """Return the covariances for the homodyne measurement sample."""
        covs = np.zeros(2 * self.size, dtype=np.float32)
        noise_q = {"p": 1 / (2 * self.delta) ** 0.5, "GKP": (self.delta / 2) ** 0.5}
        noise_p = {"p": (self.delta / 2) ** 0.5, "GKP": (self.delta / 2) ** 0.5}
        for state, indices in states:
            if self.perfect_inds:
                indices = np.array(list(set(indices).difference(self.perfect_inds)))
            if len(indices) > 0:
                covs[indices] = noise_q[state]
                covs[indices + self.size] = noise_p[state]
        return covs


class TwoStepsHomodymeSampler:
    def __init__(self, code, delta, propagate=True):
        self.delta = delta
        self.adjacency = self.egraph.adj_generator(sparse=True)
        self.size = len(code.egraph)
        self.syndrome_inds = code.all_syndrome_inds

        # Get indices of nodes we wish to designate as ideal.
        perfect_points = self.egraph.graph.get("perfect_points")
        if perfect_points is None:
            self.perfect_inds = None
        else:
            self.perfect_inds = [code.egraph.to_indices[point] for point in perfect_points]

    def measure(self, states, inds=None, rng=default_rng()):
        if inds is None:
            inds = range(self.size)
        means = self.means(states, rng)
        covs = self.covariances(inds)
        return rng.normal(means, covs)

    def means(self, states, rng=default_rng()):
        """Return the means for the homodyne measurement sample.

        Setting propagate to True applies a symplectic CZ matrix to the
        means.
        """
        means = np.zeros(2 * self.size, dtype=np.float32)
        val_funcs = {
            "p": lambda n: rng.random(size=n) * (2 * np.sqrt(np.pi)),
            "gkp": lambda n: rng.integers(0, 2, size=n) * np.sqrt(np.pi),
        }
        for state, indices in states.items():
            n_inds = len(indices)
            if n_inds > 0:
                means[indices] = val_funcs[state](n_inds)
        if self.propagate:
            means = SCZ_apply(self.adjacency, means)
        return means

    def covariances(self, inds):
        """Return the covariances for the homodyne measurement sample."""
        covs = np.zeros(2 * self.size, dtype=np.float32)
        covs = np.full(len(inds), (self.delta / 2) ** 0.5, dtype=np.float32)
        if self.perfect_inds:
            inds_to_0 = set(inds).intersection(self.perfect_inds)
            ind_arr = np.empty(len(inds_to_0), dtype=np.int64)
            for i, perfect_ind in enumerate(inds_to_0):
                ind_arr[i] = (inds == perfect_ind).nonzero()[0][0]
            covs[ind_arr] = 0
        return covs


# class CVMacroLayer(CVLayer):
#     """A class for reducing a macronode CV graph to a canonical graph.

#     Applies noise to self.egraph (assumed a macronode graph), entangles the
#     macronodes, measures the syndrome, and populates the canonical graph
#     reduced_graph with the reduced states, bit values, and error probabilities.

#     In addition to CVLayer args, the following:

#     Args:
#         bs_network (np.array, optional): the sympletic matrix corresponding to
#             the beamsplitter network entangling the macronode. By default,
#             the standard four-splitter.
#     """

#     def __init__(self, code, *, delta, bs_network=None, **kwargs):
#         # Macronize the code lattice. Pad the boundary (i.e. ensure that all
#         # macronodes have exactly 4 nodes) in the case that the boundaries
#         # are not all-periodic.
#         pad_bool = code.bound_str != "periodic"
#         macro_graph = code.graph.macronize(pad_bool)

#         # Instantiate the CVLayer parent class with the right noise model.
#         super().__init__(macro_graph, delta=delta, sampling_order="two-step", **kwargs)
#         self.reduced_graph = code.graph
#         if bs_network is None:
#             self.bs_network = splitter_symp()

#     def apply_noise(self, rng=default_rng()):
#         """Reduce the macronode code lattice to the canonical code lattice.

#         Follow the procedure in arXiv:2104.03241. Take the macronode lattice
#         macro_graph, apply noise, designate micronodes as planets and stars,
#         conduct homodyne measurements, process these measurements, and compute
#         conditional phase error probabilities.

#         This method modifies the node attributes of self.reduced_graph to include
#         effective bit values and phase error probabilities.
#         """
#         # Sample for the initial state types
#         self.populate_states(rng=rng)
#         # Obtain permuted indices, with stars (central nodes at the front)
#         self._permute_star_inds()
#         # Apply body indices to the macronode graph and effective state labels
#         # to the reduced graph.
#         self._apply_macro_and_reduced_labels()
#         # Apply symplectic matrices corresponding to CZ gates and the
#         # beamsplitter networks.
#         self._entangle_states()
#         # Measure the syndrome, corresponding to memory-mode operaiton.
#         self._measure_syndrome(rng=rng)
#         # Process homodyne outcomes to calculate effective bit values
#         # and phase error probabilities for the reduced nodes.
#         for j in range(0, self._N - 3, 4):
#             self._reduce_jth_macronode(j)

#     def _apply_macro_and_reduced_labels(self):
#         """Label body indices and types of reduced nodes.

#         In the macronode graph: for each star, associate a 'body_index' of 1,
#         and 2, 3, and 4 for the subsequent planets. In the reduced graph: if a
#         macronode contains at least one GKP state, label the reduced state as
#         'GKP' (otherwise 'p').

#         This method modifies self.egraph and self.reduced_graph.
#         """
#         for i, ind in enumerate(self.permuted_inds):
#             self.egraph.nodes[self._to_points[ind]]["body_index"] = i % 4 + 1
#         for ind in self.permuted_inds[::4]:
#             point = self._to_points[ind]
#             centre_point = tuple(round(i) for i in point)
#             self.reduced_graph.nodes[centre_point]["state"] = self.egraph.nodes[point]["state"]

#     def _entangle_states(self):
#         """Entangle the states in the macro_graph.

#         Apply CZ gates to (i.e. a symplectic CZ matrix to the quadratures of)
#         self.egraph, based on where the edges are in the graph. Then, apply the
#         four-splitter to each macronode.

#         This method sets the attributes self.permuted_quads to the permuted
#         quadratures and self.quad_permutation to the corresponding permutation
#         vector.
#         """
#         N = self._N
#         quads = self._means_sampler(propagate=True)
#         # Permute the quadrature values to align with the permuted
#         # indices in order to apply the beamsplitter network.
#         quad_permutation = np.concatenate([self.permuted_inds, N + self.permuted_inds])
#         permuted_quads = quads[quad_permutation]
#         # The beamsplitter network
#         symp_bs = self.bs_network
#         for i in range(0, N - 3, 4):
#             q_inds = np.array([i, i + 1, i + 2, i + 3])
#             p_inds = q_inds + N
#             updated_qs = symp_bs[:4, :4] @ permuted_quads[q_inds]
#             updated_ps = symp_bs[4:, 4:] @ permuted_quads[p_inds]
#             permuted_quads[q_inds] = updated_qs
#             permuted_quads[p_inds] = updated_ps
#         self.permuted_quads = permuted_quads
#         self.quad_permutation = quad_permutation

#     def _hom_outcomes(self, vertex):
#         """Measurement outcomes in the macronode containing vertex.

#         Return the values of the homodyne measurements of the macronode
#         containing vertex. Note we are only interested in q-homodyne
#         outcomes; the returned list is of the form [0, 0, q2, q3, q4].
#         If vertex is None, return a list of 0s, so that the processing
#         is unaltered by the outcomes.
#         """
#         macro_graph = self.egraph
#         if vertex is None:
#             return [0, 0, 0, 0, 0]
#         meas = np.zeros(5)
#         # The central node corresponding to the neighboring
#         # macronode.
#         central_node = tuple(round(i) for i in vertex)
#         for micro in macro_graph.macro_to_micro[central_node]:
#             index = macro_graph.nodes[micro]["body_index"]
#             # Populate meas with the q-homodyne outcomes for
#             # the planet modes.
#             if index != 1:
#                 meas[index] = macro_graph.nodes[micro]["hom_val_q"]
#         return meas

#     def _measure_syndrome(self, rng=default_rng()):
#         """Measure the syndrome of self.egraph.

#         Conduct p-homodyne measurements on the stars (central modes) of
#         self.egraph and q-homodyne measurements to the planets
#         (satellite modes). This effectively conducts an X-basis measurement
#         on the modes of the reduced lattice.

#         This method modifies self.egraph.
#         """
#         N = self._N
#         unpermuted_quads = self.permuted_quads[invert_permutation(self.quad_permutation)]
#         # Indices of stars and planets.
#         stars = self.permuted_inds[::4]
#         planets = np.delete(self.permuted_inds, np.arange(0, N, 4))
#         # Update quadrature values after CZ gate application.
#         # Measure stars in p, planets in q.
#         self.measure_hom(quad="p", inds=stars, updated_means=unpermuted_quads)
#         self.measure_hom(quad="q", inds=planets, updated_means=unpermuted_quads)

#     def _neighbor_of_micro_i_macro_j(self, i, j):
#         """Return the neighbor of the ith micronode, jth macronode in
#         self.egraph.

#         Suppose micronode i (in macronode j) is adjacent to a neighbor.
#         Return the vertex (tuple) and the body index of the neighbor to
#         help the subsequent processing rules. If there is no such
#         neighbor, return None.
#         """
#         macro_graph = self.egraph
#         # Index of ith micronode in the jth macronode.
#         ith_index = self.permuted_inds[j + i - 1]
#         ith_vertex = self._to_points[ith_index]
#         # Vertex of the neighbor of the ith micronode.
#         ith_adjacency = list(macro_graph[ith_vertex])
#         if ith_adjacency:
#             ith_neighbor = list(macro_graph[ith_vertex])[0]
#             ith_body_index = macro_graph.nodes[ith_neighbor]["body_index"]
#             return ith_neighbor, ith_body_index
#         return None

#     def _permute_star_inds(self):
#         """Obtain permuted indices for the macronode graph with stars in front.

#         For each macronode, permute indices so that the first
#         encountered GKP state comes first, designating it as the 'star'
#         ('central') mode. The first index in the resulting list and
#         every four indices thereafter correspond to star modes. The rest
#         are 'planets' ('satellite' modes).
#         """
#         inds = np.reshape(np.arange(self._N), (self._N // 4, 4))
#         self.permuted_inds = np.apply_along_axis(
#             self._permute_gkp_inds_in_macronode, 1, inds
#         ).flatten()

#     def _permute_gkp_inds_in_macronode(self, inds):
#         """Place the indices of GKP states within inds to the front."""
#         gkps = [self.egraph.nodes[self._to_points[ind]]["state"] == "GKP" for ind in inds]
#         return np.concatenate((inds[gkps], inds[[not ind for ind in gkps]]))

#     def _process_neighboring_outcomes(self, neighbor_hom_vals, neighbor_body_index):
#         """Process homodyne outcomes for a neighboring macronode.

#         Suppose some micronode is connected to a neighboring micronode,
#         i. i has a certain body index and belongs to a macronode with
#         measurement results neighbor_hom_vals. Use this information to
#         process the results (i.e. find appropriate linear combinations
#         of neighbor_hom_vals).
#         """
#         if neighbor_body_index == 1:
#             return 0
#         if neighbor_body_index == 2:
#             return neighbor_hom_vals[2] - neighbor_hom_vals[4]
#         if neighbor_body_index == 3:
#             return neighbor_hom_vals[3] - neighbor_hom_vals[4]
#         if neighbor_body_index == 4:
#             return neighbor_hom_vals[2] + neighbor_hom_vals[3]
#         return None

#     def _reduce_jth_macronode(self, j):
#         """Obtain the bit value and error probability for the jth macronode.

#         Process the measurement outcomes of the jth macronode in self.egraph to
#         populate the corresponding reduced node in reduced_graph with a bit value
#         and conditional phase error probability.

#         This method modifies self.reduced_graph.
#         """
#         macro_graph = self.egraph
#         delta = self.delta
#         star_index = self.permuted_inds[j]
#         vertex = self._to_points[star_index]

#         # Here, j corresponds to the macronode and i to to micronode.
#         # i ranges from 1 to 4, to align with manuscript.
#         verts_and_inds = [self._neighbor_of_micro_i_macro_j(i, j) for i in (1, 2, 3, 4)]
#         neighbors = [tup[0] if tup else None for tup in verts_and_inds]
#         body_indices = [tup[1] if tup else None for tup in verts_and_inds]

#         # Array of arrays of measurement outcomes in all the
#         # macronodes adjacent to j.
#         m_arr = np.array([self._hom_outcomes(neighbors[i - 1]) for i in (1, 2, 3, 4)])
#         # Array of processed q-homodyne outcomes from neighboring
#         # macronodes of the form [0, Z(1), Z(2), Z(3), Z(4)].
#         Z_arr = np.array(
#             [0]
#             + [
#                 self._process_neighboring_outcomes(m_arr[i - 1], body_indices[i - 1])
#                 for i in (1, 2, 3, 4)
#             ]
#         )
#         # p-homodyne outcome of the star node.
#         star_p_val = macro_graph.nodes[vertex]["hom_val_p"]

#         # Types of state for the four micronodes directly neighboring
#         # macronode j.
#         types = [
#             macro_graph.nodes[neighbor]["state"] if neighbor else None for neighbor in neighbors
#         ]
#         # Phase error probability and number of p-squeezed states
#         # among the four micronodes in the vicinity of macronode j
#         p_err = 0
#         num_p = 0
#         outcome = 2 * star_p_val
#         gkp_inds = []
#         for i in (1, 2, 3, 4):
#             if types[i - 1] == "p":
#                 num_p += 1
#                 outcome -= Z_arr[i]
#             if types[i - 1] == "GKP":
#                 if delta > 0:
#                     p_err += Z_err_cond(2 * delta, Z_arr[i], use_hom_val=True)
#                 gkp_inds += [i]
#         if delta > 0:
#             p_err += Z_err_cond(2 * (2 + num_p) * delta, outcome, use_hom_val=True)
#         p_err = min(p_err, 0.5)
#         p_err = max(p_err, 0)

#         bitp = GKP_binner([outcome])[0]
#         bitq = GKP_binner(Z_arr[gkp_inds].astype(np.float64)) if gkp_inds else 0

#         processed_bit_val = (bitp + np.sum(bitq)) % 2

#         # Update the reduced RHG lattice with the effective
#         # homodyne value and the phase error probability.
#         central_vert = tuple(round(i) for i in vertex)
#         self.reduced_graph.nodes[central_vert]["bit_val"] = processed_bit_val
#         self.reduced_graph.nodes[central_vert]["p_phase_cond"] = p_err
