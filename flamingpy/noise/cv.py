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
import scipy.sparse as sp

from flamingpy.codes import EGraph
from flamingpy.cv.ops import invert_permutation, SCZ_mat, SCZ_apply, splitter_symp
from flamingpy.cv.gkp import GKP_binner, Z_err_cond

# pylint: disable=too-many-instance-attributes
class CVLayer:
    """A class for applying to a code (or an EGraph) a physical layer of
    continuous-variable noise.

    Associates the nodes of an EGraph with continuous-variable quantum states,
    and its edges with continuous-variable CZ gates.

    Currently, a hybrid lattice of p-squeezed and GKP+ states is considered.

    Args:
        code (SurfaceCode or EGraph): the code object (so that code.graph is
            an EGraph), or an EGraph directly, to which the noise is applied.
        delta (float): the quadrature blurring parameter, related to the
            squeezing of the GKP states and the momentum-quadrature variance of
            the p-squeezed states.
        state (dict, optional): the dictionary of all non-GKP states and their
            indices, of the form {'state': []}. By default, all states are
            GKP states.
        p_swap (float, optional): if supplied, the probability of a node being
            a p-squeezed state. Overrides the indices given in state.
        sampling_order (str, optional): the scheme for conducting samplign
            for the homodyne measurements. Options are "initial" and "two-step".
        rng (numpy.random.Generator, optional): a random number generator
            following the NumPy API. It can be seeded for reproducibility.
            By default, numpy.random.default_rng is used without a fixed seed.

    Attributes:
        egraph (EGraph): the unerlying graph representation.
        to_points (dict): the index-to-coordinate dictionary, taken from egraph.
        N (int): the number of qubits in the lattice.
        delta (float): the delta from the Args above.
        p_swap (float): the swap-out probability from the Args above.
        _adj (sp.sparse.csr_matrix): the adjacency matrix of egraph.
        _states (dict): states along with their indices.
        _sampling_order (str): sampling order from above.
        _perfect_inds (list or NoneType): the indices of qubits to treat as
            ideal.
    """

    def __init__(self, code, **kwargs):
        if isinstance(code, EGraph):
            self.egraph = code
        else:
            self.code = code
            self.egraph = code.graph
        # Generate indices if not already done
        self.egraph.index_generator()

        # Set some egraph properties to attributes of self
        self.N = len(self.egraph)
        self.to_points = self.egraph.to_points
        self._adj = self.egraph.adj_generator(sparse=True)

        self.delta = kwargs.get("delta")
        self.p_swap = kwargs.get("p_swap")
        self._states = kwargs.get("states") or {"p": np.empty(0, dtype=int)}
        self._sampling_order = kwargs.get("sampling_order")
        self._perfect_inds = self.egraph.graph.get("perfect_inds")

    def apply_noise(self, sampling_order=None, rng=default_rng()):
        """Apply cv-level noise to the graph state.

        First, label the states; then, measure the syndrome and populate
        self.egraph with bit values.

        This method modifies self.egraph.
        """
        if sampling_order is None:
            sampling_order = self._sampling_order
        self.populate_states(rng=rng)
        self.measure_syndrome(sampling_order=sampling_order, rng=rng)
        self.inner_decoder()

    def populate_states(self, rng=default_rng()):
        """Populate the graph state with state labels.

        Assume the graph states consists of a combination of squeezed states
        and GKP+ states. A non-zero self.p_swap overrides indices specified in
        self._states and uses a binomial distribution to identify some indices
        as p-squeezed states.

        Print a message if both p_swap and p indices are supplied.

        This method modifies self.egraph.
        """
        # Generate indices of squeezed states based on swap-out
        # probability self.p_swap, if supplied.
        if self.p_swap:
            self._generate_squeezed_indices(rng)

        # Associate remaining indices with GKP states.
        self._generate_gkp_indices()

        # Give the EGraph nodes state attributes.
        self._apply_state_labels()

    def _generate_squeezed_indices(self, rng=default_rng()):
        """Generate the indices of squeezed states."""
        if len(self._states["p"]):
            print(
                "Both swap-out probability and indices of p-squeezed states supplied. "
                "Ignoring the indices."
            )
        if self.p_swap == 1:
            self._states["p"] = np.arange(self.N)
        else:
            num_p = rng.binomial(self.N, self.p_swap)
            inds = rng.choice(range(self.N), size=int(np.floor(num_p)), replace=False)
            self._states["p"] = inds

    def _generate_gkp_indices(self):
        """Associate remaining indices with GKP states."""
        used_inds = np.empty(0, dtype=int)
        for psi in self._states:
            used_inds = np.concatenate([used_inds, self._states[psi]])
        remaining_inds = list(set(range(self.N)) - set(used_inds))
        self._states["GKP"] = np.array(remaining_inds, dtype=int)

    def _apply_state_labels(self):
        """Add state labels as node attributes to the graph state."""
        for psi in self._states:
            for ind in self._states[psi]:
                self.egraph.nodes[self.to_points[ind]]["state"] = psi

    def _means_sampler(self, sampling_order="initial", propagate=False, rng=default_rng()):
        """Return the means for the homodyne measurement sample."""
        means = np.zeros(2 * self.N, dtype=np.float32)
        if sampling_order == "two-step":

            def q_val_for_p(n):
                return rng.random(size=n) * (2 * np.sqrt(np.pi))

            def q_val_for_GKP(n):
                return rng.integers(0, 2, size=n) * np.sqrt(np.pi)

            val_funcs = {"p": q_val_for_p, "GKP": q_val_for_GKP}

            for state, indices in self._states.items():
                n_inds = len(indices)
                if n_inds > 0:
                    means[indices] = val_funcs[state](n_inds)
            if propagate:
                means = SCZ_apply(self._adj, covs)
        return means

    def _covs_sampler(self, sampling_order="initial", inds=None, rng=default_rng()):
        """Return the covariances for the homodyne measurement sample."""
        delta = self.delta
        covs = np.zeros(2 * self.N, dtype=np.float32)
        if sampling_order == "initial":
            noise_q = {"p": 1 / (2 * delta) ** 0.5, "GKP": (delta / 2) ** 0.5}
            noise_p = {"p": (delta / 2) ** 0.5, "GKP": (delta / 2) ** 0.5}
            for state, inds in self._states.items():
                if self._perfect_inds:
                    inds = np.array(list(set(inds).difference(self._perfect_inds)))
                if len(inds) > 0:
                    covs[inds] = noise_q[state]
                    covs[inds + self.N] = noise_p[state]

        elif sampling_order == "two-step":
            if inds is None:
                inds = range(self.N)
            N_inds = len(inds)
            covs = np.full(N_inds, (self.delta / 2) ** 0.5, dtype=np.float32)
            if self._perfect_inds:
                inds_to_0 = set(inds).intersection(self._perfect_inds)
                ind_arr = np.empty(len(inds_to_0), dtype=np.int64)
                for i, perfect_ind in enumerate(inds_to_0):
                    ind_arr[i] = (inds == perfect_ind).nonzero()[0][0]
                covs[ind_arr] = (self.delta / 2) ** 0.5

        return covs

    def measure_hom(
        self,
        quad="p",
        inds=None,
        updated_means=None,
        updated_covs=None,
        propagate=True,
        sampling_order="initial",
        rng=default_rng(),
    ):
        """Conduct a homodyne measurement of states in the lattice.

        Simulate a homodyne measurement of quadrature quad of states at indices
        inds according to sampling order specified by self._sampling_order. Use
        the Numpy random sampling method method. If updated_quads is supplied,
        use those instead of applying an SCZ matrix to the initial quads in
        the two-step sampling.

        Args:
            rng (numpy.random.Generator, optional): a random number generator following
                NumPy API. It can be seeded for reproducibility. By default,
                numpy.random.default_rng is used without a fixed seed.
        """
        N = self.N
        if inds is None:
            inds = range(N)
        N_inds = len(inds)

        means = updated_means or self._means_sampler(
            sampling_order=sampling_order, propagate=propagate, rng=rng
        )
        covs = updated_covs or self._covs_sampler(sampling_order=sampling_order, rng=rng)
        if sampling_order == "two-step":
            if quad == "q":
                means = means[:N][inds]
            elif quad == "p":
                means = means[N:][inds]
        outcomes = rng.normal(means, covs)
        if sampling_order == "initial":
            outcomes = SCZ_apply(self._adj, covs)
            if quad == "q":
                outcomes = outcomes[:N][inds]
            elif quad == "p":
                outcomes = outcomes[N:][inds]
        for i in range(N_inds):
            self.egraph.nodes[self.to_points[inds[i]]]["hom_val_" + quad] = outcomes[i]

    def measure_syndrome(self, sampling_order="initial", rng=default_rng()):
        return self.measure_hom(
            quad="p", inds=self.code.all_syndrome_inds, sampling_order=sampling_order, rng=rng
        )

    def SCZ(self, sparse=False):
        """Return the symplectic matrix associated with CZ application.

        Returns:
            array: the symplectic matrix.
        """
        adj = self._adj
        return SCZ_mat(adj, sparse)

    def hom_outcomes(self, inds=None, quad="p"):
        """array: quad-homodyne measurement outcomes for modes inds."""
        if inds is None:
            inds = range(self.N)
        outcomes = [self.egraph.nodes[self.to_points[i]].get("hom_val_" + quad) for i in inds]
        return outcomes

    def bit_values(self, inds=None):
        """array: bit values associated with the p measurement."""
        if inds is None:
            inds = range(self.N)
        bits = [self.egraph.nodes[self.to_points[i]].get("bit_val") for i in inds]
        return bits

    @property
    def p_inds(self):
        """array: the indices of the p-squeezed states."""
        return self._states.get("p")

    @property
    def GKP_inds(self):
        """array: the indices of the GKP states."""
        return self._states.get("GKP")

    def draw(self, **kwargs):
        """Draw the CV graph state with matplotlib.

        See flamingpy.utils.viz.draw_EGraph for more details. Use the default
        colours: gold for GKP states and blue for p-squeezed
        states.
        """
        cv_opts = {"color_nodes": ("state", {"GKP": "gold", "p": "blue"})}
        updated_opts = {**cv_opts, **kwargs}
        return self.egraph.draw(**updated_opts)

    def draw_SCZ(self, **kwargs):
        """Draw the adjacency matrix of a CV graph state with matplotlib.

        See flamingpy.utils.viz.plot_mat_heat_map for more details.
        """
        from flamingpy.utils.viz import plot_mat_heat_map

        return plot_mat_heat_map(self.SCZ(), **kwargs)

    def inner_decoder(self, translator=GKP_binner):
        """Convert homodyne outcomes to bit values according to translator.

        The inner (CV) decoder, a.k.a. translator, a.k.a binning function. Set
        converted values to the bit_val attribute for nodes in G.

        Args:
            code (SurfaceCode): the qubit QEC code
            translator (func): the choice of binning function; by default, the
                standard GKP binning function that snaps to the closest
                integer multiple of sqrt(pi).
        Returns:
            None
        """
        for point in self.code.all_syndrome_coords:
            hom_val = self.code.graph.nodes[point]["hom_val_p"]
            bit_val = translator([hom_val])[0]
            self.code.graph.nodes[point]["bit_val"] = bit_val


class CVMacroLayer(CVLayer):
    """A class for reducing a macronode CV graph to a canonical graph.

    Applies noise to self.egraph (assumed a macronode graph), entangles the
    macronodes, measures the syndrome, and populates the canonical graph
    reduced_graph with the reduced states, bit values, and error probabilities.

    In addition to CVLayer args, the following --

    Args:
        reduced_graph (EGraph): the canonical (reduced) code graph.
        bs_network (np.array, optional): the sympletic matrix corresponding to
            the beamsplitter network entangling the macronode. By default,
            the standard four-splitter.
    """

    def __init__(self, code, bs_network=None, **kwargs):
        macro_graph = code.graph.macronize(pad_boundary=True)
        macro_graph.index_generator()
        macro_graph.adj_generator(sparse=True)

        perfect_points = macro_graph.graph.get("perfect_points")
        if perfect_points:
            perfect_inds = [macro_graph.to_indices[point] for point in perfect_points]
        else:
            perfect_inds = None

        noise_args = {}
        noise_args["perfect_inds"] = perfect_inds
        noise_args["sampling_order"] = "two-step"
        super().__init__(macro_graph, **kwargs)

        self.reduced_graph = code.graph
        if bs_network is None:
            self.bs_network = splitter_symp()

    def apply_noise(self, rng=default_rng()):
        """Reduce the macronode lattice macro_graph to the canonical
        reduced_graph.

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
        # Sample for the initial state
        super().populate_states(rng=rng)
        self._permute_star_inds()
        self._apply_macro_and_reduced_labels()
        self._entangle_states()
        self._measure_syndrome(rng=rng)
        # Process homodyne outcomes and calculate phase error probabilities
        for j in range(0, self.N - 3, 4):
            self._reduce_jth_macronode(j)

    def _permute_star_inds(self):
        """Obtain permuted indices for the macronode graph with stars in front.

        For each macronode, permute indices so that the first encountered GKP 
        state comes first, designating it as the 'star' ('central') mode. The 
        first index in the resulting list and every four indices thereafter 
        correspond to star modes. The rest are 'planets' ('satellite' modes).
        """
        inds = np.reshape(np.arange(self.N), (self.N // 4, 4))
        self.permuted_inds = np.apply_along_axis(self._permute_gkp_inds_in_macronode, 1, inds).flatten()

    def _permute_gkp_inds_in_macronode(self, indices):
        """Place the GKP indices in indices to the front."""
        gkps = [self.egraph.nodes[self.to_points[ind]]["state"] == "GKP" for ind in indices]
        return np.concatenate((indices[gkps], indices[[not ind for ind in gkps]]))
    
    def _apply_macro_and_reduced_labels(self):
        """Label body indices and types of reduced nodes.

        In the macronode graph: for each star, associate a 'body_index' of 1, 
        and 2, 3, and 4 for the subsequent planets. In the reduced graph: if a 
        macronode contains at least one GKP state, label the reduced state as 
        'GKP' (otherwise 'p').

        This method modifies self.egraph and self.reduced_graph.
        """
        for i, ind in enumerate(self.permuted_inds):
            self.egraph.nodes[self.to_points[ind]]["body_index"] = i % 4 + 1
        for ind in self.permuted_inds[::4]:
            point = self.to_points[ind]
            centre_point = tuple(round(i) for i in point)
            self.reduced_graph.nodes[centre_point]["state"] = self.egraph.nodes[point]["state"]

    def _entangle_states(self):
        """Entangle the states in the macro_graph.

        Apply CZ gates to (i.e. a symplectic CZ matrix to the quadratures of)
        self.egraph, based on where the edges are in the graph. Then, apply the
        four-splitter to each macronode.

        This method sets the attributes self.permuted_quads to the permuted
        quadratures and self.quad_permutation to the corresponding permutation
        vector.
        """
        quads = self._means_sampler(propagate=True)
        # Permute the quadrature values to align with the permuted
        # indices in order to apply the beamsplitter network.
        N = self.N
        quad_permutation = np.concatenate([self.permuted_inds, N + self.permuted_inds])
        permuted_quads = quads[quad_permutation]
        # The beamsplitter network
        symp_bs = self.bs_network
        for i in range(0, N - 3, 4):
            q_inds = np.array([i, i + 1, i + 2, i + 3])
            p_inds = q_inds + N
            updated_qs = symp_bs[:4, :4] @ permuted_quads[q_inds]
            updated_ps = symp_bs[4:, 4:] @ permuted_quads[p_inds]
            permuted_quads[q_inds] = updated_qs
            permuted_quads[p_inds] = updated_ps
        self.permuted_quads = permuted_quads
        self.quad_permutation = quad_permutation

    def _measure_syndrome(self, rng=default_rng()):
        """Measure the syndrome of self.egraph.

        Conduct p-homodyne measurements on the stars (central modes) of
        self.egraph and q-homodyne measurements to the planets
        (satellite modes). This effectively conducts an X-basis measurement
        on the modes of the reduced lattice.

        This method modifies self.egraph.
        """
        N = self.N
        unpermuted_quads = self.permuted_quads[invert_permutation(self.quad_permutation)]
        # Indices of stars and planets.
        stars = self.permuted_inds[::4]
        planets = np.delete(self.permuted_inds, np.arange(0, N, 4))
        # Update quadrature values after CZ gate application.
        # Measure stars in p, planets in q.
        self.measure_hom(quad="p", inds=stars, updated_means=list(unpermuted_quads))
        self.measure_hom(quad="q", inds=planets, updated_means=list(unpermuted_quads))

    def _neighbor_of_micro_i_macro_j(self, i, j):
        """Return the neighbor of the ith micronode, jth macronode in
        self.egraph.

        Suppose micronode i (in macronode j) is adjacent to a neighbor.
        Return the vertex (tuple) and the body index of the neighbor to
        help the subsequent processing rules. If there is no such
        neighbor, return None.
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
        outcomes; the returned list is of the form [0, 0, q2, q3, q4].
        If vertex is None, return a list of 0s, so that the processing
        is unaltered by the outcomes.
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

        Suppose some micronode is connected to a neighboring micronode,
        i. i has a certain body index and belongs to a macronode with
        measurement results neighbor_hom_vals. Use this information to
        process the results (i.e. find appropriate linear combinations
        of neighbor_hom_vals).
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

        Process the measurement outcomes of the jth macronode in self.egraph to
        populate the corresponding reduced node in reduced_graph with a bit value
        and conditional phase error probability.

        This method modifies self.reduced_graph.
        """
        macro_graph = self.egraph
        delta = self.delta
        star_index = self.permuted_inds[j]
        vertex = self.to_points[star_index]

        # Here, j corresponds to the macronode and i to to micronode.
        # i ranges from 1 to 4, to align with manuscript.
        verts_and_inds = [self._neighbor_of_micro_i_macro_j(i, j) for i in (1, 2, 3, 4)]
        neighbors = [tup[0] if tup else None for tup in verts_and_inds]
        body_indices = [tup[1] if tup else None for tup in verts_and_inds]

        # Array of arrays of measurement outcomes in all the
        # macronodes adjacent to j.
        m_arr = np.array([self._hom_outcomes(neighbors[i - 1]) for i in (1, 2, 3, 4)])
        # Array of processed q-homodyne outcomes from neighboring
        # macronodes of the form [0, Z(1), Z(2), Z(3), Z(4)].
        Z_arr = np.array(
            [0]
            + [
                self._process_neighboring_outcomes(m_arr[i - 1], body_indices[i - 1])
                for i in (1, 2, 3, 4)
            ]
        )
        # p-homodyne outcome of the star node.
        star_p_val = macro_graph.nodes[vertex]["hom_val_p"]

        # Types of state for the four micronodes directly neighboring
        # macronode j.
        types = [
            macro_graph.nodes[neighbor]["state"] if neighbor else None for neighbor in neighbors
        ]
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
