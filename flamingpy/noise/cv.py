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

from flamingpy.cv.ops import invert_permutation, SCZ_mat, SCZ_apply, splitter_symp
from flamingpy.cv.gkp import GKP_binner, Z_err_cond

# pylint: disable=too-many-instance-attributes
class CVLayer:
    """A class for applying to a code (or an EGraph) a physical layer of
    continuous-variable noise.

    Associates the nodes of an EGraph with continuous-variable quantum states,
    and its edges with continuous-variable CZ gates.

    For now, only a hybrid state of p-squeezed and GKP states is considered.

    Args:
        code (SurfaceCode or EGraph): the code object (so that code.graph is
            an EGraph) or an EGraph directly.
        state (dict, optional): the dictionary of all non-GKP states and their
            indices, of the form {'state': []}. By default, all states are
            GKP states.
        p_swap (float, optional): if supplied, the probability of a node being
            a p-squeezed state. Overrides the indices given in state.
        rng (numpy.random.Generator, optional): a random number generator
            following the NumPy API. It can be seeded for reproducibility.
            By default, numpy.random.default_rng is used without a fixed seed.

    Attributes:
        egraph (EGraph): the unerlying graph representation.
        _N (int): the number of qubits in the lattice.
        _adj (sp.sparse.csr_matrix): the adjacency matrix of egraph.
        _states (dict): states along with their indices.
        _delta (float): the delta from the Args above (after noise applied)
        _adj (array): adjancency matrix of the underlying graph.
        to_points (dict): pointer to self.egraph.to_points, the dictionary from
            indices to coordinates.
    """

    def __init__(self, code, states=None, p_swap=0, rng=default_rng()):
        """Initialize the CVGraph."""
        if code.__class__.__name__ == "EGraph":
            self.egraph = code
        else:
            self.egraph = code.graph
        self._N = len(self.egraph)

        self._init_quads = None
        self._noise_cov = None
        self._init_noise = None
        self._perfect_inds = None
        self._sampling_order = None
        self._delta = None
        self.to_points = None

        # Instantiate the adjacency matrix
        self._adj = self.egraph.adj_generator(sparse=True)

        self._states = states or {"p": np.empty(0, dtype=int)}

        # Generate indices of squeezed states based on swap-out
        # probability p_swap.
        if p_swap:
            self._generate_squeezed_indices(p_swap, rng)

        # Associate remaining indices with GKP states.
        self._generate_gkp_indices()

        # Give the EGraph nodes state attributes.
        self._apply_state_labels()

    def _apply_state_labels(self):
        self.egraph.index_generator()
        self.to_points = self.egraph.to_points

        for psi in self._states:
            for ind in self._states[psi]:
                self.egraph.nodes[self.to_points[ind]]["state"] = psi

    def _generate_gkp_indices(self):
        """Associate remaining indices with GKP states."""
        used_inds = np.empty(0, dtype=int)
        for psi in self._states:
            used_inds = np.concatenate([used_inds, self._states[psi]])
        remaining_inds = list(set(range(self._N)) - set(used_inds))
        self._states["GKP"] = np.array(remaining_inds, dtype=int)

    def _generate_squeezed_indices(self, p_swap, rng):
        """Use swap-out probability p_swap to hybridize the CV graph state.

        A non-zero p_swap overrides indices specified in states and uses
        a binomial distribution to associate some indices as p-squeezed
        states.

        Print a message if both p_swap and p indices are supplied.

        Args:
            p_swap (float): the swap-out probability.
            rng (numpy.random.Generator, optional): A random number generator
                following the NumPy API. It can be seeded for reproducibility.
                By default, numpy.random.default_rng is used without a fixed seed.
        """
        if len(self._states["p"]):
            print(
                "Both swap-out probability and indices of p-squeezed states supplied. "
                "Ignoring the indices."
            )
        if p_swap == 1:
            self._states["p"] = np.arange(self._N)
        else:
            num_p = rng.binomial(self._N, p_swap)
            inds = rng.choice(range(self._N), size=int(np.floor(num_p)), replace=False)
            self._states["p"] = inds

    def apply_noise(self, model=None, rng=default_rng()):
        """Apply the noise model in model with a random number generator rng.

        Args:
            model (dict, optional): the noise model dictionary of the form
                (default values displayed):

                {'model': 'grn', 'sampling_order': 'initial', 'delta': 0.01,
                 'perfect_inds': self.egraph.graph.get('perfect_inds')}

                'grn; stands for Gaussian Random Noise; sampling_order dictates
                how to simulate measurement outcomes: sample from an
                uncorrelated noise matrix initially ('initial'), a correlated
                noise matrix finally ('final'), or for ideal homodyne outcomes
                initially and from a separable noise covariance matrix finally
                ('two-step'); 'delta' is the quadrature blurring parameter,
                related to the squeezing of the GKP states and the
                momentum-quadrature variance of the p-squeezed states.
                'perfect_inds' is a list of indices to which noise should not
                be applied. By default, looks to the "perfect_inds" attribute
                of egraph.graph.

            rng (numpy.random.Generator, optional): a random number generator
                following the NumPy API. It can be seeded for reproducibility.
                By default, numpy.random.default_rng is used without a fixed
                seed.
        """
        if model is None:
            model = {}

        # Modelling the states.
        perfect_inds = self.egraph.graph.get("perfect_inds")
        default_model = {
            "noise": "grn",
            "delta": 0.01,
            "sampling_order": "initial",
            "perfect_inds": perfect_inds,
        }
        model = {**default_model, **model}
        self._delta = model["delta"]
        self._sampling_order = model["sampling_order"]
        self._perfect_inds = model["perfect_inds"]
        if model["noise"] == "grn":
            self.grn_model(rng)

    def grn_model(self, rng=default_rng()):
        """Apply Gaussian Random Noise model to the CVGraph.

        Store quadrature or noise information as attributes depending on the
        sampling order.

        Args:
            rng (numpy.random.Generator, optional): a random number generator
                following the NumPy API. It can be seeded for reproducibility.
                By default, numpy.random.default_rng is used without a fixed
                seed.
        """
        N = self._N
        delta = self._delta

        # For initial and final sampling, generate noise array depending on
        # quadrature and state.
        if self._sampling_order in ("initial", "final"):
            noise_q = {"p": 1 / (2 * delta) ** 0.5, "GKP": (delta / 2) ** 0.5}
            noise_p = {"p": (delta / 2) ** 0.5, "GKP": (delta / 2) ** 0.5}
            self._init_noise = np.zeros(2 * N, dtype=np.float32)
            for state, inds in self._states.items():
                if self._perfect_inds:
                    inds = np.array(list(set(inds).difference(self._perfect_inds)))
                if len(inds) > 0:
                    self._init_noise[inds] = noise_q[state]
                    self._init_noise[inds + N] = noise_p[state]

        # For final sampling, apply a symplectic CZ matrix to the initial noise
        # covariance.
        if self._sampling_order == "final":
            self._noise_cov = SCZ_apply(self._adj, sp.diags(self._init_noise) ** 2)

        # For two-step sampling, sample for initial (ideal) state-dependent
        # quadrature values.
        if self._sampling_order == "two-step":
            q_val_for_p = lambda n: rng.random(size=n) * (2 * np.sqrt(np.pi))
            q_val_for_GKP = lambda n: rng.integers(0, 2, size=n) * np.sqrt(np.pi)
            val_funcs = {"p": q_val_for_p, "GKP": q_val_for_GKP}
            self._init_quads = np.zeros(2 * N, dtype=np.float32)
            for state, indices in self._states.items():
                n_inds = len(indices)
                if n_inds > 0:
                    self._init_quads[indices] = val_funcs[state](n_inds)

    # pylint: disable=too-many-arguments
    def measure_hom(
        self,
        quad="p",
        inds=None,
        method="cholesky",
        updated_quads=None,
        rng=default_rng(),
    ):
        """Conduct a homodyne measurement on the lattice.

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
        N = self._N
        if inds is None:
            inds = range(N)
        N_inds = len(inds)
        if self._sampling_order == "initial":
            init_samples = rng.normal(0, self._init_noise)
            outcomes = SCZ_apply(self._adj, init_samples)
            if quad == "q":
                outcomes = outcomes[:N][inds]
            elif quad == "p":
                outcomes = outcomes[N:][inds]
        if self._sampling_order == "two-step":
            if updated_quads is not None:
                updated = updated_quads
            else:
                means = self._init_quads
                adj = self.egraph.adj_generator(sparse=True)
                updated = SCZ_apply(adj, means)
            if quad == "q":
                means = updated[:N][inds]
            elif quad == "p":
                means = updated[N:][inds]
            sigma = np.full(N_inds, (self._delta / 2) ** 0.5, dtype=np.float32)
            if self._perfect_inds:
                inds_to_0 = set(inds).intersection(self._perfect_inds)
                ind_arr = np.empty(len(inds_to_0), dtype=np.int64)
                for i, perfect_ind in enumerate(inds_to_0):
                    ind_arr[i] = (inds == perfect_ind).nonzero()[0][0]
                sigma[ind_arr] = (self._delta / 2) ** 0.5
            outcomes = rng.normal(means, sigma)
        if self._sampling_order == "final":
            cov_q = self._noise_cov[:N, :N]
            cov_p = self._noise_cov[N:, N:]
            cov_dict = {"q": cov_q, "p": cov_p}
            means = np.zeros(N_inds, dtype=bool)
            covs = cov_dict[quad][inds, :][:, inds].toarray()
            outcomes = rng.multivariate_normal(mean=means, cov=covs, method=method)
        for i in range(N_inds):
            self.egraph.nodes[self.to_points[inds[i]]]["hom_val_" + quad] = outcomes[i]

    def SCZ(self, sparse=False):
        """Return the symplectic matrix associated with CZ application.

        Returns:
            array: the symplectic matrix.
        """
        adj = self._adj
        return SCZ_mat(adj, sparse)

    def hom_outcomes(self, inds=None, quad="p"):
        """array: quad-homodyne measurement outcomes for modes inds."""
        N = self._N
        if inds is None:
            inds = range(N)
        outcomes = [self.egraph.nodes[self.to_points[i]].get("hom_val_" + quad) for i in inds]
        return outcomes

    def bit_values(self, inds=None):
        """array: bit values associated with the p measurement."""
        N = self._N
        if inds is None:
            inds = range(N)
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


four_splitter = splitter_symp()


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

    def __init__(self, *args, reduced_graph, bs_network=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.reduced_graph = reduced_graph
        if bs_network is None:
            self.bs_network = four_splitter

    def _apply_initial_noise(self, noise_model):
        """Set up the two-step noise model for macro_graph.

        Based on noise_model, populate macro_graph with states and sample for
        the initial (ideal) measurement outcomes.

        This method modifies self.egraph.
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

        This method sets the attribute self.permuted_inds to the permuted
        indices and modifies self.egraph and self.reduced_graph.
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

    def _entangle_states(self):
        """Entangle the states in the macro_graph.

        Apply CZ gates to (i.e. a symplectic CZ matrix to the quadratures of)
        self.egraph, based on where the edges are in the graph. Then, apply the
        four-splitter to each macronode.

        This method sets the attributes self.permuted_quads to the permuted
        quadratures and self.quad_permutation to the corresponding permutation
        vector.
        """
        macro_graph = self.egraph
        quads = SCZ_apply(macro_graph.adj_mat, self._init_quads)
        # Permute the quadrature values to align with the permuted
        # indices in order to apply the beamsplitter network.
        N = self._N
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

    def _measure_syndrome(self):
        """Measure the syndrome of self.egraph.

        Conduct p-homodyne measurements on the stars (central modes) of
        self.egrapg and q-homodyne measurements to the planets
        (satellite modes). This effectively conducts an X-basis measurement
        on the modes of the reduced lattice.

        This method modifies self.egraph.
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
        delta = self._delta
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

    def reduce(self, noise_model):
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
        macro_graph = self.egraph
        # Sample for the initial state
        self._apply_initial_noise(noise_model)
        self._permute_indices_and_label()
        self._entangle_states()
        self._measure_syndrome()
        # Process homodyne outcomes and calculate phase error probabilities
        for j in range(0, len(macro_graph) - 3, 4):
            self._reduce_jth_macronode(j)
