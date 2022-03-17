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

import numpy as np
from numpy.random import default_rng
import scipy.sparse as sp

# from flamingpy.cv.gkp import Z_err, Z_err_cond


def SCZ_mat(adj):
    """Return a symplectic matrix corresponding to CZ gate application.

    Give the 2N by 2N symplectic matrix for CZ gate application based on the
    adjacency matrix adj. Assumes quadrature-like convention:

        (q1, ..., qN, p_1, ..., p_N).

    Args:
        adj (array): N by N binary symmetric matrix. If modes i and j are
            linked by a CZ, then entry ij and ji is equal to the weight of the
            edge (1 by default); otherwise 0.
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
        # TODO: Specify different kind of Scipy sparse matrix?
        identity = sp.identity(N, dtype=np.int8)
        zeros = sp.csr_matrix((N, N), dtype=np.int8)
        block_func = sp.bmat
    # Construct symplectic
    symplectic = block_func([[identity, zeros], [adj, identity]])
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


class CVLayer:
    """A class for applying to an EGraph a physical layer of continuous-
    variable states.

    Has all the functionality of an EGraph, but associates its nodes with
    continuous-variable quantum states and its edges with continuous-variable
    CZ gates.

    For now, only a hybrid state of p-squeezed and GKP states is considered.

    Args:
        g (graph-type): the graph underlying the state.
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
        _sampling_order (str): the sampling order from the Args above (after
            noise applied).
        _adj (array): adjancency matrix of the underlying graph.
        to_points (dict): pointer to self.egraph.to_points, the dictionary from
            indices to coordinates.
    """

    def __init__(self, g, states={"p": np.empty(0, dtype=int)}, p_swap=0, rng=default_rng()):
        """Initialize the CVGraph."""
        self.egraph = g
        self._N = len(g)

        # Instantiate the adjacency matrix
        self._adj = self.egraph.adj_generator(sparse=True)

        if states:
            self._states = states.copy()
            # Non-zero swap-out probability overrides indices specified
            # in states and hybridizes the lattice. Print a message if
            # both supplied.
            if p_swap:
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

            # Associate remaining indices with GKP states.
            used_inds = np.empty(0, dtype=int)
            for psi in self._states:
                used_inds = np.concatenate([used_inds, self._states[psi]])
            remaining_inds = list(set(range(self._N)) - set(used_inds))
            self._states["GKP"] = np.array(remaining_inds, dtype=int)

            # Generate EGraph indices.
            self.egraph.index_generator()
            self.to_points = self.egraph.to_points

            for psi in self._states:
                for ind in self._states[psi]:
                    self.egraph.nodes[self.to_points[ind]]["state"] = psi

    def apply_noise(self, model={}, rng=default_rng()):
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

    def eval_Z_probs(self, inds=None, exact=True, cond=False):
        """Evaluate the probability of phase errors at nodes inds.

        If inds not specified, compute probabilities for all nodes.
        """
        pass

    def SCZ(self, sparse=False):
        """Return the symplectic matrix associated with CZ application.

        Returns:
            array: the symplectic matrix.
        """
        adj = self._adj
        return SCZ_mat(adj)

    def Z_probs(self, inds=None, cond=False):
        """array: the phase error probabilities of modes inds."""
        pass

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

    @property
    def noise_cov(self):
        """array: the noise covariance matrix."""
        if self._sampling_order == "final":
            return self._noise_cov
        print('Sampling order must be "final."')

    def draw(self, **kwargs):
        """Draw the CV graph state with matplotlib.

        See flamingpy.utils.viz.draw_EGraph for more details. Use the default
        colours: gold for GKP states and blue for p-squeezed
        states.
        """
        cv_opts = {"color_nodes": "state", "state_colors": {"GKP": "gold", "p": "blue"}}
        updated_opts = {**cv_opts, **kwargs}
        return self.egraph.draw(**updated_opts)
