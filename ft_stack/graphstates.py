# Copyright 2020 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Classes for representing graph states."""
import networkx as nx
import numpy as np

# TODO: Avoid Niagara errors associated with Matplotlib; e.g.:
# if __name__ != "__main__":
from numpy.random import default_rng as rng

import scipy.sparse as sp
from ft_stack.GKP import Z_err, Z_err_cond


class EGraph(nx.Graph):
    """An enhanced graph based on a NetworkX Graph.

    A class for adding some functionality to a NetworkX graph
    with some short-hand/convenience methods.

    Attributes:
        indexer (str): method for indexing the nodes; 'default' for
            Python's sorted function; 'macronodes' for rounding
            micronodes to integers, sorting those, and
            furthermore sorting the micronodes within each macronodes,
            all using Python's 'sorted'.
        to_indices (dict): if self.index_generator() has been run,
            a dictionary of the form {points: indices}
        to_points (dict): if self.index_generator() has been run,
            a dictionary of the form {indices: points}
        adj_mat (np.array): if self.adj_generator() has been run,
            the adjacency mtrix of the graph.
    """

    def __init__(self, *args, indexer="default", macronodes=False, **kwargs):
        """Initialize an EGraph (itself an NetworkX graph)."""
        super().__init__(*args, **kwargs)
        self.indexer = indexer
        self._macronodes = macronodes
        if macronodes:
            self.macro = nx.Graph()
        self.to_indices = None
        self.to_points = None
        self.adj_mat = None

    def index_generator(self):
        """Return a relabelled graph with indices as labels.

        Point tuples are stored in the 'pos' attribute of the new graph.
        Use the default sort as the index mapping.
        """
        # TODO: Let user specify index mapping.
        # TODO: SortedDict implementation.
        N = self.order()
        if self.to_indices is not None:
            return self.to_indices
        if self.indexer == "default" and not self._macronodes:
            ind_dict = dict(zip(sorted(self.nodes()), range(N)))
        elif self._macronodes:
            macro_graph = self.macro
            for node in self.nodes():
                rounded = tuple(np.round(node).astype(int))
                macro_graph.nodes[rounded]["micronodes"].append(node)
            sorted_macro = sorted(macro_graph)
            points = []
            for vertex in sorted_macro:
                points += self.macro.nodes[vertex]["micronodes"]
            ind_dict = {points[i]: i for i in range(N)}
        self.to_indices = ind_dict
        self.to_points = {index: point for point, index in ind_dict.items()}
        return ind_dict

    def adj_generator(self, sparse=True):
        """Return the adjacency matrix of the graph.

        Indices correspond to sorted nodes.
        """
        if self.adj_mat is not None:
            return self.adj_mat
        if not self.to_points:
            self.index_generator()
        # TODO: SortedDict implementation.
        sorted_nodes = [self.to_points[i] for i in range(self.order())]
        # TODO: New data type in case of fancier weights.
        if sparse:
            adj = nx.to_scipy_sparse_matrix(self, nodelist=sorted_nodes, dtype=np.int8)
        else:
            adj = nx.to_numpy_array(self, nodelist=sorted_nodes, dtype=np.int8)
        self.adj_mat = adj
        # TODO: Heat map?
        return adj

    def slice_coords(self, plane, number):
        """Obtain all the coordinates in an x, y, or z slice.

        Args:
            plane (str): 'x', 'y', or 'z', denoting the slice direction
            number (int): the index of the slice

        Returns:
            list of tuples: the coordinates of the slice.
        """
        plane_dict = {"x": 0, "y": 1, "z": 2}
        plane_ind = plane_dict[plane]
        coords = [point for point in self.nodes if point[plane_ind] == number]
        return coords


def SCZ_mat(adj):
    """Return a symplectic matrix corresponding to CZ gate application.

    Gives the 2N by 2N symplectic matrix for CZ gate application
    based on the adjacency matrix adj. Assumes quadrature-like
    convention: (q1, ..., qN, p_1, ..., p_N).

    Args:
        adj (array): N by N binary symmetric matrix. If modes i and j
            are linked by a CZ, then entry ij and ji is equal to the
            weight of the edge (1 by default); otherwise 0.
    Returns:
        np.array or sp.sparse.csr_matrix: 2N by 2N symplectic matrix.
            sparse if the adjacency matrix is sparse.
    """
    # Number of modes
    N = adj.shape[0]
    if type(adj) == np.ndarray:
        identity = np.eye(N, dtype=np.int8)
        zeros = np.zeros((N, N), dtype=np.int8)
        block_func = np.block
    else:
        # TODO: Specify kind of Scipy sparse matrix?
        identity = sp.identity(N, dtype=np.int8)
        zeros = sp.csr_matrix((N, N), dtype=np.int8)
        block_func = sp.bmat
    # Construct symplectic
    symplectic = block_func([[identity, zeros], [adj, identity]])
    return symplectic


def SCZ_apply(adj, quads, one_shot=True):
    """Apply SCZ matrix to one- or two-dimensional array quads.

    If one-shot is True, use SCZ_mat to apply a symplectic CZ matrix
    to a matrix or vector of quadratures. Otherwise, take advantage of
    the block structure of a symplectic SCZ matrix for a more memory-
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


class CVGraph:
    """A class for representing continuous-variable graph states.

    Has all the functionality of an EGraph, but associates its
    nodes with continuous-variable quantum states and its edges with
    continuous-variable CZ gates.

    For now, only a hybrid state of p-squeezed and GKP states is
    considered.

    Args:
        g (graph-type): the graph underlying the state.
        state (dict, optional): the dictionary of all non-GKP states
            and their indices, of the form {'state': []}. By default,
            all states are GKP states.
        p_swap (float, optional): if supplied, the probability of a
            node being a p-squeezed state. Overrides the indices given
            in state.

    Attributes:
        egraph (EGraph): the unerlying graph representation.
        _N (int): the number of qubits in the lattice.
        _states (dict): states along with their indices.
        _delta (float): the delta from the Args above.
        _sampling_order (str): the sampling order from the Args above.
        _adj (array): adjancency matrix of the underlying graph.
        _random_gen (Generator): numpy random number generator.
        to_points (dict): pointer to self.egraph.to_points, the
            dictionary from indices to coordinates.
    """

    def __init__(self, g, states={"p": np.empty(0, dtype=int)}, p_swap=0):
        """Initialize the CVGraph."""
        if isinstance(g, EGraph):
            self.egraph = g
        else:
            # TODO: Make sure not to confuse subsequent references to g
            # vs those to CV.egraph.
            self.egraph = EGraph(g)
        self._N = len(g)

        # Instantiate the adjacency matrix
        self._adj = self.egraph.adj_generator(sparse=True)

        # Create a generator for random numbers to be used throughout
        self._random_gen = rng()

        if states:
            self._states = states.copy()
            # Non-zero swap-out probability overrides indices specified
            # in states and hybridizes the lattice. Print a message if
            # both supplied.
            # TODO: Raise exception?
            if p_swap:
                if len(self._states["p"]):
                    print(
                        "Both swap-out probability and indices of p-squeezed states supplied. "
                        "Ignoring the indices."
                    )
                if p_swap == 1:
                    self._states["p"] = np.arange(self._N)
                else:
                    num_p = self._random_gen.binomial(self._N, p_swap)
                    inds = self._random_gen.choice(
                        range(self._N), size=int(np.floor(num_p)), replace=False
                    )
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

    def apply_noise(self, model={}):
        """Apply noise model given in model.
        
        Args:
            model (dict, optional): the noise model dictionary of the form
                (default values displayed):
    
                {'model': 'grn', 'sampling_order': 'initial', 'delta': 0.01}
    
                grn stands for Gaussian Random Noise; sampling_order
                dictates how to simulate measurement outcomes: sample from
                an uncorrelated noise matrix initially ('initial'), a
                correlated noise matrix finally ('final'), or for ideal
                homodyne outcomes initially and from a separable noise
                covariance matrix finally ('two-step'); 'delta' is the
                quadrature blurring parameter, related to the squeezing 
                of the GKP states and the momentum-quadrature variance of 
                the p-squeezed states.

        """
        # Modelling the states.
        default_model = {"noise": "grn", "delta": 0.01, "sampling_order": "initial"}
        model = {**default_model, **model}
        self._delta = model["delta"]
        self._sampling_order = model["sampling_order"]
        if model["noise"] == "grn":
            self.grn_model()

    def grn_model(self):
        """Apply Gaussian Random Noise model to the CVGraph.

        Store quadrature or noise information as attributes depnding
        on the sampling order.
        """
        N = self._N
        delta = self._delta

        if self._sampling_order == "initial":
            init_noise = np.empty(2 * N, dtype=np.float32)
            init_vals = np.empty(2 * N, dtype=np.float32)
            for state in self._states:
                indices = self._states[state]
                if state == "GKP":
                    init_noise[indices] = delta / 2
                    init_noise[indices + N] = delta / 2
                    init_vals[indices] = self._random_gen.normal(
                        0, np.sqrt(delta / 2), len(indices)
                    )
                    init_vals[indices + N] = self._random_gen.normal(
                        0, np.sqrt(delta / 2), len(indices)
                    )
                if state == "p":
                    init_noise[indices] = 1 / (2 * delta)
                    init_noise[indices + N] = delta / 2
                    init_vals[indices] = self._random_gen.normal(
                        0, np.sqrt(1 / (2 * delta)), len(indices)
                    )
                    init_vals[indices + N] = self._random_gen.normal(
                        0, np.sqrt(delta / 2), len(indices)
                    )
            self._init_noise = init_noise
            self._init_vals = init_vals

        # For final sampling, apply a symplectic CZ matrix to the
        # initial noise covariance.
        if self._sampling_order == "final":
            init_noise = np.empty(2 * N, dtype=np.float32)
            init_vals = np.empty(2 * N, dtype=np.float32)
            for state in self._states:
                indices = self._states[state]
                if state == "GKP":
                    init_noise[indices] = delta / 2
                    init_noise[indices + N] = delta / 2
                if state == "p":
                    init_noise[indices] = 1 / (2 * delta)
                    init_noise[indices + N] = delta / 2
            noise_cov_init = sp.diags(init_noise)
            self._noise_cov = SCZ_apply(self._adj, noise_cov_init)
            # TODO: Save var_p and var_q?

        # For two-step sampling, sample for initial (ideal)
        # state-dependent quadrature values.
        if self._sampling_order == "two-step":
            self._init_quads = np.zeros(2 * N, dtype=np.float32)
            for state in self._states:
                indices = self._states[state]
                for ind in indices:
                    if state == "p":
                        self._init_quads[ind] = self._random_gen.random() * (2 * np.sqrt(np.pi))
                    if state == "GKP":
                        self._init_quads[ind] = self._random_gen.integers(0, 2) * np.sqrt(np.pi)

    def measure_hom(self, quad="p", inds=[], method="cholesky", dim="single", updated_quads=[]):
        """Conduct a homodyne measurement on the lattice.

        Simulate a homodyne measurement of quadrature quad of states
        at indices inds according to sampling order specified by
        self._sampling_order. Use the Numpy random sampling method
        method. If updated_quads is supplied, use those
        instead of applying an SCZ matrix to the initial quads in
        the two-step sampling.
        """
        N = self._N
        if not len(inds):
            inds = range(N)
        N_inds = len(inds)
        if self._sampling_order == "initial":
            means = np.zeros(2 * N, dtype=bool)
            covs = self._init_noise
            if dim == "single":
                outcomes = SCZ_apply(self._adj, means + self._init_vals)
            if quad == "q":
                outcomes = outcomes[:N][inds]
            elif quad == "p":
                outcomes = outcomes[N:][inds]
        if self._sampling_order == "two-step":
            if len(updated_quads):
                updated = updated_quads
            else:
                means = self._init_quads
                adj = self.egraph.adj_generator(sparse=True)
                updated = SCZ_apply(adj, means)
            if quad == "q":
                means = updated[:N][inds]
            elif quad == "p":
                means = updated[N:][inds]
            if dim == "single":
                outcomes = np.empty(N_inds, dtype=np.float32)
                sigma = np.sqrt(self._delta / 2)
                for i in range(N_inds):
                    outcomes[i] = self._random_gen.normal(means[i], sigma)
        if self._sampling_order == "final":
            cov_q = self._noise_cov[:N, :N]
            cov_p = self._noise_cov[N:, N:]
            cov_dict = {"q": cov_q, "p": cov_p}
            means = np.zeros(N_inds, dtype=bool)
            # TODO: Is below correct?
            covs = cov_dict[quad][inds, :][:, inds].toarray()
            outcomes = self._random_gen.multivariate_normal(mean=means, cov=covs, method=method)
        for i in range(N_inds):
            self.egraph.nodes[self.to_points[inds[i]]]["hom_val_" + quad] = outcomes[i]

    def eval_Z_probs(self, inds=[], exact=True, cond=False):
        """Evaluate the probability of phase errors at nodes inds.

        If inds not specified, compute probabilities for all nodes.
        """
        N = self._N
        if not len(inds):
            inds = range(N)
        N_inds = len(inds)
        if exact:
            if self._sampling_order != "final":
                print('Sampling order must be "final"')
                raise Exception
            var_p = self._noise_cov[N:, N:][inds, :][:, inds].toarray()
            var_p = np.diag(var_p)
            if cond:
                # TODO: Fix this.
                # TODO: Account for change in hom input for two-step sampling
                hom_vals = self.hom_outcomes(inds=inds)
                errs = Z_err_cond(var_p, hom_vals)

            else:
                errs = Z_err(var_p)
            p_string = "p_phase" + "_cond" * cond
            for i in range(N_inds):
                self.egraph.nodes[self.to_points[inds[i]]][p_string] = errs[i]
        # TODO: Fix the following to account for p-squeezed states in
        # vicinity.
        # else:
        #     for i in range(N_inds):
        #         point = self.to_points[inds[i]]
        #         n_neighbors = len(self.egraph[point])
        #         delta_effective = (n_neighbors + 1) * self._delta
        #         err = Z_err([delta_effective])[0]
        #         self.egraph.nodes[point]['p_phase'] = errs[i]

    def SCZ(self, sparse=False):
        """Return the symplectic matrix associated with CZ application.

        Returns:
            array: the symplectic matrix.
        """
        adj = self._adj
        return SCZ_mat(adj)

    def Z_probs(self, inds=[], cond=False):
        """array: the phase error probabilities of modes inds."""
        N = self._N
        if not len(inds):
            inds = range(N)
        p_string = "p_phase" + "_cond" * bool(cond)
        phase_errs = [self.egraph.nodes[self.to_points[i]].get(p_string) for i in inds]
        return phase_errs

    def hom_outcomes(self, inds=[], quad="p"):
        """array: quad-homodyne measurement outcomes for modes inds."""
        N = self._N
        if not len(inds):
            inds = range(N)
        outcomes = [self.egraph.nodes[self.to_points[i]].get("hom_val_" + quad) for i in inds]
        return outcomes

    def bit_values(self, inds=[]):
        """array: bit values associated with the p measurement."""
        N = self._N
        if not len(inds):
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


