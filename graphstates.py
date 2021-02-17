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
from numpy.random import (multivariate_normal as mvn, default_rng as rng)
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import scipy.sparse as sp
from GKP import Z_err, Z_err_cond


class EGraph(nx.Graph):
    """An enhanced graph class based on a NetworkX Graph.

    A class for adding some functionality to a NetworkX graph,
    including a drawing function draw and some short-hand/conveience
    methods.

    Attributes:
        font_props (dict): graph-size-dependent font properties for use
            in the draw method and in classes that make use of EGraph.
        indexer (str): method for indexing the nodes; 'default' for
            Python's sorted function; 'macronodes' for rounding
            micronodes to integers, sorting those, and
            furthermore sorting the micronodes within each macronodes,
            all using Python's 'sorted'.
        to_indices (dict): if self.index_generator() has been run,
            a dictionary of the form {points: indices}
        to_points (dict): if self.index_generator() has been run,
            a dictionary of the form {indixes: points}
        adj_mat (np.array): if self.adj_generator() has been run,
            the adjacency mtrix of the graph.
    """

    def __init__(self, indexer='default', *args, **kwargs):
        """Initialize an EGraph (itself an NetworkX graph)."""
        super().__init__(*args, **kwargs)
        if 'dims' in self.graph:
            tot = np.sum(self.graph['dims'])
            self.font_props = {'size': 10 * tot ** (1 / 2), 'family': 'serif'}
        # TODO: If dims not specified, look at number of nodes in graph
        # to detemrine font properties.
        self.indexer = indexer
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
        if self.indexer == 'default':
            ind_dict = dict(zip(sorted(self.nodes()), range(N)))
        if self.indexer == 'macronodes':
            self.macro = nx.Graph()
            for node in self.nodes():
                rounded = tuple(np.round(node).astype(int))
                self.macro.nodes[rounded]['micronodes'].append(node)
            sorted_macro = sorted(self.macro)
            points = []
            for vertex in sorted_macro:
                points += self.macro.nodes[vertex]['micronodes']
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
        return adj

    def draw(self, color_nodes=False, color_edges=False, label=False):
        """Draw the graph.

        Args:
            color_nodes (bool): If True, color nodes based on 'color'
                attributes attached to the node. Black by default.
            color_edges (bool): If True, color edges based on 'color'
                attributes attached to the node. Grey by default.
            label (bool): if True, label the indices as per
                self.index_generator; unlabelled by default.

        Returns:
            A matplotib Axes object.
        """
        # Recommended to be viewed with IPython.
        if self.graph['dims']:
            nx, ny, nz = self.graph['dims']
        else:
            nx, ny, nz = 5, 5, 5
        fig = plt.figure(figsize=(2 * (nx + ny + nz + 2),
                                  2 * (nx + ny + nz + 2)))
        ax = fig.add_subplot(111, projection='3d')
        # Plotting points. y and z are swapped in the loops so that
        # z goes into the page; however, the axes labels are correct.
        for point in self.nodes:
            x, z, y = point

            # Color based on color attribute, if available; default if
            # unavailable; black if color_nodes is False
            color = self.nodes[point].get('color') if color_nodes else 'k'

            ax.scatter(x, y, z, s=70, c=color)
            if label:
                indices = self.index_generator()
                ax.text(x, y, z, str(indices[point]), fontdict=self.font_props,
                        color='MediumBlue', backgroundcolor='w')
        # Plotting edges.
        for edge in self.edges:

            # Color based on color attribute, if available; default if
            # unavailable; black if color_edges is False
            color = self.edges[edge].get('color') if color_edges else 'k'

            x1, z1, y1 = edge[0]
            x2, z2, y2 = edge[1]
            plt.plot([x1, x2], [y1, y2], [z1, z2], c=color)

        ax.tick_params(labelsize=self.font_props['size'])
        plt.xticks(range(0, 2 * nx + 1))
        plt.yticks(range(0, 2 * nz + 1))
        ax.set_zticks(range(0, 2 * ny + 1))
        ax.set_xlabel('x', fontdict=self.font_props, labelpad=20)
        ax.set_ylabel('z', fontdict=self.font_props, labelpad=20)
        ax.set_zlabel('y', fontdict=self.font_props, labelpad=20)
        plt.rcParams['grid.color'] = "lightgray"
        plt.tight_layout(pad=3)
        plt.draw()
        return ax


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
        # TODO: Specify kind of Scip y sparse matrix?
        identity = sp.identity(N, dtype=np.int8)
        zeros = sp.csr_matrix((N, N), dtype=np.int8)
        block_func = sp.bmat
    # Construct symplectic
    symplectic = block_func([
        [identity, zeros],
        [adj, identity]])
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
            new_quads = np.block([
                [c1, block2],
                [block3, block4]])
    return new_quads


class CVGraph:
    """A class for representing continuous-variable graph states.

    Has all the functionality of an EGraph, but associates its
    nodes with continuous-variable quantum states and its edges with
    continuous-variable CZ gates.

    For now, only a hybrid state of p-squeezed and GKP states is
    considered.

    Args:
        g (graph-type): the graph underlying the state
        model (str): the noise model; Gaussian random noise (GRN) by
            default
        p_inds (array): the indices of the p-squeezed states; empty by
            default (meaning only GKP states present)
        swap_prob (float): if supplied, the probability of a node being
            a p-squeezed state
        delta (float): the quadrature blurring parameter, related to
            the squeezing of the GKP states and the momentum-quadrature
            variance of the p-squeezed states; 0.01 by default

    Attributes:
        graph (EGraph): the unerlying graph representation
        _N (int): the number of qubits in the lattice
        _delta (float): the delta from the Args above
        _SCZ (np.array): the symplectic matrix associated with the CZ
            application at the the lattice edges
        _p_inds (array): the indices of the p-squeezed states
        _indexed_graph (nx.Graph): a graph with nodes relabelled as
            integer indices
        ind_dict(dict): a dictionary between indices and coordinates
        cov_p (array): the phase-space covariance matrix (populated
            if 'grn' model specified)
        var_ p (array): the p variances of the modes populated if
            'grn' model specified)
        Z_probs (array): if eval_z_probs has been run, the phase error
            probabilities of the modes
        Z_probs_cond (array): if eval_z_probs_cond has been run, the
            phase error probabilities of the modes conditioned on
            homodyne outcomes
        hom_outcomes (array): if measure_p has been run, the results
            of the latest p-homodyne measurement
    """

    def __init__(self, g, model='grn', p_inds=[], swap_prob=None, delta=0.01):
        """Initialize the CVGraph.

        Convert g into an EGraph and designate nodes either p-squeezed
        or GKP states depending. If p_inds is not empty, manually
        populate those nodes with p-squeezed states. If swap_prob is
        given, run the hybridize method to populate the lattice at
        random with a numer of p-squeezed states equal to swap_prob *
        (# of nodes).
        """
        if isinstance(g, EGraph):
            self.graph = g
        else:
            self.graph = EGraph(g)
        self._SCZ = SCZ_mat(self.graph.adj_mat())
        self._delta = delta
        self._N = self.graph.number_of_nodes()

        if swap_prob:
            self.hybridize(swap_prob)
        else:
            self._p_inds = p_inds

        self._indexed_graph = self.graph.index()
        idg = self._indexed_graph
        self.ind_dict = {n: idg.nodes[n]['pos'] for n in idg.nodes}
        for ind in self._p_inds:
            self.graph.nodes[self.ind_dict[ind]]['type'] = 'p'

        for ind in set(idg.nodes).difference(self._p_inds):
            self.graph.nodes[self.ind_dict[ind]]['type'] = 'GKP'

        if model == 'grn':
            self.grn_model(delta)

    def SCZ(self, heat_map=0):
        """Return the symplectic matrix associated with CZ application.

        Args:
            heat_map (bool): if True, draw a heat map of the matrix.

        Returns:
            array: the symplectic matrix.
        """
        if heat_map:
            print('The symplectic CZ matrix (dark spots 0, bright spots 1):')
            plt.figure()
            plt.matshow(self._SCZ, 0)
            plt.show()
        return self._SCZ

    def eval_Z_probs(self):
        """Evaluate the probability of phase errors at each mode."""
        errs = Z_err(self.var_p)
        for i in range(len(errs)):
            self.graph.nodes[self.ind_dict[i]]['p_phase'] = errs[i]

    def eval_Z_probs_cond(self):
        """Evaluate the conditional phase error probability of each mode."""
        try:
            hom_vals = self.hom_outcomes
        except Exception:
            return

        Z_errs = Z_err_cond(self.var_p, hom_vals)
        for i in range(len(Z_errs)):
            self.graph.nodes[self.ind_dict[i]]['p_phase_cond'] = Z_errs[i]

    def measure_p(self):
        """Conduct a p-homodyne measurement on the lattice."""
        # Randomly sample of a normal distribution centred at 0 with
        # covariances matrix given by cov_p.
        outcomes = mvn(mean=np.zeros(self._N), cov=self._cov_p)
        # bit_values = self.translator(outcomes)
        for i in range(len(outcomes)):
            self.graph.nodes[self.ind_dict[i]]['hom_val'] = outcomes[i]

    def hybridize(self, swap_prob):
        """Populate nodes with p-squeezed states at random.

        The number of p states is the sample of the binomial
        distribution with the number of trials equalling the size
        of the graph state, and a success probability swap_prob.

        Args:
            swap_prob (float): the probability of swapping out a GKP
                state for a p-squeezed state.

        Returns:
            None
        """
        num_p = rng().binomial(self._N, swap_prob)
        self._p_inds = rng().choice(range(self._N), size=int(np.floor(num_p)), replace=False)

    def grn_model(self, delta):
        """Apply Gaussian Random Noise model to the CVGraph.

        Args:
            delta (float): the sqeezing/blurring/variance parameter

        Returns:
            None
        """
        # Step 1: Construct the phase-space covariance matrix
        # Basis ordering: all N q's, then all N p's
        # Step 1a: initialize all as GKPs
        cov_phase = (delta / 2) * np.eye(2 * self._N)
        # Step 1b: replace p in relevant locations
        cov_phase[self._p_inds, self._p_inds] = 1 / (2 * delta)
        # Step 1c: apply CZ gates
        cov_phase = self._SCZ @ cov_phase @ self._SCZ.T
        # Step 1d: extract p variances
        self._cov_p = cov_phase[self._N:, self._N:]
        self.var_p = np.diag(self._cov_p)
        for i in range(self._N):
            self.graph.nodes[self.ind_dict[i]]['var_p'] = self.var_p[i]

    # Note that only the getter function has been defined for the properties
    # below because I am treating these as private. This can be changed if we
    # would like the user to be able to change the values of these attribute
    # for a given CVGraph object.
    @property
    def p_inds(self):
        """array: the indices of the p-squeezed states."""
        return self._p_inds

    @property
    def cov_p(self):
        """array: the phase-space covariance matrix."""
        return self._cov_p

    @property
    def Z_probs(self):
        """array: the phase error probabilities of the modes."""
        try:
            phase_errs = [self.graph.nodes[node]['p_phase'] for node in self.graph]
            return phase_errs
        except Exception:
            print('Z error probabilities have not yet been computed. Please '
                  'use eval_Z_probs() first.')
            return

    @property
    def Z_probs_cond(self):
        """array: the conditional phase error probabilities of the modes."""
        try:
            phase_errs = [self.graph.nodes[node]['p_phase_cond'] for node in self.graph]
            return phase_errs
        except Exception:
            print('Conditional Z error probabilities have not yet been computed. Please '
                  'use eval_Z_probs_cond() first.')
            return

    @property
    def hom_outcomes(self):
        """array: the results of the p-homodyne measurement."""
        try:
            outcomes = [self.graph.nodes[node]['hom_val'] for node in self.graph]
            return outcomes
        except Exception:
            print('A homodyne measurement has not yet been performed. Please '
                  'use measure_p() first.')
            raise

    @property
    def bit_values(self):
        """array: bit values associated with the p measurement."""
        try:
            bit_values = [self.graph.nodes[node]['bit_val'] for node in self.graph]
            return bit_values
        except Exception:
            print('A homodyne measurement or a translation has not yet'
                  'been performed.')
            return

    def sketch(self, label=None, legend=True, title=True, color_nodes=False, color_edges=False):
        """Sketch the underlying graph with CV labels.

        GKP states are black, p-squeezed states are orange. If a label
        is specified from the keys in the following dictionary, data
        corresponding to the values will be displayed:

            {'var_p': 'p Variances',
            'p_phase': 'Phase error probabilities',
            'p_phase_cond': 'Conditional phase error probabilities',
            'hom_val': 'p-homodyne outcomes',
            'bit_val': 'Bit values'}.

        Args:
            label (bool): if a string from the above dictionary, plot
                the correspond value; otherwise, no label.
            legend (bool): if True, display a legend
            title (bool): if True, display a title
            color_nodes (bool): If True, color nodes based on 'color'
                attributes attached to the node. Black by default.
            color_edges (bool): If True, color edges based on 'color'
                attributes attached to the node. Grey by default.
        Returns:
            A matplotib Axes object.
        """
        font_props = self.graph.font_props
        idg = self._indexed_graph
        p_coords = [idg.nodes[ind]['pos'] for ind in self._p_inds]
        ax = self.graph.draw(color_nodes=color_nodes, color_edges=color_edges, label=0)
        for point in p_coords:
            x, z, y = point
            ax.scatter(x, y, z, s=40, c='orange')
        orange_line = mlines.Line2D([], [], color='orange', marker='.',
                                    markersize=14, label='p')
        black_line = mlines.Line2D([], [], color='black', marker='.',
                                   markersize=14, label='GKP')
        if legend:
            ax.legend(handles=[orange_line, black_line], prop=font_props)

        title_dict = {'var_p': 'p Variances',
                      'p_phase': 'Phase error probabilities',
                      'p_phase_cond': 'Conditional phase error probabilities',
                      'hom_val': 'p-homodyne outcomes',
                      'bit_val': 'Bit values'}

        if label:
            try:
                name = title_dict[label]
            except KeyError:
                print('No such label permitted.')

            if title:
                ax.set_title(name, fontdict=font_props)

            for node in self.graph:
                try:
                    value = self.graph.nodes[node][label]
                except KeyError:
                    print(title + ' have not yet been computed.')
                    return
                x, z, y = node
                ax.text(x, y, z, '{:.0g}'.format(value),
                        fontdict=font_props, color='MediumBlue',
                        backgroundcolor='w')

        plt.draw()
        return ax


if __name__ == '__main__':
    pass
