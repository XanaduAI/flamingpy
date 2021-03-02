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
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from numpy.random import default_rng as rng

import scipy.sparse as sp
from GKP import Z_err, Z_err_cond


class EGraph(nx.Graph):
    """An enhanced graph based on a NetworkX Graph.

    A class for adding some functionality to a NetworkX graph,
    including a drawing function draw and some short-hand/convenience
    methods.

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

    def __init__(self, indexer='default', *args, **kwargs):
        """Initialize an EGraph (itself an NetworkX graph)."""
        super().__init__(*args, **kwargs)

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
        # TODO: Heat map?
        return adj

    def draw(self, color_nodes=False, color_edges=False, label_indices=False, display_axes=True):
        """Draw the graph.

        Args:
            color_nodes (bool): If True, color nodes based on 'color'
                attributes attached to the node. Black by default.
            color_edges (bool): If True, color edges based on 'color'
                attributes attached to the node. Grey by default.
            label_indices (bool): if True, label the indices as per
                self.index_generator; unlabelled by default.
            display_axes (bool): if False, turn off the axes.

        Returns:
            A matplotib Axes object.
        """
        # Recommended to be viewed with IPython.
        # Font properties
        dims = self.graph.get('dims')
        if dims:
            # TODO: Store dims as EGraph attributes, rather than a graph
            # attribute?
            xmax, ymax, zmax = dims
            font_size = 10 * sum(dims) ** (1 / 2)
        else:
            # TODO: If dims not specified find x, y, z limits of graph,
            # supposing the graph is filled. Alternatively, just change
            # the figure size?
            xmax, ymax, zmax = 5, 5, 5
            font_size = 14

        # Set plotting options
        plot_params = {'font.size': font_size, 'font.family': 'serif',
                       'axes.labelsize': font_size, 'axes.titlesize': font_size,
                       'xtick.labelsize': font_size, 'ytick.labelsize': font_size,
                       'legend.fontsize': font_size, 'grid.color': 'lightgray',
                       'lines.markersize': font_size}
        plt.rcParams.update(plot_params)

        fig = plt.figure(figsize=((2 * (sum(dims) + 2), 2 * (sum(dims) + 2))))
        ax = fig.add_subplot(111, projection='3d')

        # Plotting points. y and z are swapped in the loops so that
        # z goes into the page; however, the axes labels are correct.
        for point in self.nodes:
            x, z, y = point

            # Color nodes based on color_nodes if string, or based on
            # color attribute if True; black otherwise.
            if type(color_nodes) == str:
                color = color_nodes
            else:
                color = self.nodes[point].get('color') if color_nodes else 'k'

            ax.scatter(x, y, z, c=color, s=plt.rcParams['lines.markersize'] * 5)
            if label_indices:
                indices = self.index_generator()
                ax.text(x, y, z, str(indices[point]), backgroundcolor='w', zorder=2)

        # Plotting edges.
        for edge in self.edges:

            # Color edges based on color_edges if string, or based on
            # color ttribute if True; black otherwise.
            if type(color_edges) == str:
                color = color_edges
            else:
                color = self.edges[edge].get('color') if color_edges else 'k'

            x1, z1, y1 = edge[0]
            x2, z2, y2 = edge[1]
            plt.plot([x1, x2], [y1, y2], [z1, z2], color=color)

        plt.xticks(range(0, 2 * xmax + 1))
        plt.yticks(range(0, 2 * zmax + 1))
        ax.set_zticks(range(0, 2 * ymax + 1))
        ax.set_xlabel('x', labelpad=15)
        ax.set_ylabel('z', labelpad=15)
        ax.set_zlabel('y', labelpad=15)
        if not display_axes:
            ax.axis('off')
        plt.tight_layout(pad=5)
        plt.draw()
        return ax


def SCZ_mat(adj, heat_map=False):
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
    if heat_map:
        print('The symplectic CZ matrix (dark spots 0, bright spots 1):')
        plt.figure()
        if type(symplectic) != np.ndarray:
            symplectic = symplectic.toarray()
        plt.matshow(symplectic, 0)
        plt.show()
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
        state (dict, optional): the dictionary of all non-GKP states
            and their indices, of the form {'state': []}
        model (dict, optional): the noise model dictionary of the form
            (default values displayed):

            {'model': 'grn', 'sampling_order': 'initial', 'delta': 0.01}

            grn stands for Gaussian Random Noise; sampling_order
            dictates how to simulate measurement outcomes: sample from
            an uncorrelated noise matrix initially ('initial'), a
            correlated noise matrix finally ('final'), or for ideal
            homodyne outcomes initially and from a separable noise
            covariance matrix finally ('two-step'); 'delta' is the
            noise parameter, explained below.

        swap_prob (float, optional): if supplied, the probability of a
            node being a p-squeezed state. Overrides the indices given
            in state.
        delta (float, optional): the quadrature blurring parameter,
            related to the squeezing of the GKP states and the
            momentum-quadrature variance of the p-squeezed states;
            0.01 by default. Overrides the delta given in 'model'.

    Attributes:
        egraph (EGraph): the unerlying graph representation
        _N (int): the number of qubits in the lattice
        _states (dict): states along with their indices
        _delta (float): the delta from the Args above
        _sampling_order (str): the sampling order from the Args above
        to_points (dict): pointer to self.egraph.to_points, the
            dictionary from indices to coordinates
        Z_probs (array): if eval_z_probs has been run, the phase error
            probabilities of the modes
        Z_probs_cond (array): if eval_z_probs_cond has been run, the
            phase error probabilities of the modes conditioned on
            homodyne outcomes
        hom_outcomes (array): if measure_hom has been run, the results
            of the latest p-homodyne measurement
    """

    def __init__(self, g, states={'p': np.empty(0, dtype=int)}, model={}, p_swap=0, delta=None):
        """Initialize the CVGraph.

        Convert g into an EGraph and designate nodes either p-squeezed
        or GKP states depending. If p_inds is not empty, manually
        populate those nodes with p-squeezed states. If swap_prob is
        given, run the hybridize method to populate the lattice at
        random with a numer of p-squeezed states equal to swap_prob *
        (# of nodes).
        """
        if isinstance(g, EGraph):
            self.egraph = g
        else:
            self.egraph = EGraph(g)
        self._N = len(g)

        self._states = states
        if states:
            # Non-zero swap-out probability overrides indices specified
            # in states and hybridizes the lattice. Print a message if
            # both supplied.
            # TODO: Raise exception?
            if p_swap:
                if len(states['p']):
                    print('Both swap-out probability and indices of p-squeezed states supplied. '
                          'Ignoring the indices.')
                if p_swap == 1:
                    states['p'] = np.arange(self._N)
                else:
                    num_p = rng().binomial(self._N, p_swap)
                    inds = rng().choice(range(self._N), size=int(np.floor(num_p)), replace=False)
                    states['p'] = inds

            # Associated remaining indices with GKP states.
            used_inds = np.empty(0, dtype=int)
            for psi in states:
                used_inds = np.concatenate([used_inds, states[psi]])
            remaining_inds = list(set(range(self._N)).difference(set(used_inds)))
            states['GKP'] = np.array(remaining_inds, dtype=int)

            # Generate EGraph indices.
            self.egraph.index_generator()
            self.to_points = self.egraph.to_points

            for psi in states:
                for ind in states[psi]:
                    self.egraph.nodes[self.to_points[ind]]['state'] = psi

            # Modelling the states.
            # If both delta and delta in model specified, print a
            # message that the former will be used.
            if delta and model.get('delta'):
                print('Delta supplied twice. Using the delta given by the delta argument.')
                model['delta'] = delta
            default_model = {'noise': 'grn', 'delta': 0.01, 'sampling_order': 'initial'}
            model = {**default_model, **model}
            self._delta = model['delta']
            self._sampling_order = model['sampling_order']
            if model['noise'] == 'grn':
                self.grn_model()

    def grn_model(self):
        """Apply Gaussian Random Noise model to the CVGraph.

        Store quadrature or noise information as attributes depnding
        on the sampling order.
        """
        N = self._N
        delta = self._delta

        # For initial and final sampling, generate a sparse adjacency
        # matrix in the EGraph and state-dependent noise vectors.
        if self._sampling_order in ('initial', 'final'):
            self._adj = self.egraph.adj_generator(sparse=True)
            self._init_noise = np.empty(2 * N, dtype=np.float32)
            for state in self._states:
                indices = self._states[state]
                if state == 'GKP':
                    self._init_noise[indices] = delta / 2
                    self._init_noise[indices + N] = delta / 2
                if state == 'p':
                    self._init_noise[indices] = delta / 2
                    self._init_noise[indices + N] = 1 / (2 * delta)

        # For final sampling, apply a symplectic CZ matrix to the
        # initial noise covariance.
        if self._sampling_order == 'final':
            noise_cov_init = sp.diags(self._init_noise)
            self._noise_cov = SCZ_apply(self._adj, noise_cov_init)
            # TODO: Save var_p and var_q?

        # For two-step sampling, sample for initial (ideal)
        # state-dependent quadrature values.
        if self._sampling_order == 'two-step':
            self._init_quads = np.zeros(2 * N, dtype=np.float32)
            for state in self._states:
                indices = self._states[state]
                for ind in indices:
                    if state == 'p':
                        self._init_quads[ind] = np.random.random() * (2 * np.sqrt(np.pi))
                    if state == 'GKP':
                        self._init_quads[ind] = np.random.randint(0, 2) * np.sqrt(np.pi)

    def measure_hom(self, quad='p', inds=[], method='cholesky', dim='single', updated_quads=[]):
        """Conduct a homodyne measurement on the lattice.

        Simulate a homodyne measurement of quadrature quad of states
        at indices inds according to sampling order specified by
        self._sampling_order. Use the Numpy random sampling method
        method; by default, do many single-variable samplings where
        appropriate, otherwise set dim = 'multi' for sampling once
        from a multivariate distribution. (This might affect the speed
        of implementation). If updated_quads is supplied, use those
        instead of applying an SCZ matrix to the initial quads in
        the two-step sampling.
        """
        N = self._N
        if not len(inds):
            inds = range(N)
        N_inds = len(inds)
        if self._sampling_order == 'initial':
            means = np.zeros(2 * N, dtype=bool)
            covs = self._init_noise
            if dim == 'single':
                initial = np.empty(2 * N, dtype=np.float32)
                for i in range(2 * N):
                    initial[i] = rng().normal(means[i], np.sqrt(covs[i]))
                outcomes = SCZ_apply(self._adj, initial)
            if dim == 'multi':
                initial = rng().multivariate_normal(mean=means, cov=np.diag(covs), method=method)
                outcomes = SCZ_apply(self._adj, initial)
            if quad == 'q':
                outcomes = outcomes[:N][inds]
            elif quad == 'p':
                outcomes = outcomes[N:][inds]
        if self._sampling_order == 'two-step':
            if len(updated_quads):
                updated = updated_quads
            else:
                means = self._init_quads
                adj = self.egraph.adj_generator(sparse=True)
                updated = SCZ_apply(adj, means)
            if quad == 'q':
                means = updated[:N][inds]
            elif quad == 'p':
                means = updated[N:][inds]
            if dim == 'single':
                outcomes = np.empty(N_inds, dtype=np.float32)
                sigma = np.sqrt(self._delta / 2)
                for i in range(N_inds):
                    outcomes[i] = rng().normal(means[i], sigma)
            elif dim == 'multi':
                covs = np.eye(N_inds, dtype=np.float32) * (self._delta / 2)
                outcomes = rng().multivariate_normal(mean=means, cov=covs, method=method)
        if self._sampling_order == 'final':
            cov_q = self._noise_cov[:N, :N]
            cov_p = self._noise_cov[N:, N:]
            cov_dict = {'q': cov_q, 'p': cov_p}
            means = np.zeros(N_inds, dtype=np.bool)
            # TODO: Is below correct?
            covs = cov_dict[quad][inds, :][:, inds].toarray()
            outcomes = rng().multivariate_normal(mean=means, cov=covs, method=method)
        for i in range(N_inds):
            self.egraph.nodes[self.to_points[inds[i]]]['hom_val_' + quad] = outcomes[i]

    def eval_Z_probs(self, inds=[], exact=True, cond=False):
        """Evaluate the probability of phase errors at nodes inds.

        If inds not specified, compute probabilities for all nodes.
        """
        N = self._N
        if not len(inds):
            inds = range(N)
        N_inds = len(inds)
        if exact:
            if self._sampling_order != 'final':
                print('Sampling order must be "final"')
                raise Exception
            var_p = self._noise_cov[N:, N:][inds, :][:, inds].toarray()
            var_p = np.diag(var_p)
            if cond:
                # TODO: Fix this.
                hom_vals = self.hom_outcomes(inds=inds)
                errs = Z_err_cond(var_p, hom_vals)

            else:
                errs = Z_err(var_p)
            p_string = 'p_phase' + '_cond' * cond
            for i in range(N_inds):
                self.egraph.nodes[self.to_points[inds[i]]][p_string] = errs[i]
        # TODO: Account for change in hom input for two-step sampling
        # TODO: Fix the following to account for p-squeezed states in
        # vicinity.
        # else:
        #     for i in range(N_inds):
        #         point = self.to_points[inds[i]]
        #         n_neighbors = len(self.egraph[point])
        #         delta_effective = (n_neighbors + 1) * self._delta
        #         err = Z_err([delta_effective])[0]
        #         self.egraph.nodes[point]['p_phase'] = errs[i]

    def SCZ(self, sparse=False, heat_map=False):
        """Return the symplectic matrix associated with CZ application.

        Args:
            heat_map (bool): if True, draw a heat map of the matrix.

        Returns:
            array: the symplectic matrix.
        """
        adj = self.egraph.adj_generator(sparse=sparse)
        return SCZ_mat(adj, heat_map=heat_map)

    def Z_probs(self, inds=[], cond=False):
        """array: the phase error probabilities of modes inds."""
        N = self._N
        if not len(inds):
            inds = range(N)
        p_string = 'p_phase' + '_cond' * bool(cond)
        phase_errs = [self.egraph.nodes[self.to_points[i]].get(p_string) for i in inds]
        return phase_errs

    def hom_outcomes(self, inds=[], quad='p'):
        """array: quad-homodyne measurement outcomes for modes inds."""
        N = self._N
        if not len(inds):
            inds = range(N)
        outcomes = [self.egraph.nodes[self.to_points[i]].get('hom_val_' + quad) for i in inds]
        return outcomes

    def bit_values(self, inds=[]):
        """array: bit values associated with the p measurement."""
        N = self._N
        if not len(inds):
            inds = range(N)
        bits = [self.egraph.nodes[self.to_points[i]].get('bit_val') for i in inds]
        return bits

    @property
    def p_inds(self):
        """array: the indices of the p-squeezed states."""
        return self._states.get('p')

    @property
    def GKP_inds(self):
        """array: the indices of the GKP states."""
        return self._states.get('GKP')

    @property
    def noise_cov(self):
        """array: the noise covariance matrix."""
        if self._sampling_order == 'final':
            return self._noise_cov
        else:
            print('Sampling order must be "final."')

    def draw(self, label=None, title=True, legend=True, state_colors={}, **kwargs):
        """Draw the underlying graph with CV labels.

        Args:
            label (str, optional): next to each node, plot the value of
                the node attribute label if it exists and has been
                computed. Some possible options:

                {'p_phase': 'Phase error probabilities',
                 'p_phase_cond': 'Conditional phase error probabilities',
                 'hom_val_p': 'p-homodyne outcomes',
                 'hom_val_q': 'q-homodyne outcomes',
                 'bit_val': 'Bit values',
                 'weight': 'Weights'}

            title (bool, optional): if True, display a title
                corresponding to the label
            legend (bool, optional): if True, display a legend for
                the states
            state_color (dict, optional): dictionary of the form
                {state_name: color}; if supplied overrides the
                default color cycle for state colors.
            **kwargs: The EGraph.draw arguments. Except for color_edges,
                overwritten by the CVGraph.draw arguments.
        Returns:
            A matplotib Axes object.
        """
        ax = self.egraph.draw(**kwargs)
        handles = []
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        i = 0
        for state in self._states:
            coords = [self.to_points[ind] for ind in self._states[state]]
            color = state_colors.get(state)
            if not color:
                color = color_cycle[i]
            for point in coords:
                x, z, y = point
                ax.scatter(x, y, z, s=plt.rcParams['lines.markersize'] * 5, c=color)
            line = mlines.Line2D([], [], color=color, marker='.', label=state)
            handles += [line]
            i += 1

        if legend:
            ax.legend(handles=handles)

        title_dict = {'p_phase': 'Phase error probabilities',
                      'p_phase_cond': 'Conditional phase error probabilities',
                      'hom_val_p': 'p-homodyne outcomes',
                      'hom_val_q': 'q-homodyne outcomes',
                      'bit_val': 'Bit values',
                      'weight': 'Weights'}
        if label:
            name = title_dict.get(label) if title_dict.get(label) else label

            if title:
                ax.set_title(name)

            n_uncomputed = 0
            for node in self.egraph:
                value = self.egraph.nodes[node].get(label)
                if value:
                    x, z, y = node
                    # Raise negative sign above node.
                    sign = '^{-}' * (-int(np.sign(value)))
                    number = r'${}{:.2g}$'.format(sign, np.abs(value))
                    ax.text(x, y, z, number, color='MediumBlue', backgroundcolor='w')
                else:
                    n_uncomputed += 1
            if n_uncomputed > 0:
                message = '{} at {} node(s) have not yet been computed.'
                print(message.format(name.lower(), n_uncomputed))
        plt.draw()
        return ax


if __name__ == '__main__':
    # Bell state EGraph
    edge = [(0, 0, 0), (1, 1, 1)]
    dims = (1, 1, 1)
    bell_state = EGraph(dims=dims)
    bell_state.add_edge(*edge, color='MidnightBlue')
    # Plot the bell state
    bell_state.draw(color_nodes='magenta', color_edges=True, label_indices=True)
    bell_state.adj_generator(sparse=True)
    print('Adjacency matrix: \n', bell_state.adj_mat, '\n')

    # Noise model for CVGraph
    model = {'noise': 'grn', 'delta': 1, 'sampling_order': 'final'}
    CVbell = CVGraph(bell_state, model=model, p_swap=0)
    CVbell.measure_hom('p', [9])
    CVbell.measure_hom('q', [1])
    CVbell.eval_Z_probs(cond=False)
    CVbell.draw('hom_val_p')
    CVbell.draw('hom_val_q')
    CVbell.draw('p_phase')
    print('Nodes :', bell_state.nodes.data())
    print('Edges :', bell_state.edges.data())
    print('p indices: ', CVbell.p_inds, '\n')
    print('GKP indices: ', CVbell.GKP_inds, '\n')
    print('Symplectic CZ matrix: \n', CVbell.SCZ(heat_map=True), '\n')
