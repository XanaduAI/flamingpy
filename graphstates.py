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
'''This module provides classes for DV and CV graph states.'''
import itertools as it
import networkx as nx
import numpy as np
from numpy.random import (multivariate_normal as mvn, default_rng as rng)
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from GKP import basic_translate, Z_err, Z_err_cond

class EGraph(nx.Graph):
    '''An enhanced graph class based on networkx.Graph.'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if 'dims' in self.graph:
            tot = np.sum(self.graph['dims'])
            self.font_props = {'size': 10 * tot ** (1 / 2), 'family': 'serif'}

    def adj_mat(self):
        return nx.to_numpy_array(self, nodelist=sorted(self.nodes))


    def draw(self, color_nodes=False, color_edges=False, label=False):
        """Draw the graph.

        Args:
            color_nodes (bool): If True, color nodes based on 'color' 
                attributes attached to the node. Black by default.
            color_edges (bool): If True, color edges based on 'color'
                attributes attached to the node. Grey by default.
            label (bool): if True, label the indices; unlabelled by
                default.
            polarity (bool): if True, color edges with +1 weights blue
                and -1 weights red. All edges grey by default.

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
        # Plotting points. y and z are swapped in the loops to compare
        # the lattice diagrams from our notes, where the z axis goes
        # into the page; however the axes labels are correct.
        for point in self.nodes:
            x, z, y = point
            color_bool = int(color_nodes)
            if 'color' in self.nodes[point]:
                node_color = self.nodes[point]['color']
            ax.scatter(x, y, z, s=70, c=color_bool*node_color+(1-color_bool)*'k')
            indices = {c: n for (n, c) in enumerate(sorted(self.nodes))}
            if label:
                ax.text(x, y, z, str(indices[point]), fontdict=self.font_props,
                        color='MediumBlue', backgroundcolor='w')
        # Plotting edges.
        for edge in self.edges:
            color = int(color_edges)
            if 'color' in self.edges[edge]:
                edge_color = self.edges[edge]['color']
            x1, z1, y1 = edge[0]
            x2, z2, y2 = edge[1]
            plt.plot([x1, x2], [y1, y2], [z1, z2], c=color_bool*edge_color+(1-color_bool)*'grey')

        ax.tick_params(labelsize=self.font_props['size'])
        plt.xticks(range(0, 2*nx + 1))
        plt.yticks(range(0, 2*nz + 1))
        ax.set_zticks(range(0, 2*ny + 1))
        ax.set_xlabel('x', fontdict=self.font_props, labelpad=20)
        ax.set_ylabel('z', fontdict=self.font_props, labelpad=20)
        ax.set_zlabel('y', fontdict=self.font_props, labelpad=20)
        plt.rcParams['grid.color'] = "lightgray"
        plt.tight_layout(pad=3)
        plt.draw()
        return ax

    def index(self):
        indexed_graph = nx.convert_node_labels_to_integers(self, ordering='sorted', label_attribute='pos')
        return indexed_graph


def SCZ_mat(adj):
    """Return a symplectic matrix corresponding to CZ gate application.

    Gives the 2N by 2N symplectic matrix for CZ gate application
    based on the adjacency matrix adj.

    Args:
        adj (array): N by N binary symmetric matrix. If modes i and j
            are linked by a CZ, then entry ij and ji is 1; otherwise 0.
    Returns:
        array: 2N by 2N symplectic matrix
    """
    # Number of modes
    N = len(adj)
    # Construct symplectic
    symplectic = np.block([[np.eye(N), np.zeros((N, N))], [adj, np.eye(N)]])
    return symplectic


class CVGraph:
    '''A class for representing continuous-variable graph states.

    Has all the functionality of an EGraph, but associates its
    nodes with continuous-variable quantum states and its edges with
    continuous-variable CZ gates.

    Args:
        g (graph-type): the graph underlying the state
        model (str): the error model; gaussian random noise (GRN) by
            default
        p_inds (array): the indices of the p-squeezed states; empty by
            default (meaning only GKP states present)
        swap_prob (float): if supplied, the probability of a node being
            a p-squeezed state
        delta (float): the quadrature blurring parameter, related to
            the squeezing of the GKP states and the momentum-quadrature
            variance of the p-squeezed states; 0.01 by default

    Attributes:
        dims (int or tuple): the size of the lattice, as above
        coords (set): the coordinates of the lattice
        indices (dict): a dictionary of the form {index: coordinate},
            where the indices are derived from indfunc. Note that this
            swaps the order of the dictionary output by indfunc
        N (int): the number of qubits in the lattice
        SCZ (np.array): the symplectic matrix associated with the CZ
            application at the the lattice edges. Note that this is 
            technically a public method that can also draw SCZ.
        p_indices (array): the indices of the p-squeezed states
        cov_p (array): the phase-space covariance matrix
        p_var (array): the p variances of the modes
        Z_probs (array): if eval_z_probs has been run, the phase error
            probabilities of the modes
        Z_probs_cond (array): if eval_z_probs_cond has been run, the 
            phase  error probabilities of the modes conditioned on 
            homodyne outcomes
        hom_outcomes (array): if measure_p has been run, the results
            of the latest p-homodyne measurement
        bit_values (array): if measure_p has been run, the bit values
            associated with the latest p-homodyne measurement, assuming
            binning function translator.
    '''
    def __init__(self, g, model='grn', p_inds=[], swap_prob=None,
                 delta=0.01):
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
        """Populate nodes with p-squeezed states at random."""
        num_p = swap_prob * self._N
        self._p_inds = rng().choice(range(self._N), size=int(np.floor(num_p)), replace=False)

    def grn_model(self, delta):
        # Step 1: Construct the phase-space covariance matrix
        # Basis ordering: all N q's, then all N p's
        # Step 1a: initialize all as GKPs
        cov_phase = delta * np.eye(2*self._N)
        # Step 1b: replace p in relevant locations
        cov_phase[self._p_inds, self._p_inds] = 1/delta
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
        """Sketch the CVRHG lattice. GKP states black, p states orange.

        Args:
            label (bool): if 'v', label the p variances; if 'e', the
                phase error probabilities; if 'h', the p-homodyne
                outcomes; if 'b', the bit values; otherwise, no label.
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
