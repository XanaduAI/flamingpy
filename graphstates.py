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
from scipy.special import erf

# Style for figures
fontdict = {'fontsize': 14, 'family': 'serif'}


class EGraph(nx.Graph):
    '''An enhanced graph class based on networkx.Graph.'''

    def __init__(self):
        nx.Graph.__init__(self)
        if 'dims' not in self.graph:
            self.graph['dims'] = None

    def adj(self):
        return nx.to_numpy_array(self)

    def color(self, coord):
        return self.nodes[coord]['color']

    def draw(self, color=1, label=0):
        '''Draw the graph.
        Args:
            label (bool): if True, label the indices; unlabelled by
                default.
            color (bool): if True, color the syndrome and data qubits
                red and green, respectively; colored by default.
        Returns:
            None
        '''
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
            node_color = self.color(point)
            ax.scatter(x, y, z, s=70, c=color*node_color+(1-color)*'k')
            indices = {c: n for (n, c) in enumerate(self.nodes)}
            if label:
                ax.text(x, y, z, str(indices[point]), fontdict=fontdict,
                        color='MediumBlue', backgroundcolor='w')
        # Plotting edges.
        for edge in self.edges:
            x1, z1, y1 = edge[0]
            x2, z2, y2 = edge[1]
            plt.plot([x1, x2], [y1, y2], [z1, z2], c='grey')

        plt.xticks(range(0, 2*nx + 1))
        plt.yticks(range(0, 2*nz + 1))
        ax.set_zticks(range(0, 2*ny + 1))
        ax.set_xlabel('x', fontdict=fontdict)
        ax.set_ylabel('z', fontdict=fontdict)
        ax.set_zlabel('y', fontdict=fontdict)
        plt.rcParams['grid.color'] = "lightgray"
        plt.tight_layout(pad=3)
        plt.draw()
        return ax

    def index(self):
        indexed_graph = nx.convert_node_labels_to_integers(self, ordering='sorted', label_attribute='pos')
        return indexed_graph


def RHG_graph(dims, pol=0):
    '''Return an EGraph of the RHG lattice.'''
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims

    lattice = EGraph()
    lattice.graph['dims'] = dims
    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [(2*i+1, 2*j+1, 2*k) for (i, j, k) in
                it.product(range(nx), range(ny), range(nz+1))]
    odd_red = [(2*i, 2*j, 2*k+1) for (i, j, k) in
               it.product(range(nx+1), range(ny+1), range(nz))]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [(2*i+1, 2*j, k) for (i, j, k) in
                  it.product(range(nx), range(ny+1), range(2*nz+1))]
    odd_green = [(2*i, 2*j+1, k) for (i, j, k) in
                 it.product(range(nx+1), range(ny), range(2*nz+1))]
    all_green = set(even_green + odd_green)

    # Coordinates of all potential neighbours of red vertices.
    def red_neighbours(p):
        right = (p[0]+1, p[1], p[2])
        top = (p[0], p[1]+1, p[2])
        left = (p[0]-1, p[1], p[2])
        bottom = (p[0], p[1]-1, p[2])
        return [bottom, left, top, right]

    def green_neighbours(p):
        return {(p[0], p[1], p[2]-1), (p[0], p[1], p[2]+1)}

    for point in all_red:
        for i in range(4):
            polarity = (-1) ** (pol * (point[2] + i))
            neighbour = red_neighbours(point)[i]
            if neighbour in all_green:
                lattice.add_edge(point, neighbour, weight=polarity)
        lattice.nodes[point]['color'] = 'red'

    for point in all_green:
        polarity = (-1) ** (pol * (point[1] + 1))
        for neighbour in green_neighbours(point):
            if neighbour in all_green:
                lattice.add_edge(point, neighbour, weight=polarity)
        lattice.nodes[point]['color'] = 'green'

    return lattice


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


def basic_translate(outcomes):
    """Naively translate CV outcomes to bit values.

    The function treats values in (-sqrt(pi)/2, sqrt(pi)/2) as 0
    and those in (sqrt(pi)/2, 3sqrt(pi)/2) as 1. Values on the
    boundary are assigned a random bit value. The rest of the bins
    are defined periodically.
    Args:
        outcomes (array): the values of a p-homodyne measurement
    Retruns:
        array: the corresponding bit values.
    """
    # Bin width
    alpha = np.sqrt(np.pi)
    n = len(outcomes)
    bit_values = np.zeros(n)
    for i in range(n):
        div = np.divmod(outcomes[i], alpha)
        if div[1] < alpha/2:
            bit_values[i] = div[0] % 2
        elif div[1] > alpha / 2:
            bit_values[i] = (div[0]+1) % 2
        else:
            bit_values[i] = np.random.randint(2)
    return bit_values


def p_err(var, var_num=5, translator=basic_translate):
    """Return the probability of Z errors for a list of variances.

    Args:
        var (array): array of lattice p variances
        var_num (float): number of variances away from the origin we
            include in the integral
        translator (function): the CV-to-bit translator; the basic
            binning function by default
    Returns:
        array: probability of Z (phase flip) errors for each variance.
    """
    if translator == basic_translate:
        # Find largest bin number by finding largest variance, multiplying by
        # var_num, then rounding up to nearest integer of form 4n+1, which are
        # the left boundaries of the 0 bins mod sqrt(pi)
        n_max = int(np.ceil(var_num*np.amax(var))//2*4 + 1)
        # error = 1 - integral over the 0 bins
        # Initiate a list with length same as var
        error = np.ones(len(var))
        # Integral over 0 bins that fell within var_num*var_max away from
        # origin
        for i in range(-n_max, n_max, 4):
            error -= 0.5*(erf((i+2)*np.sqrt(np.pi)/(2*var))
                          - erf(i*np.sqrt(np.pi)/(2*var)))
    return error


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
        # indfunc (function): the function used to index the coordinates
        #     of the lattice; basic_index by default
        translator (function): the translator or binning function
            between CV outcomes and bit values; basic_translate by
            default

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
        hom_outcomes (array): if measure_p has been run, the results
            of the latest p-homodyne measurement
        bit_values (array): if measure_p has been run, the bit values
            associated with the latest p-homodyne measurement, assuming
            binning function translator.
    '''
    def __init__(self, g, model='grn', p_inds=[], swap_prob=None,
                 delta=0.01, translator=basic_translate):
        if isinstance(g, EGraph):
            self.graph = g
        else:
            self.graph = EGraph(g)
        self._translator = translator
        self._SCZ = SCZ_mat(self.graph.adj())
        self._delta = delta
        self._N = self.graph.number_of_nodes()

        if swap_prob:
            self.hybridize(swap_prob)
        else:
            self._p_inds = p_inds
        if model == 'grn':
            self.grn_model(delta)
        self._indexed_graph = self.graph.index()


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
        self._Z_probs = p_err(self._p_var, translator=self._translator)

    def measure_p(self):
        """Conduct a p-homodyne measurement on the lattice."""
        # Randomly sample of a normal distribution centred at 0 with
        # covariances matrix given by cov_p.
        outcomes = mvn(mean=np.zeros(self._N), cov=self._cov_p)
        self._hom_outcomes = outcomes
        self._bit_values = self._translator(outcomes)

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
        self._p_var = np.diag(self._cov_p)

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
    def p_var(self):
        """array: the p variances of the modes."""
        return self._p_var

    @property
    def Z_probs(self):
        """array: the phase error probabilities of the modes."""
        try:
            return self._Z_probs
        except Exception:
            print('Z error probabilities have not yet been computed. Please '
                  'use eval_Z_probs() first.')
            return

    @property
    def hom_outcomes(self):
        """array: the results of the p-homodyne measurement."""
        try:
            return self._hom_outcomes
        except Exception:
            print('A homodyne measurement has not yet been performed. Please '
                  'use measure_p() first.')
            return

    @property
    def bit_values(self):
        """array: bit values associated with the p measurement."""
        try:
            return self._bit_values
        except Exception:
            print('A homodyne measurement has not yet been performed. Please '
                  'use measure_p() first.')
            return

    def sketch(self, label=False):
        """Sketch the CVRHG lattice. GKP states black, p states orange.
        Args:
            label (bool): if 'v', label the p variances; if 'e', the
                phase error probabilities; if 'h', the p-homodyne
                outcomes; if 'b', the bit values; otherwise, no label.
        Returns:
            None
        """
        idg = self._indexed_graph
        p_coords = [idg.nodes[ind]['pos'] for ind in self._p_inds]
        ax = self.graph.draw(0, 0)
        for point in p_coords:
            x, z, y = point
            ax.scatter(x, y, z, s=40, c='orange')
        orange_line = mlines.Line2D([], [], color='orange', marker='.',
                                    markersize=14, label='p')
        black_line = mlines.Line2D([], [], color='black', marker='.',
                                   markersize=14, label='GKP')
        ax.legend(handles=[orange_line, black_line])
        if label:
            if label == 'v':
                results = self._p_var
                title = "p Variances"
            elif label == 'e':
                results = self.Z_probs
                title = "Phase Error Probabilities"
            elif label == 'h':
                results = self.hom_outcomes
                title = "p-homodyne Outcomes"
            elif label == 'b':
                results = self.bit_values
                title = "Bit Values"
            if results is None:
                return
            else:
                for i in range(self._N):
                    x, z, y = idg.nodes[i]['pos']
                    ax.set_title(title, fontsize=16, family='serif')
                    ax.text(x, y, z, '{:.0g}'.format(results[i]),
                            fontdict=fontdict, color='MediumBlue',
                            backgroundcolor='w')
        plt.draw()


if __name__ == '__main__':
    RHG = RHG_graph(1)
    RHG.draw(label=1)
    # p states at random locations.
    # Test dim X dim X dim lattice, with p% p-squeezed states randomly located,
    # and variance of of delta.
    dim = 2
    delta = 0.1
    # Number of qubits in the lattice.
    N = RHG.number_of_nodes()
    # Percent p-squeezed states.
    p = 0.3
    # p states at random locations.
    G = CVGraph(RHG, swap_prob=0.1, delta=delta)
    G.eval_Z_probs()
    G.measure_p()
    for label in {'v', 'e', 'h', 'b'}:
        G.sketch(label)