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
"""Module for the RHG code."""
import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import itertools as it
import matplotlib.pyplot as plt

from graphstates import EGraph, CVGraph

def RHG_graph(dims, boundaries='primal', pol=0):
    '''Return an EGraph of the RHG lattice.'''
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims

    # Dealing with boundaries
    if type(boundaries) == str:
        boundaries = [boundaries] * 6
    if not boundaries:
        boundaries = ['primal'] * 5 + ['dual']
        nz += 1
    type_dict = {'dual': 0, 'primal': 1}
    bound_labels = ['left', 'bottom', 'front', 'right', 'top', 'back']
    bound_labels_dict = {bound_labels[i]: boundaries[i] for i in range(6)}
    x_factor, y_factor, z_factor = [type_dict[typ] for typ in boundaries[3:]]

    # Define the EGraph of the lattice
    lattice = EGraph(dims=dims, boundaries=bound_labels_dict)

    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [(2*i+1, 2*j+1, 2*k) for (i, j, k) in
                it.product(range(nx), range(ny), range(nz+z_factor))]
    odd_red = [(2*i, 2*j, 2*k+1) for (i, j, k) in
               it.product(range(nx+x_factor), range(ny+y_factor), range(nz))]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [(2*i+1, 2*j, k) for (i, j, k) in
                  it.product(range(nx), range(ny+y_factor), range(2*nz+z_factor))]
    odd_green = [(2*i, 2*j+1, k) for (i, j, k) in
                 it.product(range(nx+x_factor), range(ny), range(2*nz+z_factor))]
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

def RHG_syndrome_coords(dims, code='primal'):
    """Return a list of lists of coordinates associated to six-body
    X stabilizer elements of the RHG lattice."""

    nx, ny, nz = dims
    combs = it.product(range(nx), range(ny), range(nz))
    if code == 'primal':
        all_syndrome6 = [[(2*i+1, 2*j+1, 2*k),
                          (2*i+1, 2*j, 2*k+1),
                          (2*i, 2*j+1, 2*k+1),
                          (2*i+1, 2*j+1, 2*k+2),
                          (2*i+1, 2*j+2, 2*k+1),
                          (2*i+2, 2*j+1, 2*k+1)
                          ] for (i, j, k) in combs]
    return all_syndrome6


def RHG_boundary_coords(dims, code='primal'):
    """Obtain the coordinates of the vertices at the centres of primal
    cubes on the boundary of the RHG lattice with dimension dims."""
    if code == 'primal':
        odds = [range(1, 2*dims[0], 2),
                range(1, 2*dims[1], 2),
                range(1, 2*dims[2], 2)]
        combs = []
        for i, j in ((0,1), (0,2), (1,2)):
            for tup in it.product(odds[i], odds[j]):
                ind = {0, 1, 2}.difference({i, j}).pop()
                l = list(tup)
                m = list(tup)
                l.insert(ind, 0)
                m.insert(ind, 2 * dims[ind])
                combs.append(tuple(l))
                combs.append(tuple(m))
    return combs


def RHG_slice_coords(RHG_lattice, plane, number):
    plane_dict = {'x': 0, 'y': 1, 'z': 2}
    plane_ind = plane_dict[plane]
    coords = [point for point in RHG_lattice.nodes if point[plane_ind] == number]
    return coords

def RHG_stabilizers(G, code='primal'):
    """Return a list of subgraphs induced by the qubits with cordinates
    from RHG_syndrome_coords."""

    syn_list = RHG_syndrome_coords(G.graph.graph['dims'], code)
    cube_list = []
    for l in syn_list:
        cube_graph = G.graph.subgraph(l)
        cube_list.append(RHGCube(CVGraph(cube_graph)))
    return cube_list


class RHGCube:
    """A class for representing an RHG latice cube. Initialized by supplying
    a CVGraph object corresponding to the subgraph induced by the six facial
    nodes of the cube."""
    def __init__(self, G):
        self.cvgraph = G
        self._parity = None

    def parity(self):
        if not self._parity is None:
            return self._parity
        G = self.cvgraph
        bit_vals = [G.graph.nodes[node]['bit_val'] for node in G.graph]
        # hom_vals = [G.graph.nodes[node]['hom_val'] for node in G.graph]
        # return basic_translate([np.sum(hom_vals)])
        par = int(np.sum(bit_vals) % 2)
        self._parity = par
        return par

    def coords(self):
        return [tup for tup in self.cvgraph.graph.nodes]

    def xlims(self):
        xs = [tup[0] for tup in self.coords()]
        xmin, xmax = np.min(xs), np.max(xs)
        return xmin, xmax

    def ylims(self):
        ys = [tup[1] for tup in self.coords()]
        ymin, ymax = np.min(ys), np.max(ys)
        return ymin, ymax

    def zlims(self):
        zs = [tup[2] for tup in self.coords()]
        zmin, zmax = np.min(zs), np.max(zs)
        return zmin, zmax

    def midpoint(self):
        return (np.average(self.xlims()), np.average(self.ylims()), np.average(self.zlims()))

if __name__ == '__main__':
    RHG = RHG_graph(1)
    RHG.draw(label=1)
    # p states at random locations.
    # Test dim X dim X dim lattice, with p% p-squeezed states randomly located,
    # and variance of of delta.
    dim = 2
    delta = 0.01
    # Number of qubits in the lattice.
    N = RHG.number_of_nodes()
    # Percent p-squeezed states.
    p = 0.3
    # p states at random locations.
    G = CVGraph(RHG, swap_prob=0.1, delta=delta)
    G.eval_Z_probs()
    G.measure_p()
    G.eval_Z_probs_cond()
    for label in {
            'var_p',
            'p_phase',
            'p_phase_cond',
            'hom_val',
            }:
        G.sketch(label)
