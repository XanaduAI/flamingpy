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

def RHG_graph(dims, boundaries='natural', pol=0):
    '''Return an EGraph of the RHG lattice.'''
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims

    # Dealing with boundaries
    if boundaries == 'natural':
        boundaries = ['periodic', 'periodic', 'primal']
    elif type(boundaries) == str:
        boundaries = [boundaries] * 3
    min_dict = {'primal': 0, 'dual': 1, 'periodic': 0}
    max_dict = {'primal': 1, 'dual': 0, 'periodic': 0}
    bound_labels = ['x', 'y', 'z']
    bound_labels_dict = {bound_labels[i]: boundaries[i] for i in range(3)}
    x_min, y_min, z_min = [min_dict[typ] for typ in boundaries]
    x_max, y_max, z_max = [max_dict[typ] for typ in boundaries]

    # Define the EGraph of the lattice
    lattice = EGraph(dims=dims, boundaries=bound_labels_dict)

    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [(2*i+1, 2*j+1, 2*k) for (i, j, k) in
                it.product(range(nx), range(ny), range(z_min, nz+z_max))]
    odd_red = [(2*i, 2*j, 2*k+1) for (i, j, k) in
               it.product(range(x_min, nx+x_max), range(y_min, ny+y_max), range(nz))]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [(2*i+1, 2*j, k) for (i, j, k) in
                  it.product(range(nx), range(y_min, ny+y_max), range(z_min, 2*nz+z_max))]
    odd_green = [(2*i, 2*j+1, k) for (i, j, k) in
                 it.product(range(x_min, nx+x_max), range(ny), range(z_min, 2*nz+z_max))]
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

    bound_arr = np.array(boundaries)
    periodic_inds = np.where(bound_arr=='periodic')[0]
    dims = (nx, ny, nz)
    for ind in periodic_inds:
        low_slice = set([point for point in lattice.nodes if point[ind] == 0])
        high_slice = set([point for point in lattice.nodes if point[ind] == 2*dims[ind]-1])
        if ind in (0, 1):
            low_reds = all_red.intersection(low_slice)
            high_reds = all_red.intersection(high_slice)
            for point in low_reds:
                polarity = (-1) ** (pol * (point[2] + ind + 1))
                high_green = list(point)
                high_green[ind] = 2*dims[ind]-1
                lattice.add_edge(point, tuple(high_green), weight=polarity)
            for point in high_reds:
                polarity = (-1) ** (pol * (point[2] + ind + 1))
                low_green = list(point)
                low_green[ind] = 0
                lattice.add_edge(point, tuple(low_green), weight=polarity)
        if ind == 2:
            low_greens = all_green.intersection(low_slice)
            for point in low_greens:
                polarity = (-1) ** (pol * (point[1] + 1))
                high_green = list(point[:])
                high_green[ind] = 2*dims[ind] - 1
                lattice.add_edge(point, tuple(high_green), weight=polarity)
    return lattice

def RHG_syndrome_coords(G):
    """Return a list of lists of coordinates associated to six-body
    X stabilizer elements of the RHG lattice."""

    dims = G.graph['dims']
    boundaries = np.array(G.graph['boundaries'])

    min_dict = {'primal': 0, 'dual': 1, 'periodic': 0}
    max_dict = {'primal': 0, 'dual': 1, 'periodic': 1}
    mins = [min_dict[typ] for typ in boundaries]
    ma = [max_dict[typ] for typ in boundaries]
    maxes = [dims[i] - ma[i] for i in range(3)]

    ranges = lambda mins, maxes: [range(mins[i], maxes[i]) for i in range(3)]
    def stabe_points(combs):
        six_body_stabes = [[
            (2*i, 2*j+1, 2*k+1),
            (2*i+1, 2*j, 2*k+1),
            (2*i+1, 2*j+1, 2*k),
            (2*i+2, 2*j+1, 2*k+1),
            (2*i+1, 2*j+2, 2*k+1),
            (2*i+1, 2*j+1, 2*k+2)] for (i, j, k) in combs]
        return six_body_stabes

    middle_combs = it.product(*ranges(mins, maxes))
    six_body_stabes = stabe_points(middle_combs)

    dual_inds = np.where(boundaries =='dual')[0]
    periodic_inds = np.where(boundaries =='periodic')[0]
    five_body_stabes = []
    periodic_six_bodies = []
    if dual_inds.size or periodic_inds.size:
        for ind in range(3):
            m2, M2 = mins[:], maxes[:]
            m2[ind] = maxes[ind]
            M2[ind] += 1
            combs_max = it.product(*ranges(m2, M2))
            stabes_max = stabe_points(combs_max)
            if ind in periodic_inds:
                for stabe in stabes_max:
                    highest_point = stabe[3+ind]
                    other_side = list(highest_point)
                    other_side[ind] = mins[ind]
                    print(highest_point, other_side)
                    stabe[3+ind] = other_side
                    periodic_six_bodies.append(stabe)
            if ind in dual_inds:
                m1, M1 = mins[:], maxes[:]
                m1[ind] -= 1
                M1[ind] = mins[ind]
                combs_min = it.product(*ranges(m1, M1))
                stabes_min = stabe_points(combs_min)
                for i in range((len(stabes_min))):
                    stabes_min[i].pop(ind)
                    stabes_max[i].pop(3+ind)
                five_body_stabes += stabes_min
                five_body_stabes += stabes_max

    all_stabes = six_body_stabes + periodic_six_bodies + five_body_stabes
    return all_stabes


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
