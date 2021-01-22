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
"""Classes for the RHG code and related functions."""
import numpy as np
import itertools as it

from graphstates import EGraph, CVGraph


def RHG_graph(dims, boundaries='natural', polarity=False):
    """Return an EGraph of a dims-dimensional RHG lattice.

    Generate an RHG lattice with dimensions given by dims (an integer
    denoting the number of stabilizer cubes in each direciton,
    or a tuple specifying all three dimensions). By default, the
    boundary conditions are periodic in 'x' and 'y' (the 'natural'
    boundary conditions) and primal in 'z'. However, any combination
    of boundary conditions can be specified. By default, all edges
    are weight-1, but if polarity is True, make edges in certain
    directions -1, for subsequent CV noise reduction.

    Args:
        dims (int or list-type): the dimensions of the lattice (the
            number of stabilizer cubes in each direction). In the case
            of dual boundaries, incomplete stabilizer cubes are treated
            as 0.5 of complete cubes
        boundaries (str or list-type): if 'natural,' assume periodic
            boundaries in x and y, and primal in z. If another string,
            assume those boundaries in every direction. If a list-type
            with three elements, assume x, y, z boundaries individually
            specified
        polarity (bool): if True, make some edges have weight -1, per
            noise-reduction strategy for the eventual CVGraph. Set
            edge 'color' attribute to blue/red for -1/+1 weights. By
            default, all edges have weight +1.
    Returns:
        EGraph: the RHG lattice.
    """
    # Dimensions of the lattice.
    # TODO: Check accuracy of dosctring description and compare with/
    # incorporate code distance.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims

    # Dealing with boundaries
    if boundaries == 'natural':
        boundaries = ['periodic', 'periodic', 'primal']
    elif type(boundaries) == str:
        boundaries = [boundaries] * 3

    # Define the EGraph of the lattice
    lattice = EGraph(dims=dims, boundaries=boundaries)

    # Constrain the ranges of the coordinates depending on the type
    # of boundaries.
    min_dict = {'primal': 0, 'dual': 1, 'periodic': 0}
    max_dict = {'primal': 1, 'dual': 0, 'periodic': 0}
    x_min, y_min, z_min = [min_dict[typ] for typ in boundaries]
    x_max, y_max, z_max = [max_dict[typ] for typ in boundaries]

    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [(2*i + 1, 2*j + 1, 2*k) for (i, j, k) in
                it.product(range(nx), range(ny), range(z_min, nz+z_max))]
    odd_red = [(2*i, 2*j, 2*k + 1) for (i, j, k) in
               it.product(range(x_min, nx+x_max), range(y_min, ny+y_max), range(nz))]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [(2*i + 1, 2*j, k) for (i, j, k) in
                  it.product(range(nx), range(y_min, ny + y_max), range(z_min, 2*nz + z_max))]
    odd_green = [(2*i, 2*j + 1, k) for (i, j, k) in
                 it.product(range(x_min, nx + x_max), range(ny), range(z_min, 2*nz + z_max))]
    all_green = set(even_green + odd_green)

    def red_neighbours(p):
        """Coordinates of all potential neighbours of red vertices."""
        right = (p[0] + 1, p[1], p[2])
        top = (p[0], p[1] + 1, p[2])
        left = (p[0] - 1, p[1], p[2])
        bottom = (p[0], p[1] - 1, p[2])
        return [bottom, left, top, right]

    def green_neighbours(p):
        """Coordinates of all potential neighbours of green vertices."""
        return {(p[0], p[1], p[2] - 1), (p[0], p[1], p[2] + 1)}

    # Polarity-dependent color function: blue for +1, red for -1.
    color = lambda pol: ((pol + 1) // 2) * 'b' + abs((pol - 1) // 2) * 'r'

    # Add edges between red points and all their neighbours.
    for point in all_red:
        is_neighbour = 0
        for i in range(4):
            pol = (-1) ** (polarity * (point[2] + i))
            neighbour = red_neighbours(point)[i]
            if neighbour in all_green:
                is_neighbour += 1
                lattice.add_edge(point, neighbour, weight=pol, color=color(pol))
        if is_neighbour:
            lattice.nodes[point]['color'] = 'red'

    # Add edges between green points and all their neighbours.
    for point in all_green:
        is_neighbour = 0
        pol = (-1) ** (polarity * (point[1] + 1))
        for neighbour in green_neighbours(point):
            if neighbour in all_green:
                is_neighbour += 1
                lattice.add_edge(point, neighbour, weight=pol, color=color(pol))
        if is_neighbour:
            lattice.nodes[point]['color'] = 'green'

    # Dealing with periodic boundary conditions.
    bound_arr = np.array(boundaries)
    periodic_inds = np.where(bound_arr == 'periodic')[0]
    for ind in periodic_inds:
        # First and last slices of the lattice in the direction
        # specified by ind.
        low_slice = set([point for point in lattice.nodes if point[ind] == 0])
        high_slice = set([point for point in lattice.nodes if point[ind] == 2*dims[ind] - 1])
        if ind in (0, 1):
            low_reds = all_red.intersection(low_slice)
            high_reds = all_red.intersection(high_slice)
            # Connect red in first slice to greens in last slice.
            for point in low_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                high_green = list(point)
                high_green[ind] = 2 * dims[ind] - 1
                lattice.add_edge(point, tuple(high_green), weight=pol, color=color(pol))
            # Connect reds last slice to greens in first slice.
            for point in high_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                low_green = list(point)
                low_green[ind] = 0
                lattice.add_edge(point, tuple(low_green), weight=pol, color=color(pol))
        # If periodic in z direction, connect greens in first slice with
        # greens in last slice.
        if ind == 2:
            low_greens = all_green.intersection(low_slice)
            for point in low_greens:
                pol = (-1) ** (polarity * (point[1] + 1))
                high_green = list(point[:])
                high_green[ind] = 2 * dims[ind] - 1
                lattice.add_edge(point, tuple(high_green), weight=pol, color=color(pol))

    return lattice

def RHG_syndrome_coords(G):
    """Return the syndrome coordinates for RHG lattice G.

    Generate a list of lists containing coordinates of six-body X
    stabilizers on non-periodic boundaries, followed by six-body X
    stabilizers on periodic boundaries, followed by five-body X
    stabilizers on dual boundaries.

    Args:
        G (EGraph): the RHG lattice

    Returns:
        list of lists of tuples: the syndrome coordinates.
    """
    # Dimensions, boundary types, max and min ranges.
    dims = list(G.graph['dims'])
    boundaries = np.array(G.graph['boundaries'])

    min_dict = {'primal': 0, 'dual': 1, 'periodic': 0}
    mins = [min_dict[typ] for typ in boundaries]
    maxes = np.array([2 * dims[i] - 1 for i in (0, 1, 2)])

    # Function for generating ranges from lists of mins and maxes.
    ranges = [range(dims[i]) for i in range(3)]
    inds = it.product(*ranges)

    # All potential six-body stabilizers
    all_six_bodies = [[
            (2*i, 2*j + 1, 2*k + 1),
            (2*i + 1, 2*j, 2*k + 1),
            (2*i + 1, 2*j + 1, 2*k),
            (2*i + 2, 2*j + 1, 2*k + 1),
            (2*i + 1, 2*j + 2, 2*k + 1),
            (2*i + 1, 2*j + 1, 2*k + 2)] for (i, j, k) in inds]

    all_stabes = []

    periodic_inds = np.where(boundaries == 'periodic')[0]
    dual_inds = np.where(boundaries == 'dual')[0]
    for stabe in all_six_bodies:
        actual_stabe = list(set(stabe).intersection(set(G)))
        if len(actual_stabe) == 6:
            all_stabes.append(actual_stabe)
        if len(actual_stabe) == 5:
            for ind in (0, 1, 2):
                lowest_point = list(stabe[ind])
                highest_point = stabe[3 + ind]
                if lowest_point[ind] == (2 * dims[ind] - 2):
                    if ind in dual_inds:
                        all_stabes += [actual_stabe]
                    if ind in periodic_inds:
                        lowest_point[ind] = 0
                        actual_stabe += [tuple(lowest_point)]
                        all_stabes.append(actual_stabe)
                if highest_point[ind] == 2 and ind in dual_inds:
                    all_stabes.append(actual_stabe)
        if len(actual_stabe) == 4:
            average_point = [sum([point[i] for point in actual_stabe]) / 4 for i in (0, 1, 2)]
            rounded_avs = np.array([round(av) for av in average_point])
            high_inds = np.where(maxes == rounded_avs)[0]
            low_inds = np.where(mins == rounded_avs)[0]
            if len(high_inds) >= 2:
                boundary_inds = high_inds
            if len(high_inds) == 1:
                boundary_inds = np.array([low_inds[0], high_inds[0]])
            if len(high_inds) == 0:
                boundary_inds = low_inds
            if boundary_inds[0] in periodic_inds:
                new_point = rounded_avs.copy()
                new_point[boundary_inds[0]] = mins[boundary_inds[0]]
                actual_stabe += [tuple(new_point)]
            if boundary_inds[1] in periodic_inds:
                new_point = rounded_avs.copy()
                new_point[boundary_inds[1]] = mins[boundary_inds[1]]
                actual_stabe += [tuple(new_point)]
            if len(boundary_inds) == 3:
                if boundary_inds[2] in periodic_inds:
                    new_point = rounded_avs.copy()
                    new_point[boundary_inds[2]] = mins[boundary_inds[2]]
                    actual_stabe += [tuple(new_point)]
            all_stabes.append(actual_stabe)

        if len(actual_stabe) == 3:
            point = [2 * dims[0] - 1] * 3
            for ind in (0, 1, 2):
                if ind in periodic_inds:
                    point_ind = point[:]
                    point_ind[ind] = 0
                    actual_stabe += [tuple(point_ind)]
            all_stabes.append(actual_stabe)

    # Dealing with six-body X stabilizers on perodic boundaries,
    # and five-body X stabilizers on dual boundaries.

    # Put all the stabilizers into a list.
    return all_stabes


def RHG_boundary_coords(G):
    """Return coordinates of syndrome qubits on primal boundaries.

    For now, omit primal boundaries in the z-direction.
    """
    # TODO: Dual boundaries, z-direction primal if necessary.
    dims = G.graph['dims']
    boundaries = np.array(G.graph['boundaries'])
    bound_inds = np.where(boundaries == 'primal')[0]
    # Odd indices, which is where primal syndrome qubits are located.
    odds = [range(1, 2 * dims[0], 2),
            range(1, 2 * dims[1], 2),
            range(1, 2 * dims[2], 2)]
    combs = []
    for ind in bound_inds:
        if ind != 2:
            for i, j in ((0, 1), (0, 2), (1, 2)):
                for tup in it.product(odds[i], odds[j]):
                    l = list(tup)
                    m = list(tup)
                    l.insert(ind, 0)
                    m.insert(ind, 2 * dims[ind])
                    combs.append(tuple(l))
                    combs.append(tuple(m))
    return combs


def RHG_slice_coords(G, plane, number, boundaries='all'):
    """Obtain all the coordinates in a slice of RHG lattice G.

    Args:
        G (EGraph): the RHG lattice
        plane (str): 'x', 'y', or 'z', denoting the slice direction
        number (int): the index of the slice

    Returns:
        list of tuples: the coordinates of the slice.
    """
    plane_dict = {'x': 0, 'y': 1, 'z': 2}
    plane_ind = plane_dict[plane]
    coords = [point for point in G.nodes if point[plane_ind] == number]
    dims = G.graph['dims']
    if boundaries in ('primal', 'dual'):
        remaining_inds = {0, 1, 2}.difference({plane_ind})
        undesirables = []
        for i in range(len(coords)):
            for ind in remaining_inds:
                if boundaries == 'dual':
                    if coords[i][ind] == 0:
                        undesirables.append(i)
                if boundaries == 'primal':
                    if coords[i][ind] == (2*dims[ind] - 1):
                        undesirables.append(i)
        coords = [coords[i] for i in range(len(coords)) if i not in undesirables]
        return coords

    return coords


class RHGCube:
    """A class for representing an RHG latice stabilizer cube.

    Arguments:
        G (CVGraph): a CVGraph whose EGraph is the subgraph induced
            by the (five or six) stabilizer vertices.

    Attributes:
        cvgraph (CVGraph): the corresponding CVGraph
        _parity (int): the total parity of the cube, if self.parity has
            been run.
        type (str): the type of stabilizer ('six-body' or 'five-body').
    """

    def __init__(self, G):
        """Initialize the RHGCube with its associated cvgraph."""
        self.cvgraph = G
        self._parity = None
        self.type = 'six-body'

    def parity(self):
        """Compute total parity of the cube.

        For now, add the bit_values found in the 'bit_val' attribute
        of the nodes of the cvgraph.

        Returns:
            int: the total parity
        """
        # If parity already known, set self._parity to it.
        if self._parity is not None:
            return self._parity
        G = self.cvgraph
        bit_vals = [G.graph.nodes[node]['bit_val'] for node in G.graph]

        # TODO: Process the homodyne outcomes first, as below.
        # hom_vals = [G.graph.nodes[node]['hom_val'] for node in G.graph]
        # return basic_translate([np.sum(hom_vals)])

        par = int(np.sum(bit_vals) % 2)
        self._parity = par
        return par

    def coords(self):
        """Return the coordinates of the syndrome vertices."""
        return [tup for tup in self.cvgraph.graph.nodes]

    def xlims(self):
        """Return the smallest and largest x vals of the coordinates."""
        xs = [tup[0] for tup in self.coords()]
        xmin, xmax = np.min(xs), np.max(xs)
        return xmin, xmax

    def ylims(self):
        """Return the smallest and largest y vals of the coordinates."""
        ys = [tup[1] for tup in self.coords()]
        ymin, ymax = np.min(ys), np.max(ys)
        return ymin, ymax

    def zlims(self):
        """Return the smallest and largest z vals of the coordinates."""
        zs = [tup[2] for tup in self.coords()]
        zmin, zmax = np.min(zs), np.max(zs)
        return zmin, zmax

    def midpoint(self):
        """Return the midpoint of the cube."""
        # TODO: Perhaps choose a different point for a five-body
        # X stabilizer, for plotting purposes.
        return (np.average(self.xlims()), np.average(self.ylims()), np.average(self.zlims()))


def RHG_stabilizers(G):
    """Return a list of RHGCubes corresponding to RHG stabilizers.

    For each stabilizer in G, generate an RHGCube and add it to the
    list. For five-body X stabilizers, change the type attribute to
    'five-body'.

    Args:
        G (EGraph): the RHG lattice

    Returns:
        list: the RHGCubes.
    """
    syn_list = RHG_syndrome_coords(G.graph)
    cube_list = []
    for stabe in syn_list:
        cube_graph = G.graph.subgraph(stabe)
        cube = RHGCube(CVGraph(cube_graph))
        if len(cube_graph) == 5:
            cube.type = 'five-body'
        cube_list.append(cube)
    return cube_list


if __name__ == '__main__':
    boundaries = 'natural'
    # # boundaries = 'primal'
    # # boundaries = 'periodic'
    # # boundaries = ['primal', 'dual', 'primal']
    RHG = RHG_graph(2, boundaries=boundaries, polarity=True)
    ax = RHG.draw(color_nodes=True, color_edges=True, label=True)

    # Test dim X dim X dim lattice, with p% p-squeezed states randomly
    # located, and variance of of delta.
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
            'hom_val'
            }:
        G.sketch(label)
