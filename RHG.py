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
from matplotlib import pyplot as plt
import itertools as it
from graphstates import EGraph, CVGraph


class RHGCode:
    """RHG code class."""

    def __init__(self, distance, boundaries='finite', polarity=False):
        """Initialize the RHG code."""
        # TODO: Check code distance convention.
        if boundaries == 'finite':
            self.boundaries = ['primal', 'dual', 'dual']
        elif type(boundaries) == str:
            self.boundaries = [boundaries] * 3
        self.dims = (distance, distance, distance)
        self.polarity = polarity

        self.graph = self.graph_generator(self.dims, self.boundaries, self.polarity)
        self.stabilizers = self.identify_syndrome()
        self.boundary_coords = self.identify_boundary()

    def graph_generator(self, dims, boundaries='finite', polarity=False):
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
        nx, ny, nz = dims
        # Dealing with boundaries

        # Define the EGraph of the lattice
        lattice = EGraph(dims=dims)

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

    def identify_stabilizers(self):
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
        G = self.graph
        # Dimensions, boundary types, max and min ranges.
        dims = list(self.dims)
        boundaries = np.array(self.boundaries)
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

        all_cubes = []

        periodic_inds = np.where(boundaries == 'periodic')[0]
        dual_inds = np.where(boundaries == 'dual')[0]
        for stabe in all_six_bodies:
            actual_stabe = list(set(stabe).intersection(set(G)))
            if len(actual_stabe) == 6:
                cube = RHGCube(G.subgraph(actual_stabe))
                cube.physical = stabe
                all_cubes.append(cube)
            if len(actual_stabe) == 5:
                for ind in (0, 1, 2):
                    lowest_point = list(stabe[ind])
                    highest_point = stabe[3 + ind]
                    if lowest_point[ind] == (2 * dims[ind] - 2):
                        if ind in dual_inds:
                            cube = RHGCube(G.subgraph(actual_stabe))
                            cube.physical = stabe
                            all_cubes.append(cube)
                        if ind in periodic_inds:
                            lowest_point[ind] = 0
                            virtual_point = tuple(lowest_point)
                            actual_stabe += [virtual_point]
                            cube = RHGCube(G.subgraph(actual_stabe))
                            cube.physical = stabe
                            all_cubes.append(cube)
                    if highest_point[ind] == 2 and ind in dual_inds:
                        cube = RHGCube(G.subgraph(actual_stabe))
                        cube.physical = stabe
                        all_cubes.append(cube)
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
                    virtual_point = tuple(new_point)
                    actual_stabe += [virtual_point]
                    cube = RHGCube(G.subgraph(actual_stabe))
                    cube.physical = stabe
                elif boundary_inds[0] in dual_inds:
                    cube = RHGCube(G.subgraph(actual_stabe))
                    cube.physical = stabe
                if boundary_inds[1] in periodic_inds:
                    new_point = rounded_avs.copy()
                    new_point[boundary_inds[1]] = mins[boundary_inds[1]]
                    virtual_point = tuple(new_point)
                    actual_stabe += [virtual_point]
                    cube = RHGCube(G.subgraph(actual_stabe))
                    cube.physical = stabe
                elif boundary_inds[1] in dual_inds:
                    cube = RHGCube(G.subgraph(actual_stabe))
                    cube.physical = stabe
                if len(boundary_inds) == 3:
                    if boundary_inds[2] in periodic_inds:
                        new_point = rounded_avs.copy()
                        new_point[boundary_inds[2]] = mins[boundary_inds[2]]
                        virtual_point = tuple(new_point)
                        actual_stabe += [virtual_point]
                        cube = RHGCube(G.subgraph(actual_stabe))
                        cube.physical = stabe
                    elif boundary_inds[2] in dual_inds:
                        cube = RHGCube(G.subgraph(actual_stabe))
                        cube.physical = stabe
                all_cubes.append(cube)

            if len(actual_stabe) == 3:
                point = [2 * dims[0] - 1] * 3
                for ind in (0, 1, 2):
                    if ind in periodic_inds:
                        point_ind = point[:]
                        point_ind[ind] = 0
                        virtual_point = tuple(point_ind)
                        actual_stabe += [virtual_point]
                cube = RHGCube(G.subgraph(actual_stabe))
                cube.physical = stabe
                all_cubes.append(cube)
        # Dealing with six-body X stabilizers on perodic boundaries,
        # and five-body X stabilizers on dual boundaries.
        return all_cubes

    def identify_boundary(self):
        """Return coordinates of syndrome qubits on primal boundaries.

        For now, omit primal boundaries in the z-direction.
        """
        # TODO: Dual boundaries, z-direction primal if necessary.
        dims = self.dims
        boundaries = np.array(self.boundaries)
        bound_inds = np.where(boundaries == 'primal')[0]
        # Odd indices, which is where primal syndrome qubits are located.
        odds = [range(1, 2 * dims[0], 2),
                range(1, 2 * dims[1], 2),
                range(1, 2 * dims[2], 2)]
        low = []
        high = []
        for ind in bound_inds:
            for i, j in ((0, 1), (0, 2), (1, 2)):
                for tup in it.product(odds[i], odds[j]):
                    l = list(tup)
                    m = list(tup)
                    l.insert(ind, 0)
                    m.insert(ind, 2 * dims[ind])
                    low.append(tuple(l))
                    high.append(tuple(m))
        return list(set(low)) + list(set(high))

    # TODO: slice_coords function that constructs rather than iterates,
    # like the EGraph.


class RHGCube:
    """A class for representing an RHG latice stabilizer cube.

    Arguments:
        G (EGraph): the subgraph of the RHG lattice EGraph induced by
            the (three to six) stabilizer vertices
        physical (NoneType): the list of physical coordinates
            associated with the cube, for plotting purposes. Specified
            during syndrome identification.

    Attributes:
        egraph (EGraph): the corresponding EGraph
        type (str): the type of stabilizer ('six-body' or 'five-body').
    """

    def __init__(self, G):
        """Initialize the RHGCube with its associated egraph."""
        self.egraph = G
        self.physical = None

    @property
    def parity(self):
        """Compute total parity of the cube.

        For now, add the bit_values found in the 'bit_val' attribute
        of the nodes of the EGraph.

        Returns:
            int: the total parity
        """
        G = self.egraph
        bit_vals = [G.nodes[node]['bit_val'] for node in G]
        # TODO: Option for processing the homodyne outcomes first, as
        # below.
        # hom_vals = [G.graph.nodes[node]['hom_val'] for node in G.graph]
        # return basic_translate([np.sum(hom_vals)])
        par = int(np.sum(bit_vals) % 2)
        return par

    def coords(self):
        """Return the coordinates of the syndrome vertices."""
        return list(self.egraph.nodes)

    def xlims(self):
        """Return the smallest and largest x vals of the coordinates."""
        xs = [tup[0] for tup in self.physical]
        xmin, xmax = min(xs), max(xs)
        return xmin, xmax

    def ylims(self):
        """Return the smallest and largest y vals of the coordinates."""
        ys = [tup[1] for tup in self.physical]
        ymin, ymax = min(ys), max(ys)
        return ymin, ymax

    def zlims(self):
        """Return the smallest and largest z vals of the coordinates."""
        zs = [tup[2] for tup in self.physical]
        zmin, zmax = min(zs), max(zs)
        return zmin, zmax

    def midpoint(self):
        """Return the midpoint of the cube."""
        # TODO: Perhaps choose a different point for a five-body
        # X stabilizer, for plotting purposes.
        return (np.average(self.xlims()), np.average(self.ylims()), np.average(self.zlims()))


if __name__ == '__main__':
    # Simple tests
    # Instantiate an RHG latice of a certian distance, with certain
    # boundaries. Draw the EGraph.
    d = 2
    boundaries = 'periodic'
    # boundaries = 'primal'
    # boundaries = 'periodic'
    # for boundaries in it.product(['primal', 'dual', 'periodic'], repeat=3):
    RHG = RHGCode(d, boundaries=boundaries, polarity=True)
    RHG_lattice = RHG.graph
    ax = RHG_lattice.draw(color_nodes=False, color_edges=False, label_indices=False)

    # Check edges between boundaries for periodic boundary conditions.
    all_boundaries = []
    for plane in ('x', 'y', 'z'):
        for i in (0, 2 * d - 1):
            all_boundaries += RHG.graph.slice_coords(plane, i)
    RHG_subgraph = RHG_lattice.subgraph(all_boundaries)
    RHG_subgraph.draw(color_edges=True)

    # Check stabilizer coordinates
    syndrome = RHG.stabilizers
    print('6-body stabilizers :', len(syndrome))
    for i in range(len(syndrome)):
        cube = syndrome[i]
        color = np.random.rand(3)
        for point in cube.egraph:
            x, z, y = point
            ax.scatter(x, z, y, color=color, s=200)
    ax.set_title(str(boundaries).capitalize() + ' boundaries')

    # Check sampling
    delta = 0.001
    # Percent p-squeezed states.
    p_swap = 0
    for sampling_order in ['initial', 'final', 'two-step']:
        model = {'noise': 'grn', 'delta': delta, 'sampling_order': sampling_order}
        CVRHG = CVGraph(RHG_lattice, model=model, p_swap=p_swap)
        CVRHG.measure_hom('p')
        outcomes = CVRHG.hom_outcomes()
        plt.figure(figsize=(16, 9))
        plt.hist(outcomes, bins=100)
        # CVRHG.draw('hom_val_p')
