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


def RHG_graph(dims, boundaries='finite', macronodes=False, polarity=False):
    """Create an EGraph of a dims-dimensional RHG lattice.

    Generate an RHG lattice with dimensions given by dims (an integer
    denoting the number of stabilizer cubes in each direciton,
    or a tuple specifying all three dimensions). By default, a useful
    set of finite boundary conditions assumed, but any combination can
    be specified.

    Args:
        dims (int or list): the dimensions of the lattice (the
            number of stabilizer cubes in each direction, complete or
            incomplete).
        boundaries (str or list-type, optional): the boundaries in each
            direction. We use primal/smooth and dual/rough
            interchangeably. If a string, 'primal', 'dual', 'periodic',
            assumes those boundaries in all three directions; or,
            accepts a list that specifies which boundary in which
            direction. Set  to 'finite' or ['primal', 'dual', 'dual']
            by default.
        polarity (bool): if True, make some edges have weight -1 to be
            helpful for some noise models, e.g. in the CV context. Set
            edge 'color' attribute to blue/red for -1/+1 weights. By
            default, all edges have weight +1.

    Returns:
        EGraph: the RHG lattice.
    """
    # TODO: Compactify construction by identifying syndrome qubits
    # first.
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims
    # Dealing with boundaries
    if boundaries == 'finite':
        boundaries = ['primal', 'dual', 'dual']
    elif type(boundaries) == str:
        boundaries = [boundaries] * 3
    # Define the EGraph of the lattice
    if macronodes:
        indexer = 'macronodes'
    else:
        indexer = 'default'
    lattice = EGraph(dims=dims, boundaries=boundaries, indexer=indexer)

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

    def red_neighbours(p, displace=1):
        """Coordinates of all potential neighbours of red vertices."""
        right = (p[0] + displace, p[1], p[2])
        top = (p[0], p[1] + displace, p[2])
        left = (p[0] - displace, p[1], p[2])
        bottom = (p[0], p[1] - displace, p[2])
        return [bottom, left, top, right]

    def green_neighbours(p, displace=1):
        """Coordinates of all potential neighbours of green vertices."""
        return [(p[0], p[1], p[2] - displace), (p[0], p[1], p[2] + displace)]

    # Polarity-dependent color function: blue for +1, red for -1.
    color = lambda pol: ((pol + 1) // 2) * 'b' + abs((pol - 1) // 2) * 'r'

    if macronodes:
        macro_graph = lattice.macro

    # Add edges between red points and all their neighbours.
    for point in all_red:
        neighbours = red_neighbours(point)
        if macronodes:
            planetary_bodies = red_neighbours(point, displace=0.1)
            neighbouring_bodies = red_neighbours(point, displace=0.9)
            macro_graph.add_node(point, micronodes=[])
        for i in range(4):
            pol = (-1) ** (polarity * (point[2] + i))
            neighbour = neighbours[i]
            if neighbour in all_green:
                if macronodes:
                    body = planetary_bodies[i]
                    nearby_body = neighbouring_bodies[i]
                    nearby_body = tuple([round(num, 1) for num in nearby_body])
                    lattice.add_edge(body, nearby_body, weight=pol, color=color(pol))
                    lattice.nodes[body]['color'] = 'red'
                    lattice.nodes[nearby_body]['color'] = 'green'
                    # macro_graph.nodes[point]['micronodes'].append(body)
                    # macro_graph.nodes[neighbour]['micronodes'].append(nearby_body)
                else:
                    lattice.add_edge(point, neighbour, weight=pol, color=color(pol))
                    lattice.nodes[point]['color'] = 'red'

    # Add edges between green points and all their neighbours.
    for point in all_green:
        neighbours = green_neighbours(point)
        if macronodes:
            planetary_bodies = green_neighbours(point, displace=0.1)
            neighbouring_bodies = green_neighbours(point, displace=0.9)
            macro_graph.add_node(point, micronodes=[])
        pol = (-1) ** (polarity * (point[1] + 1))
        for i in (0, 1):
            if neighbours[i] in all_green:
                if macronodes:
                    body = planetary_bodies[i]
                    nearby_body = neighbouring_bodies[i]
                    nearby_body = tuple([round(num, 1) for num in nearby_body])
                    lattice.add_edge(body, nearby_body, weight=pol, color=color(pol))
                    lattice.nodes[body]['color'] = 'green'
                    lattice.nodes[nearby_body]['color'] = 'green'
                    # macro_graph.nodes[point]['micronodes'].append(body)
                else:
                    lattice.add_edge(point, neighbours[i], weight=pol, color=color(pol))
                    lattice.nodes[point]['color'] = 'green'

    # Dealing with periodic boundary conditions.
    bound_arr = np.array(boundaries)
    periodic_inds = np.where(bound_arr == 'periodic')[0]
    for ind in periodic_inds:
        # First and last slices of the lattice in the direction
        # specified by ind.
        if macronodes:
            integer_vertices = [point for point in macro_graph.nodes]
        else:
            integer_vertices = [point for point in lattice.nodes]
        low_slice = set([point for point in integer_vertices if point[ind] == 0])
        high_slice = set([point for point in integer_vertices if point[ind] == 2*dims[ind] - 1])
        if ind in (0, 1):
            low_reds = all_red.intersection(low_slice)
            high_reds = all_red.intersection(high_slice)
            # Connect red in first slice to greens in last slice.
            for point in low_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                high_green = list(point)
                high_green[ind] = 2 * dims[ind] - 1
                high_green = tuple(high_green)
                if macronodes:
                    low_macro = macro_graph.nodes[point]['micronodes']
                    high_macro = macro_graph.nodes[high_green]['micronodes']
                    point = red_neighbours(point, displace=0.1)[1 - ind]
                    high_green = red_neighbours(high_green, displace=0.1)[-1 - ind]
                    # low_macro.add_node(point)
                    # high_macro.add_node(high_green)
                lattice.add_edge(point, high_green, weight=pol, color=color(pol))
                lattice.nodes[point]['color'] = 'red'
                lattice.nodes[high_green]['color'] = 'green'
            # Connect reds in last slice to greens in first slice.
            for point in high_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                low_green = list(point)
                low_green[ind] = 0
                low_green = tuple(low_green)
                if macronodes:
                    high_macro = macro_graph.nodes[point]['micronodes']
                    low_macro = macro_graph.nodes[low_green]['micronodes']
                    point = red_neighbours(point, displace=0.1)[-1 - ind]
                    low_green = red_neighbours(low_green, displace=0.1)[1 - ind]
                    # high_macro.add_node(point)
                    # low_macro.add_node(low_green)
                lattice.add_edge(point, low_green, weight=pol, color=color(pol))
                lattice.nodes[point]['color'] = 'red'
                lattice.nodes[low_green]['color'] = 'green'
        # If periodic in z direction, connect greens in first slice with
        # greens in last slice.
        if ind == 2:
            low_greens = all_green.intersection(low_slice)
            for point in low_greens:
                pol = (-1) ** (polarity * (point[1] + 1))
                high_green = list(point[:])
                high_green[ind] = 2 * dims[ind] - 1
                high_green = tuple(high_green)
                if macronodes:
                    low_macro = macro_graph.nodes[point]['micronodes']
                    high_macro = macro_graph.nodes[high_green]['micronodes']
                    point = green_neighbours(point, displace=0.1)[0]
                    high_green = green_neighbours(high_green, displace=0.1)[1]
                    # low_macro.add_node(point)
                    # high_macro.add_node(high_green)
                lattice.add_edge(point, high_green, weight=pol, color=color(pol))
                lattice.nodes[point]['color'] = 'green'
                lattice.nodes[high_green]['color'] = 'green'

    return lattice


class RHGCode:
    """A class for representing the RHG code.

    Represent the Raussendorf-Harrington-Goyal (RHG) code, which can
    be viewed as a foliated surface code, in its measurement-based
    description. By specifying the distance and choice of boundaries,
    store the graph state corresponding to the code as an EGraph, the
    set of stabilizers (RHGCube objects), and the boundary vertices.

    Attributes:
        distance (int): the code distance. Corresponds to the number
            of stabilizer cubes (complete or incomplete) in each
            x, y, z direction.
        dims (tup): a tuple of the spatial extent in x, y, z.
        complex (str): the error complex (primal or dual). For now
            only primal implemented.
        boundaries (list or str): the boundaries in x, y, z. We use the
            identification primal = smooth and dual = rough, to align
            with surface code terminology. Available choices in the
            order x, y, z are 'finite' (primal, dual, dual)

                'finite': primal, dual, dual
                'b': 'b', 'b', 'b'
                ['b1', 'b2', 'b3']: b1, b2, b3.

        _polarity (bool): if True, lattice is constructed with two
            edges of +1 weights and perpendicular edges with -1 weights
            for noise cancellation in a subsequent CV lattice.
        graph (EGraph): the EGraph correspond to the code.
        stabilizers (list of RHGCubes): the stabilizer elements of the
            code according to the error complex.
        syndrome_coords (list of tup): the coordinates of the syndrome
            vertices, according to the error complex.
        boundary_coords (list of tup): the coordinates of the boundary
            according to the error complex.
    """

    def __init__(self,
                 distance,
                 error_complex='primal',
                 boundaries='finite',
                 polarity=False):
        """Initialize the RHG code."""
        # TODO: Check code distance convention.
        self.distance = distance
        self.dims = (distance, distance, distance)
        self.complex = error_complex
        if boundaries == 'finite':
            self.boundaries = ['primal', 'dual', 'dual']
        elif type(boundaries) == str:
            self.boundaries = [boundaries] * 3
        self._polarity = polarity

        self.graph = RHG_graph(self.dims, boundaries=self.boundaries, polarity=polarity)
        self.graph.index_generator()
        # The following line also defines the self.syndrome_coords
        # attribute.
        self.stabilizers = self.identify_stabilizers(self.complex)
        self.syndrome_inds = [self.graph.to_indices[point] for point in self.syndrome_coords]
        self.boundary_coords = self.identify_boundary(self.complex)

    def identify_stabilizers(self, error_complex='primal'):
        """Return the syndrome coordinates for the RHG lattice G.

        Generate a list of RHGCube objects containing coordinates of
        all the stabilizer elements according to error_complex, starting
        with six-body X stabilizers, followed by five, four, and three,
        if required by choice of boundary.
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

        # TODO: Implement dual error complex.
        if error_complex == 'primal':
            # All potential six-body stabilizers
            all_six_bodies = [[
                (2*i, 2*j + 1, 2*k + 1),
                (2*i + 1, 2*j, 2*k + 1),
                (2*i + 1, 2*j + 1, 2*k),
                (2*i + 2, 2*j + 1, 2*k + 1),
                (2*i + 1, 2*j + 2, 2*k + 1),
                (2*i + 1, 2*j + 1, 2*k + 2)] for (i, j, k) in inds]

        all_cubes = []
        syndrome_coords = []

        periodic_inds = np.where(boundaries == 'periodic')[0]
        dual_inds = np.where(boundaries == 'dual')[0]
        for stabe in all_six_bodies:
            actual_stabe = list(set(stabe).intersection(set(G)))
            if len(actual_stabe) == 6:
                cube = RHGCube(G.subgraph(actual_stabe))
                cube.physical = stabe
                syndrome_coords += actual_stabe
                all_cubes.append(cube)
            if len(actual_stabe) == 5:
                for ind in (0, 1, 2):
                    lowest_point = list(stabe[ind])
                    highest_point = stabe[3 + ind]
                    if lowest_point[ind] == (2 * dims[ind] - 2):
                        if ind in dual_inds:
                            cube = RHGCube(G.subgraph(actual_stabe))
                            cube.physical = stabe
                            syndrome_coords += actual_stabe
                            all_cubes.append(cube)
                        if ind in periodic_inds:
                            lowest_point[ind] = 0
                            virtual_point = tuple(lowest_point)
                            actual_stabe += [virtual_point]
                            cube = RHGCube(G.subgraph(actual_stabe))
                            cube.physical = stabe
                            syndrome_coords += actual_stabe
                            all_cubes.append(cube)
                    if highest_point[ind] == 2 and ind in dual_inds:
                        cube = RHGCube(G.subgraph(actual_stabe))
                        cube.physical = stabe
                        syndrome_coords += actual_stabe
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
                syndrome_coords += actual_stabe
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
                syndrome_coords += actual_stabe
                all_cubes.append(cube)
        # Dealing with six-body X stabilizers on perodic boundaries,
        # and five-body X stabilizers on dual boundaries.
        self.syndrome_coords = syndrome_coords
        return all_cubes

    def identify_boundary(self, error_complex='primal'):
        """Obtain coordinates of syndrome qubits on the boundary.

        The relevant boundary is determined by the error_complex string.
        """
        # TODO: Dual boundaries.
        dims = self.dims
        boundaries = np.array(self.boundaries)
        if error_complex == 'primal':
            # Odd indices, which is where primal syndrome qubits are located.
            odds = [range(1, 2 * dims[0], 2),
                    range(1, 2 * dims[1], 2),
                    range(1, 2 * dims[2], 2)]
        low = []
        high = []
        bound_inds = np.where(boundaries == error_complex)[0]
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
    """A class for representing an RHG lattice stabilizer cube.

    Arguments:
        G (EGraph): the subgraph of the RHG lattice EGraph induced by
            the (three to six) stabilizer vertices
        physical (NoneType): the list of physical coordinates
            associated with the cube, for plotting purposes. Specified
            during syndrome identification.

    Attributes:
        egraph (EGraph): the corresponding EGraph.
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
    # Check maronode lattice
    # RHG_lattice = RHG_graph(d, boundaries=boundaries, macronodes=True)
    ax = RHG_lattice.draw(color_nodes=False, color_edges=False, label='index')

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
    CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)
    for sampling_order in ['initial', 'final', 'two-step']:
        model = {'noise': 'grn', 'delta': delta, 'sampling_order': sampling_order}
        CVRHG.apply_noise(model)
        CVRHG.measure_hom('p')
        outcomes = CVRHG.hom_outcomes()
        plt.figure(figsize=(16, 9))
        plt.hist(outcomes, bins=100)
        CVRHG.draw(label='hom_val_p')
