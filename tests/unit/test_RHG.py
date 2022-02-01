# Copyright 2022 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
""""Unit tests for RHG classes and methods in RHG.py."""
from flamingpy.RHG import RHG_graph, RHGCode, RHGCube
from flamingpy.graphstates import EGraph
from numpy.random import default_rng as rng
from networkx import fast_gnp_random_graph
from networkx.algorithms.operators import difference

import numpy as np
import itertools as it
import pytest


def RHG_graph_old(dims, boundaries="finite", macronodes=False, polarity=False):
    """Create an EGraph of a dims-dimensional RHG lattice.

    Generate an RHG lattice with dimensions given by dims (an integer
    denoting the number of stabilizer cubes in each direciton,
    or a tuple specifying all three dimensions). By default, a useful
    set of finite boundary conditions assumed, but any combination can
    be specified.

    Included in this test file for a cross-check with new approaches
    of constructing the lattice.

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
    # TODO: Compactify construction by identifying syndrome qubits;
    # potentially use NetworkX lattice generators.
    # first.
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims
    # Dealing with boundaries
    if boundaries == "finite":
        boundaries = ["primal", "dual", "dual"]
    elif type(boundaries) == str:
        boundaries = [boundaries] * 3
    # Define the EGraph of the lattice
    if macronodes:
        indexer = "macronodes"
    else:
        indexer = "default"
    lattice = EGraph(dims=dims, indexer=indexer, macronodes=macronodes)

    # Constrain the ranges of the coordinates depending on the type
    # of boundaries.
    min_dict = {"primal": 0, "dual": 1, "periodic": 0}
    max_dict = {"primal": 1, "dual": 0, "periodic": 0}
    x_min, y_min, z_min = [min_dict[typ] for typ in boundaries]
    x_max, y_max, z_max = [max_dict[typ] for typ in boundaries]

    # TODO: Change notation to data / ancillae notation
    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [
        (2 * i + 1, 2 * j + 1, 2 * k)
        for (i, j, k) in it.product(range(nx), range(ny), range(z_min, nz + z_max))
    ]
    odd_red = [
        (2 * i, 2 * j, 2 * k + 1)
        for (i, j, k) in it.product(range(x_min, nx + x_max), range(y_min, ny + y_max), range(nz))
    ]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [
        (2 * i + 1, 2 * j, k)
        for (i, j, k) in it.product(
            range(nx), range(y_min, ny + y_max), range(z_min, 2 * nz + z_max)
        )
    ]
    odd_green = [
        (2 * i, 2 * j + 1, k)
        for (i, j, k) in it.product(
            range(x_min, nx + x_max), range(ny), range(z_min, 2 * nz + z_max)
        )
    ]
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
    color = lambda pol: ((pol + 1) // 2) * "b" + abs((pol - 1) // 2) * "r"

    if macronodes:
        macro_graph = lattice.macro

    # Add edges between red points and all their neighbours.
    for point in all_red:
        neighbours = red_neighbours(point)
        if macronodes:
            planetary_bodies = red_neighbours(point, displace=0.1)
            neighbouring_bodies = red_neighbours(point, displace=0.9)
            macro_graph.add_node(point, micronodes=[])
            # TODO: Append micronodes to macro_graph so indexer does
            # not have to be run.
        for i in range(4):
            pol = (-1) ** (polarity * (point[2] + i))
            neighbour = neighbours[i]
            if neighbour in all_green:
                if macronodes:
                    body = planetary_bodies[i]
                    nearby_body = neighbouring_bodies[i]
                    nearby_body = tuple([round(num, 1) for num in nearby_body])
                    lattice.add_edge(body, nearby_body, weight=pol, color=color(pol))
                    lattice.nodes[body]["color"] = "red"
                    lattice.nodes[nearby_body]["color"] = "green"
                else:
                    lattice.add_edge(point, neighbour, weight=pol, color=color(pol))
                    lattice.nodes[point]["color"] = "red"

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
                    lattice.nodes[body]["color"] = "green"
                    lattice.nodes[nearby_body]["color"] = "green"
                else:
                    lattice.add_edge(point, neighbours[i], weight=pol, color=color(pol))
                    lattice.nodes[point]["color"] = "green"

    # Dealing with periodic boundary conditions.
    bound_arr = np.array(boundaries)
    periodic_inds = np.where(bound_arr == "periodic")[0]
    for ind in periodic_inds:
        # First and last slices of the lattice in the direction
        # specified by ind.
        if macronodes:
            integer_vertices = macro_graph.nodes
        else:
            integer_vertices = lattice.nodes
        low_slice = {point for point in integer_vertices if point[ind] == 0}
        high_slice = {point for point in integer_vertices if point[ind] == 2 * dims[ind] - 1}
        if ind in (0, 1):
            low_reds = all_red & low_slice
            high_reds = all_red & high_slice
            # Connect red in first slice to greens in last slice.
            for point in low_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                high_green = list(point)
                high_green[ind] = 2 * dims[ind] - 1
                high_green = tuple(high_green)
                if macronodes:
                    point = red_neighbours(point, displace=0.1)[1 - ind]
                    high_green = red_neighbours(high_green, displace=0.1)[-1 - ind]
                lattice.add_edge(point, high_green, weight=pol, color=color(pol))
                lattice.nodes[point]["color"] = "red"
                lattice.nodes[high_green]["color"] = "green"
            # Connect reds in last slice to greens in first slice.
            for point in high_reds:
                pol = (-1) ** (polarity * (point[2] + ind + 1))
                low_green = list(point)
                low_green[ind] = 0
                low_green = tuple(low_green)
                if macronodes:
                    point = red_neighbours(point, displace=0.1)[-1 - ind]
                    low_green = red_neighbours(low_green, displace=0.1)[1 - ind]
                lattice.add_edge(point, low_green, weight=pol, color=color(pol))
                lattice.nodes[point]["color"] = "red"
                lattice.nodes[low_green]["color"] = "green"
        # If periodic in z direction, connect greens in first slice with
        # greens in last slice.
        if ind == 2:
            low_greens = all_green & low_slice
            for point in low_greens:
                pol = (-1) ** (polarity * (point[1] + 1))
                high_green = list(point[:])
                high_green[ind] = 2 * dims[ind] - 1
                high_green = tuple(high_green)
                if macronodes:
                    point = green_neighbours(point, displace=0.1)[0]
                    high_green = green_neighbours(high_green, displace=0.1)[1]
                lattice.add_edge(point, high_green, weight=pol, color=color(pol))
                lattice.nodes[point]["color"] = "green"
                lattice.nodes[high_green]["color"] = "green"

    return lattice


@pytest.mark.parametrize("d", range(2, 5))
class TestRHGGraph:
    """Test the RHG_graph function."""

    def test_boundary_combinations(self, d):
        # Check that different boundary conditions produce a nonempty
        # lattice, and equal the lattice constructed via an alternative
        # approach.
        for boundaries in it.product(["primal", "dual", "periodic"], repeat=3):
            RHG_lattice = RHG_graph(d, boundaries)
            assert len(RHG_lattice)
            assert not difference(RHG_lattice, RHG_graph_old(d, boundaries)).edges
        RHG_lattice_finite = RHG_graph(d, "finite")
        assert len(RHG_lattice_finite)
        assert not difference(RHG_lattice_finite, RHG_graph_old(d, "finite")).edges

    def test_periodic_boundaries(self, d):
        RHG_lattice = RHG_graph(d, "periodic")
        assert len(RHG_lattice) == 6 * (d**3)

        all_boundaries = []
        for plane in ("x", "y", "z"):
            for i in (0, 2 * d - 1):
                all_boundaries += RHG_lattice.slice_coords(plane, i)
            # For macronode lattice, change to (-0.1, 2 * d - 0.9)
        # Check that every boundary coordinate has four neighbours.
        for point in all_boundaries:
            assert len(RHG_lattice[point]) == 4

        assert not difference(RHG_lattice, RHG_graph_old(d, "periodic")).edges

    def test_polarity(self, d):
        pass

    def test_macronodes(self, d):
        for boundaries in {"finite", "periodic"}:
            macronode_lattice = RHG_graph(d, boundaries, macronodes=True)
            macronode_lattice.index_generator()
            RHG_macronodes = macronode_lattice.macro
            RHG_reduced = RHG_graph(d, boundaries)
            # Check that the reduced lattice has the same node set
            # as the nodes stored in the macro attribute of the
            # macronode lattice.
            assert not set(RHG_macronodes) - set(RHG_reduced)
            # Check that the macronode lattice in the case of periodic
            # boundaries is exactly four times larger than the reduced
            # lattice.
            if boundaries == "periodic":
                assert len(macronode_lattice) == 4 * len(RHG_reduced)
            assert not difference(
                macronode_lattice, RHG_graph_old(d, boundaries, macronodes=True)
            ).edges


code_params = it.product(range(2, 5), ["primal", "dual", "finite", "periodic"])


@pytest.fixture(scope="module", params=code_params)
def RHG_code(request):
    distance, boundaries = request.param
    return RHGCode(distance, error_complex="primal", boundaries=boundaries)


class TestRHGCode:
    """ "Test the RHGCode class."""

    def test_init(self, RHG_code):
        # Check proper initialization
        distance = RHG_code.distance
        assert RHG_code.dims == (distance, distance, distance)
        assert RHG_code.complex
        assert RHG_code.boundaries
        assert RHG_code.syndrome_inds

    def test_stabilizers(self, RHG_code):
        cubes = RHG_code.stabilizers
        # Check that there are d^3 stabilizers
        assert len(cubes) == RHG_code.distance**3
        # Check that each stabilizer has 6 corresponing physical
        # vertices, no matter if it's a 3, 4, 5, or 6 body stabilizer.
        for cube in cubes:
            assert len(cube.physical) == 6
        # For all primal or periodic boundaries, check that there are
        # only 6-body stabilizers.
        if RHG_code.boundaries in [["primal"] * 3, ["periodic"] * 3]:
            assert len(cube.egraph) == 6
        # TODO: "finite" boundary test

    def test_syndrome(self, RHG_code):
        syndrome_coords = RHG_code.syndrome_coords
        syndrome_from_stabes = [set(cube.egraph) for cube in RHG_code.stabilizers]
        flattened_syndrome = [a for b in syndrome_from_stabes for a in b]
        # Check that syndrome coordinates in the code class match the
        # coordinates in the stabilizers.
        assert not set(syndrome_coords) - set(flattened_syndrome)

    def test_boundary(self, RHG_code):
        if "primal" in RHG_code.boundaries:
            assert RHG_code.boundary_coords
        else:
            assert not RHG_code.boundary_coords


class TestRHGCube:
    """Test the RHGCube class."""

    def test_parity(self):
        # Check that parity is properly computed.
        n = rng().integers(10)
        random_graph = fast_gnp_random_graph(n, 0.5)
        bit_list = []
        for node in random_graph:
            random_bit = rng().integers(2)
            bit_list += [random_bit]
            random_graph.nodes[node]["bit_val"] = random_bit
        random_cube = RHGCube(random_graph)
        assert random_cube.parity == sum(bit_list) % 2
