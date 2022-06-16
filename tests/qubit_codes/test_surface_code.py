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
""""Unit tests for classes and methods in the surface_code module."""

# pylint: disable=too-many-statements

import itertools as it

import networkx as nx
from networkx import fast_gnp_random_graph
from networkx.algorithms.operators import difference
import numpy as np
from numpy.random import default_rng as rng
import pytest

from flamingpy.codes.graphs import EGraph
from flamingpy.codes import RHG_graph, Stabilizer, SurfaceCode, alternating_polarity

# All possible boundary combinations.
all_bound_combs = it.product(["primal", "dual", "periodic"], repeat=3)
all_bound_combs = [np.array(bound) for bound in all_bound_combs]


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
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims
    # Dealing with boundaries
    if boundaries == "finite":
        boundaries = ["primal", "dual", "dual"]
    elif isinstance(boundaries, str):
        boundaries = [boundaries] * 3
    # Define the EGraph of the lattice
    lattice = EGraph(dims=dims, macronodes=macronodes)

    # Constrain the ranges of the coordinates depending on the type
    # of boundaries.
    min_dict = {"primal": 0, "dual": 1, "periodic": 0}
    max_dict = {"primal": 1, "dual": 0, "periodic": 0}
    x_min, y_min, z_min = [min_dict[typ] for typ in boundaries]
    x_max, y_max, z_max = [max_dict[typ] for typ in boundaries]

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

    def color(pol):
        """Polarity-dependent color function: blue for +1, red for -1."""
        return ((pol + 1) // 2) * "b" + abs((pol - 1) // 2) * "r"

    if macronodes:
        lattice.macro_to_micro = {}
        macro_dict = lattice.macro_to_micro

    # Add edges between red points and all their neighbours.
    for point in all_red:
        neighbours = red_neighbours(point)
        if macronodes:
            planetary_bodies = red_neighbours(point, displace=0.1)
            neighbouring_bodies = red_neighbours(point, displace=0.9)
            macro_dict[point] = []
        for i in range(4):
            pol = (-1) ** (polarity * (point[2] + i))
            neighbour = neighbours[i]
            if neighbour in all_green:
                if macronodes:
                    body = planetary_bodies[i]
                    nearby_body = neighbouring_bodies[i]
                    nearby_body = tuple(round(num, 1) for num in nearby_body)
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
            macro_dict[point] = []
        pol = (-1) ** (polarity * (point[1] + 1))
        for i in (0, 1):
            if neighbours[i] in all_green:
                if macronodes:
                    body = planetary_bodies[i]
                    nearby_body = neighbouring_bodies[i]
                    nearby_body = tuple(round(num, 1) for num in nearby_body)
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
            integer_vertices = macro_dict
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
    if macronodes:
        for node in lattice.nodes():
            rounded = tuple(np.round(node).astype(int))
            if not lattice.macro_to_micro.get(rounded):
                lattice.macro_to_micro[rounded] = []
            lattice.macro_to_micro[rounded].append(node)
    return lattice


@pytest.mark.parametrize("d", range(2, 5))
class TestRHGGraph:
    """Test the RHG_graph function."""

    def test_boundary_combinations(self, d):
        """Check that different boundary conditions produce a nonempty
        lattice."""
        for boundaries in it.product(["primal", "dual", "periodic"], repeat=3):
            boundaries = np.array(boundaries)
            RHG_lattice = RHG_graph(d, boundaries)
            assert len(RHG_lattice)
            # The following test requires a modification to RHG_graph_old,
            # since the new method creates a lattice with a different shape.
            # assert not set(RHG_lattice.edges) - set(RHG_graph_old(d, boundaries).edges)
        for boundaries in ["primal", "dual", "periodic", "open_primal", "open_dual"]:
            RHG_lattice = RHG_graph(d, boundaries)
            assert len(RHG_lattice)
            # The following test requires a modification to RHG_graph_old,
            # since the new method creates a lattice with a different shape.
            # assert not set(RHG_lattice_finite.edges) - set(RHG_graph_old(d, "finite").edges)

    def test_periodic_boundaries(self, d):
        """Test whether periodic boundary conditions produce a lattice with the
        expected size."""
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

        assert not set(RHG_lattice.edges) - set(RHG_graph(d, "periodic").edges)

    @pytest.mark.parametrize("boundaries", all_bound_combs)
    def test_polarity(self, d, boundaries):
        """Test whether lattice polarity behaves as expected."""
        RHG_reduced = RHG_graph(d, boundaries)
        RHG_macro = RHG_reduced.macronize()
        for graph in (RHG_reduced, RHG_macro):
            weights = list(nx.get_edge_attributes(graph, "weight").values())
            assert np.all(np.array(weights) == 1)
        # Check alternating polarity.
        RHG_reduced = RHG_graph(d, boundaries, polarity=alternating_polarity)
        RHG_macro = RHG_reduced.macronize()
        for graph in (RHG_reduced, RHG_macro):
            weights = np.array(list(nx.get_edge_attributes(graph, "weight").values()))
            # Check that all weights are either 1 or -1
            assert np.all((weights == 1) | (weights == -1))
            # For all-periodic boundaries, check that each node in the
            # reduced lattice has two +1 and two -1 edges; for
            # the macronode lattice, that each macronode contains
            # four edges, two with weight +1 and two with weight -1.
            if tuple(boundaries) == ("periodic", "periodic", "periodic"):
                if graph == RHG_reduced:
                    adj_mat = graph.adj_generator(sparse=False)
                    for row in adj_mat:
                        assert np.count_nonzero(row == 1) == 2
                        assert np.count_nonzero(row == -1) == 2
                if graph == RHG_macro:
                    for macronode in RHG_macro.macro_to_micro:
                        micronodes = RHG_macro.macro_to_micro[macronode]
                        macro_edges = RHG_macro.edges(micronodes)
                        weights = np.array(
                            [RHG_macro.edges[edge]["weight"] for edge in macro_edges]
                        )
                        assert np.count_nonzero(weights == 1) == 2
                        assert np.count_nonzero(weights == -1) == 2

    def test_macronodes(self, d):
        """Check that the macronode lattice was generated as expected."""
        # Tests for periodic boundaries.
        for boundaries in ["periodic"]:
            RHG_reduced = RHG_graph(d, boundaries)
            macronode_lattice = RHG_reduced.macronize()
            RHG_macronodes = macronode_lattice.macro_to_micro
            # Check that the reduced lattice has the same node set
            # as the nodes stored in the macro attribute of the
            # macronode lattice.
            assert not set(RHG_macronodes) - set(RHG_reduced)
            # Check that the macronode lattice in the case of periodic
            # boundaries is exactly four times larger than the reduced
            # lattice.
            assert len(macronode_lattice) == 4 * len(RHG_reduced)
            # Check that there is no difference in the edge set
            # between the macronode_lattice and the macronode lattice
            # create by the RHG_graph_old function.
            assert not difference(
                macronode_lattice, RHG_graph_old(d, boundaries, macronodes=True)
            ).edges
        # Tests for all boundaries.
        for boundaries in it.product(["primal", "dual", "periodic"], repeat=3):
            boundaries = np.array(boundaries)
            RHG_reduced = RHG_graph(d, boundaries)
            macronode_lattice = RHG_reduced.macronize(pad_boundary=True)
            # Check that the macronode lattice in the case of padded
            # boundaries is exactly four times larger than the reduced
            # lattice.
            assert len(macronode_lattice) == 4 * len(RHG_reduced)


code_params = it.product(range(2, 5), ["primal", "dual"], ["open", "periodic"])


@pytest.fixture(scope="module", params=code_params)
def surface_code(request):
    """A handy function to define the surface code for use in this module."""
    distance, ec, boundaries = request.param
    return SurfaceCode(distance, ec, boundaries)


class TestSurfaceCode:
    """Test the SurfaceCode class."""

    def test_init(self, surface_code):
        """Check the proper initialization of SurfaceCode."""
        distance = surface_code.distance
        assert surface_code.dims == (distance, distance, distance)
        assert surface_code.ec
        assert surface_code.bound_str
        assert list(surface_code.boundaries)
        for ec in surface_code.ec:
            attributes = [
                "stabilizers",
                "syndrome_inds",
                "syndrome_coords",
                "bound_points",
                "stab_graph",
            ]
            for att in attributes:
                getattr(surface_code, ec + "_" + att)
            for att in ["all_syndrome_inds", "all_syndrome_coords"]:
                getattr(surface_code, att)

    def test_stabilizers(self, surface_code):
        """Check whether all stabilizers were generated as expected."""
        for ec in surface_code.ec:
            d = surface_code.distance
            cubes = getattr(surface_code, ec + "_stabilizers")
            # Check that there are a correct number of stabilizer elements
            # depending on the boundary conditions.
            if surface_code.bound_str.startswith("open"):
                if len(surface_code.ec) == 2 and ec == "dual":
                    assert len(cubes) == (d - 1) ** 2 * d
                else:
                    assert len(cubes) == d**2 * (d - 1)
            elif surface_code.bound_str == "periodic":
                assert len(cubes) == d**3
            # Check that each stabilizer has 6 corresponing physical
            # vertices, even if it's an n-body stabilizer with n < 6.
            for cube in cubes:
                assert len(cube.physical) == 6
                # For periodic boundaries, check that there are
                # only 6-body stabilizers.
                if surface_code.bound_str == "periodic":
                    assert len(cube.egraph) == 6

    def test_syndrome(self, surface_code):
        """Check that the syndrome coordinates have been appropriately
        identified."""
        # Check that syndrome coordinates in the code class match the
        # coordinates in the stabilizers.
        for ec in surface_code.ec:
            syndrome_coords = getattr(surface_code, ec + "_syndrome_coords")
            stabilizers = getattr(surface_code, ec + "_stabilizers")
            syndrome_from_stabes = [set(cube.egraph) for cube in stabilizers]
            flattened_syndrome = [a for b in syndrome_from_stabes for a in b]
            assert not set(syndrome_coords) - set(flattened_syndrome)
            # Check that the total number of syndrome coordinates is correct.
            d = surface_code.distance
            if surface_code.bound_str.startswith("open"):
                if ec == "dual" and len(surface_code.ec) == 2:
                    assert len(syndrome_coords) == (d - 1) ** 3 + 2 * (d - 1) * d**2
                else:
                    assert len(syndrome_coords) == d**3 + 2 * d * (d - 1) ** 2
            elif surface_code.bound_str == "periodic":
                assert len(syndrome_coords) == 3 * d**3

    def test_boundary(self, surface_code):
        """Check that the length of boundary coordinates is correct."""
        for ec in surface_code.ec:
            bound_points = getattr(surface_code, ec + "_bound_points")
            if surface_code.bound_str.startswith("open"):
                d = surface_code.distance
                if ec == "dual" and len(surface_code.ec) == 2:
                    assert len(bound_points) == 2 * d * (d - 1)
                else:
                    assert len(bound_points) == 2 * d**2
            elif surface_code.bound_str == "periodic":
                assert not bound_points


code_params_rectangular = it.product(
    rng().integers(low=2, high=7, size=3), ["primal", "dual"], ["open", "periodic"]
)


class TestRectangularSurfaceCode:
    """Test the SurfaceCode class for non-cubic lattices."""

    @pytest.fixture(scope="module", params=code_params_rectangular)
    def test_init_rectangular(self, request):
        """Check the proper initialization of SurfaceCode."""
        # initialize a rectangular surface code
        dx, dy, dz, ec, boundaries = request.param
        rect_sc = SurfaceCode((dx, dy, dz), ec, boundaries)

        # assert that the surface code is initialized correctly
        assert rect_sc.dims == (dx, dy, dz)
        assert rect_sc.ec
        assert rect_sc.bound_str
        assert list(rect_sc.boundaries)
        for ec in rect_sc.ec:
            attributes = [
                "stabilizers",
                "syndrome_inds",
                "syndrome_coords",
                "bound_points",
                "stab_graph",
            ]
            for att in attributes:
                getattr(rect_sc, ec + "_" + att)
            for att in ["all_syndrome_inds", "all_syndrome_coords"]:
                getattr(rect_sc, att)


class TestStabilizer:
    """Test the Stabilizer class."""

    def test_parity(self):
        """Check that parity is properly computed."""
        n = rng().integers(10)
        random_graph = fast_gnp_random_graph(n, 0.5)
        bit_list = []
        for node in random_graph:
            random_bit = rng().integers(2)
            bit_list += [random_bit]
            random_graph.nodes[node]["bit_val"] = random_bit
        random_cube = Stabilizer(random_graph)
        assert random_cube.parity == sum(bit_list) % 2
