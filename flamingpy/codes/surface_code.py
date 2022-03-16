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
"""Class for the measurement-based surface code and related functions."""

import itertools as it

import numpy as np

from flamingpy.codes import Stabilizer
from flamingpy.codes.graphs import EGraph
from flamingpy.codes.graphs import NxStabilizerGraph, RxStabilizerGraph


def alternating_polarity(edge):
    """Return +1 or -1 depending on the vertices that form edge edge.

    Help with the assignment of edge weights (i.e. polarity) for
    the RHG graph. This particular alternating pattern ensures that
    every vertex has two +1 and two -1 weight edges incident on it. The
    precise combination depends on the direction of the edge and the
    coordinates of the vertices. This pattern may be helpful to reduce, e.g.
    CV noise.

    Args:
        edge (list-type): a pair of tuples, denoting lattice vertices.

    Returns:
        int: +1 or -1.
    """
    point1, point2 = np.array(edge[0]), np.array(edge[1])
    direction = np.where(point2 - point1)[0][0]
    if direction == 0:
        pol = (-1) ** point1[1]
    elif direction == 1:
        pol = (-1) ** point1[2]
    elif direction == 2:
        pol = (-1) ** point1[0]
    else:
        print("Vertices must be separated by one unit on the integer lattice.")
        return 1
    return pol


def dual_neighbours(p, displace=1):
    """All the dual neighbours of primal vertex p in the RHG lattice.

    A helper function for RHG_graph. Given a primal vertex p, returns
    the coordinates of all the dual neighbours. Assumes each neighbour
    is 1 unit away by default.

    Args:
        p (tuple): the coordinates of the primal vertex.
        displace (float): how much to displace the neighbour by. Useful
            to change when creating maronodes.

    Returns:
        List[Tuple]: the coordinates of the four neighbours.
    """
    x, y, z = p[0], p[1], p[2]
    top = (x, y + displace, z)
    bottom = (x, y - displace, z)
    left = (x - displace, y, z)
    right = (x + displace, y, z)
    if z % 2:
        front = (x, y, z + displace)
        back = (x, y, z - displace)
        if x % 2:
            return [back, left, front, right]
        else:
            return [back, top, front, bottom]
    else:
        return [bottom, left, top, right]


def str_to_bound(bound_name):
    """Return a list of x-y-z boundaries corresponding to bound_name.

    The options are:

        'open_primal': [primal, dual, dual]
        'open_dual': [primal, dual, primal]
        '{b}': [b, b, b], where b can be 'primal', 'dual', or 'periodic'.
    """
    if bound_name == "open_primal":
        boundaries = ["primal", "dual", "dual"]
    elif bound_name == "open_dual":
        boundaries = ["primal", "dual", "primal"]
    elif bound_name in ("primal", "dual", "periodic"):
        boundaries = [bound_name] * 3
    elif not isinstance(bound_name, str):
        print("Boundary type must be string.")
        raise Exception
    return np.array(boundaries)


def RHG_graph(
    dims,
    boundaries="primal",
    polarity=None,
):
    """Return an EGraph of a dims-dimensional RHG lattice.

    Generate a Raussendorf-Harrington-Goyal (RHG) lattice, which can
    be viewed as the measurement-based version or foliation of
    the surface code, with specified dimensions and boundary types.

    Args:
        dims (int or list-type): the dimensions of the lattice. If int,
            generates a cube corresponding to a code of distance dims.
            If a three-element list [dx, dy, dz], assumes distances
            dx, dy, dz in x, y, z directions, respectively.
        boundaries (str or list-type, optional): the boundary types
            in x, y, z. We use the identification primal = smooth and
            dual = rough, to align with surface code terminology.
            Available choices in the order x, y, z are:

                'open_primal': primal, dual, dual
                'open_dual':,  primal, dual, primal
                '{b}': b, b, b,
                ['{b1}', '{b2}', '{b3}']: b1, b2, b3,

            where each b above can be 'primal', 'dual', or 'periodic'.
            By default, 'primal'  is used (i.e. ['primal', 'primal', 'primal']).
        polarity (func): a function that specifies edge weights. It
            must be of the following form:

                polarity(edge) = weight.

            If not supplied, assumes all edges have weight 1.
    Returns:
        EGraph: the RHG lattice.
    """
    # Create an EGraph with the graph attribute 'dims' (used for
    # plotting purposes.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    G = EGraph(dims=dims)

    # Dealing with boundaries.
    if isinstance(boundaries, str):
        boundaries = str_to_bound(boundaries)

    # Locations of all primal vertices.
    max_dict = {"primal": 1, "dual": 0, "periodic": 0}
    range_max = dims - np.array([max_dict[typ] for typ in boundaries])
    ranges = [range(range_max[i]) for i in (0, 1, 2)]
    inds = it.product(*ranges)
    # TODO: Possibly change z-direction extent of dual complex for 'both'
    # option in SurfaceCode.
    # Primal vertices are combined into lists of six to be later usable
    # by the syndrome indentification in SurfaceCode.
    all_six_bodies = [
        [
            (2 * i, 2 * j + 1, 2 * k + 1),
            (2 * i + 1, 2 * j, 2 * k + 1),
            (2 * i + 1, 2 * j + 1, 2 * k),
            (2 * i + 2, 2 * j + 1, 2 * k + 1),
            (2 * i + 1, 2 * j + 2, 2 * k + 1),
            (2 * i + 1, 2 * j + 1, 2 * k + 2),
        ]
        for (i, j, k) in inds
    ]
    denested_six_bodies = set([a for b in all_six_bodies for a in b])

    # Tuple indices corresponding to dual and periodic boundaries.
    dual_inds = set((boundaries == "dual").nonzero()[0])
    periodic_inds = set((boundaries == "periodic").nonzero()[0])
    for vertex in denested_six_bodies:
        where_vertex_0, where_vertex_max = set(), set()
        # Ensure no vertices are included if they extend beyond
        # requested boundary conditions. Note this can also be achieved
        # by changing 'inds'.
        for i in range(3):
            if vertex[i] == 0:
                where_vertex_0 ^= {i}
            elif vertex[i] == 2 * dims[i]:
                where_vertex_max ^= {i}
        if not (
            where_vertex_max & dual_inds
            or where_vertex_0 & dual_inds
            or where_vertex_max & periodic_inds
        ):
            for neighbor in dual_neighbours(vertex):
                # Ensure no neighbours are included if they extend beyond
                # requested boundary conditions.
                where_neighbor_0, where_neighbor_max = set(), set()
                for i in range(3):
                    if neighbor[i] == 0:
                        where_neighbor_0 ^= {i}
                    elif neighbor[i] == 2 * dims[i]:
                        where_neighbor_max ^= {i}
                if not (
                    where_neighbor_max & dual_inds
                    or where_neighbor_0 & dual_inds
                    or where_neighbor_max & periodic_inds
                ):
                    edge = (vertex, neighbor)
                    weight = polarity(edge) if polarity else 1
                    G.add_node(vertex, type="primal")
                    G.add_node(neighbor, type="dual")
                    G.add_edge(vertex, neighbor, weight=weight)

                    # Additional edges for periodic boundaries.
                    for ind in where_neighbor_0 & periodic_inds:
                        max_coord = 2 * dims[ind] - 1
                        high_primal_vertex = list(neighbor)
                        high_primal_vertex[ind] = max_coord
                        neighbor_other_side = tuple(high_primal_vertex)
                        edge = (neighbor, neighbor_other_side)
                        weight = polarity(edge) if polarity else 1

                        G.add_node(neighbor_other_side, type="primal")
                        G.add_edge(neighbor, neighbor_other_side, weight=weight, periodic=True)

    # Store coordinates of primal cubes for later use.
    G.graph["primal_cubes"] = all_six_bodies
    return G


class SurfaceCode:
    """A class for representing the surface code.

    Represent the surface code in its measurement-based description. By
    specifying the distance, error complex, choice of boundaries, and polarity,
    store the graph state corresponding to the code, the set of stabilizer
    elements and the stabilizer graph, as well as the syndrome and boundary
    vertices.

    Attributes:
        distance (int): the code distance.
        dims (tup): a tuple of the spatial extent in x, y, z.
        ec (str): the error complex ('primal', 'dual', or 'both').
        boundaries (str): the boundary conditions. The options are:

            'open': ['primal', 'dual', 'dual'] for 'primal' or 'both' EC
                    ['primal', 'dual', 'primal'] for 'dual' EC
            'periodic': 'periodic' in all three directions.

        polarity (func): a function that specifies edge weights. It
            must be of the following form:

                polarity(edge) = weight.

            If not supplied, assumes all edges have weight 1.
        backend (string): The backend to use for the stabilizer graph.
            Can be "networkx" (the default) or "retworkx".
            The retworkx backend should be used when speed is a concern.

        graph (EGraph): the EGraph corresponding to the code, representing the
            graph state.
        '{b}'_stab_graph (StabilizerGraph): the stabilizer graph combining
            stabilizers from error complex b ('primal' or 'dual'). The
            particular implementation depends on backend.

        '{b}'_stabilizers (List[Stabilizer]): the stabilizer elements of the
            code according to the error complex b ('primal'/'dual').
        '{b}'_syndrome_coords (List[Tuple]): the coordinates of the syndrome
            vertices according to the error complex b ('primal'/'dual'/'all').
        '{b}'_syndrome_inds (List[Int]): the integer indices of the syndrome
            vertices according to the error comple b ('primal'/'dual'/'all').
        '{b}'_boundary_coords (list of tup): the coordinates of the boundary
            according to the error complex b ('primal'/'dual').
    """

    # TODO: Allow for codes with different aspect ratios.
    # TODO: Check distance convention for periodic boundaries.
    # TODO: Add x-y-but-not-z periodic boundaries.
    def __init__(
        self,
        distance,
        ec="primal",
        boundaries="open",
        polarity=None,
        backend="networkx",
    ):
        self.distance = distance
        self.dims = (distance, distance, distance)
        self.ec = ["primal", "dual"] if ec == "both" else [ec]

        if boundaries == "open":
            self.bound_str = "open_primal" if ec in ("primal", "both") else "open_dual"
        else:
            self.bound_str = boundaries
        self.boundaries = str_to_bound(self.bound_str)

        self.polarity = polarity

        self.graph = RHG_graph(self.dims, self.boundaries, polarity=polarity)
        self.graph.index_generator()
        # The following line defines the stabilizer, syndrome coordinate,
        # and syndrome index attributes.
        self.identify_stabilizers()
        # The following line defines the boundary points attribute.
        self.identify_boundary()

        if ec == "both":
            # For both error complexes, designate certain qubits as perfect
            # so that the correction check proceeds as expected. In particular
            # the qubits on the first and last temporal (z-direction) slice
            # are made perfect.
            perfect_qubits = self.graph.slice_coords("z", 1) + self.graph.slice_coords(
                "z", 2 * self.dims[2] - 1
            )
            self.graph.graph["perfect_points"] = perfect_qubits
            self.graph.graph["perfect_inds"] = [
                self.graph.to_indices[point] for point in perfect_qubits
            ]

        for ec in self.ec:
            if backend == "networkx":
                stabilizer_graph = NxStabilizerGraph(ec, self)
            elif backend == "retworkx":
                stabilizer_graph = RxStabilizerGraph(ec, self)
            else:
                raise ValueError("Invalid backend; options are 'networkx' and 'retworkx'.")
            setattr(self, ec + "_stab_graph", stabilizer_graph)

    def identify_stabilizers(self):
        """Set the stabilizer and syndrome coordinates of self.

        Generate a list of Stabilizer objects containing coordinates of
        all the stabilizer elements according to error complex ec.
        Furthermore, generate a list of all the relevant syndrome
        coordinates.

        In the end, the {ec}_syndrome_coords, {ec}_syndrome_inds,
        and {ec}_stabilizers attributes (where ec can be 'primal' or
        'dual') as well as all_syndrome_inds and all_syndrome_coords are set.
        """
        rhg_lattice = self.graph
        # Dimensions, boundary types, max and min ranges.
        dims = np.array(self.dims)

        all_six_bodies = {}

        if "primal" in self.ec:
            all_six_bodies["primal"] = self.graph.graph["primal_cubes"]
        if "dual" in self.ec:
            min_dict = {"primal": -1, "dual": 0, "periodic": -1}
            max_dict = {"primal": 1, "dual": 1, "periodic": 1}
            range_min = np.array([min_dict[typ] for typ in self.boundaries])
            range_max = dims - np.array([max_dict[typ] for typ in self.boundaries])
            ranges = [range(range_min[i], range_max[i]) for i in (0, 1, 2)]
            inds = it.product(*ranges)
            # All potential six-body stabilizers
            stabes = [
                [
                    (2 * i + 1, 2 * j + 2, 2 * k + 2),
                    (2 * i + 2, 2 * j + 1, 2 * k + 2),
                    (2 * i + 2, 2 * j + 2, 2 * k + 1),
                    (2 * i + 3, 2 * j + 2, 2 * k + 2),
                    (2 * i + 2, 2 * j + 3, 2 * k + 2),
                    (2 * i + 2, 2 * j + 2, 2 * k + 3),
                ]
                for (i, j, k) in inds
            ]
            all_six_bodies["dual"] = stabes

        periodic_inds = np.where(self.boundaries == "periodic")[0]

        for ec in self.ec:
            all_cubes = []
            syndrome_coords = []
            for stabe in all_six_bodies[ec]:
                actual_stabe = list(set(stabe).intersection(rhg_lattice))
                # Dealing with stabilizers at periodic boundaries
                if len(actual_stabe) < 6:
                    for ind in periodic_inds:
                        if ec == "dual":
                            highest_point = list(stabe[3 + ind])
                            if highest_point[ind] == 1:
                                highest_point[ind] = 2 * dims[ind] - 1
                                virtual_point = tuple(highest_point)
                                actual_stabe += [virtual_point]
                        else:
                            lowest_point = list(stabe[ind])
                            if lowest_point[ind] == 2 * dims[ind] - 2:
                                lowest_point[ind] = 0
                                virtual_point = tuple(lowest_point)
                                actual_stabe += [virtual_point]
                cube = Stabilizer(rhg_lattice.subgraph(actual_stabe))
                cube.physical = stabe
                syndrome_coords += actual_stabe
                all_cubes.append(cube)
            setattr(self, ec + "_syndrome_coords", list(set(syndrome_coords)))
            setattr(
                self,
                ec + "_syndrome_inds",
                [self.graph.to_indices[point] for point in syndrome_coords],
            )
            setattr(self, ec + "_stabilizers", all_cubes)
        for att in ["_syndrome_inds", "_syndrome_coords"]:
            new_attr = sum((getattr(self, ec + att) for ec in self.ec), start=[])
            setattr(self, "all" + att, new_attr)

    def identify_boundary(self):
        """Obtain coordinates of syndrome qubits on the boundary.

        The relevant boundaries are determined by the ec string. In the
        end, the attributes {b}_bound_points are set, where b can be
        'primal' or 'dual'.
        """
        for ec in self.ec:
            if self.bound_str == "periodic":
                setattr(self, ec + "_bound_points", [])
            else:
                dims = self.dims
                syndrome_coords = getattr(self, ec + "_syndrome_coords")
                bound_ind = np.where(self.boundaries == ec)[0][0]
                plane_dict = {0: "x", 1: "y", 2: "z"}

                low_index = 0 if ec == "primal" else 1
                low_bound_points = self.graph.slice_coords(plane_dict[bound_ind], low_index)
                final_low_set = set(low_bound_points).intersection(syndrome_coords)

                high_index = 2 * dims[bound_ind] - 2 if ec == "primal" else 2 * dims[bound_ind] - 1
                high_bound_points = self.graph.slice_coords(plane_dict[bound_ind], high_index)
                final_high_set = set(high_bound_points).intersection(syndrome_coords)

                setattr(self, ec + "_bound_points", list(final_low_set) + list(final_high_set))

    # TODO: tailored slice_coords function that constructs rather than iterates,
    # improving over EGraph method.

    def draw(self, **kwargs):
        """Draw the cluster state with matplotlib.

        See flamingpy.utils.viz.draw_EGraph for mor details. Use the
        default colour options: black for primal nodes, grey for dual
        nodes; blue for weight +1 edges, red for weight -1 edges.
        """
        edge_colors = {1: "b", -1: "r"} if self.polarity == alternating_polarity else "grey"
        default_opts = {
            "color_nodes": {"primal": "k", "dual": "grey"},
            "color_edges": edge_colors,
        }
        updated_opts = {**default_opts, **kwargs}
        return self.graph.draw(**updated_opts)
