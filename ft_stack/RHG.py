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
import itertools as it
import numpy as np

from matplotlib import pyplot as plt
from ft_stack.graphstates import EGraph, CVGraph
import networkx as nx


def alternating_polarity(edge):
    """Return +1 or -1 depending on the vertices that form edge edge.

    Help with the assignment of edge weights (i.e. polarity) for 
    the RHG graph. This particular alternating pattern ensures that 
    every vertex has two +1 and two -1 weight edges incident on it. The 
    precise combination depends on the direction of the edge and the 
    coordinates of the vertices. This pattern may be helpful to reduce 
    CV-level noise.
    
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
        list of tuple: the coordinates of the four neighbours.
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


def RHG_graph(
    dims,
    boundaries="finite",
    macronodes=False,
    polarity=None,
    node_color={"primal": "k", "dual": "grey"},
    edge_color={1: "b", -1: "r"},
    memory_saver=False,
):
    """Create an EGraph of a dims-dimensional RHG lattice.
    
    Generate an RHG lattice with dimensions given by dims (an integer
    denoting the number of stabilizer cubes in each direction,
    or a tuple specifying all three dimensions). By default, a useful
    set of finite boundary conditions is assumed, but any combination 
    can be specified.

    Args:
        dims (int or list): the dimensions of the lattice (the
            number of stabilizer cubes in each direction, complete or
            incomplete).
        boundaries (str or list-type, optional): the boundaries in each
            direction. We use primal/smooth and dual/rough
            interchangeably. If a string, 'primal', 'dual', 'periodic',
            assumes those boundaries in all three directions; or,
            accepts a list that specifies which boundary in which
            direction. Set  to 'finite' == ['primal', 'dual', 'dual']
            by default.
        macronodes (bool): if True, generates a macronode version of 
            the lattice, where each pair of vertices connected by
            an edge is replaced with a dumbbell .-., causing each
            vertex to be replaced by four (in bulk).
        polarity (func): a function that specifies edge weights. The
            input to the function should be an edge (i.e. list of two
            vertices) and the output should be the edge weight.
        node_color (dict): a dictionary of the form 
            {vertex_attribute: color}, where vertex_attribute can
            be 'primal' or 'dual'.
        edge_color (dict): a dictionary of the form {weight: color}
        memory_saver (bool): if True, omits string attributes in nodes 
            and edges.
            
    Returns:
        EGraph: the RHG lattice.
    """
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    G = EGraph(dims=dims, macronodes=macronodes)
    # Dealing with boundaries
    if boundaries == "finite":
        boundaries = ["primal", "dual", "dual"]
    elif type(boundaries) == str:
        boundaries = [boundaries] * 3
    # Locations of all primal vertices.
    inds = it.product(range(dims[0]), range(dims[1]), range(dims[2]))
    # Primal vertices are combined into lists of six to be later usable
    # by the syndrome indentification in the RHGCode.
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

    dual_inds = set(np.where(np.array(boundaries) == "dual")[0])
    periodic_inds = set(np.where(np.array(boundaries) == "periodic")[0])
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
            if macronodes:
                G.macro.add_node(vertex, micronodes=[])
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
                    color = edge_color.get(weight)
                    if macronodes:
                        central_vertex, central_neighbor = (
                            np.array(vertex),
                            np.array(neighbor),
                        )
                        direction_vec = 0.1 * (central_neighbor - central_vertex)
                        displaced_vertex = np.round(central_vertex + direction_vec, 1)
                        displaced_vertex = tuple(displaced_vertex)
                        displaced_neighbor = np.round(central_neighbor - direction_vec, 1)
                        displaced_neighbor = tuple(displaced_neighbor)
                        G.add_edge(
                            displaced_vertex, displaced_neighbor, weight=weight, color=color,
                        )
                        G.macro.add_node(neighbor, micronodes=[])
                    else:
                        if memory_saver:
                            G.add_edge(vertex, neighbor, weight=weight)
                        else:
                            G.add_node(vertex, type="primal", color=node_color["primal"])
                            G.add_node(neighbor, type="dual", color=node_color["dual"])
                            G.add_edge(vertex, neighbor, color=color, weight=weight)

                    # Additional edges for periodic boundaries.
                    for ind in where_neighbor_0 & periodic_inds:
                        max_coord = 2 * dims[ind] - 1
                        high_primal_vertex = list(neighbor)
                        high_primal_vertex[ind] = max_coord
                        neighbor_other_side = tuple(high_primal_vertex)
                        edge = (neighbor, neighbor_other_side)
                        weight = polarity(edge) if polarity else 1
                        color = edge_color.get(weight)
                        if macronodes:
                            central_neighbor, other_side = (
                                np.array(neighbor),
                                high_primal_vertex,
                            )
                            direction_vec = other_side - central_neighbor
                            shortened_vec = -0.1 * direction_vec / np.linalg.norm(direction_vec)
                            displaced_neighbor = np.round(central_neighbor + shortened_vec, 1)
                            displaced_neighbor = tuple(displaced_neighbor)
                            displaced_other_side = np.round(other_side - shortened_vec, 1)
                            displaced_other_side = tuple(displaced_other_side)
                            G.add_edge(
                                displaced_neighbor,
                                displaced_other_side,
                                weight=weight,
                                color=color,
                            )
                        else:
                            if memory_saver:
                                G.add_edge(neighbor, neighbor_other_side, weight=weight)
                            else:
                                G.add_node(
                                    neighbor_other_side, type="primal", color=node_color["primal"],
                                )
                                G.add_edge(
                                    neighbor, neighbor_other_side, color=color, weight=weight,
                                )

    G.graph["primal_cubes"] = all_six_bodies
    return G


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
        decoding_graph (Graph): decoding graph constructed from the stabilizers
            without edge weights populated
        _decoder_mapping (dict): mapping between nodes and indices for the
            decoding graph
    """

    def __init__(
        self,
        distance,
        error_complex="primal",
        boundaries="finite",
        polarity=None,
        decoding_graph=True,
    ):
        """Initialize the RHG code."""
        # TODO: Check code distance convention.
        self.distance = distance
        self.dims = (distance, distance, distance)
        self.complex = error_complex
        if boundaries == "finite":
            self.boundaries = ["primal", "dual", "dual"]
        elif type(boundaries) == str:
            self.boundaries = [boundaries] * 3
        else:
            self.boundaries = boundaries
        self._polarity = polarity

        self.graph = RHG_graph(self.dims, boundaries=self.boundaries, polarity=polarity)
        self.graph.index_generator()
        # The following line also defines the self.syndrome_coords
        # attribute.
        self.stabilizers = self.identify_stabilizers(self.complex)
        self.syndrome_inds = [self.graph.to_indices[point] for point in self.syndrome_coords]
        self.boundary_coords = self.identify_boundary(self.complex)

        if decoding_graph:
            self.decoding_graph, self._decoder_mapping = self.construct_decoding_graph()

    def identify_stabilizers(self, error_complex="primal"):
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
        min_dict = {"primal": 0, "dual": 1, "periodic": 0}
        mins = [min_dict[typ] for typ in boundaries]
        maxes = np.array([2 * dims[i] - 1 for i in (0, 1, 2)])

        # Function for generating ranges from lists of mins and maxes.
        ranges = [range(dims[i]) for i in range(3)]
        inds = it.product(*ranges)

        # TODO: Implement dual error complex.
        if error_complex == "primal":
            # All potential six-body stabilizers
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

        all_cubes = []
        syndrome_coords = []

        periodic_inds = np.where(boundaries == "periodic")[0]
        dual_inds = np.where(boundaries == "dual")[0]
        for stabe in all_six_bodies:
            actual_stabe = list(set(stabe) & set(G))
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

    def identify_boundary(self, error_complex="primal"):
        """Obtain coordinates of syndrome qubits on the boundary.

        The relevant boundary is determined by the error_complex string.
        """
        # TODO: Dual boundaries.
        dims = self.dims
        boundaries = np.array(self.boundaries)
        if error_complex == "primal":
            # Odd indices, which is where primal syndrome qubits are located.
            odds = [
                range(1, 2 * dims[0], 2),
                range(1, 2 * dims[1], 2),
                range(1, 2 * dims[2], 2),
            ]
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

    def construct_decoding_graph(self):
        """Create a decoding graph from the RHG lattice G.

        The decoding graph has as its nodes every stabilizer in G and a
        every boundary point (for now coming uniquely from a primal
        boundary). The output graph has stabilizer nodes relabelled to 
        integer indices, but still points to the original stabilizer with 
        the help of the 'stabilizer' attribute. Common vertices are stored 
        under the 'common_vertex' edge attribute.
        """
        # An empty decoding graph.
        G_dec = nx.Graph(title="Decoding Graph")
        cubes = self.stabilizers
        G_dec.add_nodes_from(cubes)
        # For stabilizer cubes sharing a vertex, define an edge between
        # them with weight equal to the weight assigned to the vertex.
        for (cube1, cube2) in it.combinations(G_dec, 2):
            common_vertex = set(cube1.coords()) & set(cube2.coords())
            if common_vertex:
                coordinate = common_vertex.pop()
                G_dec.add_edge(cube1, cube2, common_vertex=coordinate, weight=None)

        # Include primal boundary vertices in the case of non-periodic
        # boundary conditions. (The following list will be empty absent a
        # primal boundary).
        bound_points = self.boundary_coords

        # For boundary points sharing a vertex with a stabilizer cube,
        # add an edge between them with weight equal to the weight assigned
        # to the vertex.
        for (cube, point) in it.product(G_dec, bound_points):
            if point in cube.coords():
                G_dec.add_edge(cube, point, common_vertex=point, weight=None)

        # Relabel the nodes of the decoding graph to integers and define
        # the mapping between nodes and indices.
        G_relabelled = nx.convert_node_labels_to_integers(G_dec, label_attribute="stabilizer")
        mapping = dict(zip(G_dec.nodes(), range(0, G_dec.order())))

        # For non-periodic boundary conditions, include boundary vertices.
        # Get the indices of the boundary vertices from the decoding
        # graph.
        bound_inds = [mapping[point] for point in bound_points]
        G_relabelled.graph["boundary_points"] = bound_inds[:]
        bound_points = G_relabelled.graph["boundary_points"]
        # Connect two helper points, one per primal boundary, with weight 0
        # to all the boundary points. This will speed up Dijkstra.
        n_b = len(bound_points) // 2
        low_bound_points, high_bound_points = bound_points[:n_b], bound_points[n_b:]
        v_low = "low"
        v_high = "high"
        for point in low_bound_points:
            G_relabelled.add_edge(v_low, point, weight=0)
        for point in high_bound_points:
            G_relabelled.add_edge(v_high, point, weight=0)
        
        # Identify the real points so "high" and "low" can be skipped in some of 
        # the decoding functions down the line
        G_relabelled.graph["real_points"] = set(G_relabelled.nodes)-set(('high','low'))
        
        return G_relabelled, mapping


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
        bit_vals = [G.nodes[node]["bit_val"] for node in G]
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
        return (
            np.average(self.xlims()),
            np.average(self.ylims()),
            np.average(self.zlims()),
        )


if __name__ == "__main__":
    # Simple tests
    # Instantiate an RHG latice of a certian distance, with certain
    # boundaries. Draw the EGraph.
    d = 2
    # boundaries = "finite"
    # boundaries = "primal"
    boundaries = "periodic"
    # boundaries = ["primal", "dual", "periodic"]
    # For iterating through boundaries
    # for boundaries in it.product(['primal', 'dual', 'periodic'], repeat=3):

    RHG = RHGCode(d, boundaries=boundaries, polarity=alternating_polarity)
    RHG_lattice = RHG.graph
    # Check maronode lattice
    # RHG_lattice = RHG_graph(d, boundaries=boundaries, macronodes=True)
    ax = RHG_lattice.draw(color_nodes="MidnightBlue", color_edges=False, label="index")

    # # Check edges between boundaries for periodic boundary conditions.
    # all_boundaries = []
    # for plane in ("x", "y", "z"):
    #     for i in (0, 2 * d - 1):
    #         all_boundaries += RHG.graph.slice_coords(plane, i)
    #     # For macronode lattice, change to (-0.1, 2 * d - 0.9)
    # RHG_subgraph = RHG_lattice.subgraph(all_boundaries)
    # RHG_subgraph.draw()

    # # Check stabilizer coordinates
    # syndrome = RHG.stabilizers
    # print("Number of six-body stabilizers :", len(syndrome))
    # for i in range(len(syndrome)):
    #     cube = syndrome[i]
    #     color = np.random.rand(3)
    #     for point in cube.egraph:
    #         x, z, y = point
    #         ax.scatter(x, z, y, color=color, s=200)
    # ax.set_title(str(boundaries).capitalize() + " boundaries")

    # # Check sampling
    # delta = 0.001
    # # Percent p-squeezed states.
    # p_swap = 0
    # CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)
    # for sampling_order in ["initial", "final", "two-step"]:
    #     model = {"noise": "grn", "delta": delta, "sampling_order": sampling_order}
    #     CVRHG.apply_noise(model)
    #     CVRHG.measure_hom("p")
    #     outcomes = CVRHG.hom_outcomes()
    #     plt.figure(figsize=(16, 9))
    #     plt.hist(outcomes, bins=100)
    #     CVRHG.draw(label="hom_val_p")
    #     CVRHG.draw(label="hom_val_p")
