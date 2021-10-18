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
"""Decoders for measurement-based codes."""
import sys
import itertools as it
import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ft_stack.graphstates import CVGraph
from ft_stack.GKP import GKP_binner, Z_err_cond
from ft_stack.RHG import alternating_polarity, RHGCube, RHGCode

# Smallest and largest numbers representable.
smallest_number = sys.float_info.min
largest_number = sys.float_info.max


def graph_drawer(G, label_edges=True):
    """Draw decoding and matching graphs G with a color legend."""
    title = G.graph["title"]
    plt.figure()
    plt.title(title, family="serif", size=10)
    # NetworkX drawing function for circular embedding of graphs.
    nx.draw_circular(
        G,
        edgelist=[],
        with_labels=True,
        node_color="k",
        font_size=7,
        font_color="w",
        font_family="serif",
    )
    # Color edges based on weight, and draw a colobar.
    weight_list = [G.edges[edge]["weight"] for edge in G.edges]
    weight_dict = {edge: "{:.2f}".format(G.edges[edge]["weight"]) for edge in G.edges}
    if label_edges:
        nx.draw_networkx_edge_labels(G, nx.circular_layout(G), edge_labels=weight_dict, font_size=7)
    r = nx.draw_networkx_edges(G, nx.circular_layout(G), edge_color=weight_list)
    plt.colorbar(r)


def syndrome_plot(code, G_dec, index_dict=None, drawing_opts=None):
    """Draw the syndrome plot for a CVGraph G.

    A comprehensive graphing tool for drawing the error syndrome of
    a CVGraph G. Labelling options are specified with the help of
    drawing_opts, and can include:

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per CVGraph.draw
                'label_cubes' -> indices of the stabilizers
                'label_boundary'-> indices of the boundary points
                'legend' -> legends for the nodes and cubes

    Cubes are shown as transparent voxels, green for even parity and
    red for odd. For now, cubes on periodic boundaries are not shown,
    and stabilizers on dual boundaries occupy are shown to occupy
    the full space of a six-body cube.

    Args:
        G_dec (networkx.Graph): the decoding graph
        G (CVGraph): the encoded CVGraph
        index_dict (dict): the stabiizer-to-index mapping
        drawing_opts (dict): a dictionary of drawing options, with
            all possibilities described above.

    Returns:
        matplotlib.pyplot.axes: the 'axes' object
    """
    # Font properties
    # TODO: Make consistent with EGraph fontprops
    font_size = 10 * sum(code.dims) ** (1 / 2)
    # Set plotting options
    plot_params = {
        "font.size": font_size,
        "font.family": "serif",
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "grid.color": "lightgray",
        "lines.markersize": font_size,
    }
    plt.rcParams.update(plot_params)

    bound_inds = G_dec.graph["boundary_points"]
    points = [G_dec.nodes[i]["stabilizer"] for i in bound_inds]

    cubes = code.stabilizers
    # Default drawing options.
    draw_dict = {
        "show_nodes": False,
        "color_nodes": "state",
        "label": None,
        "legend": True,
        "title": True,
        "state_colors": {"p": None, "GKP": None},
        "display_axes": True,
        "label_cubes": True,
        "label_boundary": False,
    }
    # Combine default dictionary with supplied dictionary, duplicates
    # favor supplied dictionary.
    if drawing_opts is None:
        drawing_opts = {}
    drawing_opts = {**draw_dict, **drawing_opts}

    # Shape and font properties from the original graph.
    shape = np.array(code.dims)
    # font_props = state.font_props
    # If show_nodes is True, get the axes object and legend from
    # CVGraph.sketch (this also plots the graph in the console).
    if drawing_opts["show_nodes"]:
        # TODO: If draw method moved out of CVGraph and into EGraph,
        # the state argument would be unnecessary here.
        egraph_args = [
            "color_nodes",
            "label",
            "legend",
            "title",
            "state_colors",
            "display_axes",
        ]
        egraph_opts = {k: drawing_opts[k] for k in egraph_args}
        ax = code.graph.draw(**egraph_opts)
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        # TODO: Initialize axes based on empty ax object from state.draw()
        # but prevent from state.draw() from plotting.
        fig = plt.figure(figsize=(2 * (np.sum(shape) + 2), 2 * (np.sum(shape) + 2)))
        ax = fig.gca(projection="3d")
        # ax.tick_params(labelsize=font_props['size'])
        plt.xticks(range(0, 2 * shape[0] + 1))
        plt.yticks(range(0, 2 * shape[1] + 1))
        ax.set_zticks(range(0, 2 * shape[2] + 1))
        ax.set_xlabel(
            "x",
            # fontdict=font_props,
            labelpad=20,
        )
        ax.set_ylabel(
            "z",
            # fontdict=font_props,
            labelpad=20,
        )
        ax.set_zlabel(
            "y",
            # fontdict=font_props,
            labelpad=20,
        )
        plt.rcParams["grid.color"] = "lightgray"
        leg = None
    # Illustrate stabilizers with voxels colored green for even
    # parity and red for odd pariy.
    filled = np.zeros(shape, dtype=object)
    for cube in cubes:

        # TODO: Deal appropriately with cubes on periodic and dual
        # boundaries.

        # Obtain smallest, largest, and middle coordinates for each
        # cube. Divided by 2 becaues voxels are 1X1X1.
        xmin, xmax = np.array(cube.xlims(), dtype=int) // 2
        ymin, ymax = np.array(cube.ylims(), dtype=int) // 2
        zmin, zmax = np.array(cube.zlims(), dtype=int) // 2
        xmid, ymid, zmid = np.array(cube.midpoint())
        # Fill in the color arrays depending on parity.
        if cube.parity:
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = "#FF000015"
        else:
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = "#00FF0015"
        if drawing_opts["label_cubes"] and index_dict:
            if cube in index_dict:
                ax.text(
                    xmid,
                    ymid,
                    zmid,
                    index_dict[cube],
                    # fontdict=font_props
                )

    # This portion adapted from a Matplotlib official example to fix
    # an issue with filling in the insides of voxels: the code
    # expands the indices and creates small gaps between the voxels.

    def explode(data):
        size = np.array(data.shape) * 2
        data_e = np.zeros(size - 1, dtype=data.dtype)
        data_e[::2, ::2, ::2] = data
        return data_e

    # upscale the above voxel image, leaving gaps
    filled_e = explode(filled)
    # Shrink the gaps
    x, y, z = np.indices(np.array(filled_e.shape) + 1, dtype=float)
    x[0::2, :, :] += 0.05
    y[:, 0::2, :] += 0.05
    z[:, :, 0::2] += 0.05
    x[1::2, :, :] += 0.95
    y[:, 1::2, :] += 0.95
    z[:, :, 1::2] += 0.95
    ax.voxels(x, y, z, filled_e, facecolors=filled_e)

    if drawing_opts["label_boundary"]:
        for point in points:
            ax.scatter(point[0], point[1], point[2], s=70, c="k")
            ax.text(
                point[0],
                point[1],
                point[2],
                index_dict[point],
                # fontdict=font_props
            )

    # Define a legend for red/green cubes.
    legend_elements = [
        Patch(facecolor="#00FF0050", label="even parity"),
        Patch(facecolor="#FF000050", label="odd parity"),
    ]
    if drawing_opts["legend"]:
        ax.legend(
            handles=legend_elements,
            # prop=font_props,
            loc="upper left",
        )
    # Since CVGraph.sketch() legend has been overwritten, readd
    # it to the plot.
    if leg:
        ax.add_artist(leg)
    return ax


def assign_weights(code, **kwargs):
    """Assign weights to qubits in a hybrid CV graph state CVG.

    Args:
        code (code class): the qubit code
        state (CVGraph): the CVGraph whose syndrome has been measured
        method (str, optional): the method for weight assignment. By
            default, 'unit', denoting weight 1 everyoewhere. For
            heuristic and analog weight assignment from blueprint, use
            'blueprint'
        integer (bool, optional): whether to convert weights to
            integers using Python's round function; False by default
        multiplier (int, optional): multiply the weight by multiplier
            before rounding; 1 by default.

    Returns:
        None
    """
    default_options = {"method": "unit", "integer": False, "multiplier": 1}
    weight_options = {**default_options, **kwargs}
    G = code.graph
    # Get and denest the syndrome coordinates.
    syndrome_coords = code.syndrome_coords
    # Blueprint weight assignment dependent on type of neighbour.
    if weight_options.get("method") == "blueprint":
        for node in syndrome_coords:
            neighbors = G[node]
            # List and number of p-squeezed states in neighborhood of node.
            p_list = [G.nodes[v]["state"] for v in neighbors if G.nodes[v]["state"] == "p"]
            p_count = len(p_list)
            if p_count in (0, 1):
                if weight_options.get("prob_precomputed"):
                    err_prob = G.nodes[node]["p_phase_cond"]
                else:
                    delta_effective = (len(neighbors) + 1) * weight_options.get("delta")
                    hom_val = G.nodes[node]["hom_val_p"]
                    err_prob = Z_err_cond(delta_effective, hom_val)
                # Allow for taking log of 0.
                # TODO: Is this the best way to do it? Or can I just choose
                # an arbitrary small number?
                if err_prob > 0.5:
                    err_prob = 0.5
                if err_prob == 0:
                    err_prob = smallest_number
                if weight_options.get("integer"):
                    multiplier = weight_options.get("multiplier")
                    weight = round(-multiplier * np.log(err_prob))
                else:
                    weight = -np.log(err_prob)
                G.nodes[node]["weight"] = weight
            else:
                # Dictionary of the form number of swapouts: error probability.
                weight_dict = {2: 1 / 4, 3: 1 / 3, 4: 2 / 5}
                err_prob = weight_dict[p_count]
                if weight_options.get("integer"):
                    multiplier = weight_options.get("multiplier")
                    weight = round(-multiplier * np.log(err_prob))
                else:
                    weight = -np.log(err_prob)
                G.nodes[node]["weight"] = weight
        return
    # Naive weight assignment, unity weights.
    if weight_options.get("method") == "unit":
        for node in syndrome_coords:
            G.nodes[node]["weight"] = 1


# TODO: General functions for applying noise and measuring syndrome.


def CV_decoder(code, translator=GKP_binner):
    """Convert homodyne outcomes to bit values according to translate.

    The inner (CV) decoder, aka translator, aka binning function. Set
    converted values to the bit_val attribute for nodes in G.

    Args:
        state: the CVGraph with homodyne outcomes computed.
        translator: the choice of binning function; by default, the
            standard GKP binning function that snaps to the closest
            integer multiple of sqrt(pi).
    Returns:
        None
    """
    syndrome = code.syndrome_coords
    # TODO: Generalize to nonlocal translators
    # TODO: Vectorize?
    for point in syndrome:
        hom_val = code.graph.nodes[point]["hom_val_p"]
        bit_val = translator([hom_val])[0]
        code.graph.nodes[point]["bit_val"] = bit_val


def decoding_graph(code, draw=False, drawing_opts=None, label_edges=False):
    """Populate the edge weights of the decoding graph from the RHG lattice G,
    and determine the parity of the stabilizer cubes.

    The decoding graph has as its nodes every stabilizer in G and a
    every boundary point (for now coming uniquely from a primal
    boundary). Two stabilizers (or a stabilizer and a boundary point)
    sharing a vertex are connected by an edge whose weight is equal to
    the weight assigned to that vertex in G. The output graph has
    stabilizer nodes relabelled to integer indices, but still points
    to the original stabilizer with the help of the 'stabilizer'
    attribute. The output graph furthermore stores the indices of
    odd-parity cubes (under an 'odd_cubes' graph attribute) and of
    the boundary points (under 'boundary_points'). Common vertices
    are stored under the 'common_vertex' edge attribute. Note that one
    ought first compute the phase error probabilities, conduct a
    homodyne measurement, and translate the outcomes on G.

    Args:
        G (CVGraph): the CVGraph to decode
        draw (bool): if True, draw the decoding graph and syndrome plot
        drawing_opts (dict): the drawing options, as in syndrome_plot:

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per CVGraph.draw
                'label_cubes' -> indices of the stabilizers
                'label_boundary'-> indices of the boundary points
                'legend' -> legends for the nodes and cubes

    Returns:
        networkx.Graph: the decoding graph.
    """
    cubes = code.stabilizers

    decoding_graph = code.decoding_graph
    mapping = code._decoder_mapping

    # Assign edge weights based on the common vertices between stabilizers
    for edge in decoding_graph.edges:
        if "common_vertex" in decoding_graph.edges[edge]:
            common = decoding_graph.edges[edge]["common_vertex"]
            decoding_graph.edges[edge]["weight"] = code.graph.nodes[common]["weight"]

    # Indices of odd parity cubes and boundary vertices; add these
    # index lists to the graph attributes dictionary, to be used by the
    # matching graph.
    odd_parity_cubes = [cube for cube in cubes if cube.parity]
    odd_parity_inds = [mapping[cube] for cube in odd_parity_cubes]
    decoding_graph.graph["odd_cubes"] = odd_parity_inds[:]

    # Draw the lattice and the abstract decoding graph.
    if draw:
        graph_drawer(decoding_graph, label_edges=label_edges)
        syndrome_plot(code, decoding_graph, index_dict=mapping, drawing_opts=drawing_opts)
    return decoding_graph


def matching_graph(G_dec, alg="dijkstra", draw=False, label_edges=False):
    """Create a matching graph from the decoding graph G.

    Generate a matching graph from the decoding graph G according to
    algorithm alg. By default, this is the NetworkX Dijkstra shortest-
    path algorithm. This graph will be fed into a subsequent minimum-
    weight-perfect-matching algorithm. The matching graph has as half
    of its nodes the odd-parity stabilizers. The edge connecting two
    nodes corresponds to the weight of the minimum-weight-path between
    the nodes in the decoding graph. Additionally, each unsatisfied
    stabilizer is connected to a unique boundary point (for now from
    a primal bundary) located at the shortest weighted distance from
    the stabilizer. Between each other, the boundary points are
    connected by an edge of weight 0. The output graph stores the
    indices of the used boundary points under the 'used_boundary_point'
    attribute. Paths are stored under the 'paths' attribute of edges,
    and 'inverse_weights' are also stored, for the benefit of maximum-
    weight-matching algorithms

    Args:
        G (networkx.Graph): the decoding graph, storing information
            about indices of odd-parity-cubes (under 'odd_cubes' graph
            attribute) and boundary points (under 'boundary_points').
        alg (str): the algorithm for shortest-path finding. By default,
            uses variations of Dijkstra functions from NetworkX
        draw (bool): if True, draws the matching graph

    Returns:
        networkx.Graph: the matching graph.
    """
    # An empty matching graph.
    G_match = nx.Graph(title="Matching Graph")

    # Get the indices of the odd parity cubes from teh decoding graph.
    odd_parity_inds = G_dec.graph["odd_cubes"]

    # Give shorter names to the Dijkstra shortest path algorithms.
    if alg == "dijkstra":
        alg = sp.single_source_dijkstra

    # Combinations of odd-parity cubes.
    odd_ind_dict = {i: [] for i in odd_parity_inds[:-1]}
    odd_combs = it.combinations(odd_parity_inds, 2)
    for pair in odd_combs:
        odd_ind_dict[pair[0]] += [pair[1]]
    # Find the shortest paths between odd-parity cubes.
    for cube1 in odd_parity_inds[:-1]:
        lengths, paths = alg(G_dec, cube1)
        for cube2 in odd_ind_dict[cube1]:
            length = lengths[cube2]
            path = paths[cube2]
            # Add edge to the matching graph between the cubes, with weight
            # equal to the length of the shortest path.
            # TODO: Is the behavior correct for negative weights, or do I
            # want 1/weight or max_num - weight?
            G_match.add_edge(cube1, cube2, weight=length, inverse_weight=-length, path=path)

    virtual_points = []
    if G_dec.graph["boundary_points"]:

        i = 0
        low_lengths, low_paths = alg(G_dec, "low")
        high_lengths, high_paths = alg(G_dec, "high")
        for cube in odd_parity_inds:
            distances = (low_lengths[cube], high_lengths[cube])
            where_shortest = np.argmin(distances)
            if where_shortest == 0:
                length = low_lengths[cube]
                full_path = low_paths[cube]
            if where_shortest == 1:
                length = high_lengths[cube]
                full_path = high_paths[cube]
            point = full_path[1]
            virtual_point = (point, i)
            path = full_path[1:]
            # Add edge to the matching graph between the cube and
            # the virtual excitation corresponding to the boundary
            # vertex, with weight equal to the length of the shortest
            # path.
            G_match.add_edge(cube, virtual_point, weight=length, inverse_weight=-length, path=path)
            i += 1
            virtual_points += [virtual_point]
        # Add edge with weight 0 between any two virtual excitations.
        for (point1, point2) in it.combinations(virtual_points, 2):
            G_match.add_edge(point1, point2, weight=0, inverse_weight=0)

    G_match.graph["virtual_points"] = virtual_points

    if draw:
        graph_drawer(G_match, label_edges=label_edges)

    return G_match


def MWPM(G_match, G_dec, alg="blossom_nx", draw=False, label_edges=False):
    """Run minimum-weight-perfect matching on matching graph G_match.

    Run a minimum-weight-perfect-matching (MWPM) algorithm (the
    BlossomV/Edmunds aglorithm as implemented in NetworkX by default)
    on the matching graph G_match. Under the hood, runs maximum weight
    matching with maximum cardinality on inverse weights. The matching
    combines pairs of nodes in the matching graph in a way that
    minimizes the total weight. Perfect matching combine all pairs of
    nodes.

    Args:
        G_match (networkx.Graph): the matching graph, storing inverse
            weights under the 'inverse_weight' edge attribute ad paths
            under 'path'
        G_dec (networkx.Graph): the decoding graph, pointing to the
            stabilizer cubes under the 'stabilizer' node attribute
        alg (string): the matching algorithm; by default, NetworkX
            implementation of BlossomV/Edmunds maximum-weight matching
        draw (bool): if True, visualize the matching on top of the
            syndrome plot (for now requires the syndrome plot from the
            decoding graph to be drawn immediately prior)

    Return:
        set of tuples: pairs of all matched nodes.
    """
    if alg == "blossom_nx":
        alg = nx.max_weight_matching
    matching = alg(G_match, maxcardinality=True, weight="inverse_weight")
    # TODO: Drop the requirement of the syndrome plot from having
    # to be plotted immediately prior to the matching.
    if draw:
        virtual_points = G_match.graph["virtual_points"]
        for pair in matching:
            if pair not in it.product(virtual_points, virtual_points):
                xlist, ylist, zlist = [], [], []
                path = G_match.edges[pair]["path"]
                for node in path:
                    stabe = G_dec.nodes[node]["stabilizer"]
                    if isinstance(stabe, RHGCube):
                        x, y, z = stabe.midpoint()
                    else:
                        x, y, z = stabe
                        plt.plot(x, y, z, marker="2", ms=50, c="k")
                    xlist += [x]
                    ylist += [y]
                    zlist += [z]
                plt.title("Minimum-weight perfect matching", family="serif", size=20)
                plt.plot(xlist, ylist, zlist, "o-k", ms=20, linewidth=5, c=np.random.rand(3))
        graph_drawer(G_match, label_edges=label_edges)
    return matching


def recovery(code, G_match, G_dec, matching, sanity_check=False):
    """Run the recovery operation on graph G.

    Fip the bit values of all the vertices in the path connecting each
    pair of stabilizers according to the matching. If check, verify
    that there are no odd-parity cubes remaining, or display their
    indices of there are.

    Args:
        G_match (networkx.Graph): the matching graph
        G_dec (networkx.Graph): the decoding graph
        G (CVGraph): the CVGraph to correct
        matching (set of tuples): the minimum-weight perfect matching
        check (bool): if True, check if the recovery has succeeded.

    Returns:
        None or bool: if check, False if the recovery failed, True if
            it succeeded. Otherwise None.
    """
    virtual_points = G_match.graph["virtual_points"]
    for pair in matching:
        if pair not in it.product(virtual_points, virtual_points):
            path = G_match.edges[pair]["path"]
            pairs = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for pair in pairs:
                if "common_vertex" in G_dec.edges[pair]:
                    common_vertex = G_dec.edges[pair]["common_vertex"]
                    code.graph.nodes[common_vertex]["bit_val"] ^= 1

    if sanity_check:
        G_dec_new = decoding_graph(code, draw=False)
        odd_cubes = G_dec_new.graph["odd_cubes"]
        if odd_cubes:
            print("Unsatisfied stabilizers:", odd_cubes)
        print("Recovery succeeded - no unsatisfied stabilizers.")


# TODO: Rename to logical_error_check or someting like that. Clarify
# what correlation surface is, per RHG paper.
def check_correction(code, plane=None, sheet=0, sanity_check=False):
    """Perform a correlation-surface check.

    Check the total parity of a correlation surface specified by
    direction plane and index sheet to check if the error correction
    procedure has succeeded. By default, checks the x = 0 plane.

    Args:
        G (CVGraph): the recovered graph
        plane (str): 'x', 'y', or 'z', determining the direction of the
            correlation surface; if None, select plane to check
            depending on the boundary conditions
        sheet (int): the sheet index (from 0 to one less than the
            largest coordinate in the direction plane).
        sanity_check (bool): if True, display the total parity of
            all correlation surfaces in the direction plane to verify
            if parity is conserved.

    Returns:
        None
    """
    dims = np.array(code.dims)
    dir_dict = {"x": 0, "y": 1, "z": 2}
    truth_dict = {"x": [], "y": [], "z": []}
    boundaries = code.boundaries

    planes_to_check = []
    if plane:
        planes_to_check += [plane]
    elif boundaries == ["periodic"] * 3:
        planes_to_check = ["x", "y", "z"]
    elif tuple(boundaries) in set(it.permutations(["primal", "dual", "dual"])):
        where_primal = np.where(np.array(boundaries) == "primal")[0][0]
        planes_to_check = [["x", "y", "z"][where_primal]]

    minimum, maximum = sheet, sheet + 2
    for plane in planes_to_check:
        if sanity_check:
            minimum, maximum = 0, 2 * dims[dir_dict[plane]]
        for sheet in range(minimum, maximum, 2):
            slice_verts = code.graph.slice_coords(plane, sheet)
            syndrome_verts = code.syndrome_coords
            only_primal = set(slice_verts) & set(syndrome_verts)
            parity = 0
            for node in only_primal:
                parity ^= code.graph.nodes[node]["bit_val"]
            truth_dict[plane].append(bool(1 - parity))

    if sanity_check:
        print(truth_dict)

    all_surfaces = np.array([truth_dict[i][0] for i in planes_to_check])
    return np.all(all_surfaces)


def correct(
    code, decoder, weight_options=None, draw=False, drawing_opts=None, sanity_check=False,
):
    """Run through all the error-correction steps.

    Combines weight assignment, decoding and matching graph creation,
    minimum-weight-perfect matching, recovery, and correctness check.

    Args:
        code (code): the code class to decode and correct
        inner_decoder (str): the CV decoder; GKP_binner by default
        outer_decoder (str): the DV decoder; MWPM by default
        weight_options (dict): how to assign weights; options are

                'method': 'unit' or 'blueprint'
                'integer': True (for rounding) or False (for not)
                'multiplier': integer denoting multiplicative factor
                    before rounding

        draw (bool): if True, draw the decoding graph, matching graph,
            syndrome plot, and minimum-weight matching
        drawing_opts (dict): the drawing options, as in syndrome_plot
            and decoding_graph

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per CVGraph.draw
                'label_cubes' -> indices of the stabilizers
                'label_boundary'-> indices of the boundary points
                'legend' -> legends for the nodes and cubes

        sanity_check (bool): if True, check that the recovery
            operation succeeded and verify that parity is conserved
            among all correlation surfaces.
    Returns:
        bool: True if error correction succeeded, False if not.
    """
    inner_dict = {"basic": GKP_binner}
    outer_dict = {"MWPM": "MWPM"}

    inner_decoder = decoder.get("inner")
    outer_decoder = decoder.get("outer")
    if inner_decoder:
        CV_decoder(code, translator=inner_dict[inner_decoder])
    if outer_dict[outer_decoder] == "MWPM":
        if drawing_opts is None:
            drawing_opts = {}
        label_edges = drawing_opts.get("label_edges")
        if weight_options is None:
            weight_options = {}
        assign_weights(code, **weight_options)
        G_dec = decoding_graph(code, draw=draw, drawing_opts=drawing_opts, label_edges=label_edges)
        G_match = matching_graph(G_dec)
        matching = MWPM(G_match, G_dec, draw=draw, label_edges=label_edges)
        recovery(code, G_match, G_dec, matching, sanity_check=sanity_check)
    result = check_correction(code, sanity_check=sanity_check)
    return result


if __name__ == "__main__":
    # DV (outer) code
    distance = 3
    boundaries = "periodic"
    RHG_code = RHGCode(distance=distance, boundaries=boundaries, polarity=alternating_polarity)
    RHG_lattice = RHG_code.graph
    # CV (inner) code/state
    p_swap = 0
    CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)

    # Noise model
    delta = 0.1
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

    # Decoding options
    decoder = {"inner": "basic", "outer": "MWPM"}
    weight_options = {
        "method": "blueprint",
        "integer": True,
        "multiplier": 100,
        "delta": delta,
    }

    # Drawing options
    dw = {
        "show_nodes": True,
        "color_nodes": "state",
        "label": "bit_val",
        "legend": True,
        "title": True,
        "display_axes": True,
        "label_edges": True,
        "label_cubes": True,
        "label_boundary": False,
    }

    trials = 1
    success = 0
    for trial in range(trials):
        # Apply noise
        CVRHG.apply_noise(cv_noise)
        # Measure syndrome
        CVRHG.measure_hom("p", RHG_code.syndrome_inds)
        c = correct(
            code=RHG_code,
            decoder=decoder,
            weight_options=weight_options,
            draw=True,
            drawing_opts=dw,
            sanity_check=True,
        )
        success += c
    error = (trials - success) / trials
    print(error)
