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
import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import itertools as it
import matplotlib.pyplot as plt

from graphstates import CVGraph
from GKP import basic_translate
import RHG

# Smallest and largest numbers representable.
smallest_number = sys.float_info.min
largest_number = sys.float_info.max


def graph_drawer(G, label_edges=True):
    """Draw decoding and matching graphs with a color legend."""
    title = G.graph['title']
    plt.figure()
    plt.title(title, family='serif', size=10)
    # NetworkX drawing function for circular embedding of graphs.
    nx.draw_circular(G,
                     edgelist=[],
                     with_labels=True,
                     node_color='k',
                     font_size=7,
                     font_color='w',
                     font_family='serif')
    # Color edges based on weight, and draw a colobar.
    weight_list = [G.edges[edge]['weight'] for edge in G.edges]
    weight_dict = {edge: '{:.2f}'.format(G.edges[edge]['weight']) for edge in G.edges}
    if label_edges:
        nx.draw_networkx_edge_labels(G, nx.circular_layout(G), edge_labels=weight_dict,
                                     font_size=7)
    r = nx.draw_networkx_edges(G, nx.circular_layout(G), edge_color=weight_list)
    plt.colorbar(r)


def syndrome_plot(G_dec, G, index_dict=None, drawing_opts={}):
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
    bound_inds = G_dec.graph['boundary_points']
    points = [G_dec.nodes[i]['stabilizer'] for i in bound_inds]

    cube_inds = set(G_dec.nodes).difference(bound_inds)
    cubes = [G_dec.nodes[i]['stabilizer'] for i in cube_inds]
    # Default drawing options.
    draw_dict = {'show_nodes': False,
                 'label_nodes': None,
                 'label_cubes': True,
                 'label_boundary': False,
                 'legend': True}
    # Combine default dictionary with supplied dictionary, duplicates
    # favor supplied dictionary.
    drawing_opts = {**draw_dict, **drawing_opts}

    # Shape and font properties from the original graph.
    shape = np.array(G.graph.graph['dims'])
    font_props = G.graph.font_props
    # If show_nodes is True, get the axes object and legend from
    # CVGraph.sketch (this also plots the graph in the console).
    if drawing_opts['show_nodes']:
        ax = G.sketch(drawing_opts['label_nodes'])
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        # TODO: Initialize axes based on empty ax object from G.sketch()
        # but prevent from G.sketch() from plotting.
        fig = plt.figure(figsize=(2 * (np.sum(shape) + 2), 2 * (np.sum(shape) + 2)))
        ax = fig.gca(projection='3d')
        ax.tick_params(labelsize=font_props['size'])
        plt.xticks(range(0, 2 * shape[0] + 1))
        plt.yticks(range(0, 2 * shape[1] + 1))
        ax.set_zticks(range(0, 2 * shape[2] + 1))
        ax.set_xlabel('x', fontdict=font_props, labelpad=20)
        ax.set_ylabel('z', fontdict=font_props, labelpad=20)
        ax.set_zlabel('y', fontdict=font_props, labelpad=20)
        plt.rcParams['grid.color'] = "lightgray"
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
        if cube.parity():
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = '#FF000015'
        else:
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = '#00FF0015'
        if drawing_opts['label_cubes'] and index_dict:
            if cube in index_dict:
                ax.text(xmid, ymid, zmid, index_dict[cube], fontdict=font_props)

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

    if drawing_opts['label_boundary']:
        for point in points:
            ax.scatter(point[0], point[1], point[2], s=70, c='k')
            ax.text(point[0], point[1], point[2], index_dict[point], fontdict=font_props)

    # Define a legend for red/green cubes.
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#00FF0050', label='even parity'),
                       Patch(facecolor='#FF000050', label='odd parity')]
    if drawing_opts['legend']:
        ax.legend(handles=legend_elements, prop=font_props, loc='upper left')
    # Since CVGraph.sketch() legend has been overwritten, readd
    # it to the plot.
    if leg:
        ax.add_artist(leg)
    return ax


def assign_weights(CVG, method='unit'):
    """Assign weights to qubits in a hybrid CV graph state CVG.

    Args:
        CVG (CVGraph): the CVGraph whose error probabilities have been
            computed
        method (str): the method for weight assignment. By default,
            'unit', denoting weight 1 everyoewhere. For heuristic
            weight assignment from blueprint, use 'blueprint'

    Returns:
        None
    """
    G = CVG.graph
    # Get and denest the syndrome coordinates.
    syndrome_coords = [a for b in RHG.RHG_syndrome_coords(CVG.graph) for a in b]
    error_printed = False
    # Blueprint weight assignment dependent on type of neighbour.
    if method =='blueprint':
        for node in syndrome_coords:
            neighbors = G.subgraph(G[node]).nodes
            # List and number of p-squeezed states in neighborhood of node.
            p_list = [neighbors[v]['type'] for v in neighbors if neighbors[v]['type'] == 'p']
            p_count = len(p_list)
            # If conditional phase error information available, use it; otherwise set
            # error probability to 0 and print a message.
            if 'p_phase_cond' in G.nodes[node]:
                err_prob = G.nodes[node]['p_phase_cond']
            else:
                err_prob = 0
                if not error_printed:
                    print('Conditional Z error probabilities have not yet been computed. Please '
                          'use eval_Z_probs_cond() first. Setting error probability in 0 '
                          'and 1 swap-out case to 0. \n')
                    error_printed = True
            # Allow for taking log of 0.
            # TODO: Is this the best way to do it? Or can I just choose
            # an arbitrary small number?
            if err_prob == 0:
                err_prob = smallest_number
            # Dictionary of the form number of swapouts: error probability.
            weight_dict = {0: err_prob, 1: err_prob, 2: 1/4, 3: 1/3, 4: 2/5}
            # Do we a multiplicative factor here, followed by rounding to get
            # integral weights?
            G.nodes[node]['weight'] = - np.log(weight_dict[p_count])
        return
    # Naive weight assignment, unity weights.
    if method == 'unit':
        for node in syndrome_coords:
            G.nodes[node]['weight'] = 1


def CV_decoder(G, translator=basic_translate):
    """Convert homodyne outcomes to bit values according to translate.

    The inner (CV) decoder, aka translator, aka binning function. Set
    converted values to the bit_val attribute for nodes in G.

    Args:
        G: the CVGraph with homodyne outcomes computed.
        translator: the choice of binning function; basic_translate
            by default, which is the GKP-sqrt(pi) grid snapper.
    Returns:
        None
    """
    try:
        cv_values = G.hom_outcomes
    except Exception:
        print('A homodyne measurement has not yet been performed. Please '
              'use measure_p() first.')
        return
    bit_values = translator(cv_values)
    for i in range(len(bit_values)):
        G.graph.nodes[G.ind_dict[i]]['bit_val'] = bit_values[i]


def decoding_graph(G, draw=False, drawing_opts={}):
    """Create a decoding graph from the RHG lattice G.

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
    # An empty decoding graph.
    G_dec = nx.Graph(title='Decoding Graph')
    cubes = RHG.RHG_stabilizers(G)
    G_dec.add_nodes_from(cubes)
    # For stabilizer cubes sharing a vertex, define an edge between
    # them with weight equal to the weight assigned to the vertex.
    for (cube1, cube2) in it.combinations(G_dec, 2):
        common_vertex = set(cube1.coords()).intersection(cube2.coords())
        if common_vertex:
            coordinate = common_vertex.pop()
            weight = G.graph.nodes[coordinate]['weight']
            G_dec.add_edge(cube1, cube2, weight=weight, common_vertex=coordinate)

    # Include boundary vertices in the case of non-periodic boundary
    # conditions. (The following list will be empty absent a primal
    # boundary).
    # TODO: Deal with dual boundaries.
    bound_points = RHG.RHG_boundary_coords(G.graph)
    # For boundary points sharing a vertex with a stabilizer cube,
    # add an edge between them with weight equal to the weight assigned
    # to the vertex.
    for (cube, point) in it.product(G_dec, bound_points):
        if point in cube.coords():
            weight = G.graph.nodes[point]['weight']
            G_dec.add_edge(cube, point, weight=weight, common_vertex=point)

    # Relabel the nodes of the decoding graph to integers and define
    # the mapping between nodes and indices.
    G_relabelled = nx.convert_node_labels_to_integers(G_dec, label_attribute='stabilizer')
    mapping = dict(zip(G_dec.nodes(), range(0, G_dec.order())))

    # Indices of odd parity cubes and boundary vertices; add these
    # index lists to the graph attributes dictionary, to be used by the
    # matching graph.
    odd_parity_cubes = [cube for cube in cubes if cube.parity()]
    odd_parity_inds = [mapping[cube] for cube in odd_parity_cubes]
    G_relabelled.graph['odd_cubes'] = odd_parity_inds[:]
    bound_inds = [mapping[point] for point in bound_points]
    G_relabelled.graph['boundary_points'] = bound_inds[:]

    # Draw the lattice and the abstract decoding graph.
    if draw:
        graph_drawer(G_relabelled)
        syndrome_plot(G_relabelled,
                      G,
                      index_dict=mapping,
                      drawing_opts=drawing_opts)
    return G_relabelled


def matching_graph(G, alg='dijkstra', draw=False):
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
    G_match = nx.Graph(title='Matching Graph')

    # Get the indices of the odd parity cubes from teh decoding graph.
    odd_parity_inds = G.graph['odd_cubes']

    # Give shorter names to the Dijkstra shortest path algorithms.
    if alg == 'dijkstra':
        alg1 = sp.single_source_dijkstra
        alg2 = sp.multi_source_dijkstra

    # Combinations of odd-parity cubes.
    odd_ind_dict = {i: [] for i in odd_parity_inds[:-1]}
    odd_combs = it.combinations(odd_parity_inds, 2)
    for pair in odd_combs:
        odd_ind_dict[pair[0]] += [pair[1]]
    # Find the shortest paths between odd-parity cubes.
    for cube1 in odd_parity_inds[:-1]:
        lengths, paths = alg1(G, cube1)
        for cube2 in odd_ind_dict[cube1]:
            length = lengths[cube2]
            path = paths[cube2]
            # Add edge to the matching graph between the cubes, with weight
            # equal to the length of the shortest path.
            # TODO: Is the behavior correct for negative weights, or do I
            # want 1/weight or max_num - weight?
            G_match.add_edge(cube1, cube2, weight=length, inverse_weight=-length, path=path)

    # For non-periodic boundary conditions, include boundary vertices.
    # Get the indices of the boundary vertices from the decoding
    # graph.
    boundary_inds = G.graph['boundary_points']
    # Keep track of the boundary points that have been connected
    # to a cube.
    used_boundary_points = []
    if boundary_inds:
        for cube in odd_parity_inds:
            remaining_boundary = list(set(boundary_inds).difference(set(used_boundary_points)))
            # Find the shortest path from the any of the boundary
            # vertices to the cube. Note that we might wish to change
            # the sources to unused boundary vertices so that each
            # cube is connected to a unique vertex.
            point_paths = alg2(G, sources=remaining_boundary, target=cube)
            path = point_paths[1]
            point = path[0]
            length = point_paths[0]
            # Add edge to the matching graph between the cube and
            # the boundary vertex, with weight equal to the length
            # of the shortest path.
            G_match.add_edge(cube, point, weight=length, inverse_weight=-length, path=path)
            # Add to the list of used boundary vertices.
            used_boundary_points.append(point)

            # Add edge with weight 0 between any two boundary points.
            for (point1, point2) in it.combinations(used_boundary_points, 2):
                G_match.add_edge(point1, point2, weight=0, inverse_weight=0)

    # Add indices of used boundary points as a graph attribute.
    G_match.graph['used_boundary_points'] = used_boundary_points[:]

    if draw:
        graph_drawer(G_match)

    return G_match


def MWPM(G_match, G_dec, alg='blossom_nx', draw=False):
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
    if alg == 'blossom_nx':
        alg = nx.max_weight_matching
    matching = alg(G_match, maxcardinality=True, weight="inverse_weight")
    # TODO: Drop the requirement of the syndrome plot from having
    # to be plotted immediately prior to the matching.
    if draw:
        boundary = G_match.graph['used_boundary_points']
        for pair in matching:
            if pair not in it.product(boundary, boundary):
                xlist, ylist, zlist = [], [], []
                path = G_match.edges[pair]['path']
                for node in path:
                    stabe = G_dec.nodes[node]['stabilizer']
                    if isinstance(stabe, RHG.RHGCube):
                        x, y, z = stabe.midpoint()
                    else:
                        x, y, z = stabe
                        plt.plot(x, y, z, marker='2', ms=50, c='k')
                    xlist += [x]
                    ylist += [y]
                    zlist += [z]
                plt.title('Minimum-weight perfect matching', family='serif', size=20)
                plt.plot(xlist, ylist, zlist, 'o-k', ms=20, linewidth=5, c=np.random.rand(3))
        graph_drawer(G_match)
    return matching


def recovery(G_match, G_dec, G, matching, sanity_check=False):
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
    boundary = G_match.graph['used_boundary_points']
    for pair in matching:
        if pair not in it.product(boundary, boundary):
            path = G_match.edges[pair]['path']
            pairs = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for pair in pairs:
                common_vertex = G_dec.edges[pair]['common_vertex']
                G.graph.nodes[common_vertex]['bit_val'] ^= 1

    if sanity_check:
        G_dec_new = decoding_graph(G, draw=False)
        odd_cubes = G_dec_new.graph['odd_cubes']
        if odd_cubes:
            print('Unsatisfied stabilizers:', odd_cubes)
            return False
        else:
            print('Recovery succeeded!')
            return True
    # if check:
    #     parity = 0
    #     for stabe in G_dec:
    #         if isinstance(stabe, RHG.RHGCube):
    #             stabe._parity = None
    #             parity ^= stabe.parity()
    #     return (G, bool(1-parity))


def check_correction(G, plane=None, sheet=0, sanity_check=False):
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
    dims = np.array(G.graph.graph['dims'])
    dir_dict = {'x': 0, 'y': 1, 'z': 2}
    truth_dict = {'x': [], 'y': [], 'z': []}
    boundaries = G.graph.graph['boundaries']

    planes_to_check = []
    if plane:
        planes_to_check += [plane]
    else:
        for plane in ('x', 'y', 'z'):
            if boundaries[dir_dict[plane]] == 'periodic':
                planes_to_check.append(plane)

    minimum, maximum = sheet, sheet + 2
    for plane in planes_to_check:
        if sanity_check:
            minimum, maximum = 0, 2 * dims[dir_dict[plane]]
        for sheet in range(minimum, maximum, 2):
            slice_verts = RHG.RHG_slice_coords(G.graph, plane, sheet, boundaries='dual')
            syndrome_verts = [a for b in RHG.RHG_syndrome_coords(G.graph) for a in b]
            only_primal = set(slice_verts).intersection(set(syndrome_verts))
            parity = 0
            for node in only_primal:
                parity ^= G.graph.nodes[node]['bit_val']
            truth_dict[plane].append(bool(1 - parity))

    if sanity_check:
        print(truth_dict)

    all_surfaces = np.array([truth_dict[i][0] for i in planes_to_check])
    return np.all(all_surfaces)


def correct(G,
            inner='basic',
            outer='MWPM',
            weights='unit',
            draw=False,
            drawing_opts={},
            sanity_check=False):
    """Run through all the error-correction steps.

    Combines weight assignment, decoding and matching graph creation,
    minimum-weight-perfect matching, recovery, and correctness check.

    Args:
        G (CVGraph): the graph to decode and correct
        inner (str): the inner decoder; basic_translate by default
        outer (str): the outer decoder; MWPM by default
        weights (str): method for weight assignment; unit by default
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
    inner_dict = {'basic': basic_translate}
    outer_dict = {'MWPM': 'MWPM'}

    if inner:
        CV_decoder(G, translator=inner_dict[inner])

    if outer_dict[outer] == 'MWPM':
        assign_weights(G, method=weights)
        G_dec = decoding_graph(G, draw=draw, drawing_opts=drawing_opts)
        G_match = matching_graph(G_dec)
        matching = MWPM(G_match, G_dec, draw=draw)
        recovery(G_match, G_dec, G, matching, sanity_check=sanity_check)
    result = check_correction(G, sanity_check=sanity_check)
    return result


if __name__ == '__main__':
    boundaries = 'periodic'
    RHG_lattice = RHG.RHG_graph(2, boundaries=boundaries, polarity=1)

    swap_prob = 0
    delta = 0.1

    G = CVGraph(RHG_lattice, swap_prob=swap_prob, delta=delta)
    G.measure_p()
    G.eval_Z_probs_cond()

    dw = {'show_nodes': False, 'label_nodes': '', 'label_cubes': True,
          'label_boundary': False, 'legend': True}

    correct(G, inner='basic', outer='MWPM', weights='blueprint', draw=True, drawing_opts=dw)

