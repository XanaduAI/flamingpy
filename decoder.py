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
"""The decoder module."""

import sys
import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import itertools as it
import matplotlib.pyplot as plt

from graphstates import EGraph, CVGraph, basic_translate
import RHG

smallest_number = sys.float_info.min
largest_number = sys.float_info.max

def graph_drawer(G):
    """Convenience function for drawing decoding and matching graphs."""
    title = G.graph['title']
    fig = plt.figure()
    plt.title(title, family='serif', size=10)
    weight_list = [edge[2]['weight'] for edge in G.edges.data()]
    drawing = nx.draw_circular(G,
                               edgelist=[],
                               with_labels=True,
                               node_color='k',
                               font_size=7,
                               font_color='w',
                               font_family='serif')
    r = nx.draw_networkx_edges(G,
                               nx.circular_layout(G),
                               edge_color=weight_list)
    plt.colorbar(r)


def syndrome_plot(G_dec, G, index_dict=None, code='primal', drawing_opts={}, bc='periodic'):
    """Convenience function for drawing the syndrome plot given the underlying
    graph G and a list of cubes cubes."""

    if bc != 'periodic':
        bound_inds = G_dec.graph['boundary_points']
        points = [G_dec.nodes[i]['stabilizer'] for i in bound_inds]
    else:
        bound_inds = {}

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
    drawing_opts = dict(draw_dict, **drawing_opts)

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
        fig = plt.figure(figsize=(2 * (np.sum(shape)+2),  2 * (np.sum(shape)+2)))
        ax = fig.gca(projection='3d')
        ax.tick_params(labelsize=font_props['size'])
        plt.xticks(range(0, 2*shape[0] + 1))
        plt.yticks(range(0, 2*shape[1] + 1))
        ax.set_zticks(range(0, 2*shape[2] + 1))
        ax.set_xlabel('x', fontdict=font_props, labelpad=20)
        ax.set_ylabel('z', fontdict=font_props, labelpad=20)
        ax.set_zlabel('y', fontdict=font_props, labelpad=20)
        plt.rcParams['grid.color'] = "lightgray"
        leg = None
    # Illustrate stabilizers with voxels colored green for even
    # parity and red for odd pariy.
    filled = np.zeros(shape, dtype=object)
    for cube in cubes:
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
        if draw_dict['label_cubes'] and index_dict:
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

    if bc != 'periodic' and drawing_opts['label_boundary']:
        print(len(points))
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


def assign_weights(CVG, method='naive', code='primal'):
    """Assign weights to qubits in a hybrid CV graph state.

    By default, use unity weights; otherwise use the heristic weights
    from the blueprint."""

    G = CVG.graph
    dims = G.graph['dims']
    # Get and denest the syndrome coordinates.
    syndrome_coords = [a for b in RHG.RHG_syndrome_coords(dims, code) for a in b]
    error_printed = False
    for node in syndrome_coords:
        # Naive weight assignment, unity weights.
        if method=='naive':
            G.nodes[node]['weight'] = 1
        # Blueprint weight assignment dependent on type of neighbour.
        if method=='blueprint':
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


def CV_decoder(G, translator=basic_translate):
    """The inner (CV) decoder, aka translator, aka binning function.

    Convert homodyne outcomes to bit values according to translator."""
    try:
        cv_values = G.hom_outcomes
    except Exception:
        print('A homodyne measurement has not yet been performed. Please '
              'use measure_p() first.')
        return
    bit_values = translator(cv_values)
    for i in range(len(bit_values)):
        G.graph.nodes[G.ind_dict[i]]['bit_val'] = bit_values[i]


def decoding_graph(G, code='primal', bc='periodic', draw=False, drawing_opts={}):
    """Create a decoding graph from the RHG lattice G. 

    Note that one ought first compute the phase error probabilities, 
    conduct a homodyne measurement, and translate the outcomes on G."""

    # An empty decoding graph.
    G_dec = nx.Graph(title='Decoding Graph')
    cubes = RHG.RHG_stabilizers(G, code)
    G_dec.add_nodes_from(cubes)

    # For stabilizer cubes sharing a vertex, define an edge between
    # them with weight equal to the weight assigned to the vertex.
    for (cube1, cube2) in it.combinations(G_dec, 2):
        common_vertex = set(cube1.coords()).intersection(cube2.coords())
        if common_vertex:
            weight = G.graph.nodes[common_vertex.pop()]['weight']
            G_dec.add_edge(cube1, cube2, weight=weight)

    # Include boundary vertices in the case of non-periodic boundary
    # conditions. Note that info about boundary conditions will
    # probably be located later in the graph G.
    if bc != 'periodic':
        # For boundary points sharing a vertex with a stabilizer cube,
        # add an edge between them with weight equal to the weight assigned
        # to the vertex.
        bound_points = RHG.RHG_boundary_coords(np.array(G.graph.graph['dims']))
        for (cube, point) in it.product(G_dec, bound_points):
                if point in cube.coords():
                    weight = G.graph.nodes[point]['weight']
                    G_dec.add_edge(cube, point, weight=weight)

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
    if bc != 'periodic':
        bound_inds = [mapping[point] for point in bound_points]
        G_relabelled.graph['boundary_points'] = bound_inds[:]

    # Draw the lattice and the abstract decoding graph.
    if draw:
        graph_drawer(G_relabelled)
        ax = syndrome_plot(G_relabelled,
                           G,
                           index_dict=mapping,
                           code=code,
                           drawing_opts=drawing_opts,
                           bc=bc)
    return G_relabelled


def matching_graph(G, bc='periodic', alg='dijkstra', draw=False, drawing_opts={}):
    """Create a matching graph from the decoding graph G.

    Create graph according to algorithm alg. If draw, draw the matching
    graph with a color bar."""

    # An empty matching graph.
    G_match = nx.Graph(title='Matching Graph')

    # Get the indices of the odd parity cubes from teh decoding graph.
    odd_parity_inds = G.graph['odd_cubes']

    # Give shorter names to the Dijkstra shortest path algorithms.
    if alg == 'dijkstra':
        alg1 = sp.single_source_dijkstra
        alg2 = sp.multi_source_dijkstra

    # Find the shortest paths between odd-parity cubes.
    for (cube1, cube2) in it.combinations(odd_parity_inds, 2):
        length, path = alg1(G, cube1, cube2)
        # Add edge to the matching graph between the cubes, with weight
        # equal to the length of the shortest path.
        if length == 0:
            length = smallest_number
        G_match.add_edge(cube1, cube2, weight=length, inverse_weight=1/length, path=path)
    # For non-periodic boundary conditions, include boundary vertices.
    used_boundary_points = []
    if bc != 'periodic':
        # Get the indics of the boundary vertices from the decoding 
        # graph.
        boundary_inds = G.graph['boundary_points']
        # Keep track of the boundary points that have been connected
        # to a cube.
        for cube in odd_parity_inds:
            # Find the shortest path from the any of the boundary
            # vertices to the cube. Note that we might wish to change
            # the sources to unused boundary vertices so that each
            # cube is connected to a unique vertex.
            point_paths = alg2(G, sources=boundary_inds, target=cube)
            path = point_paths[1]
            point = path[0]
            length = point_paths[0]
            if length == 0:
                length = smallest_number
            # Add edge to the matching graph between the cube and
            # the boundary vertex, with weight equal to the length
            # of the shortest path.
            G_match.add_edge(cube, point, weight=length, inverse_weight=1/length, path=path)
            # Add to the list of used boundary vertices.
            used_boundary_points.append(point)

        # Add edge with weight 0 between any two boundary points.
        for (point1, point2) in it.combinations(used_boundary_points, 2):
            G_match.add_edge(point1, point2, weight=smallest_number, inverse_weight=largest_number)

    # Add indices of used boundary points as a graph attribute.
    G_match.graph['used_boundary_points'] = used_boundary_points[:]

    if draw:
        graph_drawer(G_match)

    return G_match


def MWPM(G_match, G_dec, alg='blossom_nx', bc='periodic', draw=False):
    """Minimum-weight-perfect matching on matching graph G."""

    if alg == 'blossom_nx':
        alg = nx.max_weight_matching
    matching = alg(G_match, maxcardinality=True, weight="inverse_weight")
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

if __name__ == '__main__':
    RHG_lattice = RHG.RHG_graph(2, pol=1)

    swap_prob = 0.2
    delta = 0.01

    G = CVGraph(RHG_lattice, swap_prob=swap_prob, delta=delta)
    G.measure_p()
    G.eval_Z_probs_cond()
    CV_decoder(G)
    assign_weights(G, method='blueprint')

    dw = {'show_nodes': False, 'label_nodes': '', 'label_cubes': True,
          'label_boundary': False, 'legend':False}
    G_dec = decoding_graph(G, bc='', drawing_opts=dw, draw=True)
    G_match = matching_graph(G_dec, bc='', draw=False)
    matching = MWPM(G_match, G_dec, draw=True)