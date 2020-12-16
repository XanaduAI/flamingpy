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

import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import itertools as it
import matplotlib.pyplot as plt

from graphstates import EGraph, CVGraph, basic_translate
import RHG


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
            # If phase error information available, use it; otherwise set
            # error probability to 0 and print a message.
            if 'p_phase' in G.nodes[node]:
                err_prob = G.nodes[node]['p_phase']
            else:
                err_prob = 0
                if not error_printed:
                    print('Z error probabilities have not yet been computed. Please '
                          'use eval_Z_probs() first. Setting error probability in 0 '
                          'and 1 swap-out case to 0. \n')
                    error_printed = True
            # Allow for taking log of 0.
            if err_prob == 0:
                err_prob = 1e-10
            # Dictionary of the form number of swapouts: error probability.
            weight_dict = {0: err_prob, 1: err_prob, 2: 1/4, 3: 1/3, 4: 2/5}
            # Multiplicative factor.
            mult = 100
            G.nodes[node]['weight'] = - mult * np.log(weight_dict[p_count])
    return


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
        # Default drawing options.
        draw_dict = {'show_nodes': False, 'label_nodes': None, 'label_cubes': True}
        # Combine default dictionary with supplied dictionary, duplicates
        # favor supplied dictionary.
        drawing_opts = dict(draw_dict, **drawing_opts)

        shape = 2 * np.array(G.graph.graph['dims'])
        # If show_nodes is True, get the axes object and legend from
        # CVGraph.sketch (this also plots the graph in the console).
        if drawing_opts['show_nodes']:
            ax = G.sketch(drawing_opts['label_nodes'])
            leg = ax.get_legend()
        # If show_nodes is False, create a new figure with size
        # determined by the dimensions of the lattice.
        else:
            fig = plt.figure(figsize=(np.sum(shape)+2, np.sum(shape)+2))
            ax = fig.gca(projection='3d')
            leg = None
        # Illustrate stabilizers with voxels colored green for even
        # parity and red for odd pariy.
        voxels = np.zeros(shape)
        colors = np.zeros(shape, dtype=object)
        for cube in cubes:
            # Obtain smallest, largest, and middle coordinates for each
            # cube.
            xmin, xmid, xmax = cube.xlims()
            ymin, ymid, ymax = cube.ylims()
            zmin, zmid, zmax = cube.zlims()
            # Fill in the color arrays depending on parity.
            if cube.parity():
                colors[xmin:xmax, ymin:ymax, zmin:zmax] = '#FF000050'
            else:
                colors[xmin:xmax, ymin:ymax, zmin:zmax] = '#00FF0050'
            if cube in mapping:
                ax.text(xmid, ymid, zmid, mapping[cube], fontsize=30)
        # Use colors array in voxel parameter because it is 0 in the
        # right places.
        ax.voxels(colors, facecolors=colors)
        # Define a legend for red/green cubes.
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='#00FF0050', label='+1 parity'),
                           Patch(facecolor='#FF000050', label='-1 parity')]
        font_props = G.graph.font_props
        ax.legend(handles=legend_elements, prop=font_props, loc='upper left')
        # Since CGGraph.sketch() legend has been overwritten, readd
        # it to the plot.
        if leg:
            ax.add_artist(leg)
        # Draw the abstract graph representation.
        graph_drawer(G_relabelled)
    return G_relabelled

def matching_graph(G, bc='periodic', alg='dijkstra', draw=False):
    """Create a matching graph from the decoding graph G.

    Create graph according to algorithm alg. If draw, draw the matching
    graph with a color bar."""

    # An empty matching graph.
    G_match = nx.Graph(title='Matching Graph')

    # Get the indices of the odd parity cubes from teh decoding graph.
    odd_parity_inds = G.graph['odd_cubes']

    # Give shorter names to the Dijkstra shortest path algorithms.
    if alg == 'dijkstra':
        alg1 = sp.dijkstra_path_length
        alg2 = sp.multi_source_dijkstra

    # Find the shortest paths between odd-parity cubes.
    for (cube1, cube2) in it.combinations(odd_parity_inds, 2):
        path_length = alg1(G, cube1, cube2)
        # Add edge to the matching graph between the cubes, with weight
        # equal to the length of the shortest path.
        G_match.add_edge(cube1, cube2, weight=path_length)

    # For non-periodic boundary conditions, include boundary vertices.
    if bc != 'periodic':
        # Get the indics of the boundary vertices from the decoding 
        # graph.
        boundary_inds = G.graph['boundary_points']
        # Keep track of the boundary points that have been connected
        # to a cube.
        used_boundary_points = []
        for cube in odd_parity_inds:
            # Find the shortest path from the any of the boundary
            # vertices to the cube. Note that we might wish to change
            # the sources to unused boundary vertices so that each
            # cube is connected to a unique vertex.
            point_paths = alg2(G, sources=boundary_inds, target=cube)
            point = point_paths[1][0]
            # Add edge to the matching graph between the cube and
            # the boundary vertex, with weight equal to the length
            # of the shortest path.
            G_match.add_edge(cube, point, weight=point_paths[0])
            # ADd to the list of used boundary vertices.
            used_boundary_points.append(point)

        # Add edge with weight 0 between any two boundary points.
        for (point1, point2) in it.combinations(used_boundary_points, 2):
            G_match.add_edge(point1, point2, weight=0)

    if draw:
        graph_drawer(G_match)

    return G_match


if __name__ == '__main__':
    RHG_lattice = RHG.RHG_graph(2, pol=1)

    swap_prob = 0.2
    delta = 0.01

    G = CVGraph(RHG_lattice, swap_prob=swap_prob, delta=delta)
    G.eval_Z_probs()
    G.measure_p()
    G.translate_outcomes()
    assign_weights(G, method='blueprint')

    dw = {'show_nodes': True, 'label_nodes': False, 'label_cubes': False}
    G_dec = decoding_graph(G, bc='non-periodic', draw=True, drawing_opts=dw)
    G_match = matching_graph(G_dec, bc='non-periodic', draw=True)