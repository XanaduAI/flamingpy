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


def decoding_graph(G, code='primal', draw=False, drawing_opts={}):
    """Create a decoding graph from the RHG lattice G. Note that one 
    must first compute the phase error probabilities, conduct a 
    homodyne |measurement, and translate the outcomes on G."""

    draw_dict = {'show_nodes': False, 'label_nodes': None, 'label_cubes': True}
    # Combine default dictionary with supplied dictionary, duplicates
    # favor supplied dictionary.
    drawing_opts = dict(draw_dict, **drawing_opts)

    G_dec = nx.Graph(title='Decoding Graph')
    cubes = RHG.RHG_stabilizers(G, code)
    G_dec.add_nodes_from(cubes)
    for (cube1, cube2) in it.combinations(G_dec, 2):
        common_vertex = set(cube1.coords()).intersection(cube2.coords())
        if common_vertex:
            weight = G.graph.nodes[common_vertex.pop()]['weight']
            G_dec.add_edge(cube1, cube2, weight=weight)
    G_relabelled = nx.convert_node_labels_to_integers(G_dec, label_attribute='stabilizer')
    mapping = dict(zip(G_dec.nodes(), range(0, G_dec.order())))
    if draw:
        shape = 2 * np.array(G.graph.graph['dims'])
        if drawing_opts['show_nodes']:
            ax = G.sketch(drawing_opts['label_nodes'])
            leg = ax.get_legend()
        else:
            fig = plt.figure(figsize=(np.sum(shape)+2, np.sum(shape)+2))
            ax = fig.gca(projection='3d')
            leg = None
        # ax = G.graph.draw(0, 0)
        voxels = np.zeros(shape)
        colors = np.zeros(shape, dtype=object)
        for cube in cubes:
            xmin, xmed, xmax = cube.xlims()
            ymin, ymed, ymax = cube.ylims()
            zmin, zmed, zmax = cube.zlims()
            voxels[xmin:xmax, ymin:ymax, zmin:zmax] = True
            if cube.parity() == 0:
                colors[xmin:xmax, ymin:ymax, zmin:zmax] = '#00FF0050'
            else:
                colors[xmin:xmax, ymin:ymax, zmin:zmax] = '#FF000050'
            if cube in mapping:
                ax.text(xmed, ymed, zmed, mapping[cube], fontsize=30)
        ax.voxels(voxels, facecolors=colors)
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='#00FF0050', label='+1 parity'),
                           Patch(facecolor='#FF000050', label='-1 parity')]
        font_props = G.graph.font_props
        ax.legend(handles=legend_elements, prop=font_props, loc='upper left')
        if leg:
            ax.add_artist(leg)
        graph_drawer(G_relabelled)
    return G_relabelled


def matching_graph(G, alg='dijkstra', draw=False):
    """Create a matching graph from the decoding graph G according to
    algorithm alg. If draw, draw the matching graph with a color bar."""

    if alg == 'dijkstra':
        d_alg = sp.all_pairs_dijkstra_path_length
    path_lengths = d_alg(G)
    G_match = nx.Graph(title='Matching Graph')
    for tup in path_lengths:
        for cube in tup[1]:
            parity1 = G.nodes[cube]['stabilizer'].parity()
            parity2 = G.nodes[tup[0]]['stabilizer'].parity()
            if parity1 and parity2:
                if tup[0] != cube:
                    w = tup[1][cube]
                    G_match.add_edge(tup[0], cube, weight=w)
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
    assign_weights(G)

    dw = {'show_nodes': True, 'label_nodes': False, 'label_cubes': False}
    G_dec = decoding_graph(G, draw=True, drawing_opts=dw)
    G_match = matching_graph(G_dec, draw=True)
