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
'''The decoder module'''

import numpy as np
import networkx as nx
import networkx.algorithms.shortest_paths as sp
import itertools as it
import matplotlib.pyplot as plt

from graphstates import EGraph, CVGraph, RHG_graph


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
                               font_size=8,
                               font_color='w',
                               font_family='serif')
    r = nx.draw_networkx_edges(G,
                               nx.circular_layout(G),
                               edge_color=weight_list)
    plt.colorbar(r)


def RHG_syndrome_coords(dims, code='primal'):
    """Return a list of lists of coordinates associated to six-body
    X stabilizer elements of the RHG lattice."""

    nx, ny, nz = dims
    combs = it.product(range(nx), range(ny), range(nz))
    if code == 'primal':
        all_syndrome6 = [[(2*i+1, 2*j+1, 2*k),
                          (2*i+1, 2*j, 2*k+1),
                          (2*i, 2*j+1, 2*k+1),
                          (2*i+1, 2*j+1, 2*k+2),
                          (2*i+1, 2*j+2, 2*k+1),
                          (2*i+2, 2*j+1, 2*k+1)
                          ] for (i, j, k) in combs]
    return all_syndrome6


def RHG_stabilizers(G, code='primal'):
    """Return a list of subgraphs induced by the qubits with cordinates
    from RHG_syndrome_coords."""

    syn_list = RHG_syndrome_coords(G.graph.graph['dims'], code)
    cube_list = []
    for l in syn_list:
        cube_graph = G.graph.subgraph(l)
        cube_list.append(CVGraph(cube_graph))
    return cube_list


def parity(cube):
    """Obtain the parity of the bit values following a measurement of
    the a six-body X stabilizer element represented in cube."""

    bit_vals = [cube.graph.nodes[node]['bit_val'] for node in cube.graph]
    return np.sum(bit_vals) % 2


def assign_weights(CVG):
    """Assign weights to qubits in a hybrid p-squeezed/GKP CV graph 
    state per the heuristic procedure in the blueprint."""

    G = CVG.graph
    for node in G:
        neighbors = G.subgraph(G[node]).nodes
        p_list = [neighbors[v]['type'] for v in neighbors if neighbors[v]['type'] == 'p']
        p_count = len(p_list)
        try:
            err_prob = G.nodes[node]['p_phase']
        except Exception:
            print('Z error probabilities have not yet been computed. Please '
                  'use eval_Z_probs() first.')
            return
        weight_dict = {0: err_prob, 1: err_prob, 2: 1/4, 3: 1/3, 4: 2/5}
        G.nodes[node]['weight'] = -np.log(weight_dict[p_count])
    return


def decoding_graph(G, code='primal', draw=False):
    """Create a decoding graph from the RHG lattice G. Note that one 
    must first compute the phase error probabilities, conduct a 
    homodyne |measurement, and translate the outcomes on G."""

    G_dec = nx.Graph(title='Decoding Graph')
    stabes = RHG_stabilizers(G, code)
    for stabe in stabes:
        if parity(stabe) == 1:
            G_dec.add_node(stabe)
    for (stabe1, stabe2) in it.combinations(G_dec, 2):
        common_vertex = set(stabe1.graph.nodes).intersection(stabe2.graph.nodes)
        if common_vertex:
            weight = G.graph.nodes[common_vertex.pop()]['weight']
            G_dec.add_edge(stabe1, stabe2, weight=weight)
    G_relabelled = nx.convert_node_labels_to_integers(G_dec, label_attribute='stabilizer')
    if draw:
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
        for stabe in tup[1]:
            if tup[0] != stabe:
                w = tup[1][stabe]
                G_match.add_edge(tup[0], stabe, weight=w)
    if draw:
        graph_drawer(G_match)
    return G_match


if __name__ == '__main__':
    RHG = RHG_graph(2)

    swap_prob = 0.2
    delta = 0.1

    G = CVGraph(RHG, swap_prob=swap_prob, delta=delta)
    G.eval_Z_probs()
    G.measure_p()
    G.translate_outcomes()
    assign_weights(G)

    G.sketch('bit_val')
    G_dec = decoding_graph(G, draw=True)
    G_match = matching_graph(G_dec, draw=True)
