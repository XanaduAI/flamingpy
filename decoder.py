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


class RHGCube:
    """A class for representing an RHG latice cube. Initialized by supplying
    a CVGraph object corresponding to the subgraph induced by the six facial
    nodes of the cube."""
    def __init__(self, G):
        self.cvgraph = G

    def parity(self):
        G = self.cvgraph
        bit_vals = [G.graph.nodes[node]['bit_val'] for node in G.graph]
        return int(np.sum(bit_vals) % 2)

    def coords(self):
        return [tup for tup in self.cvgraph.graph.nodes]

    def xlims(self):
        xs = [tup[0] for tup in self.coords()]
        xmin, xmax = np.min(xs), np.max(xs)
        xmed = (xmin + xmax) / 2
        return (xmin, xmed, xmax)

    def ylims(self):
        ys = [tup[1] for tup in self.coords()]
        ymin, ymax = np.min(ys), np.max(ys)
        ymed = (ymin + ymax) / 2
        return (ymin, ymed, ymax)

    def zlims(self):
        zs = [tup[2] for tup in self.coords()]
        zmin, zmax = np.min(zs), np.max(zs)
        zmed = (zmin + zmax) / 2
        return (zmin, zmed, zmax)


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
        if err_prob == 0:
            err_prob = 1e-10
        weight_dict = {0: err_prob, 1: err_prob, 2: 1/4, 3: 1/3, 4: 2/5}
        G.nodes[node]['weight'] = -np.log(weight_dict[p_count])
    return


def decoding_graph(G, code='primal', draw=False):
    """Create a decoding graph from the RHG lattice G. Note that one 
    must first compute the phase error probabilities, conduct a 
    homodyne |measurement, and translate the outcomes on G."""

    G_dec = nx.Graph(title='Decoding Graph')
    stabes = RHG_stabilizers(G, code)
    cubes = [RHGCube(stabe) for stabe in stabes]
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
        fig = plt.figure(figsize=(np.sum(shape)+2, np.sum(shape)+2))
        ax = fig.gca(projection='3d')
        # ax = G.sketch()
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
        ax.legend(handles=legend_elements, fontsize=30)
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
