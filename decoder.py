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
# import itertools as it
# import networkx as nx
# import numpy as np
# from numpy.random import (multivariate_normal as mvn, default_rng as rng)
# import matplotlib.pyplot as plt
# import matplotlib.lines as mlines
# from scipy.special import erf

from graphstates import EGraph, CVGraph, RHG_graph


def decoding_graph(G):
    return G


def assign_weights(CVG):
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
        G.nodes[node]['weight'] = weight_dict[p_count]
    return

if __name__ == '__main__':
    RHG = RHG_graph(1)
    # RHG.draw(label=1)

    dim = 2
    delta = 0.1

    N = RHG.number_of_nodes()

    p = 0.3
    G = CVGraph(RHG, swap_prob=0.5, delta=delta)
    G.eval_Z_probs()
    # G.measure_p()
    # for label in {'v', 'e', 'h', 'b'}:
    #     G.sketch(label)
    assign_weights(G)