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
"""Module for Monte Carlo simulations for estimating FT thresholds."""

#%%
import numpy as np
from decoder import *
#%%

def monte_carlo(iterations,L,swap_prob,delta):
    corrected = 0
    for i in range(iterations):
        RHG_lattice = RHG.RHG_graph(L, pol=1)
        G = CVGraph(RHG_lattice, swap_prob=swap_prob, delta=delta)
        G.eval_Z_probs()
        G.measure_p()
        G.translate_outcomes()
        assign_weights(G, method='blueprint')
        
        dw = {'show_nodes': False, 'label_nodes': '', 'label_cubes': True,
              'label_boundary': False, 'legend':False}
        G_dec = decoding_graph(G, bc='', drawing_opts=dw, draw=False)
        G_match = matching_graph(G_dec, bc='', draw=False)
        matching = MWPM(G_match, G_dec, draw=False)
        
        #if correction_consitions_satisfied:
            #corrected = corrected+1
        
        
    return corrected/iterations

#%%
iterations = 10000
L = 2
swap_prob = 0.2
delta = 0.01

monte_carlo(iterations,L,swap_prob,delta)