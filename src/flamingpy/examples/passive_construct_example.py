# Copyright 2022 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Example of simulating Xanadu's passive architecture."""
import time
from flamingpy import RHG
from flamingpy.decoder import correct
from flamingpy.graphstates import CVGraph
from flamingpy.passive_construct import BS_network, reduce_macro_and_simulate


total = 1
d = 2
delta = 0.01
p_swap = 0.9
boundaries = "periodic"
# The lattice with macronodes.
RHG_macro = RHG.RHG_graph(d, boundaries=boundaries, macronodes=True, polarity=False)
RHG_macro.index_generator()
RHG_macro.adj_generator(sparse=True)
# The reduced lattice.
RHG_code = RHG.RHGCode(d, boundaries=boundaries)
RHG_reduced = RHG_code.graph
RHG_reduced.index_generator()
# The empty CV state, uninitiated with any error model.
CVRHG_reduced = CVGraph(RHG_reduced)
# Define the 4X4 beamsplitter network for a given macronode.
# star at index 0, planets at indices 1-3.
bs_network = BS_network(4)

start = time.time()
correction_time = 0
successes = 0
for trial in range(total):
    # The empty CV state, uninitiated with any error model.
    reduce_macro_and_simulate(RHG_macro, RHG_reduced, CVRHG_reduced, bs_network, p_swap, delta)
    dw = {
        "show_nodes": True,
        "label": "bit_val",
        "label_cubes": False,
        "label_boundary": False,
        "legend": True,
    }
    weight_options = {
        "method": "blueprint",
        "prob_precomputed": True,
    }
    decoder = {"outer": "MWPM"}
    correct_start = time.time()
    c = correct(
        code=RHG_code,
        decoder=decoder,
        weight_options=weight_options,
    )
    successes += int(c)
    correct_end = time.time()
    correction_time += correct_end - correct_start
end = time.time()
error = (total - successes) / total
print("Error rate: ", error)
print("Prep time: ", end - start - correction_time)
print("Correction time: ", correction_time)
print("Total time: ", end - start)
