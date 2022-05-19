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
"""Example of simulating Xanadu's passive and static architecture."""

from flamingpy.codes import SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import reduce_macronode_graph, splitter_symp
from flamingpy.decoders.decoder import correct

# Number of trials
total = 100
# Code parameters
d = 3
boundaries = "periodic"
ec = "primal"
# Noise parameters
delta = 0.1
p_swap = 0.25

# The reduced lattice.
RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)
RHG_reduced = RHG_code.graph
RHG_reduced.index_generator()
# The lattice with macronodes
pad_bool = boundaries != "periodic"
RHG_macro = RHG_code.graph.macronize(pad_boundary=pad_bool)
RHG_macro.index_generator()
RHG_macro.adj_generator(sparse=True)

# Define the 4X4 beamsplitter network for a given macronode.
# star at index 0, planets at indices 1-3.
bs_network = splitter_symp(4)

noise_model = {"noise": "grn", "delta": delta, "p_swap": p_swap, "sampling_order": "two-step"}

successes = 0
for trial in range(total):
    # The empty CV state, uninitiated with any error model.
    reduce_macronode_graph(RHG_macro, RHG_reduced, CVLayer, noise_model, bs_network)
    decoder = {"outer": "MWPM"}
    decoder_opts = {"backend": "networkx"}
    if decoder["outer"] == "MWPM":
        weight_options = {
            "method": "blueprint",
            "prob_precomputed": True,
        }
    else:
        weight_options = None
    c = correct(
        code=RHG_code, decoder=decoder, weight_options=weight_options, decoder_opts=decoder_opts
    )
    successes += int(c)

error = (total - successes) / total
print("Error rate: ", error)
