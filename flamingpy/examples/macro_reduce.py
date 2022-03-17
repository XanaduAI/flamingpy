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
from flamingpy.cv.macro_reduce import BS_network, reduce_macro_and_simulate
from flamingpy.decoders.decoder import correct

# Number of trials
total = 1
# Code parameters
d = 3
boundaries = "open"
ec = "primal"
# Noise parameters
delta = 0.1
p_swap = 0.8

# The reduced lattice.
RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)
RHG_reduced = RHG_code.graph
RHG_reduced.index_generator()
# The lattice with macronodes
pad_bool = True if boundaries == "open" else False
RHG_macro = RHG_code.graph.macronize(pad_boundary=pad_bool)
RHG_macro.index_generator()
RHG_macro.adj_generator(sparse=True)

# The empty CV layer, uninitiated with any error model.
CVRHG_reduced = CVLayer(RHG_reduced)
# Define the 4X4 beamsplitter network for a given macronode.
# star at index 0, planets at indices 1-3.
bs_network = BS_network(4)

successes = 0
for trial in range(total):
    # The empty CV state, uninitiated with any error model.
    reduce_macro_and_simulate(RHG_macro, RHG_reduced, CVRHG_reduced, bs_network, p_swap, delta)
    weight_options = {
        "method": "blueprint",
        "prob_precomputed": True,
    }
    decoder = {"outer": "MWPM"}
    c = correct(code=RHG_code, decoder=decoder, weight_options=weight_options)
    successes += int(c)

error = (total - successes) / total
print("Error rate: ", error)
