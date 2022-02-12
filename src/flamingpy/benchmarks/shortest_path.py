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
"""Benchmark shortest-path-finding algorithms between networkx and retworkx."""
import time

import matplotlib.pyplot as plt

from flamingpy import decoder as dec
from flamingpy.graphstates import CVGraph
from flamingpy.RHG import RHGCode, alternating_polarity

# How many simulations to do for each algorithm.
num_trials = 10

# DV (outer) code
distance = 5
boundaries = "periodic"

# Noise model
p_swap = 0.2
delta = 0.1
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Decoding options
decoder = {"inner": "basic", "outer": "MWPM"}
weight_options = {
    "method": "blueprint",
    "integer": True,
    "multiplier": 100,
    "delta": delta,
}

times = {
    "networkx": [],
    "retworkx": [],
}

for backend in ["networkx", "retworkx"]:
    print(f"* {backend}")
    # Build code
    RHG_code = RHGCode(
        distance=distance, boundaries=boundaries, polarity=alternating_polarity, backend=backend
    )
    RHG_lattice = RHG_code.graph
    for i in range(num_trials):
        print(f"-- {i} --")
        CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)
        # Apply noise
        CVRHG.apply_noise(cv_noise)
        # Measure syndrome
        CVRHG.measure_hom("p", RHG_code.syndrome_inds)
        # Inner decoder
        dec.CV_decoder(RHG_code, translator=dec.GKP_binner)

        before = time.time()
        matching_graph = dec.build_match_graph(RHG_code, weight_options, backend)
        after = time.time()
        times[backend].append(after - before)

plt.figure()
for backend in ["networkx", "retworkx"]:
    plt.hist(times[backend], label=backend)
plt.legend()
plt.xlabel("Times [seconds]")
plt.ylabel("Count")
plt.title(f"Building matching graph for code distance {distance}")
plt.savefig(f"benchmark_shortest_path_distance_{distance}.pdf")
