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
"""Benchmark minimum-weight matching via NetworkX, retworkx, and lemon
backends."""
import time

import matplotlib.pyplot as plt

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.decoders import decoder as dec
from flamingpy.decoders.mwpm.algos import build_match_graph
from flamingpy.noise import CVLayer

# How many simulations to do for each algorithm
num_trials = 10

# DV (outer) code
distance = 3
boundaries = "periodic"
RHG_code = SurfaceCode(distance=distance, boundaries=boundaries, polarity=alternating_polarity)

# Instantiate the noise layer
p_swap = 0.2
delta = 0.1
CVRHG = CVLayer(RHG_code, delta=delta, p_swap=p_swap)

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
    "lemon": [],
    "retworkx": [],
}

for alg in ["networkx", "lemon", "retworkx"]:
    print(f"\n* {alg}")
    for i in range(num_trials):
        print(f"-- {i} --")
        # Apply noise
        CVRHG.apply_noise()
        # Manually decode so as to benchmark just the matching portion
        dec.assign_weights(RHG_code, "MWPM", **weight_options)
        G_match = build_match_graph(RHG_code, "primal", alg)
        before = time.time()
        matching = G_match.min_weight_perfect_matching()
        after = time.time()
        times[alg].append(after - before)

plt.figure()
# bins = np.logspace(-3, 0, 30)
for alg in ["networkx", "retworkx", "lemon"]:
    plt.hist(times[alg], bins=100, label=alg)
plt.legend()
plt.xscale("log")
plt.xlabel("Times [seconds]")
plt.ylabel("Count")
plt.title(f"Matching for code distance {distance}")

if __name__ == "__main__":
    plt.savefig(f"benchmark_matching_distance_{distance}.pdf")
