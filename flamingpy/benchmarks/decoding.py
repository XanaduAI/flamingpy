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
"""Benchmark for the complete decoding procedure comparing NetworkX and
retworkx."""
import time

import matplotlib.pyplot as plt

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.decoders import decoder as dec
from flamingpy.noise import CVLayer

# How many simulations to do for each algorithm
num_trials = 10

# DV (outer) code parameters
distance = 5

# Boundaries ("open", "toric" or "periodic")
boundaries = "periodic"

# Noise parameters
p_swap = 0.1
delta = 0.1

# Decoding options
decoder = "MWPM"
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
    RHG_code = SurfaceCode(
        distance=distance, boundaries=boundaries, polarity=alternating_polarity, backend=backend
    )
    # Instantiate the noise layer
    CVRHG = CVLayer(RHG_code, delta=delta, p_swap=p_swap)
    for i in range(num_trials):
        print(f"-- {i} --")
        # Apply noise
        CVRHG.apply_noise()
        # Inner decoder
        before = time.time()
        dec.correct(
            RHG_code,
            decoder,
            weight_options,
            decoder_opts={"backend": backend},
        )
        after = time.time()
        times[backend].append(after - before)

plt.figure()
for backend in ["networkx", "retworkx"]:
    plt.hist(times[backend], label=backend)
plt.legend()
plt.xlabel("Times [seconds]")
plt.ylabel("Count")
plt.title(f"Decoding for code distance {distance}")

if __name__ == "__main__":
    plt.savefig(f"benchmark_decoding_distance_{distance}.pdf")
