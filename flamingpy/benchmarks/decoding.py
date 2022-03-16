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

# pylint: disable=no-member

import time

import matplotlib.pyplot as plt

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders import decoder as dec

# How many simulations to do for each algorithm
num_trials = 10

# DV (outer) code parameters
distance = 5
boundaries = "periodic"

# Noise model
p_swap = 0.1
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
    RHG_code = SurfaceCode(
        distance=distance, boundaries=boundaries, polarity=alternating_polarity, backend=backend
    )
    RHG_lattice = RHG_code.graph
    for i in range(num_trials):
        print(f"-- {i} --")
        CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
        # Apply noise
        CVRHG.apply_noise(cv_noise)
        # Measure syndrome
        CVRHG.measure_hom("p", RHG_code.all_syndrome_inds)
        # Inner decoder
        before = time.time()
        dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
        dec.correct(
            RHG_code,
            decoder,
            weight_options,
            matching_backend=backend,
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
