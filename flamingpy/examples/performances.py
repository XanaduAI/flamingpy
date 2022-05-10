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
"""Example of instantiating, applying noise, decoding, recovering, and
visualizing this procedure for the measurement-based surface code."""

# pylint: disable=no-member

import matplotlib.pyplot as plt
from numpy.random import default_rng

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders import decoder as dec
from flamingpy.noise import IidNoise
from flamingpy.utils import viz

# QEC code parameters
distance = 5
# Boundaries ("open" or "periodic")
boundaries = "open"
# Error complex ("primal", "dual", or "both")
ec = "primal"

# Code and code lattice (cluster state)
code = SurfaceCode(
    distance=distance,
    ec=ec,
    boundaries=boundaries,
    polarity=alternating_polarity,
    backend="retworkx",
)

# CV (inner) code / state preparation
p_swap = 0.05  # probability of having squeezed states (the rest are GKPs)
CVRHG = CVLayer(code.graph, p_swap=p_swap)
# Noise model
delta = 0.1  # GKP squeezing parameter
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Decoding options
weight_options = {
    "method": "blueprint",
    "integer": True,
    "multiplier": 100,
    "delta": delta,
}

decoder = {"inner": "basic", "outer": "MWPM"}

print("Number of successes for 50 samples")
for backend in ["lemon", "retworkx"]:
    num_successes = 0
    rng = default_rng(123)
    for i in range(50):
        CVRHG.apply_noise(cv_noise, rng)
        CVRHG.measure_hom("p", code.all_syndrome_inds, rng)
        dec.CV_decoder(code, translator=dec.GKP_binner)
        dec.assign_weights(code, **weight_options)
        # Automatic decoding
        num_successes += dec.correct(
            code=code,
            decoder=decoder,
            weight_options=weight_options,
            sanity_check=False,
            matching_backend=backend,
        )
    print(f"{backend}: {num_successes}")
