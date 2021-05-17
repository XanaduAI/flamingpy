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
""""Unit tests for decoder funcions in decoder.py."""

import numpy as np
from graphstates import CVGraph
from decoder import (assign_weights,
                     CV_decoder,
                     decoding_graph,
                     matching_graph,
                     MWPM,
                     recovery,
                     check_correction,
                     correct)
from RHG import RHGCode

distance = 2
boundaries = "periodic"
RHG_code = RHGCode(distance=distance, boundaries=boundaries, polarity=True)
RHG_lattice = RHG_code.graph
# CV (inner) code/state
p_swap = 0
CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)

# Noise model
delta = 0.01
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Decoding options
decoder = {"inner": "basic", "outer": "MWPM"}
# Apply noise
CVRHG.apply_noise(cv_noise)
# Measure syndrome
CVRHG.measure_hom("p", RHG_code.syndrome_inds)

RHG_lattice.draw()


class TestAssignWeights:
    "Test the weight assignment in decoder.py."
    weight_options = {
        "method": "blueprint",
        "integer": False,
        "multiplier": 1,
        "delta": delta,
        }
    assign_weights(RHG_code, **weight_options)
    RHG_lattice.draw(label='weight')

    def test_unit_weights(self):
        assign_weights(RHG_code, method="unit")
        for point in RHG_code.syndrome_coords:
            assert RHG_lattice.nodes[point]['weight'] == 1

    def test_blueprint_weights(self):
        weight_options = {
            "method": "blueprint",
            "integer": True,
            "multiplier": 100,
            "delta": delta,
            }
        assign_weights(RHG_code, **weight_options)
        for point in RHG_code.syndrome_coords:
            weight = RHG_lattice.nodes[point]['weight']
            assert isinstance(weight, int)
            assert weight >= 0


class TestCVDecover:
    pass


class TestDecodingGraph:
    pass


class TestMatchingGraph:
    pass


class TestMWPM:
    pass


class TestRecovery:
    pass


class TestCorrectionCheck:
    pass


class TestCorrect:
    pass
