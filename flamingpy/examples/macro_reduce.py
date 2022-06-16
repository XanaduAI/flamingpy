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
from flamingpy.decoders.decoder import correct
from flamingpy.noise.cv import CVMacroLayer

# Number of trials
total = 10
# Code parameters and definition
d = 3
boundaries = "open"
ec = "primal"
RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)

# Noise parameters and layer definition
delta = 0.1
p_swap = 0.25
CV_macro = CVMacroLayer(RHG_code, delta=delta, p_swap=p_swap)

successes = 0
for trial in range(total):
    # The CV macronode noise layer and reduction
    CV_macro.apply_noise()
    decoder = "MWPM"
    decoder_opts = {"backend": "networkx"}
    if decoder == "MWPM":
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
