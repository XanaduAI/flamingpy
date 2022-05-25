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
"""Cython-based Monte Carlo simulations for estimating FT thresholds."""
from flamingpy.decoders.decoder import correct
from flamingpy.noise import CVLayer, CVMacroLayer

cpdef int cpp_mc_loop(object coder, int trials, dict decoder, dict weight_options, object macro_graph, double p_swap, double delta, dict cv_noise):
    """Exclusively run loop section of Monte Carlo simulations of error-correction on code=code."""
    noise_model = {"noise": "grn", "delta": delta}
    cdef int successes = 0
    cdef int result, errors, ii
    for ii in range(trials):
        if macro_graph:
            CV_macro = CVMacroLayer(macro_graph, p_swap=p_swap, reduced_graph=code.graph)
            CV_macro.reduce(noise_model)
        else:
            # Apply noise
            CVRHG = CVLayer(code, p_swap=p_swap)
            CVRHG.apply_noise(cv_noise)
            # Measure syndrome
            CVRHG.measure_hom("p", code.syndrome_inds)
        result = correct(code=code, decoder=decoder, weight_options=weight_options)
        successes += result
    errors = trials - successes
    return errors
