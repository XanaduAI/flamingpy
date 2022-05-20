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
from flamingpy.cv.macro_reduce import reduce_macronode_graph
from flamingpy.noise import CVLayer

cpdef int cpp_mc_loop(object code, int trials, dict decoder, dict weight_options, object passive_objects, double p_swap, double delta, dict cv_noise):
    """Exclusively run loop section of Monte Carlo simulations of error-correction on code=code."""
    cdef int successes = 0
    cdef int result, errors, ii
    for ii in range(trials):
        if passive_objects:
            reduce_macronode_graph(*passive_objects, p_swap, delta)
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
