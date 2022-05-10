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
""""Unit tests for MatchingGraph classes in matching.py.

The networkx implementation is used as a reference.
"""
import itertools as it
from copy import deepcopy

import numpy as np
import pytest
import networkx as nx

from flamingpy.decoders.decoder import correct, assign_weights
from flamingpy.codes.surface_code import SurfaceCode
from flamingpy.noise.iid import IidNoise


@pytest.fixture(scope="module", params=it.product(["networkx", "lemon", "retworkx"]))
def backends(request):
    return request.param

def test_compare_success_to_nx(backends):
    code = SurfaceCode(distance=3, boundaries="open")    
    noise = IidNoise(code, error_probability=0.1)

    rng = np.random.default_rng(123)
    nx_rng = np.random.default_rng(123)

    for _ in range(5):
        noise.apply_noise(rng)
        assign_weights(code, method="unit")
        result, matching = correct(code, {"outer": "MWPM"}, matching_backend=backends[0])

        noise.apply_noise(nx_rng)
        assign_weights(code, method="unit")
        nx_result, nx_matching = correct(code, {"outer": "MWPM"}, matching_backend="networkx")

        # assert matching == nx_matching
        assert result == nx_result
