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
"""Unit tests for iid noise generation and decoding."""


from flamingpy.codes import SurfaceCode
from flamingpy.noise.iid import IidNoise
from flamingpy.decoders.decoder import correct


def test_zero_noise():
    """Check that when the error probability is 0,
    then all bit values are 0.
    """
    code = SurfaceCode(3)
    noise = IidNoise(code, 0.0)
    noise.apply_noise()
    for _, node_data in code.graph.nodes(data=True):
        assert node_data["bit_val"] == 0


def test_full_noise():
    """Check that when the error probability is 1
    then all bit values are 1.
    """
    code = SurfaceCode(3)
    noise = IidNoise(code, 1.0)
    noise.apply_noise()
    for _, node_data in code.graph.nodes(data=True):
        assert node_data["bit_val"] == 1


def test_finite_prob_noise():
    """Check that when the error probability is 0.1
    then all bit values are either 0 or 1.
    """
    code = SurfaceCode(3)
    noise = IidNoise(code, 0.1)
    noise.apply_noise()
    for _, node_data in code.graph.nodes(data=True):
        assert node_data["bit_val"] in [0, 1]


def test_decoding(request):
    """Check that we can use the correct function to decode the code
    after applying iid noise. """
    code = SurfaceCode(3)
    noise = IidNoise(code, 0.1)
    noise.apply_noise()
    assert correct(code, {"outer": "MWPM"},
                   weight_options="uniform") in [True, False]
