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
"""Check that we can run all the example and benchmark files without showing
the plots."""

# pylint: disable=import-outside-toplevel,unused-import

import pytest
from flamingpy.codes import alternating_polarity


@pytest.mark.parametrize("noise", ["cv", "dv"])
@pytest.mark.parametrize("decoder", ["MWPM", "UF"])
def test_decoder_example(noise, decoder):
    """Simple test for the decoding module in flamingpy.examples."""
    from flamingpy.examples.decoding import decode_surface_code

    distance = 3
    boundaries = "open"
    ec = "primal"

    result = decode_surface_code(distance, boundaries, ec, noise, decoder, draw=True)
    assert result.__class__.__name__ == "bool_"


def test_gkp_example():
    """Simple test for the gkp module in flamingpy.examples."""
    from flamingpy.examples import gkp


def test_graphstates_example():
    """Simple test for the graphstates module in flamingpy.examples."""
    from flamingpy.examples import graphstates


def test_macro_reduce_example():
    """Simple test for the macro_reduce module in flamingpy.examples."""
    from flamingpy.examples import macro_reduce


@pytest.mark.parametrize("boundaries", ["periodic", "open"])
def test_surface_code_example(boundaries):
    """Simple test for the surface_code module in flamingpy.examples."""
    from flamingpy.examples.surface_code import illustrate_surface_code

    d = 2
    err = "primal"
    polarity = None
    illustrate_surface_code(d, boundaries, err, polarity, show=False)


def test_lemon_benchmark():
    """Simple test for the lemon module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import lemon


def test_matching_benchmark():
    """Simple test for the matching module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import matching


def test_shortest_path_benchmark():
    """Simple test for the shortest_path module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import shortest_path
