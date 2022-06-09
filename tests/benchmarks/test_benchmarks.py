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


def test_decoding_benchmark():
    """Simple test for the decoding module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import decoding as dc_benchmarks


def test_lemon_benchmark():
    """Simple test for the lemon module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import lemon


def test_matching_benchmark():
    """Simple test for the matching module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import matching


def test_shortest_path_benchmark():
    """Simple test for the shortest_path module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import shortest_path
