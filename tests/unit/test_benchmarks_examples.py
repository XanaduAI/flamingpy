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

# pylint: disable=import-outside-toplevel


def test_decoder_example():
    """Simple test for the decoding module in flamingpy.examples."""
    from flamingpy.examples import decoding as dc_examples


def test_decoding_benchmark():
    """Simple test for the decoding module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import decoding as dc_benchmarks


def test_gkp_example():
    """Simple test for the gkp module in flamingpy.examples."""
    from flamingpy.examples import gkp


def test_graphstates_example():
    """Simple test for the graphstates module in flamingpy.examples."""
    from flamingpy.examples import graphstates


def test_macro_reduce_example():
    """Simple test for the macro_reduce module in flamingpy.examples."""
    from flamingpy.examples import macro_reduce


def test_surface_code_example():
    """Simple test for the surface_code module in flamingpy.examples."""
    from flamingpy.examples import surface_code


def test_performance_example():
    """Simple test for the decoding module in flamingpy.examples."""
    from flamingpy.examples import performances


def test_lemon_benchmark():
    """Simple test for the lemon module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import lemon


def test_matching_benchmark():
    """Simple test for the matching module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import matching


def test_shortest_path_benchmark():
    """Simple test for the shortest_path module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import shortest_path


def test_simulations_benchmark():
    """Simple test for the simulations module in flamingpy.benchmarks."""
    from flamingpy.benchmarks import simulations
