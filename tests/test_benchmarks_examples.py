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
""" Check that we can run all the example and benchmark files
    without showing the plots. """


def test_decoder_example():
    from flamingpy.examples import decoding


def test_gkp_example():
    from flamingpy.examples import gkp


def test_graphstates_example():
    from flamingpy.examples import graphstates


def test_macro_reduce_example():
    from flamingpy.examples import macro_reduce


def test_surface_code_example():
    from flamingpy.examples import surface_code


def test_decoding_benchmark():
    from flamingpy.benchmarks import decoding


def test_lemon_benchmark():
    from flamingpy.benchmarks import lemon


def test_matching_benchmark():
    from flamingpy.benchmarks import matching


def test_shortest_path_benchmark():
    from flamingpy.benchmarks import shortest_path


def test_simulations_benchmark():
    from flamingpy.benchmarks import simulations
