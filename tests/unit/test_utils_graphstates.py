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
""""Unit tests for the utils.graph_states module"""

from flamingpy.utils import graph_states
import pytest


def test_star_graph():
    """Check that star state with n qubits has n-1 edges"""
    n = 10
    star_state = graph_states.star(n)
    assert star_state.number_of_edges() == (n - 1)

    n = 5
    star_state = graph_states.star(n)
    assert star_state.number_of_edges() == (n - 1)


def test_ghz_state():
    """Check that ghz state with n qubits has n(n-1)/2 edges"""
    n = 10
    ghz = graph_states.ghz(n)
    assert ghz.number_of_edges() == (n * (n - 1) / 2)

    n = 5
    ghz = graph_states.ghz(n)
    assert ghz.number_of_edges() == (n * (n - 1) / 2)


def test_ring_state():
    """Check that a ring state with n qubits has n edges"""
    n = 10
    ring = graph_states.ring(n)
    assert ring.number_of_edges() == n

    n = 5
    ring = graph_states.ring(n)
    assert ring.number_of_edges() == n


def test_linear_state():
    """Check that a linear state with n qubits has n-1 edges"""
    n = 10
    linear = graph_states.linear(n)
    assert linear.number_of_edges() == n - 1

    n = 5
    linear = graph_states.linear(n)
    assert linear.number_of_edges() == n - 1


def test_bell_state():
    """Check that the bell state generates"""
    bell = graph_states.bell()
    assert bell.number_of_edges() == 1
