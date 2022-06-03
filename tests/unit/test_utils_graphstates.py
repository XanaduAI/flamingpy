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
import pytest
from flamingpy.utils import graph_states


def test_star_graph():
    """Check that star state with n qubits has n-1 edges"""
    n = 10
    star_state = graph_states.star_graph(n)
    assert star_state.number_of_edges() == (n - 1)

    n = 5
    star_state = graph_states.star_graph(n)
    assert star_state.number_of_edges() == (n - 1)


def test_star_graph_inputs():
    """Checks that invalid inputs in star_graph raise a warning"""
    n = 1
    with pytest.raises(Exception) as e:
        graph_states.star_graph(n)
    assert e.type == ValueError
    assert "Input n should be 2 or larger." in str(e.value)

    n = "hello"
    with pytest.raises(Exception) as e:
        graph_states.star_graph(n)
    assert e.type == ValueError
    assert "Input n should be an integer." in str(e.value)


def test_complete_graph():
    """Check that complete graph state with n qubits has n(n-1)/2 edges"""
    n = 15
    complete_graph = graph_states.complete_graph(n)
    assert complete_graph.number_of_edges() == (n * (n - 1) / 2)
    assert complete_graph.number_of_nodes() == n

    n = 5
    complete_graph = graph_states.complete_graph(n)
    assert complete_graph.number_of_edges() == (n * (n - 1) / 2)
    assert complete_graph.number_of_nodes() == n


def test_complete_graph_inputs():
    """Checks that invalid inputs in complete_graph raise a warning"""
    n = 2
    with pytest.raises(Exception) as e:
        graph_states.complete_graph(n)
    assert e.type == ValueError
    assert "Input n should be 3 or larger." in str(e.value)

    n = "hello"
    with pytest.raises(Exception) as e:
        graph_states.complete_graph(n)
    assert e.type == ValueError
    assert "Input n should be an integer." in str(e.value)


def test_ring_state():
    """Check that a ring state with n qubits has n edges."""
    n = 15
    ring_graph = graph_states.ring_graph(n)
    assert ring_graph.number_of_edges() == n
    assert ring_graph.number_of_nodes() == n

    n = 5
    ring_graph = graph_states.ring_graph(n)
    assert ring_graph.number_of_edges() == n
    assert ring_graph.number_of_nodes() == n


def test_ring_state_inputs():
    """Checks that invalid inputs in ring_state raise a warning"""
    n = 2
    with pytest.raises(Exception) as e:
        graph_states.ring_graph(n)
    assert e.type == ValueError
    assert "Input n should be 3 or larger." in str(e.value)

    n = "hello"
    with pytest.raises(Exception) as e:
        graph_states.ring_graph(n)
    assert e.type == ValueError
    assert "Input n should be an integer." in str(e.value)


def test_linear_state():
    """Check that a linear state with n qubits has n-1 edges"""
    n = 15
    linear_cluster = graph_states.linear_cluster(n)
    assert linear_cluster.number_of_edges() == n - 1
    assert linear_cluster.number_of_nodes() == n

    n = 5
    linear_cluster = graph_states.linear_cluster(n)
    assert linear_cluster.number_of_edges() == n - 1
    assert linear_cluster.number_of_nodes() == n


def test_linear_state_inputs():
    """Checks that invalid inputs in linear_state raise a warning"""
    n = 1
    with pytest.raises(Exception) as e:
        graph_states.linear_cluster(n)
    assert e.type == ValueError
    assert "Input n should be 2 or larger." in str(e.value)

    n = "hello"
    with pytest.raises(Exception) as e:
        graph_states.linear_cluster(n)
    assert e.type == ValueError
    assert "Input n should be an integer." in str(e.value)


def test_bell_state():
    """Check that the bell state generates"""
    bell = graph_states.bell()
    assert bell.number_of_edges() == 1
    assert bell.number_of_nodes() == 2
