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
"""A module with common cluster states."""

import numpy as np
from flamingpy.codes.graphs import EGraph


def star_graph(n):
    """EGraph of a star graph state with n nodes.

    Args:
        n (int): Number of qubits

    Returns:
        EGraph: the star graph state.
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 2, f"Input n should be 2 or larger. Current value is {n}"
    star_graph_state = EGraph()
    for i in range(n - 1):
        degs = 2 * np.pi * i / (n - 1)
        x, y = np.cos(degs), np.sin(degs)
        edge = [(0, 0, 0), (x, y, 0)]
        star_graph_state.add_edge(*edge, color="MidnightBlue")
    return star_graph_state


def complete_graph(n):
    """EGraph of a complete graph state with n nodes.

    Args:
        n (int): Number of qubits

    Returns:
        EGraph: the complete graph state.
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 3, f"Input n should be 3 or larger. Current value is {n}"
    complete_graph = EGraph()
    for i in range(n):
        degs = 2 * np.pi * i / n
        x, y = np.cos(degs), np.sin(degs)
        for j in range(1, n - i):
            degs_adj = 2 * np.pi * j / n + degs
            x_adj, y_adj = np.cos(degs_adj), np.sin(degs_adj)
            edge = [(x, y, 0), (x_adj, y_adj, 0)]
            complete_graph.add_edge(*edge, color="MidnightBlue")
    return complete_graph


def linear_cluster(n):
    """EGraph of a linear cluster state with n nodes.

    Args:
        n (int): Number of qubits

    Returns:
        linear_state (EGraph)
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 2, f"Input n should be 2 or larger. Current value is {n}"
    linear_state = EGraph()
    for i in range(n - 1):
        edge = [(i, 0, 0), (i + 1, 0, 0)]
        linear_state.add_edge(*edge, color="MidnightBlue")
    return linear_state


def ring_graph(n):
    """EGraph of a ring graph state with n nodes.

    Args:
        n (int): Number of qubits

    Returns:
        ring_graph_state (EGraph)
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 3, f"Input n should be 3 or larger. Current value is {n}"
    ring_graph_state = EGraph()
    for i in range(n):
        degs = 2 * np.pi * i / (n)
        degs_next = 2 * np.pi * (i + 1) / (n)
        x, y = np.cos(degs), np.sin(degs)
        x_next, y_next = np.cos(degs_next), np.sin(degs_next)
        edge = [(x, y, 0), (x_next, y_next, 0)]
        ring_graph_state.add_edge(*edge, color="MidnightBlue")
    return ring_graph_state


def bell():
    """EGraph of a two-qubit bell state.

    Returns:
        linear(2) (EGraph)
    """
    return linear(2)
