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


def star_graph_state(n):
    """
    Returns the EGraph of a star graph state with n nodes.
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


def ghz_state(n):
    """
    Returns the EGraph of a GHZ state with n nodes.
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 3, f"Input n should be 3 or larger. Current value is {n}"
    ghz_state = EGraph()
    for i in range(n):
        degs = 2 * np.pi * i / n
        x, y = np.cos(degs), np.sin(degs)
        for j in range(1, n):
            degs_adj = 2 * np.pi * j / n + degs
            x_adj, y_adj = np.cos(degs_adj), np.sin(degs_adj)
            edge = [(x, y, 0), (x_adj, y_adj, 0)]
            ghz_state.add_edge(*edge, color="MidnightBlue")
    return ghz_state


def linear_state(n):
    """
    Returns the EGraph of a linear cluster state with n nodes.
    """
    assert isinstance(n, int), f"Input n should be an integer. Current type is {type(n)}"
    assert n >= 2, f"Input n should be 2 or larger. Current value is {n}"
    linear_state = EGraph()
    for i in range(n - 1):
        edge = [(i, 0, 0), (i + 1, 0, 0)]
        linear_state.add_edge(*edge, color="MidnightBlue")
    return linear_state


def ring_state(n):
    """
    Returns the EGraph of a ring graph state with n nodes.
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


def bell_state(n=2):
    """
    Returns the EGraph of a two-qubit bell state.
    """
    assert (
        n == 2
    ), "Bell states are defined for two qubits. For larger n, you might want to consider GHZ states"
    return linear_state(n)
