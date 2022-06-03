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


def _generic_graph_state(n):
    """Return an EGraph when n is an integer and n>=1."""
    if not isinstance(n, int):
        raise ValueError(f"Input n should be an integer. Current type is {type(n)}")
    if n < 1:
        raise ValueError(f"Input n should be 1 or larger. Current value is {n}")
    return EGraph()


def star_graph(n):
    """Return an EGraph of a star graph state with n nodes."""
    star_graph_state = _generic_graph_state(n)
    if n == 1:
        star_graph_state.add_node((0, 0, 0))
    else:
        for i in range(n - 1):
            degs = 2 * np.pi * i / (n - 1)
            x, y = np.cos(degs), np.sin(degs)
            edge = [(0, 0, 0), (x, y, 0)]
            star_graph_state.add_edge(*edge)
    return star_graph_state


def complete_graph(n):
    """Return an EGraph of a complete graph state with n nodes."""
    complete_graph = _generic_graph_state(n)
    if n == 1:
        complete_graph.add_node((0, 0, 0))
    else:
        for i in range(n):
            degs = 2 * np.pi * i / n
            x, y = np.cos(degs), np.sin(degs)
            for j in range(1, n - i):
                degs_adj = 2 * np.pi * j / n + degs
                x_adj, y_adj = np.cos(degs_adj), np.sin(degs_adj)
                edge = [
                    (float(round(x, 5)), float(round(y, 5)), 0),
                    (float(round(x_adj, 5)), float(round(y_adj, 5)), 0),
                ]
                complete_graph.add_edge(*edge)
    return complete_graph


def linear_cluster(n):
    """Return an EGraph of a linear cluster state with n nodes."""
    linear_state = _generic_graph_state(n)
    if n == 1:
        linear_state.add_node((0, 0, 0))
    else:
        for i in range(n - 1):
            edge = [(i, 0, 0), (i + 1, 0, 0)]
            linear_state.add_edge(*edge)
    return linear_state


def ring_graph(n):
    """Return an EGraph of a ring graph state with n nodes."""
    ring_graph_state = _generic_graph_state(n)
    if n == 1:
        ring_graph_state.add_node((0, 0, 0))
    else:
        x, y = np.cos(0), np.sin(0)
        for i in range(1, n + 1):
            degs = 2 * np.pi * i / (n)
            x_next, y_next = np.cos(degs), np.sin(degs)
            edge = [
                (float(round(x, 5)), float(round(y, 5)), 0),
                (float(round(x_next, 5)), float(round(y_next, 5)), 0),
            ]
            ring_graph_state.add_edge(*edge)
            x, y = x_next, y_next
    return ring_graph_state


def bell():
    """Return an EGraph of the two-qubit Bell pair."""
    return linear_cluster(2)
