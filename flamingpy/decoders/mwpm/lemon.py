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
"""Define max_weight_matching based on the lemonpy module."""
import warnings

import networkx as nx

try:
    import flamingpy.cpp.lemonpy as lp
except ImportError:  # pragma: no cover
    warnings.warn("Failed to import flamingpy.cpp.lemonpy library.", ImportWarning)


def min_weight_matching(G_match):
    """Compute the minimum-weight perfect matching using lemon.

    Assumptions:
        1. Symmetric adjacency matrix.
        2. Adjacency matrix has zeros along diagonal.
    """
    adjacency = nx.to_numpy_array(G_match, weight="inverse_weight", nonedge=-(2**64) + 1)
    lemon_matching = lp.mwpm(adjacency)

    # lemon uses different node ids, so we convert back to networkx node ids
    nx_map = list(G_match.nodes())
    matching = set()
    for i in lemon_matching:
        matching.add((nx_map[i[0]], nx_map[i[1]]))
    return matching
