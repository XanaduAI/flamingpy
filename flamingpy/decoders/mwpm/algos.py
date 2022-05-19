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
"""The minimum-weight perfect matching decoding algorithms."""

# pylint: disable=import-outside-toplevel

import itertools as it

from flamingpy.decoders.mwpm.matching import NxMatchingGraph, RxMatchingGraph, LemonMatchingGraph


def build_match_graph(code, ec, matching_backend="retworkx"):
    """Build the matching graph for the given code.

    Args:
        code (SurfaceCode): the code class to decode
        ec (string): the error complex ("primal" or "dual")
        matching_backend (str or flamingpy.matching.MatchingGraph, optional):
            The type of matching graph to build. If providing a string,
            it must be either "networkx", "retworkx" or "lemon" to pick one
            of the already implemented backends. Else, the provided type should
            inherit from the MatchingGraph abstract base class and have an empty init.
            The default is retworkx.

    Returns:
        MatchingGraph: The matching graph.
    """
    default_backends = {
        "networkx": NxMatchingGraph,
        "retworkx": RxMatchingGraph,
        "lemon": LemonMatchingGraph,
    }
    if matching_backend in default_backends:
        matching_backend = default_backends[matching_backend]

    return matching_backend(ec, code)


def mwpm_decoder(code, ec, backend="retworkx", draw=False, drawing_opts=None):
    """Run the minimum-weight perfect matching decoder on code.

    Args:
        code (SurfaceCode): the code class to decode and correct
        ec (string): the error complex ("primal" or "dual")
        backend (str or flamingpy.matching.MatchingGraph, optional):
            The type of matching graph to build. If providing a string,
            it must be either "networkx", "retworkx" or "lemon" to pick one
            of the already implemented backends. Else, the provided type should
            inherit from the MatchingGraph abstract base class and have an empty init.
            The default is the retworkx.
        draw (bool): whether or not to draw the MWPM decoding process:
                stabilizer graph, matching graph, syndrome plot, and
                matching.
        drawing_opts (dict): the drawing options, to be fed into
            viz.draw_mwpm_decoding (see that function for more details).

    Returns:
        set[tuples]: the nodes (representing qubits) be fed into the recovery
            (i.e. whose bit values must be flipped).
    """
    stab_graph = getattr(code, ec + "_stab_graph")
    matching_graph = build_match_graph(code, ec, backend)
    matching = matching_graph.min_weight_perfect_matching()
    virtual_points = matching_graph.virtual_points

    # Draw the stabilizer graph, matching graph, and syndrome, if
    # desired.
    if draw:
        from flamingpy.utils.viz import draw_decoding

        dec_objects = {"matching_graph": matching_graph, "matching": matching}
        draw_decoding(code, ec, dec_objects, drawing_opts)

    qubits_to_flip = set()
    for match in matching:
        if match not in it.product(virtual_points, virtual_points):
            path = matching_graph.edge_path(match)
            pairs = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for pair in pairs:
                common_vertex = stab_graph.edge_data(*pair)["common_vertex"]
                qubits_to_flip ^= {common_vertex}
    return qubits_to_flip
