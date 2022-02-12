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
"""Decoders for measurement-based codes."""
import sys
import itertools as it
import numpy as np

from flamingpy.matching import NxMatchingGraph, RxMatchingGraph, LemonMatchingGraph
from flamingpy.GKP import GKP_binner, Z_err_cond

# Smallest and largest numbers representable.
smallest_number = sys.float_info.min
largest_number = sys.float_info.max


def assign_weights(code, **kwargs):
    """Assign weights to qubits in a hybrid CV graph state CVG.

    Args:
        code (code class): the qubit code
        state (CVGraph): the CVGraph whose syndrome has been measured
        method (str, optional): the method for weight assignment. By
            default, 'unit', denoting weight 1 everyoewhere. For
            heuristic and analog weight assignment from blueprint, use
            'blueprint'
        integer (bool, optional): whether to convert weights to
            integers using Python's round function; False by default
        multiplier (int, optional): multiply the weight by multiplier
            before rounding; 1 by default.

    Returns:
        None
    """
    default_options = {"method": "unit", "integer": False, "multiplier": 1}
    weight_options = {**default_options, **kwargs}
    G = code.graph
    # Get and denest the syndrome coordinates.
    syndrome_coords = code.syndrome_coords
    # Blueprint weight assignment dependent on type of neighbour.
    if weight_options.get("method") == "blueprint":
        for node in syndrome_coords:
            neighbors = G[node]
            # List and number of p-squeezed states in neighborhood of node.
            p_list = [G.nodes[v]["state"] for v in neighbors if G.nodes[v]["state"] == "p"]
            p_count = len(p_list)
            if p_count in (0, 1):
                if weight_options.get("prob_precomputed"):
                    err_prob = G.nodes[node]["p_phase_cond"]
                else:
                    delta_effective = (len(neighbors) + 1) * weight_options.get("delta")
                    hom_val = G.nodes[node]["hom_val_p"]
                    err_prob = Z_err_cond(delta_effective, hom_val)
                # Allow for taking log of 0.
                # TODO: Is this the best way to do it? Or can I just choose
                # an arbitrary small number?
                if err_prob > 0.5:
                    err_prob = 0.5
                if err_prob == 0:
                    err_prob = smallest_number
                if weight_options.get("integer"):
                    multiplier = weight_options.get("multiplier")
                    weight = round(-multiplier * np.log(err_prob))
                else:
                    weight = -np.log(err_prob)
                G.nodes[node]["weight"] = weight
            else:
                # Dictionary of the form number of swapouts: error probability.
                weight_dict = {2: 1 / 4, 3: 1 / 3, 4: 2 / 5}
                err_prob = weight_dict[p_count]
                if weight_options.get("integer"):
                    multiplier = weight_options.get("multiplier")
                    weight = round(-multiplier * np.log(err_prob))
                else:
                    weight = -np.log(err_prob)
                G.nodes[node]["weight"] = weight
        return
    # Naive weight assignment, unity weights.
    if weight_options.get("method") == "unit":
        for node in syndrome_coords:
            G.nodes[node]["weight"] = 1


# TODO: General functions for applying noise and measuring syndrome.


def CV_decoder(code, translator=GKP_binner):
    """Convert homodyne outcomes to bit values according to translate.

    The inner (CV) decoder, aka translator, aka binning function. Set
    converted values to the bit_val attribute for nodes in G.

    Args:
        state: the CVGraph with homodyne outcomes computed.
        translator: the choice of binning function; by default, the
            standard GKP binning function that snaps to the closest
            integer multiple of sqrt(pi).
    Returns:
        None
    """
    syndrome = code.syndrome_coords
    # TODO: Generalize to nonlocal translators
    # TODO: Vectorize?
    for point in syndrome:
        hom_val = code.graph.nodes[point]["hom_val_p"]
        bit_val = translator([hom_val])[0]
        code.graph.nodes[point]["bit_val"] = bit_val


def recovery(code, G_match, matching, sanity_check=False):
    """Run the recovery operation on graph G.

    Fip the bit values of all the vertices in the path connecting each
    pair of stabilizers according to the matching. If check, verify
    that there are no odd-parity cubes remaining, or display their
    indices of there are.

    Args:
        code (RHGCode): the code
        G_match (networkx.Graph): the matching graph
        matching (set of tuples): the minimum-weight perfect matching
        check (bool): if True, check if the recovery has succeeded.

    Returns:
        None or bool: if check, False if the recovery failed, True if
            it succeeded. Otherwise None.
    """
    virtual_points = G_match.virtual_points
    for match in matching:
        if match not in it.product(virtual_points, virtual_points):
            path = G_match.edge_path(match)
            pairs = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for pair in pairs:
                common_vertex = code.stab_graph.edge_data(*pair)["common_vertex"]
                code.graph.nodes[common_vertex]["bit_val"] ^= 1

    if sanity_check:
        odd_cubes = list(code.stab_graph.odd_parity_stabilizers())
        if odd_cubes:
            print("Unsatisfied stabilizers:", odd_cubes)
        print("Recovery succeeded - no unsatisfied stabilizers.")


# TODO: Rename to logical_error_check or someting like that. Clarify
# what correlation surface is, per RHG paper.
def check_correction(code, plane=None, sheet=0, sanity_check=False):
    """Perform a correlation-surface check.

    Check the total parity of a correlation surface specified by
    direction plane and index sheet to check if the error correction
    procedure has succeeded. By default, checks the x = 0 plane.

    Args:
        G (CVGraph): the recovered graph
        plane (str): 'x', 'y', or 'z', determining the direction of the
            correlation surface; if None, select plane to check
            depending on the boundary conditions
        sheet (int): the sheet index (from 0 to one less than the
            largest coordinate in the direction plane).
        sanity_check (bool): if True, display the total parity of
            all correlation surfaces in the direction plane to verify
            if parity is conserved.

    Returns:
        None
    """
    dims = np.array(code.dims)
    dir_dict = {"x": 0, "y": 1, "z": 2}
    truth_dict = {"x": [], "y": [], "z": []}
    boundaries = code.boundaries

    planes_to_check = []
    if plane:
        planes_to_check += [plane]
    elif boundaries == ["periodic"] * 3:
        planes_to_check = ["x", "y", "z"]
    elif tuple(boundaries) in set(it.permutations(["primal", "dual", "dual"])):
        where_primal = np.where(np.array(boundaries) == "primal")[0][0]
        planes_to_check = [["x", "y", "z"][where_primal]]

    minimum, maximum = sheet, sheet + 2
    for plane in planes_to_check:
        if sanity_check:
            minimum, maximum = 0, 2 * dims[dir_dict[plane]]
        for sheet in range(minimum, maximum, 2):
            slice_verts = code.graph.slice_coords(plane, sheet)
            syndrome_verts = code.syndrome_coords
            only_primal = set(slice_verts) & set(syndrome_verts)
            parity = 0
            for node in only_primal:
                parity ^= code.graph.nodes[node]["bit_val"]
            truth_dict[plane].append(bool(1 - parity))

    if sanity_check:
        print(truth_dict)

    all_surfaces = np.array([truth_dict[i][0] for i in planes_to_check])
    return np.all(all_surfaces)


def build_match_graph(code, weight_options, matching_backend="networkx"):
    """
    Build the matching graph for the given code.

    Combines weight assignment and matching graph creation.

    Args:
        code (code): the code class to decode and correct
        weight_options (dict): how to assign weights; options are
                'method': 'unit' or 'blueprint'
                'integer': True (for rounding) or False (for not)
                'multiplier': integer denoting multiplicative factor
                    before rounding
        matching_backend (str or flamingpy.matching.MatchingGraph, optional):
            The type of matching graph to build. If providing a string,
            it most be either "networkx", "retworkx" or "lemon" to pick one
            of the already implemented backends. Else, the provided type should
            inherit from the MatchingGraph abstract base class and have an empty init.
            The default is the networkx backend since it is the reference implementation.
            However, both retworkx and lemon and orders of magnitude faster.
    Returns:
        MatchingGraph: The matching graph.
    """
    if weight_options is None:
        weight_options = {}
    assign_weights(code, **weight_options)
    default_backends = {
        "networkx": NxMatchingGraph,
        "retworkx": RxMatchingGraph,
        "lemon": LemonMatchingGraph,
    }
    if matching_backend in default_backends:
        matching_backend = default_backends[matching_backend]
    return matching_backend(code)


def correct(
    code,
    decoder,
    weight_options=None,
    sanity_check=False,
    matching_backend="networkx",
):
    """Run through all the error-correction steps.

    Combines weight assignment matching graph creation,
    minimum-weight-perfect matching, recovery, and correctness check.

    Args:
        code (code): the code class to decode and correct
        inner_decoder (str): the CV decoder; GKP_binner by default
        outer_decoder (str): the DV decoder; MWPM by default
        weight_options (dict): how to assign weights; options are
            'method': 'unit' or 'blueprint'
            'integer': True (for rounding) or False (for not)
            'multiplier': integer denoting multiplicative factor
            before rounding
        sanity_check (bool): if True, check that the recovery
            operation succeeded and verify that parity is conserved
            among all correlation surfaces.
        matching_backend (str or flamingpy.matching.MatchingGraph, optional):
            The backend to generate the matching graph. See build_dec_and_match_graphs
            for more details.
    Returns:
        bool: True if error correction succeeded, False if not.
    """

    inner_dict = {"basic": GKP_binner}
    outer_dict = {"MWPM": "MWPM"}

    inner_decoder = decoder.get("inner")
    outer_decoder = decoder.get("outer")
    if inner_decoder:
        CV_decoder(code, translator=inner_dict[inner_decoder])
    if outer_dict[outer_decoder] == "MWPM":
        matching_graph = build_match_graph(code, weight_options, matching_backend)
        matching = matching_graph.min_weight_perfect_matching()
        recovery(code, matching_graph, matching, sanity_check=sanity_check)
    result = check_correction(code, sanity_check=sanity_check)
    return result
