# Copyright 2020 Xanadu Quantum Technologies Inc.

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
import networkx as nx
import networkx.algorithms.shortest_paths as sp

import ft_stack.lemon as lemon
from ft_stack.GKP import GKP_binner, Z_err_cond

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
            p_list = [
                G.nodes[v]["state"] for v in neighbors if G.nodes[v]["state"] == "p"
            ]
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


def decoding_graph(code):
    """Populate the edge weights of the decoding graph from the RHG lattice G,
    and determine the parity of the stabilizer cubes.

    The decoding graph has as its nodes every stabilizer in G and a
    every boundary point (for now coming uniquely from a primal
    boundary). Two stabilizers (or a stabilizer and a boundary point)
    sharing a vertex are connected by an edge whose weight is equal to
    the weight assigned to that vertex in G. The output graph has
    stabilizer nodes relabelled to integer indices, but still points
    to the original stabilizer with the help of the 'stabilizer'
    attribute. The output graph furthermore stores the indices of
    odd-parity cubes (under an 'odd_cubes' graph attribute) and of
    the boundary points (under 'boundary_points'). Common vertices
    are stored under the 'common_vertex' edge attribute. Note that one
    ought first compute the phase error probabilities, conduct a
    homodyne measurement, and translate the outcomes on G.

    Args:
        G (CVGraph): the CVGraph to decode
    Returns:
        networkx.Graph: the decoding graph.
    """
    cubes = code.stabilizers

    decoding_graph = code.decoding_graph
    mapping = code._decoder_mapping

    # Assign edge weights based on the common vertices between stabilizers,
    # and not to edges with the 'high' and 'low points
    real_points = decoding_graph.graph["real_points"]
    for edge in decoding_graph.subgraph(real_points).edges:
        common = decoding_graph.edges[edge]["common_vertex"]
        decoding_graph.edges[edge]["weight"] = code.graph.nodes[common]["weight"]

    # Indices of odd parity cubes and boundary vertices; add these
    # index lists to the graph attributes dictionary, to be used by the
    # matching graph.
    odd_parity_cubes = [cube for cube in cubes if cube.parity]
    odd_parity_inds = [mapping[cube] for cube in odd_parity_cubes]
    decoding_graph.graph["odd_cubes"] = odd_parity_inds[:]

    return decoding_graph


def matching_graph(G_dec, alg="dijkstra"):
    """Create a matching graph from the decoding graph G.

    Generate a matching graph from the decoding graph G according to
    algorithm alg. By default, this is the NetworkX Dijkstra shortest-
    path algorithm. This graph will be fed into a subsequent minimum-
    weight-perfect-matching algorithm. The matching graph has as half
    of its nodes the odd-parity stabilizers. The edge connecting two
    nodes corresponds to the weight of the minimum-weight-path between
    the nodes in the decoding graph. Additionally, each unsatisfied
    stabilizer is connected to a unique boundary point (for now from
    a primal bundary) located at the shortest weighted distance from
    the stabilizer. Between each other, the boundary points are
    connected by an edge of weight 0. The output graph stores the
    indices of the used boundary points under the 'used_boundary_point'
    attribute. Paths are stored under the 'paths' attribute of edges,
    and 'inverse_weights' are also stored, for the benefit of maximum-
    weight-matching algorithms

    Args:
        G (networkx.Graph): the decoding graph, storing information
            about indices of odd-parity-cubes (under 'odd_cubes' graph
            attribute) and boundary points (under 'boundary_points').
        alg (str): the algorithm for shortest-path finding. By default,
            uses variations of Dijkstra functions from NetworkX
    Returns:
        networkx.Graph: the matching graph.
    """
    # An empty matching graph.
    G_match = nx.Graph(title="Matching Graph")

    # Get the indices of the odd parity cubes from the decoding graph.
    odd_parity_inds = G_dec.graph["odd_cubes"]

    # Give shorter names to the Dijkstra shortest path algorithms.
    if alg == "dijkstra":
        alg = sp.single_source_dijkstra

    # Run the matching algorithm first without the 'high' and 'low points
    real_points = G_dec.graph["real_points"]

    # Combinations of odd-parity cubes.
    odd_ind_dict = {i: [] for i in odd_parity_inds[:-1]}
    odd_combs = it.combinations(odd_parity_inds, 2)
    for pair in odd_combs:
        odd_ind_dict[pair[0]] += [pair[1]]
    # Find the shortest paths between odd-parity cubes.
    for cube1 in odd_parity_inds[:-1]:
        lengths, paths = alg(G_dec.subgraph(real_points), cube1)
        for cube2 in odd_ind_dict[cube1]:
            length = lengths[cube2]
            path = paths[cube2]
            # Add edge to the matching graph between the cubes, with weight
            # equal to the length of the shortest path.
            # TODO: Is the behavior correct for negative weights, or do I
            # want 1/weight or max_num - weight?
            G_match.add_edge(
                cube1, cube2, weight=length, inverse_weight=-length, path=path
            )

    virtual_points = []
    if G_dec.graph["boundary_points"]:

        i = 0
        low_lengths, low_paths = alg(G_dec, "low")
        high_lengths, high_paths = alg(G_dec, "high")
        for cube in odd_parity_inds:
            distances = (low_lengths[cube], high_lengths[cube])
            where_shortest = np.argmin(distances)
            if where_shortest == 0:
                length = low_lengths[cube]
                full_path = low_paths[cube]
            if where_shortest == 1:
                length = high_lengths[cube]
                full_path = high_paths[cube]
            point = full_path[1]
            virtual_point = (point, i)
            path = full_path[1:]
            # Add edge to the matching graph between the cube and
            # the virtual excitation corresponding to the boundary
            # vertex, with weight equal to the length of the shortest
            # path.
            G_match.add_edge(
                cube, virtual_point, weight=length, inverse_weight=-length, path=path
            )
            i += 1
            virtual_points += [virtual_point]
        # Add edge with weight 0 between any two virtual excitations.
        for (point1, point2) in it.combinations(virtual_points, 2):
            G_match.add_edge(point1, point2, weight=0, inverse_weight=0)

    G_match.graph["virtual_points"] = virtual_points

    return G_match


def MWPM(G_match, G_dec, alg="blossom_nx"):
    """Run minimum-weight-perfect matching on matching graph G_match.

    Run a minimum-weight-perfect-matching (MWPM) algorithm (the
    BlossomV/Edmunds aglorithm as implemented in NetworkX by default)
    on the matching graph G_match. Under the hood, runs maximum weight
    matching with maximum cardinality on inverse weights. The matching
    combines pairs of nodes in the matching graph in a way that
    minimizes the total weight. Perfect matching combine all pairs of
    nodes.

    Args:
        G_match (networkx.Graph): the matching graph, storing inverse
            weights under the 'inverse_weight' edge attribute ad paths
            under 'path'
        G_dec (networkx.Graph): the decoding graph, pointing to the
            stabilizer cubes under the 'stabilizer' node attribute
        alg (string): the matching algorithm; by default, NetworkX
            implementation of BlossomV/Edmunds maximum-weight matching

    Return:
        set of tuples: pairs of all matched nodes.
    """
    if alg == "blossom_nx":
        matching = nx.max_weight_matching(G_match, maxcardinality=True, weight="inverse_weight")
    elif alg == "lemon":
        matching = lemon.max_weight_matching(G_match, weight="inverse_weight")
    return matching


def recovery(code, G_match, G_dec, matching, sanity_check=False):
    """Run the recovery operation on graph G.

    Fip the bit values of all the vertices in the path connecting each
    pair of stabilizers according to the matching. If check, verify
    that there are no odd-parity cubes remaining, or display their
    indices of there are.

    Args:
        G_match (networkx.Graph): the matching graph
        G_dec (networkx.Graph): the decoding graph
        G (CVGraph): the CVGraph to correct
        matching (set of tuples): the minimum-weight perfect matching
        check (bool): if True, check if the recovery has succeeded.

    Returns:
        None or bool: if check, False if the recovery failed, True if
            it succeeded. Otherwise None.
    """
    virtual_points = G_match.graph["virtual_points"]
    for pair in matching:
        if pair not in it.product(virtual_points, virtual_points):
            path = G_match.edges[pair]["path"]
            pairs = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for pair in pairs:
                common_vertex = G_dec.edges[pair]["common_vertex"]
                code.graph.nodes[common_vertex]["bit_val"] ^= 1

    if sanity_check:
        G_dec_new = decoding_graph(code)
        odd_cubes = G_dec_new.graph["odd_cubes"]
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


def build_dec_and_match_graphs(code, weight_options):
    """Build the decoding and matching graphs.

    Combines weight assignment, decoding and matching graph creation.

    Args:
        code (code): the code class to decode and correct
        weight_options (dict): how to assign weights; options are

                'method': 'unit' or 'blueprint'
                'integer': True (for rounding) or False (for not)
                'multiplier': integer denoting multiplicative factor
                    before rounding
    Returns:
        (EGraph, EGraph): The decoding and matching graphs.
    """
    if weight_options is None:
        weight_options = {}
    assign_weights(code, **weight_options)
    G_dec = decoding_graph(code)
    G_match = matching_graph(G_dec)
    return G_dec, G_match


def correct(
    code,
    decoder,
    weight_options=None,
    sanity_check=False,
):
    """Run through all the error-correction steps.

    Combines weight assignment, decoding and matching graph creation,
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
        G_dec, G_match = build_dec_and_match_graphs(code, weight_options)
        matching = MWPM(G_match, G_dec)
        recovery(code, G_match, G_dec, matching, sanity_check=sanity_check)
    result = check_correction(code, sanity_check=sanity_check)
    return result
