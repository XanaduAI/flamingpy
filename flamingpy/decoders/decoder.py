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
"""Decoding and recovery functions."""

# pylint: disable=import-outside-toplevel,too-many-locals,too-many-arguments

import sys
import numpy as np

from flamingpy.cv.gkp import GKP_binner, Z_err_cond

from flamingpy.decoders.mwpm import mwpm_decoder
from flamingpy.decoders.unionfind import uf_decoder

# Smallest and largest numbers representable.
smallest_number = sys.float_info.min
largest_number = sys.float_info.max


def assign_weights(code, decoder, **kwargs):
    """Assign weights, reflecting error probabilities, to qubits in code.

    Args:
        code (SurfaceCode): the qubit QEC code
        decoder (str): the decoder used
        method (str, optional): the method for weight assignment. By
            default, 'uniform', denoting equal weights everywhere. For
            heuristic and analog weight assignment from Xanadu's blueprint, use
            'blueprint' (compatible with the MWPM decoder).
        integer (bool, optional): whether to convert weights to
            integers using Python's round function; False by default
        multiplier (int, optional): multiply the weight by multiplier
            before rounding; 1 by default
        delta (float, optional): the CV noise parameter.

    Returns:
        None
    """
    default_options = {"method": "uniform", "integer": False, "multiplier": 1}
    weight_options = {**default_options, **kwargs}
    G = code.graph
    # Get the set of qubits that are used for parity check measurements.
    qubit_coords = set(code.all_syndrome_coords)

    # Blueprint weight assignment or weighted-union-find weight assignment
    # dependent on the type of neighbours.
    if weight_options.get("method") == "blueprint" and decoder == "MWPM":
        for node in qubit_coords:
            neighbors = G[node]
            # Obtain the list and the number of p-squeezed states in the neighborhood of the node.
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
                err_prob = min(err_prob, 0.5)
                # TODO: Can I just choose an arbitrary small number?
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
        if decoder == "UF":
            raise Exception("Incompatible decoder/weight options.")
    # Naive weight assignment, unity weights.
    elif weight_options.get("method") == "uniform":
        if decoder == "UF":
            # Weight assignment for Union-Find decoder
            for node in qubit_coords:
                neighbors = G[node]
                # List and number of p-squeezed states in neighborhood of node.
                p_count = 0
                for neighbor in neighbors:
                    state_name = G.nodes[neighbor].get("state")
                    if state_name == "p":
                        p_count += 1
                if p_count in (0, 1):
                    # We consider weight-2 edges as we need half edges
                    G.nodes[node]["weight"] = 2
                else:
                    # These edges correspond to the erased edges fed toUnion-Find
                    G.nodes[node]["weight"] = -1
        else:
            for node in qubit_coords:
                G.nodes[node]["weight"] = 1
    # Also assign the weights to the stabilizer graph edges.
    for ec in code.ec:
        getattr(code, f"{ec}_stab_graph").assign_weights(code)


def CV_decoder(code, translator=GKP_binner):
    """Convert homodyne outcomes to bit values according to translator.

    The inner (CV) decoder, a.k.a. translator, a.k.a binning function. Set
    converted values to the bit_val attribute for nodes in G.

    Args:
        code (SurfaceCode): the qubit QEC code
        translator (func): the choice of binning function; by default, the
            standard GKP binning function that snaps to the closest
            integer multiple of sqrt(pi).
    Returns:
        None
    """
    for point in code.all_syndrome_coords:
        hom_val = code.graph.nodes[point]["hom_val_p"]
        bit_val = translator([hom_val])[0]
        code.graph.nodes[point]["bit_val"] = bit_val


def recovery(qubits_to_flip, code, ec, sanity_check=False):
    """Run recovery operations based on code based on qubits_to_flip.

    Fip the bit values of all the vertices in qubits_to_flip.
    If sanity_check is True, verify that there are no odd-parity cubes
    remaining, or display their indices of there are.

    Args:
        qubits_to_flip (iterable): the set of nodes of the code EGraph
            whose bit values oughtto be flipped
        code (SurfaceCode): the qubit QEC code
        ec (string): the error complex ('primal' or 'dual')
        sanity_check (bool): if True, check if the recovery has succeeded
            and print a message.

    Returns:
        None
    """
    for qubit in qubits_to_flip:
        code.graph.nodes[qubit]["bit_val"] ^= 1
    if sanity_check:
        stab_graph = getattr(code, ec + "_stab_graph")
        odd_cubes = list(stab_graph.odd_parity_stabilizers())
        if odd_cubes:
            print("Unsatisfied " + ec + " stabilizers:", odd_cubes)
        else:
            print(ec.capitalize() + " recovery succeeded - no unsatisfied stabilizers.")


def check_correction(code, sanity_check=False):
    """Check whether the error correction has succeded or failed.

    Verify that no logical operator has been applied at the end of the
    recovery. This is achieved by checking that the total parity on the
    appropriate correlation or gauge surface(s) is even. Such as surface
    consists of the set of all syndrome qubits (of the kind in code.ec),
    in a plane of the correct kind ("primal" or "dual") along a
    specified direction.

    For all-periodic boundaries, total parities of planes along all
    three directions are compared. For open boundaries, a plane parallel
    to the primal or dual spatial (x or y) boundary (if the error complex
    is primal or dual, respectively) is checked. In every case, the first
    plane along the specified direction of the right type is checked
    (although this doesn't matter for the following reason).

    Since a logical operator cuts through a set of parallel planes at
    a single point per plane, we expect parities along all parallel
    planes to be the same: this is the idea behind sanity_check.

    Args:
        code (SurfaceCode): the qubit QEC code. At the stage this functions is
            run, the recovery has already been applied, so that code.graph
            represents the error-corrected graph state.
        sanity_check (bool): if True, display the total parity of
            all parallel correlation surfaces to verify if parity
            is conserved.

    Returns:
        list or (list, list): a list of bools indicating whether error
            correction succeeded for each complex. If sanity_check is set to
            True, also output a dictionary between planes and results of
            the parallel-plane-parity sanity check.
    """
    dims = np.array(code.dims)
    dir_dict = {"x": 0, "y": 1, "z": 2}

    if sanity_check:
        print()

    ec_checks = []
    truth_dicts = []
    for ec in code.ec:
        planes_to_check = []
        truth_dict = {"x": [], "y": [], "z": []}
        if code.bound_str == "periodic":
            planes_to_check = ["x", "y", "z"]
        elif code.bound_str.startswith("open"):
            planes_to_check = ["x"] if ec == "primal" else ["y"]

        minimum = 0 if ec == "primal" else 1
        for plane_str in planes_to_check:
            maximum = 2 * dims[dir_dict[plane_str]] if sanity_check else minimum + 2
            for sheet in range(minimum, maximum, 2):
                slice_verts = code.graph.slice_coords(plane_str, sheet)
                syndrome_verts = getattr(code, ec + "_syndrome_coords")
                only_syndrome = set(slice_verts) & set(syndrome_verts)
                parity = 0
                for node in only_syndrome:
                    parity ^= code.graph.nodes[node]["bit_val"]
                truth_dict[plane_str].append(bool(1 - parity))

        if sanity_check:
            print(ec.capitalize() + " error correction check --", truth_dict)

        all_surfaces = [truth_dict[i][0] for i in planes_to_check]
        ec_checks += [np.all(all_surfaces)]
        truth_dicts += [truth_dict]
    if sanity_check:
        return ec_checks, truth_dicts
    return ec_checks


def correct(
    code,
    decoder=None,
    weight_options=None,
    sanity_check=False,
    decoder_opts=None,
    draw=False,
    drawing_opts=None,
):
    """Run through all the error-correction steps.

    Combines weight assignment, inner decoding and outer decoding. The last
    of these includes matching graph creation, minimum-weight-perfect matching
    or union find decoding, recovery, and correctness check.

    Args:
        code (code): the code class to decode and correct
        decoder (dict, optional): a dictionary of the form

            {"inner": f, "outer": s},

            where f is the string indicating inner/CV decoder and s is the
            string indicating the outer/DV decoder to use.

            Currently available inner decoders are:

                "basic": standard GKP binning function

            and outer decoders:

                "MWPM": minimum-weight perfect matching (the default)
                "UF": Union-Find.
        weight_options (dict, optional): how to assign weights for the outer
            decoder; options are

            'method': 'uniform' or 'blueprint' (latter for MWPM decoder)
            'integer': True (for rounding) or False (for not)
            'multiplier': integer denoting multiplicative factor
            before rounding

            Uniform weights by default.
        sanity_check (bool, optional): if True, check that the recovery
            operation has succeeded and verify that parity is conserved
            among all correlation surfaces
        decoder_opts (dict, optional): a dictionary of decoding options,
            including the backend ("networkx" or "retworkx" for "MWPM")
        draw (bool, optional): set to True to illustrate the decoding
            procedure, including the stabilizer graph, syndrome plot, and recovery.
        drawing_opts (dict, optional): the drawing options to be fed into
            viz.draw_decoding (see that function for more details).
    Returns:
        result (bool): True if error correction succeded, False if not.
    """
    if decoder is None:
        decoder = {}

    default_decoder = {"inner": None, "outer": "MWPM"}
    updated_decoder = {**default_decoder, **decoder}
    outer_dict = {"MWPM": mwpm_decoder, "UF": uf_decoder}

    if updated_decoder["inner"] == "basic":
        CV_decoder(code, translator=GKP_binner)

    if weight_options is None:
        weight_options = {}

    outer_decoder_str = updated_decoder["outer"]
    assign_weights(code, outer_decoder_str, **weight_options)

    if decoder_opts is None:
        decoder_opts = {}
    default_decoder_opts = {"backend": "retworkx", "draw": draw, "drawing_opts": drawing_opts}
    updated_decoder_opts = {**default_decoder_opts, **decoder_opts}
    outer_decoder = outer_dict[outer_decoder_str]

    for ec in code.ec:
        recovery_set = outer_decoder(code, ec, **updated_decoder_opts)
        recovery(recovery_set, code, ec, sanity_check=sanity_check)
    result = check_correction(code, sanity_check=sanity_check)
    if sanity_check:
        return np.all(result[0])
    return np.all(result)
