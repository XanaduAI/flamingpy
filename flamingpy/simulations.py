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
"""Monte Carlo simulations for estimating FT thresholds."""

# pylint: disable=too-many-locals,too-many-arguments

import argparse
import csv
import os
import sys

from datetime import datetime
from time import perf_counter

from flamingpy.codes import SurfaceCode
from flamingpy.decoders.decoder import correct
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import BS_network, reduce_macro_and_simulate


def ec_monte_carlo(
    code, trials, delta, p_swap, decoder="MWPM", passive_objects=None, return_decoding_time=False
):
    """Run Monte Carlo simulations of error-correction for the given code.

    Given a code object code, a noise parameter delta, and a
    swap-out probably p_swap, run a number of Monte Carlo
    simulations equal to trials of the complete error-corection
    procedure.

    Args:
        code (code object): the abstract code.
        trials (int): the number of trials.
        delta (float): the noise/squeezing/width parameter.
        p_swap (float): the probability of replacing a GKP state
            with a p-squeezed state in the lattice
        decoder (str): the decoding algorithm ('MWPM' or 'UF')
        passive_objects (NoneType or list, optional): the arguments for
            reduce_macro_and_simulate for passive architecture simulations.
        return_decoding_time (bool, optional): total decoding time is returned when set to True

    Returns:
        errors (integer): the number of errors.
        decoding_time (float): the total number of seconds taken by the decoder. This parameter is
            returned only if return_decoding_time is set to True
    """
    if passive_objects is not None:
        decoder = {"outer": decoder}
        if decoder["outer"] == "MWPM":
            weight_options = {"method": "blueprint", "prob_precomputed": True}
        else:
            weight_options = None
    else:
        # Noise model
        cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
        # Decoding options
        decoder = {"inner": "basic", "outer": decoder}
        if decoder["outer"] == "MWPM":
            weight_options = {
                "method": "blueprint",
                "integer": False,
                "multiplier": 1,
                "delta": delta,
            }
        else:
            weight_options = None

    successes = 0
    if return_decoding_time:
        decoding_time = 0
    for _ in range(trials):
        if passive_objects is not None:
            reduce_macro_and_simulate(*passive_objects, p_swap, delta)
        else:
            # Apply noise
            CVRHG = CVLayer(code, p_swap=p_swap)
            # Measure syndrome
            CVRHG.apply_noise(cv_noise)
            CVRHG.measure_hom("p", code.all_syndrome_inds)

        if return_decoding_time:
            decoding_start_time = perf_counter()
        result = correct(code=code, decoder=decoder, weight_options=weight_options)
        if return_decoding_time:
            decoding_stop_time = perf_counter()
            decoding_time += decoding_stop_time - decoding_start_time
        successes += result
    errors = trials - successes
    if return_decoding_time:
        return errors, decoding_time
    else:
        return errors


# pylint: disable=too-many-arguments
def run_ec_simulation(
    distance, ec, boundaries, delta, p_swap, trials, passive, decoder="MWPM", fname=None
):
    """Run full Monte Carlo error-correction simulations for the surface
    code."""

    # The Monte Carlo simulations

    # The qubit code
    RHG_code = SurfaceCode(distance, ec, boundaries, backend="retworkx")
    RHG_lattice = RHG_code.graph
    RHG_lattice.index_generator()
    if passive:
        # The lattice with macronodes.
        pad_bool = boundaries != "periodic"
        RHG_macro = RHG_lattice.macronize(pad_boundary=pad_bool)
        RHG_macro.index_generator()
        RHG_macro.adj_generator(sparse=True)
        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVLayer(RHG_lattice)
        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = BS_network(4)
        passive_objects = [RHG_macro, RHG_lattice, CVRHG_reduced, bs_network]
    else:
        passive_objects = None

    # Perform the simulation
    simulation_start_time = perf_counter()
    errors, decoding_time = ec_monte_carlo(
        RHG_code, trials, delta, p_swap, decoder, passive_objects, True
    )
    simulation_stop_time = perf_counter()

    # Store results in the provided file-path or by default in
    # a sims_data directory in the file simulations_results.csv.
    file_name = fname or "./flamingpy/sims_data/sims_results.csv"

    # Create a CSV file if it doesn't already exist.
    # pylint: disable=consider-using-with
    try:
        file = open(file_name, "x", newline="", encoding="utf8")
        writer = csv.writer(file)
        writer.writerow(
            [
                "distance",
                "passive",
                "ec",
                "boundaries",
                "delta",
                "p_swap",
                "decoder",
                "errors_py",
                "trials",
                "current_time",
                "decoding_time",
                "simulation_time",
            ]
        )
    # Open the file for appending if it already exists.
    except FileExistsError:
        file = open(file_name, "a", newline="", encoding="utf8")
        writer = csv.writer(file)
    current_time = datetime.now().time().strftime("%H:%M:%S")
    writer.writerow(
        [
            distance,
            passive,
            ec,
            boundaries,
            delta,
            p_swap,
            decoder,
            errors,
            trials,
            current_time,
            decoding_time,
            (simulation_stop_time - simulation_start_time),
        ]
    )
    file.close()


if __name__ == "__main__":
    if len(sys.argv) != 1:
        # Parsing input parameters
        parser = argparse.ArgumentParser(description="Arguments for Monte Carlo FT simulations.")
        parser.add_argument("-distance", type=int)
        parser.add_argument("-ec", type=str)
        parser.add_argument("-boundaries", type=str)
        parser.add_argument("-delta", type=float)
        parser.add_argument("-p_swap", type=float)
        parser.add_argument("-trials", type=int)
        parser.add_argument("-passive", type=lambda s: s == "True")
        parser.add_argument("-decoder", type=str)
        parser.add_argument(
            "-dir", type=str, help="The directory where the result file should be stored"
        )

        args = parser.parse_args()
        params = {
            "distance": args.distance,
            "ec": args.ec,
            "boundaries": args.boundaries,
            "delta": args.delta,
            "p_swap": args.p_swap,
            "trials": args.trials,
            "passive": args.passive,
            "decoder": args.decoder,
        }

    else:
        # User-specified values, if not using command line.
        params = {
            "distance": 2,
            "ec": "primal",
            "boundaries": "open",
            "delta": 0.04,
            "p_swap": 0.5,
            "trials": 100,
            "passive": True,
            "decoder": "MWPM",
        }
    # Checking that a valid decoder choice is provided
    if params["decoder"].lower() in ["unionfind", "uf", "union-find", "union find"]:
        params["decoder"] = "UF"
    elif params["decoder"].lower() in [
        "mwpm",
        "minimum-weight-perfect-matching",
        "minimum weight perfect matching",
    ]:
        params["decoder"] = "MWPM"
    else:
        raise ValueError(f"Decoder {params['decoder']} is either invalid or not yet implemented.")

    # The Monte Carlo simulations
    run_ec_simulation(**params)
