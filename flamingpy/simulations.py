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
from flamingpy.cv.ops import splitter_symp
from flamingpy.decoders.decoder import correct
from flamingpy.noise import CVLayer, CVMacroLayer, IidNoise


def ec_monte_carlo(
    code,
    noise,
    noise_args,
    decoder,
    decoder_args,
    trials,
    return_decoding_time=False,
):
    """Run Monte Carlo simulations of error-correction for the given code.

    Given a code object code, a noise parameter delta, and a
    swap-out probably p_swap, run a number of Monte Carlo
    simulations equal to trials of the complete error-corection
    procedure.

    Args:
        code (instance code object): the initialized qubit code
        noise (noise object): the noise layer (CVLayer, CVMacroLayer, or IidNoise)
        noise_args (dict): the arguments to the noise layer
        decoder (str): the decoding algorithm ("MWPM" or "UF")
        trials (int): the number of trials.
        return_decoding_time (bool, optional): total decoding time is returned when set to True

    Returns:
        errors (integer): the number of errors.
        decoding_time (float): the total number of seconds taken by the decoder. This parameter is
            returned only if return_decoding_time is set to True
    """
    weight_opts = decoder_args["weight_opts"]
    decoder = {"outer": decoder}
    if noise in (CVLayer, CVMacroLayer):
        delta, p_swap = noise_args["delta"], noise_args["p_swap"]
        noise_model = {"noise": "grn", "delta": delta}
    if noise == CVLayer:
        noise_model["sampling_order"] = "initial"
        decoder["inner"] = "basic"
    elif noise == CVMacroLayer:
        macro_graph, bs_network = noise_args["macro_graph"], noise_args["bs_network"]
    elif noise == IidNoise:
        p_err = noise_args["p_err"]

    successes = 0
    if return_decoding_time:
        decoding_time = 0
    for _ in range(trials):
        if noise == CVLayer:
            CVRHG = CVLayer(code, p_swap=p_swap)
            CVRHG.apply_noise(noise_model)
            CVRHG.measure_hom("p", code.all_syndrome_inds)
        elif noise == CVMacroLayer:
            CV_macro = CVMacroLayer(macro_graph, p_swap=p_swap, reduced_graph=code.graph)
            CV_macro.reduce(noise_model, bs_network)
        elif noise == IidNoise:
            IidNoise(code, p_err).apply_noise()

        if return_decoding_time:
            decoding_start_time = perf_counter()

        result = correct(code=code, decoder=decoder, weight_options=weight_opts)

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
    trials, code, code_args, noise, noise_args, decoder, decoder_args=None, fname=None
):
    """Run full Monte Carlo error-correction simulations."""
    # Set up the objects common to all trials.
    if decoder_args is None:
        decoder_args = {}

    # Instance of the qubit QEC code
    code_instance = code(**code_args)
    code_instance.graph.index_generator()

    # The noise model
    # For the blueprint
    if noise == CVLayer:
        if decoder == "MWPM":
            weight_opts = decoder_args.get("weight_opts")
            if weight_opts is None:
                weight_opts = {}
            default_weight_opts = {
                "method": "blueprint",
                "integer": False,
                "multiplier": 1,
                "delta": noise_args.get("delta"),
            }
            weight_opts = {**default_weight_opts, **weight_opts}
        else:
            weight_opts = None

    # For the passive architecture
    elif noise == CVMacroLayer:
        pad_bool = code_args["boundaries"] != "periodic"
        # Instantiate macronode graph and beamsplitter network
        macro_graph = code_instance.graph.macronize(pad_boundary=pad_bool)
        macro_graph.index_generator()
        macro_graph.adj_generator(sparse=True)
        bs_network = splitter_symp()
        noise_args.update({"bs_network": bs_network, "macro_graph": macro_graph})
        if decoder == "MWPM":
            weight_opts = {"method": "blueprint", "prob_precomputed": True}
        else:
            weight_opts = None

    # For iid Z errors
    elif noise == IidNoise:
        weight_opts = {"method": "uniform"}
    decoder_args.update({"weight_opts": weight_opts})

    # Perform and time the simulation
    simulation_start_time = perf_counter()
    errors, decoding_time = ec_monte_carlo(
        code_instance, noise, noise_args, decoder, decoder_args, trials, return_decoding_time=True
    )
    simulation_stop_time = perf_counter()

    # Store results in the provided file-path or by default in
    # a sims_data directory in the file simulations_results.csv.
    file_name = fname or "./sims_data/sims_results.csv"

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
            code_args["distance"],
            code_args["ec"],
            code_args["boundaries"],
            noise_args.get("delta"),
            noise_args.get("p_swap"),
            noise_args.get("p_err"),
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
        parser.add_argument("-noise", type=str)
        parser.add_argument("-distance", type=int)
        parser.add_argument("-ec", type=str)
        parser.add_argument("-boundaries", type=str)
        parser.add_argument("-delta", type=float)
        parser.add_argument("-p_swap", type=float)
        parser.add_argument("-trials", type=int)
        parser.add_argument("-decoder", type=str)
        parser.add_argument(
            "-dir", type=str, help="The directory where the result file should be stored"
        )

        args = parser.parse_args()
        params = {
            "type": args.type,
            "distance": args.distance,
            "ec": args.ec,
            "boundaries": args.boundaries,
            "delta": args.delta,
            "p_swap": args.p_swap,
            "trials": args.trials,
            "decoder": args.decoder,
        }

    else:
        # User-specified values, if not using command line.
        params = {
            "noise": "passive",
            "distance": 3,
            "ec": "primal",
            "boundaries": "open",
            "delta": 0.09,
            "p_swap": 0.25,
            "p_err": 0.1,
            "trials": 100,
            "decoder": "MWPM",
        }
    # Checking that a valid decoder choice is provided
    if params["decoder"].lower() in ["unionfind", "uf", "union-find", "union find"]:
        params["decoder"] = "UF"
    elif params["decoder"].lower() in ["mwpm"]:
        params["decoder"] = "MWPM"
    else:
        raise ValueError(f"Decoder {params['decoder']} is either invalid or not yet implemented.")

    # The Monte Carlo simulations
    code = SurfaceCode
    code_args = {key: params[key] for key in ["distance", "ec", "boundaries"]}

    noise_dict = {"blueprint": CVLayer, "passive": CVMacroLayer, "iid": IidNoise}
    noise = noise_dict[params["noise"]]
    noise_args = {key: params[key] for key in ["delta", "p_swap", "p_err"]}

    decoder = params["decoder"]
    args = [params["trials"], code, code_args, noise, noise_args, decoder]
    run_ec_simulation(*args)
