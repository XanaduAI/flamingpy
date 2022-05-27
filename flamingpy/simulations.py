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

# pylint: disable=too-many-locals,too-many-arguments,wrong-import-position

import argparse
import csv
import sys
import warnings

from datetime import datetime
from time import perf_counter

int_time = int(str(datetime.now().timestamp()).replace(".", ""))

try:
    import mpi4py.rc

    mpi4py.rc.threaded = False
    from mpi4py import MPI
except ImportError:  # pragma: no cover
    warnings.warn("Failed to import mpi4py libraries.", ImportWarning)

import numpy as np
from numpy.random import default_rng

from flamingpy.codes import SurfaceCode
from flamingpy.decoders.decoder import correct
from flamingpy.noise import CVLayer, CVMacroLayer, IidNoise


noise_dict = {"blueprint": CVLayer, "passive": CVMacroLayer, "iid": IidNoise}
reverse_dict = {b: a for a, b in noise_dict.items()}


def ec_mc_trial(
    bs_network,
    p_swap,
    p_err,
    noise_model,
    code,
    decoder,
    macro_graph,
    weight_options,
    rng=default_rng(),
):
    """Runs a single trial of Monte Carlo simulations of error-correction for the given code."""
    if noise == CVLayer:
        CVRHG = CVLayer(code, p_swap=p_swap, rng=rng)
        CVRHG.apply_noise(noise_model, rng=rng)
        CVRHG.measure_hom("p", code.all_syndrome_inds, rng=rng)
    elif noise == CVMacroLayer:
        CV_macro = CVMacroLayer(
            macro_graph, p_swap=p_swap, reduced_graph=code.graph, bs_network=bs_network, rng=rng
        )
        CV_macro.reduce(noise_model)
    elif noise == IidNoise:
        IidNoise(code, p_err, rng=rng).apply_noise()

    decoding_start_time = perf_counter()

    result = correct(code=code, decoder=decoder, weight_options=weight_options)

    decoding_stop_time = perf_counter()
    decoding_time = decoding_stop_time - decoding_start_time

    return result, decoding_time


def ec_monte_carlo(
    trials,
    code,
    noise,
    noise_args,
    decoder,
    decoder_args,
    return_decoding_time=False,
    world_comm=None,
    mpi_rank=0,
    mpi_size=1,
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
        decoding_time_total (float): the total time is seconds taken by the decoder steps. This
            parameter is returned only if return_decoding_time is set to True
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
        macro_graph, bs_network = noise_args["macro_graph"], noise_args.get("bs_network")
    elif noise == IidNoise:
        p_err = noise_args["p_err"]

    successes = np.zeros(1)
    local_successes = np.zeros(1)

    rng = np.random.default_rng(mpi_rank + int_time)

    if return_decoding_time:
        decoding_time_total = 0

    for i in range(trials):
        if i % mpi_size == mpi_rank:
            result, decoding_time = ec_mc_trial(
                bs_network,
                p_swap,
                p_err,
                noise_model,
                code,
                decoder,
                macro_graph,
                weight_opts,
                rng,
            )
            if return_decoding_time:
                decoding_time_total += decoding_time
            local_successes[0] += result

    if "MPI" in globals():
        world_comm.Reduce(local_successes, successes, op=MPI.SUM, root=0)

    errors = int(trials - successes[0])

    if return_decoding_time:
        return errors, decoding_time_total

    return errors


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
        noise_args.update({"bs_network": None, "macro_graph": macro_graph})
        if decoder == "MWPM":
            weight_opts = {"method": "blueprint", "prob_precomputed": True}
        else:
            weight_opts = None

    # For iid Z errors
    elif noise == IidNoise:
        weight_opts = {"method": "uniform"}
    decoder_args.update({"weight_opts": weight_opts})

    if "MPI" in globals():
        world_comm = MPI.COMM_WORLD
        mpi_size = world_comm.Get_size()
        mpi_rank = world_comm.Get_rank()
    else:
        world_comm = None
        mpi_size = 1
        mpi_rank = 0

    # Perform and time the simulation
    simulation_start_time = perf_counter()
    errors, decoding_time = ec_monte_carlo(
        trials,
        code,
        noise,
        noise_args,
        decoder,
        decoder_args,
        True,
        world_comm,
        mpi_rank,
        mpi_size,
    )
    simulation_stop_time = perf_counter()

    if mpi_rank == 0:
        # Store results in the provided file-path or by default in
        # a sims_data directory in the file simulations_results.csv.
        file_name = fname or ".flamingpy/sims_data/sims_results.csv"

        # Create a CSV file if it doesn't already exist.
        # pylint: disable=consider-using-with
        try:
            file = open(file_name, "x", newline="", encoding="utf8")
            writer = csv.writer(file)
            writer.writerow(
                [
                    "noise",
                    "distance",
                    "ec",
                    "boundaries",
                    "delta",
                    "p_swap",
                    "p_err",
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
                reverse_dict[noise],
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
            "noise": args.noise,
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
            "distance": 2,
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

    noise = noise_dict[params["noise"]]
    noise_args = {key: params[key] for key in ["delta", "p_swap", "p_err"]}

    decoder = params["decoder"]
    args = [params["trials"], code, code_args, noise, noise_args, decoder]
    run_ec_simulation(*args)
