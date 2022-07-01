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

# pylint: disable=wrong-import-position,consider-using-with

import argparse
import csv
import sys
import warnings
import logging

from datetime import datetime
from time import perf_counter

int_time = int(str(datetime.now().timestamp()).replace(".", ""))
logging.info("the following seed was used for random number generation: %i", int_time)

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
reverse_noise_dict = {b: a for a, b in noise_dict.items()}


def ec_mc_trial(
    code_instance,
    noise_instance,
    decoder,
    weight_options,
    rng=default_rng(),
):
    """Runs a single trial of Monte Carlo simulations of error-correction for
    the given code."""
    noise_instance.apply_noise(rng)

    result = correct(code=code_instance, decoder=decoder, weight_options=weight_options)

    return result


def ec_monte_carlo(
    trials,
    code_instance,
    noise_instance,
    decoder,
    decoder_args,
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
        trials (int): the number of trials.
        code_instance (code object): the initialized qubit code
        noise_instance (noise object): the initialized noise layer
            (CVLayer, CVMacroLayer, or IidNoise)
        decoder (str): the decoding algorithm ("MWPM" or "UF")
        deocder_args (dict): arguments for the decoder (such as weight options)
        world_comm, mpi_rank, mpi_size: arguments for the MPI library.

    Returns:
        errors (integer): the number of errors.
    """
    weight_opts = decoder_args["weight_opts"]

    successes = np.zeros(1)
    local_successes = np.zeros(1)

    rng = np.random.default_rng(mpi_rank + int_time)

    for i in range(trials):
        if i % mpi_size == mpi_rank:
            result = ec_mc_trial(
                code_instance,
                noise_instance,
                decoder,
                weight_opts,
                rng,
            )
            local_successes[0] += result

    if "MPI" in globals():
        world_comm.Reduce(local_successes, successes, op=MPI.SUM, root=0)
    else:
        successes[0] = local_successes[0]

    errors = int(trials - successes[0])

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
    noise_instance = noise(code_instance, **noise_args)

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
    errors = ec_monte_carlo(
        trials,
        code_instance,
        noise_instance,
        decoder,
        decoder_args,
        world_comm,
        mpi_rank,
        mpi_size,
    )
    simulation_stop_time = perf_counter()

    if mpi_rank == 0:
        # Store results in the provided file-path or by default in
        # a .sims_data directory in the file simulations_results.csv.
        file_name = fname or "./flamingpy/.sims_data/sims_results.csv"

        # Create a CSV file if it doesn't already exist.
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
                    "error_probability",
                    "decoder",
                    "errors",
                    "trials",
                    "current_time",
                    "simulation_time",
                    "mpi_size",
                ]
            )
        # Open the file for appending if it already exists.
        except FileExistsError:
            file = open(file_name, "a", newline="", encoding="utf8")
            writer = csv.writer(file)
        current_time = datetime.now().time().strftime("%H:%M:%S")
        for key in ["delta", "p_swap", "error_probability"]:
            if key not in noise_args:
                noise_args.update({key: "None"})
        writer.writerow(
            [
                reverse_noise_dict[noise],
                code_args["distance"],
                code_args["ec"],
                code_args["boundaries"],
                noise_args.get("delta"),
                noise_args.get("p_swap"),
                noise_args.get("error_probability"),
                decoder,
                errors,
                trials,
                current_time,
                (simulation_stop_time - simulation_start_time),
                mpi_size,
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
        parser.add_argument("-pswap", type=float)
        parser.add_argument("-errprob", type=float)
        parser.add_argument("-trials", type=int)
        parser.add_argument("-decoder", type=str)

        args = parser.parse_args()
        params = {
            "noise": args.noise,
            "distance": args.distance,
            "ec": args.ec,
            "boundaries": args.boundaries,
            "delta": args.delta,
            "p_swap": args.pswap,
            "error_probability": args.errprob,
            "trials": args.trials,
            "decoder": args.decoder,
        }

    else:
        # User can specify values here, if not using command line.
        params = {
            "noise": "passive",
            "distance": 3,
            "ec": "primal",
            "boundaries": "open",
            "delta": 0.09,
            "p_swap": 0.25,
            "error_probability": 0.1,
            "trials": 100,
            "decoder": "MWPM",
        }
    # Checking that a valid decoder choice is provided
    if params["decoder"].lower() in ["unionfind", "uf", "union-find", "union find"]:
        params["decoder"] = "UF"
    elif params["decoder"].lower() in ["mwpm", "minimum weight perfect matching"]:
        params["decoder"] = "MWPM"
    else:
        raise ValueError(f"Decoder {params['decoder']} is either invalid or not yet implemented.")

    # The Monte Carlo simulations
    code = SurfaceCode
    code_args = {key: params[key] for key in ["distance", "ec", "boundaries"]}

    noise = noise_dict[params["noise"]]
    if params.get("noise") == "iid":
        if params.get("error_probability") is None:
            raise ValueError("No argument `err_prob` found for `iid` noise.")
        noise_args = {"error_probability": params.get("error_probability")}
    else:
        noise_args = {key: params[key] for key in ["delta", "p_swap"]}

    decoder = params["decoder"]
    args = {
        "trials": params["trials"],
        "code": code,
        "code_args": code_args,
        "noise": noise,
        "noise_args": noise_args,
        "decoder": decoder,
    }
    run_ec_simulation(**args)
