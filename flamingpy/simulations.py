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

import csv
import warnings
import logging
import cProfile
import pstats

from datetime import datetime
from time import perf_counter
from pathlib import Path
import argh

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

DEFAULT_PATH_DATA = Path("flamingpy/.sims_data/")
DEFAULT_PATH_PROFILING = Path("flamingpy/.profiling/")
DEFAULT_PATH_PROFILES = DEFAULT_PATH_PROFILING / "profiles"


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
    trials,
    code,
    code_args,
    noise,
    noise_args,
    decoder="MWPM",
    decoder_args=None,
    fname=None,
):
    """Run full Monte Carlo error-correction simulations."""
    # Start timing the simulation
    simulation_start_time = perf_counter()

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

    # Perform the simulation
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
        file_name = fname or Path("./flamingpy/.sims_data/sims_results.csv")

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
                    "backend_code",
                    "backend_decoder",
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
                code_args["backend"],
                decoder_args["backend"],
            ]
        )
        file.close()


def run_ec_simulation_with_profiler(
    trials,
    code,
    code_args,
    noise,
    noise_args,
    decoder="MWPM",
    decoder_args=None,
    fname=None,
    profilename=None,
):

    # Running simulation with profiler
    with cProfile.Profile() as profiler:
        run_ec_simulation(
            trials,
            code,
            code_args,
            noise,
            noise_args,
            decoder=decoder,
            decoder_args=decoder_args,
            fname=fname,
        )

    # Saving profiling results
    stats = pstats.Stats(profiler).sort_stats("cumtime")

    stats.dump_stats(profilename)


def _simulation_name(
    decoder, boundaries, noise, ec, error_probability=None, delta=None, p_swap=None
):
    """Generate a unique name for a simulation.

    Args:
        decoder (str): The decoder used in the simulation.
        boundaries (str): The boundaries used in the simulation.
        noise (str): The noise used in the simulation.
        error_probability (float): The error probability used in the simulation.
        delta (float): The delta used in the simulation.
        p_swap (float): The probability of a swap used in the simulation.

    Returns:
        str: The name of the simulation formatted as
            f"{decoder}_{code_args['boundaries']}_{code_args['ec']}_{noise_params}
    """
    noise_params = f"{noise}"
    if noise == "iid":
        noise_params += f"_perr-{error_probability}"
    elif noise in ("passive", "blueprint"):
        noise_params += f"_delta-{delta}_pswap-{p_swap}"
    noise_params = noise_params.replace(".", "")

    return f"{decoder}_{boundaries}_{ec}_{noise_params}"


def simulations(
    noise="passive",
    distance=3,
    ec="primal",
    boundaries="open",
    delta=0.09,
    p_swap=0.25,
    error_probability=0.1,
    trials=100,
    decoder="MWPM",
    backend="retworkx",
    profile=False,
    data_folder=DEFAULT_PATH_DATA,
    profile_folder=DEFAULT_PATH_PROFILES,
):
    """Run simulations for the given parameters.

    Args:
        noise (str): The noise used in the simulation. Options are "iid", "passive", and
            "blueprint".
        distance (int): Code distance of the surface code.
        ec (str): The error complex used in the simulation. Options are "primal" or "dual".
        boundaries (str): The boundaries of the surface code. Options are "open", "periodic", or
            "toric".
        delta (float): The finite-energy parameter delta in case of "blueprint" or "passive" noise.
        p_swap (float): The probability of a swap-out in the hybrid CV cluster state in case of
            "blueprint" or "passive" noise.
        error_probability (float): The error probability used in the simulation in case of "iid"
            noise.
        trials (int): The number of trials used in the Monte-Carlo simulation.
        decoder (str): The decoder used. Options are "MWPM" (Minimum Weight Perfect Matching
            decoder) and "UF" (Union Find decoder)
        backend (str): The backend used for the code and decoder
        profile (bool): Whether to run the simulation with profiler.
        data_folder (str): The path to the folder where the data will be saved.
        profile_folder (str): The path to the folder where the profiling data will be saved if
            `profiling == True`.

    Example for running on the command line::
        python flamingpy/simulations.py --distance 5 --noise passive --ec primal --boundaries open
            --delta 0.08 --p-swap 0.25 --decoder MWPM --trials 1000 --backend retworkx --profile

    """

    # Checking that a valid decoder choice is provided
    if decoder.lower() in ["unionfind", "uf", "union-find", "union find"]:
        decoder = "UF"
    elif decoder.lower() in ["mwpm", "minimum weight perfect matching"]:
        decoder = "MWPM"
    else:
        raise ValueError(f"Decoder {decoder} is either invalid or not yet implemented.")

    noise_class = noise_dict[noise]
    noise_args = {"delta": delta, "p_swap": p_swap}
    # check that arg error_probability is provided for iid noise
    if noise == "iid":
        if error_probability is None:
            raise ValueError(f"No argument `error_probability` found for iid noise.")

        # set to None unused args and update noise_args
        delta, p_swap = None, None
        noise_args = {"error_probability": error_probability}

    # The Monte Carlo simulations
    code = SurfaceCode

    # Code args
    code_args = {"distance": distance, "ec": ec, "boundaries": boundaries, "backend": backend}

    # Decoder args
    decoder_args = {"backend": backend}

    # Collective args
    args = [trials, code, code_args, noise_class, noise_args, decoder, decoder_args]

    # Setting up file paths for data storage
    simname = _simulation_name(
        decoder,
        boundaries,
        noise,
        ec,
        error_probability=error_probability,
        delta=delta,
        p_swap=p_swap,
    )
    fname = data_folder / f"sims_{simname}.csv"

    # Running the simulations
    if not profile:
        # Simulation without profiler
        run_ec_simulation(*args, fname=fname)
    else:
        # Simulation with profiler
        profilename = profile_folder / f"profile_distance_{code_args['distance']}_{simname}"
        run_ec_simulation_with_profiler(
            *args,
            fname=fname,
            profilename=profilename,
        )


if __name__ == "__main__":
    argh.dispatch_command(simulations)
