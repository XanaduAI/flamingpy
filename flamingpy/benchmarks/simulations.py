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
"""Benchmark for the Monte Carlo simulations estimating FT thresholds and
comparing python and cpp loops."""
import csv

from time import process_time
from datetime import datetime

from flamingpy.codes import SurfaceCode, alternating_polarity
from flamingpy.decoders.decoder import correct
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import BS_network, reduce_macro_and_simulate
import flamingpy.cpp.cpp_mc_loop as cmc


def ec_monte_carlo(code, trials, delta, p_swap, passive_objects=None, backend="cpp"):
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
            with a p-squeezed state in the lattice.
        backend (string): choose from "cpp" Monte Carlo sampling
            loops or pure "python" version
        passive_objects (NoneType or list): the arguments for
            reduce_macro_and_simulate for passive architecture simulations.

    Returns:
        errors (integer): the number of errors.
    """
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
    if passive_objects is not None:
        RHG_macro, RHG_reduced, CVRHG_reduced, bs_network = passive_objects
        decoder = {"outer": "MWPM"}
        weight_options = {"method": "blueprint", "prob_precomputed": True}
    else:
        code_lattice = code.graph
        # Noise model
        # Decoding options
        decoder = {"inner": "basic", "outer": "MWPM"}
        weight_options = {
            "method": "blueprint",
            "integer": False,
            "multiplier": 1,
            "delta": delta,
        }
    if backend == "python":
        successes = 0
        for _ in range(trials):
            if passive_objects is not None:
                reduce_macro_and_simulate(*passive_objects, p_swap, delta)
            else:
                # Apply noise
                CVRHG = CVLayer(code_lattice, p_swap=p_swap)
                CVRHG.apply_noise(cv_noise)
                # Measure syndrome
                CVRHG.measure_hom("p", code.all_syndrome_inds)
            result = correct(code=code, decoder=decoder, weight_options=weight_options)
            successes += result
        errors = trials - successes
    elif backend == "cpp":
        errors = cmc.cpp_mc_loop(
            code, trials, decoder, weight_options, passive_objects, p_swap, delta, cv_noise
        )

    return errors


# The Monte Carlo simulations

# Arguments
distance = 2
ec = "primal"
boundaries = "open"
delta = 0.01
p_swap = 0.5
trials = 100
passive = True

# The qubit code
RHG_code = SurfaceCode(distance, ec, boundaries, alternating_polarity)
RHG_lattice = RHG_code.graph
RHG_lattice.index_generator()
if passive:
    # The lattice with macronodes.
    pad_bool = False if boundaries == "periodic" else True
    RHG_macro = RHG_lattice.macronize(pad_boundary=pad_bool)
    RHG_macro.index_generator()
    RHG_macro.adj_generator(sparse=True)
    # The empty CV state, uninitiated with any error model.
    CVRHG_reduced = CVLayer(RHG_lattice)
    # Define the 4X4 beamsplitter network for a given macronode.
    # star at index 0, planets at indices 1-3.
    bs_network = BS_network(4)
    passive_objects = [RHG_macro, RHG_lattice, CVRHG_reduced, bs_network]

tic = process_time()
errors = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects, backend="python")
walltime_py = process_time() - tic

tic = process_time()
ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects, backend="cpp")
walltime_cpp = process_time() - tic

# Store results in a sims_data directory in the file simulations_results.csv.
file_name = "./flamingpy/sims_data/sims_benchmark_results.csv"
# Create a CSV file if it doesn't already exist.
try:
    file = open(file_name, "x")
    writer = csv.writer(file)
    writer.writerow(
        [
            "distance",
            "ec",
            "boundaries",
            "delta",
            "p_swap",
            "errors_py",
            "trials",
            "current_time",
            "cpp_to_py_speedup",
        ]
    )
# Open the file for appending if it already exists.
except FileExistsError:
    file = open(file_name, "a", newline="")
    writer = csv.writer(file)
current_time = datetime.now().time().strftime("%H:%M:%S")
writer.writerow(
    [
        distance,
        ec,
        boundaries,
        delta,
        p_swap,
        errors,
        trials,
        current_time,
        walltime_py / walltime_cpp,
    ]
)
file.close()
