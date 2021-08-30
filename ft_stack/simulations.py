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
"""Monte Carlo simulations for estimating FT thresholds."""
import argparse
import csv
import sys
from datetime import datetime
from ft_stack.decoder import correct
from ft_stack.graphstates import CVGraph
from ft_stack.RHG import RHG_graph, RHGCode
from ft_stack.passive_construct import BS_network, reduce_macro_and_simulate


def ec_monte_carlo(code, trials, delta, p_swap, passive_objects):
    """Run Monte Carlo simulations of error-correction on code code.

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

    Returns:
        float: the failure probability.
    """
    if passive_objects:
        RHG_macro, RHG_reduced, CVRHG_reduced, bs_network = passive_objects
        decoder = {"outer": "MWPM"}
        weight_options = {"method": "blueprint", "prob_precomputed": True}
    else:
        # TODO: Noise model input.
        code_lattice = code.graph
        # Noise model
        cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
        # Decoding options
        decoder = {"inner": "basic", "outer": "MWPM"}
        weight_options = {
            "method": "blueprint",
            "integer": False,
            "multiplier": 1,
            "delta": delta,
        }

    successes = 0
    for _ in range(trials):
        if passive_objects:
            reduce_macro_and_simulate(*passive_objects, p_swap, delta)
        else:
            # Apply noise
            CVRHG = CVGraph(code_lattice, p_swap=p_swap)
            CVRHG.apply_noise(cv_noise)
            # Measure syndrome
            CVRHG.measure_hom("p", code.syndrome_inds)

        result = correct(code=code, decoder=decoder, weight_options=weight_options)
        successes += result
        errors = trials - successes

    return errors


if __name__ == "__main__":
    # TODO: Intention of below is to allow not to use the command line
    # if desired. Is this appropriate?
    if len(sys.argv) != 1:
        print(sys.argv)
        # Parsing input parameters
        parser = argparse.ArgumentParser(
            description="Arguments for Monte Carlo FT simulations."
        )
        parser.add_argument("distance", type=int)
        parser.add_argument("delta", type=float)
        parser.add_argument("p_swap", type=float)
        parser.add_argument("trials", type=int)
        parser.add_argument("passive", type=bool)

        args = parser.parse_args()
        distance, delta, p_swap, trials = (
            args.distance,
            args.delta,
            args.p_swap,
            args.trials,
        )

    else:
        # User-specified values, if not using command line.
        distance, delta, p_swap, trials, passive = 2, 0.01, 0.5, 100, False

    # The Monte Carlo simulations
    if passive:
        # The lattice with macronodes.
        RHG_macro = RHG_graph(
            distance, boundaries="periodic", macronodes=True, polarity=False
        )
        RHG_macro.index_generator()
        RHG_macro.adj_generator(sparse=True)
        # The reduced lattice.
        RHG_code = RHGCode(distance, boundaries="periodic")
        RHG_reduced = RHG_code.graph
        RHG_reduced.index_generator()
        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVGraph(RHG_reduced)

        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = BS_network(4)
        passive_objects = [RHG_macro, RHG_reduced, CVRHG_reduced, bs_network]
    else:
        boundaries = "finite"
        RHG_code = RHGCode(distance, boundaries=boundaries, polarity=True)
        passive_objects = None

    errors = ec_monte_carlo(RHG_code, trials, delta, p_swap, passive_objects)

    # Store results in the data directory in the file results.csv.
    file_name = "data/results.csv"
    # Create a CSV file if it doesn't already exist.
    try:
        file = open(file_name, "x")
        writer = csv.writer(file)
        writer.writerow(["distance", "delta", "p_swap", "errors", "trials", "time"])
    # Open the file for appending if it already exists.
    except FileExistsError:
        file = open(file_name, "a", newline="")
        writer = csv.writer(file)
    # TODO: Do we need to record time?
    current_time = datetime.now().time().strftime("%H:%M:%S")
    writer.writerow([distance, delta, p_swap, errors, trials, current_time])
    file.close()
