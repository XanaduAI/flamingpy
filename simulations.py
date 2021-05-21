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
from decoder import correct
from graphstates import CVGraph
from RHG import RHGCode


def ec_monte_carlo(code, trials, delta, p_swap):
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
    # TODO: Noise model input.
    code_lattice = code.graph
    # Noise model
    CVRHG = CVGraph(code_lattice, p_swap=p_swap)
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
        # Apply noise
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
    boundaries = "finite"
    RHG_code = RHGCode(distance, boundaries=boundaries, polarity=True)
    errors = ec_monte_carlo(RHG_code, trials, delta, p_swap)

    # Store results in the data directory in the file results.csv.
    file_name = "./data/results.csv"
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
