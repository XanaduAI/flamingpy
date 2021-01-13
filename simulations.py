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

import numpy as np
import pandas as pd
import itertools as it
from matplotlib import pyplot as plt
from datetime import date, datetime
import csv

from decoder import correct
from graphstates import CVGraph
from RHG import RHG_graph


def monte_carlo(code_lattice, trials, delta, swap_prob):
    """Run Monte Carlo simulations of error-correction on code_lattice.

    Given a code_lattice EGraph, a noise parameter delta, and a
    swap-out probably swap_prob, run a number of Monte Carlo
    simulations equal to trials of the complete error-corection
    procedure.

    Args:
        code_lattice (EGraph): the abstract code lattice
        trials (int): the number of trials
        delta (float): the noise/squeezing/width parameter
        swap_prob (float): the probability of replacing a GKP state
            with a p-squeezed state in the lattice

    Returns:
        float, float: the failure probability and the associated error.
    """
    successes = 0
    for i in range(trials):
        G = CVGraph(code_lattice, model='grn', swap_prob=swap_prob, delta=delta)
        G.measure_p()
        G.eval_Z_probs_cond()
        result = correct(G, inner='basic', outer='MWPM')[0]
        successes += result
        # successes += np.random.randint(2)
    p_fail = (trials - successes) / trials
    err = np.sqrt((p_fail * (1 - p_fail)) / trials)
    return p_fail, err


if __name__ == '__main__':
    # File stored in data/ with today's date.
    todays_date = date.today().strftime("%d-%m-%Y")
    file_name = 'data/' + todays_date + '.csv'
    # Create a CSV file if it doesn't already exist.
    try:
        file = open(file_name, 'x')
        writer = csv.writer(file)
        writer.writerow(['time', 'distance', 'delta', 'p_swap', 'p_fail', 'error', 'trials'])

    # Open the file for appending if it already exists.
    except Exception:
        file = open(file_name, 'a')
        writer = csv.writer(file)

    # The Monte Carlo simulations
    trials = 1
    distances = [2]
    deltas = [0.01, 0.1]
    probs = [0, 0.25]
    for distance in distances:
        L = distance
        boundaries = ['primal', 'dual', 'primal']
        RHG_lattice = RHG_graph(L, boundaries=boundaries, polarity=1)
        for (d, p) in it.product(deltas, probs):
            p_fail, err = monte_carlo(RHG_lattice, trials, d, p)
            current_time = datetime.now().time().strftime("%H:%M:%S")
            writer.writerow([current_time, distance, d, p, p_fail, err, trials])
    file.close()

    # Read the file and plot.
    table = pd.read_csv(file_name)
    fig, ax = plt.subplots()
    table.plot.scatter(x='delta', y='p_fail', ax=ax)
    ax.set_xlabel(r'$\delta$')
    ax.set_ylabel('$p_{fail}$')
    plt.show()
