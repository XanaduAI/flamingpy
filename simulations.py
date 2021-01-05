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
"""Module for Monte Carlo simulations for estimating FT thresholds."""

import numpy as np
import pandas as pd
import itertools as it
from matplotlib import pyplot as plt
from datetime import date, datetime
import csv

from decoder import correct
from graphstates import CVGraph
from RHG import RHG_graph


def monte_carlo(trials, distance, delta, swap_prob, save=False):
    # L = distance
    # RHG_lattice = RHG_graph(L, pol=1)

    successes = 0
    for i in range(trials):
        # G = CVGraph(RHG_lattice, model='grn', swap_prob=swap_prob, delta=delta)
        # G.measure_p()
        # G.eval_Z_probs_cond()
        # result = correct(G, inner='basic', outer='MWPM', bc='non-periodic')
        # successes += result
        successes += np.random.randint(2)
    p_fail = (trials - successes) / trials
    err = np.sqrt((p_fail * (1 - p_fail)) / trials)
    return p_fail, err

if __name__ == '__main__':
    todays_date = date.today().strftime("%d-%m-%Y")
    file_name = 'data/' + todays_date + '.csv'
    try:
        file = open(file_name, 'x')
        # writer = csv.DictWriter(file, ['time', 'delta', 'p_swap', 'p_fail', 'error', 'trials'])
        # writer.writeheader()
        writer = csv.writer(file)
        writer.writerow(['time', 'distance', 'delta', 'p_swap', 'p_fail', 'error', 'trials'])

    except Exception:
        file = open(file_name, 'a')
        writer = csv.writer(file)

    trials = 100
    distances = [1, 2]
    deltas = [0.01, 0.1]
    probs = [0, 0.25]
    for (l, d, p) in it.product(distances, deltas, probs):
        p_fail, err = monte_carlo(trials, l, d, p, True)
        current_time = datetime.now().time().strftime("%H:%M:%S")
        writer.writerow([current_time, l, d, p, p_fail, err, trials])

    file.close()
    table = pd.read_csv(file_name)

    fig, ax = plt.subplots()
    table.plot.scatter(x='delta', y='p_fail', ax=ax)

    ax.set_xlabel('$\delta$')
    ax.set_ylabel('$p_{fail}$')
    plt.show()

    print(table)

