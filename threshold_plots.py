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
"""Plot data from FT simulations."""
import argparse
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12


def process_results(file_name, unit=None, save=True):
    """Process the threshold estimation results in file file_name.

    Add the numbers of trials per distance per delta together, compute
    the failure rate, and create a Pandas DataFrame with columns
    ['p_swap', 'distance', 'delta', 'p_fail', 'error_bar']. Convert
    x-axis to dBs if unit specified, and save if save=True.
    """
    df = pd.read_csv(file_name)
    p_swap_set = set(df.p_swap)
    distance_set = set(df.distance)
    delta_set = list(set(df.delta))
    delta_set.sort(reverse=bool(unit))
    list_of_rows = []
    for p_swap in p_swap_set:
        df_p = df[df.p_swap == p_swap]
        for distance in distance_set:
            df_d = df_p[df.distance == distance]
            for delta in delta_set:
                df_d_del = df_d[df_d.delta == delta]
                total_errors = np.sum(df_d_del.errors)
                total_trials = np.sum(df_d_del.trials)
                if total_trials != 0:
                    p_fail = total_errors / total_trials
                    error_bar = np.sqrt(p_fail * (1-p_fail) / total_trials)
                    if unit == 'dB':
                        delta = -10 * np.log10(delta)
                    row = [p_swap, distance, delta, p_fail, error_bar]
                    list_of_rows += [row]
    columns = ['p_swap', 'distance', 'delta', 'p_fail', 'error_bar']
    df_new = pd.DataFrame(list_of_rows, columns=columns)
    if save:
        split_dir = file_name.split('/')
        split_file_name = split_dir[-1].split('.')
        new_file_name = ('/'.join(split_dir[:-1]) + '/'
                         + split_file_name[0] + '_processed.' + split_file_name[1])
        df_new.to_csv(new_file_name)
    return df_new


def plot_results(data, p_swap=0,
                 file_name=None, unit=None, show=False, rescale=None, threshold=None):
    """Plot processed threshold estimation data.

    Create a plate for swap-out probability p_swap. If save_file
    specified, save to a file with this file name. If unit is dB,
    adjust x-axis label accordingly. If show, display the plot.
    """
    df = data[data.p_swap == p_swap]
    ds = set(df.distance)
    if rescale:
        sigma_t, mu = rescale
        y_err, logy = None, False
    else:
        y_err, logy = 'error_bar', True

    axs = plt.subplot()
    for distance in ds:
        new = df[df.distance == distance].copy()
        if rescale:
            new.delta = (new.delta - sigma_t) * (distance ** (1/mu))
        label = '$d = {}$'.format(int(distance))
        new.plot(x='delta', y='p_fail', style='.', marker='.',
                 markersize=10, ax=axs, logy=logy, yerr=y_err, label=label)
        plt.legend()

    delta_str = r'\Delta_{\mathrm{dB}}' * bool(unit) + r'\delta' * (1 - bool(unit))
    if rescale:
        x_str = r'$({0} - {0}^t)L^{{-\mu}}$'.format(delta_str)
    else:
        x_str = r'${}{}$'.format(delta_str, bool(unit) * r'=-10\log_{{10}}\delta')
    box_props = dict(boxstyle='square', facecolor='w')
    if threshold:
        plt.axvline(x=threshold[0], ls=':', c='dimgrey')
        plt.axhline(y=threshold[1], ls=':', c='dimgrey')
        x_offset = 12 * bool(unit) - 60 * (1 - bool(unit))
        text = r'${}^t = {:.3g}$'.format(delta_str, threshold[0])
        axs.annotate(text, (threshold[0], threshold[1]), xytext=(x_offset, 12),
                     textcoords='offset points', bbox=box_props)
    plt.xlabel(x_str)
    plt.ylabel(r'$p_\mathrm{fail}$')
    p_str = r"$p_{{\mathrm{{swap}}}} = {:.2g}$".format(p_swap)
    axs.annotate(p_str, (0.5, 0.15), xycoords='figure fraction', bbox=box_props)
    plt.style.use('default')
    plt.tight_layout()
    if file_name:
        plt.savefig(file_name)
    if show:
        plt.show()
    return axs


def find_threshold(data, p_swap, unit='dB', file_name=None, plot=True):
    """Estimate the threshold from processed data data.

    Estimate the fault-tolerant threshold from processed data data,
    for swap probability set by p_swap. Optionally, plot the quadratic
    fit for the failure probability over the rescaled delta.
    """
    df = data[data.p_swap == p_swap][data.p_fail < 0.5].copy()
    df = df[df.p_fail > 0.01]

    def delta_rescaled(delta_distance, delta_t, mu):
        result = (delta_distance[0] - delta_t) * delta_distance[1] ** (1/mu)
        return result
    quad_fit = lambda x, a0, a1, a2: a0 + a1 * x + a2 * (x ** 2)

    def fit_func(delta_distance, delta_t, mu, a0, a1, a2):
        result = quad_fit(delta_rescaled(delta_distance, delta_t, mu), a0, a1, a2)
        return result
    try:
        fitted = curve_fit(fit_func, (df.delta, df.distance), df.p_fail)
    except RuntimeError:
        print('No threshold found.')
        return
    delta_t, mu, a0, a1, a2 = fitted[0]
    if plot:
        plot_results(df, p_swap, rescale=(delta_t, mu), unit=unit)
        if unit == 'dB':
            x = np.arange(-3.5, 4.5, 0.01)
        else:
            x = np.arange(-0.1, 0.08, 0.001)
        plt.plot(x, quad_fit(x, a0, a1, a2), ':', c='dimgrey')
        if file_name:
            split_dir = file_name.split('/')
            split_file_name = split_dir[-1].split('.')
            new_file_name = ('/'.join(split_dir[:-1]) + '/'
                         + split_file_name[0] + '_threshold_fit.' + split_file_name[1])
            plt.savefig(new_file_name)
        plt.show()
    print_str = ('The threshold is estimated at delta = {:.3g}{}. '
                 'Here, the failure probability is {:.2g}.'.format(delta_t, ' dB' * bool(unit), a0)
                 )
    print(print_str)
    return delta_t, a0


if __name__ == '__main__':
    # Change options here.
    options = {'input_file': './data/test.csv',
               'save_file': './data/results.png',
               'p_swap': 0,
               'unit': 'dB'
               }

    # For users who wish to use command line.
    parser = argparse.ArgumentParser()
    # Input file
    parser.add_argument('-i', type=str, default=options['input_file'])
    # Save file
    parser.add_argument('-s', type=str, default=options['save_file'])
    # Swap-out probability
    parser.add_argument('-p', type=float, default=options['p_swap'])
    # Unit
    parser.add_argument("-u", type=str, default=options['unit'])
    args = parser.parse_args()

    # Read and process simulation results, and save collected
    results = process_results(args.i, args.u, save=True)
    # Find the threshold and plot the results
    d_p_t = find_threshold(results, args.p, args.u, file_name=args.s)
    plot_results(results, args.p, args.s, args.u, threshold=d_p_t)
