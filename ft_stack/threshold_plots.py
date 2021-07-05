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

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "DejaVu Serif"
plt.rcParams["font.size"] = 11


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
                    error_bar = np.sqrt(p_fail * (1 - p_fail) / total_trials)
                    if unit == "dB":
                        delta = -10 * np.log10(delta)
                    elif unit == "loss":
                        delta = 1 / (1 + delta)
                    row = [p_swap, distance, delta, p_fail, error_bar]
                    list_of_rows += [row]
    columns = ["p_swap", "distance", "delta", "p_fail", "error_bar"]
    df_new = pd.DataFrame(list_of_rows, columns=columns)
    if save:
        split_dir = file_name.split("/")
        split_file_name = split_dir[-1].split(".")
        new_file_name = (
            "/".join(split_dir[:-1])
            + "/"
            + split_file_name[0]
            + "_processed."
            + split_file_name[1]
        )
        df_new.to_csv(new_file_name)
    return df_new


def plot_results(
    data,
    p_swap=0,
    distances="all",
    unit=None,
    threshold=None,
    rescale=None,
    show=False,
    file_name=None,
    swap_tol_plot=None,
    inset=False,
    break_axis=True,
):
    """Plot processed threshold estimation data.

    Create a plot of error rate vs. delta for a given swap-out
    probability p_swap.

    Args:
        data (pd.DataFrame): the simulations data with columns
            ['p_swap', 'distance', 'delta', 'p_fail', 'error_bar']
            (not necessarily in this order).
        p_swap (float, optional): the swap-out probability.
        unit (str, optional): 'dB' for dBs
        threshold (float, float): the error threshold and logical
            error rate and threshold. If supplied, identifies and
            labels the threshold on the plot.
        rescale (tup): threshold and fitting parameter mu for plots
            with rescaled physical errors (see find_threshold).
        show (bool): if True, display the plot.
        file_name (str): if supplied, where to save the plot.

    Returns:
        mpl.Axes: the Matplotlib axes object.
    """
    delta_str = r"\epsilon_{\mathrm{dB}}" * bool(unit) + r"\delta" * (1 - bool(unit))
    # delta_str = r"\eta"
    box_props = dict(boxstyle="square", facecolor="w")

    if swap_tol_plot:
        swap_tol_data = swap_tol_plot
        zipped = list(zip(*swap_tol_data))
        deltas, ps = zipped[0], zipped[1]
        if break_axis:
            factor = 15
            fig, (main_axs, main_axs2) = plt.subplots(
                1,
                2,
                sharey=True,
                figsize=(6, 4.5),
                gridspec_kw={"width_ratios": [factor, 1]},
            )
            main_axs2.plot(deltas, ps, ".-", markersize=12, color="xkcd:navy")
            main_axs2.fill_between(deltas, ps, color="whitesmoke")
            main_axs2.set_xlim(58, 62)
            main_axs2.spines["left"].set_visible(False)
            main_axs.spines["right"].set_visible(False)
            main_axs2.yaxis.tick_right()
            main_axs2.tick_params(right=False, labelright="off")
        else:
            _, main_axs = plt.subplots(figsize=(6, 4.5))

        main_axs.plot(deltas, ps, ".-", markersize=12, color="xkcd:navy")
        main_axs.fill_between(deltas, ps, color="whitesmoke")

        main_axs.annotate("correctable region", (0.6, 0.75), xycoords="axes fraction")

        if break_axis:
            main_axs.set_xlim(10, 25)
            xticks = main_axs2.get_xticks()
            updated = [""] * (len(xticks) - 2) + [r"$\infty$"] + [""]
            main_axs2.set_xticks(xticks)
            main_axs2.set_xticklabels(updated)
            fig.subplots_adjust(wspace=0)

            d = 0.015
            kwargs = dict(transform=main_axs.transAxes, color="k", clip_on=False)
            main_axs.plot((1 - d, 1 + d), (-d, d), **kwargs)
            main_axs.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

            kwargs.update(transform=main_axs2.transAxes)  # switch to the bottom axes
            main_axs2.plot((-d * factor, d * factor), (1 - d, 1 + d), **kwargs)
            main_axs2.plot((-d * factor, d * factor), (-d, d), **kwargs)
        else:
            xticks = main_axs.get_xticks()[1:-1]
            updated = ["{:.2g}".format(i) for i in xticks[:-1]] + [r"$\infty$"]
            # updated = ['{:.2g}'.format(i) for i in xticks[:-1]] + ['1']
            plt.xticks(xticks, updated)

    if (swap_tol_plot and inset) or (not swap_tol_plot):
        df = data[data.p_swap == p_swap][data.delta < 10.5]
        ds = set(df.distance) if distances == "all" else distances
        if rescale:
            sigma_t, mu = rescale
            y_err, logy = None, False
        else:
            y_err, logy = "error_bar", True

        if inset:
            inset_x, inset_y = 0.45, 0.2
            axs = main_axs.inset_axes([inset_x, inset_y, 0.48, 0.48])
            plt.rcParams["font.size"] = 8
            if threshold:
                delta_inset = threshold[0]
            else:
                delta_inset_index = np.where(np.array(zipped[1]) == p_swap)[0][0]
                delta_inset = zipped[0][delta_inset_index]
            main_axs.annotate(
                "",
                xy=(delta_inset + 0.1, p_swap + 0.005),
                # xy=(0.915, p_swap),
                xytext=(inset_x, inset_y),
                textcoords="axes fraction",
                arrowprops=dict(arrowstyle="->"),
            )
        else:
            _, axs = plt.subplots()

        for distance in ds:
            new = df[df.distance == distance].copy()
            if rescale:
                new.delta = (new.delta - sigma_t) * (distance ** (1 / mu))
            label = "$d = {}$".format(int(distance))
            new.plot(
                x="delta",
                y="p_fail",
                style=".",
                marker=".",
                markersize=10,
                ax=axs,
                logy=logy,
                yerr=y_err,
                label=label,
            )
            axs.legend()

        if threshold:
            axs.axvline(x=threshold[0], ls=":", c="dimgrey")
            axs.axhline(y=threshold[1], ls=":", c="dimgrey")
            x_offset = 12 * bool(unit) - 60 * (1 - bool(unit))
            text = r"${}^t = {:.3g}$".format(delta_str, threshold[0])
            if not swap_tol_plot:
                axs.annotate(
                    text,
                    (threshold[0], threshold[1]),
                    xytext=(x_offset - 75, 12 - 30),
                    textcoords="offset points",
                    bbox=box_props,
                )

        if inset:
            axs.set_xlabel("")
            axs.set_ylabel(r"$p_{\mathrm{fail}}$")
        # if not swap_tol_plot:
        #     p_str = r"$p_{{\mathrm{{swap}}}} = {:.2g}$".format(p_swap)
        #     axs.annotate(p_str, (0.5, 0.15), xycoords="figure fraction", bbox=box_props)

    if rescale:
        x_str = r"$({0} - {0}^t)L^{{-\mu}}$".format(delta_str)
    else:
        x_str = r"${}{}$".format(delta_str, bool(unit) * r"=-10\log_{{10}}2\epsilon")
    y_str = r"$p_\mathrm{swap}$" if swap_tol_plot else r"$p_\mathrm{fail}$"
    if swap_tol_plot:
        ax = main_axs
    else:
        ax = axs
    ax.set_xlabel(x_str, fontsize=13)
    ax.set_ylabel(y_str, fontsize=13)
    # plt.style.use("default")
    plt.tight_layout()
    # fig.patch.set_alpha(0)
    if file_name:
        plt.savefig(file_name, dpi=300)
    if show:
        plt.show()
    return ax


def find_threshold(data, p_swap, unit="dB", file_name=None, plot=True):
    """Estimate the threshold from processed data data.

    Estimate the fault-tolerant threshold from processed data data
    for swap probability set by p_swap. Optionally, plot the quadratic
    fit for the failure probability over the rescaled delta. Uses the
    fitting procedure from https://doi.org/10.7907/AHMQ-EG82
    (Harrington 2014).
    """
    df = data[data.p_swap == p_swap][data.p_fail < 0.5].copy()
    df = df[df.p_fail > 0.01]

    def delta_rescaled(delta_distance, delta_t, mu):
        result = (delta_distance[0] - delta_t) * delta_distance[1] ** (1 / mu)
        return result

    quad_fit = lambda x, a0, a1, a2: a0 + a1 * x + a2 * (x ** 2)

    def fit_func(delta_distance, delta_t, mu, a0, a1, a2):
        result = quad_fit(delta_rescaled(delta_distance, delta_t, mu), a0, a1, a2)
        return result

    try:
        fitted = curve_fit(fit_func, (df.delta, df.distance), df.p_fail)
    except RuntimeError:
        print("No threshold found.")
        return
    delta_t, mu, a0, a1, a2 = fitted[0]
    if plot:
        plot_results(df, p_swap, rescale=(delta_t, mu), unit=unit)
        if unit == "dB":
            x = np.arange(-3.5, 4.5, 0.01)
        else:
            x = np.arange(-0.04, 0.04, 0.001)
        plt.plot(x, quad_fit(x, a0, a1, a2), ":", c="dimgrey")
        if file_name:
            split_dir = file_name.split("/")
            split_file_name = split_dir[-1].split(".")
            new_file_name = (
                "/".join(split_dir[:-1])
                + "/"
                + split_file_name[0]
                + "_threshold_fit."
                + split_file_name[1]
            )
            plt.savefig(new_file_name)
        plt.show()
    print_str = (
        "The threshold is estimated at delta = {:.3g}{}. "
        "Here, the failure probability is {:.3g}.".format(
            delta_t, " dB" * bool(unit), a0
        )
    )
    print(print_str)
    return delta_t - 0.02, a0


if __name__ == "__main__":
    # Change options here.
    options = {
        "input_file": "./data/test_processed.csv",
        "save_file": "./data/results.pdf",
        "p_swap": 0,
        "unit": "dB",
    }

    # For users who wish to use command line.
    parser = argparse.ArgumentParser()
    # Input file
    parser.add_argument("-i", type=str, default=options["input_file"])
    # Save file
    parser.add_argument("-s", type=str, default=options["save_file"])
    # Swap-out probability
    parser.add_argument("-p", type=float, default=options["p_swap"])
    # Unit
    parser.add_argument("-u", type=str, default=options["unit"])
    args = parser.parse_args()

    # Read and process simulation results, and save collected
    # results = process_results(args.i, args.u, save=True)
    results = pd.read_csv(args.i)
    # Find the threshold and plot the results
    d_p_t = find_threshold(results, args.p, args.u, file_name=args.s)

    p_swap = np.array(
        [0, 0.06, 0.12, 0.18, 0.24, 0.30, 0.36, 0.42, 0.48, 0.54, 0.6, 0.66, 0.70, 0.71]
    )
    sqz = np.array(
        [10.1, 10.4, 10.7, 11.1, 11.5, 11.8, 12.4, 12.9, 13.7, 14.9, 17, 24, 35, 60]
    )

    swap_tol_data = zip(sqz, p_swap)
    plot = plot_results(
        results,
        p_swap=args.p,
        distances=[7, 9, 11, 13],
        unit=args.u,
        threshold=d_p_t,
        file_name=args.s,
        show=False,
        swap_tol_plot=swap_tol_data,
        inset=True,
        break_axis=True,
    )
