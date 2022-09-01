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
"""Helper functions to visualize real-numbers-mode-alpha and GKP states.

Plots are configured via the ``plot_params`` dictionary. These parameters
are associated with Matplolib's rc settings and are modified within the
plotting functions using the ``rc_context`` context manager. This approach
avoids having to modify the global Matplotlib ``rc_params``.

To modify the plot parameters use, for example,

.. code-block:: python

    from flamingpy.viz.GKP_plotter import plot_params as fp_plot_params
    fp_plot_params["font.size"] = 20
"""

# pylint: disable=too-many-statements,singleton-comparison,too-many-lines

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from flamingpy.cv import gkp

plot_params = {
    "font.size": 10,
    "font.family": "serif",
    "axes.labelsize": 11,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "grid.color": "lightgray",
    "lines.markersize": 5,
    "lines.linewidth": 4,
    "figure.figsize": (8, 6),
}


@mpl.rc_context(plot_params)
def plot_integer_part(xs, ns, alpha, show=True):
    """Plot the integer part of real numbers mod alpha."""
    fig = plt.figure()
    ax = plt.gca()

    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    ax.xaxis.set_major_formatter(PiFormatter())

    plt.plot(xs, ns, ".")
    plt.title("Integer Part")
    plt.xlabel("$x$")
    plt.xticks(newxticks)
    plt.ylabel(r"$\mathrm{int}(x)$")

    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def plot_fractional_part(xs, fs, alpha, show=True):
    """Plot the fractional part of real numbers mod alpha."""
    fig = plt.figure()
    ax = plt.gca()

    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
    ax.xaxis.set_major_formatter(PiFormatter())
    ax.yaxis.set_major_formatter(PiFormatter())

    plt.plot(xs, fs, ".")
    plt.title("Fractional Part")
    plt.xticks(newxticks)
    plt.xlabel("$x$")
    plt.yticks(newyticks)
    plt.ylabel(r"$\mathrm{frac}(x)$")

    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def plot_GKP_bins(outcomes, bit_values, alpha, show=True):
    """Plot binned real numbers mod alpha."""
    fig = plt.figure()
    ax = plt.gca()

    xmin, xmax = alpha * (outcomes[0] // alpha), alpha * (outcomes[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    ax.xaxis.set_major_formatter(PiFormatter())

    plt.plot(outcomes, bit_values, ".")
    plt.title("Binned values")
    plt.xticks(newxticks)
    plt.xlabel("Outcomes")
    plt.yticks([0, 1], [0, 1])
    plt.ylabel("Bit values")

    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def plot_Z_err_cond(hom_val, error, alpha, use_hom_val, show=True):
    """Plot conditional phase probabilities for GKP states."""
    fig = plt.figure()
    ax = plt.gca()

    _, frac = gkp.GKP_binner(hom_val, return_fraction=True)
    val = hom_val if use_hom_val else frac
    # bounds for the plot
    if use_hom_val:
        xmin, xmax = alpha * np.array([hom_val[0] // alpha, hom_val[-1] // alpha + 1])
    else:
        xmin, xmax = -alpha / 2, alpha / 2

    print(xmin, xmax, min(val), max(val))

    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    ax.xaxis.set_major_formatter(PiFormatter())

    plt.plot(val, error, ".")
    plt.xticks(newxticks)
    plt.xlabel("Homodyne value")
    plt.ylabel("Error")
    plt.title(
        "Conditional phase probabilities: "
        + ("Full homodyne value" if use_hom_val else "Central peak")
    )

    if show:
        plt.show()

    return fig, ax
