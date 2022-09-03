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
"""Functions to plot the heat maps of input matrices.

Plots are configured via the ``plot_params`` dictionary. These parameters
are associated with Matplolib's rc settings and are modified within the
plotting functions using the ``rc_context`` context manager. This approach
avoids having to modify the global Matplotlib ``rc_params``.

To modify the plot parameters use, for example,

.. code-block:: python

    from flamingpy.viz.plot_mat_heat_map import plot_params as fp_plot_params
    fp_plot_params["font.size"] = 20
"""

# pylint: disable=too-many-statements,singleton-comparison, too-many-lines

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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
def plot_mat_heat_map(mat, show=True, title=None):
    """Plot the heat map of a matrix."""
    fig = plt.figure()
    ax = plt.gca()
    if not isinstance(mat, np.ndarray):
        mat = mat.toarray()
    plt.matshow(mat, 0)
    if title:
        plt.title(title)
    cbar = plt.colorbar()
    cbar.set_label(
        "value", rotation=270, fontsize=plot_params.get("axes.labelsize", 10) * 1.2, labelpad=20
    )
    if show:
        plt.show()

    axs = [ax, cbar.ax]
    return fig, axs
