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
"""A series of functions to visualize Error Correction process and results.

In most cases, plots are configured via the ``plot_params`` dictionary.
These parameters are associated with Matplolib's rc settings and are
modified within the plotting functions using the ``rc_context`` context
manager. This approach avoids having to modify the global Matplotlib
``rc_params``.

To modify the plot parameters use, for example,

.. code-block:: python

    from flamingpy.viz.syndrome_plot import plot_params as fp_plot_params
    fp_plot_params["font.size"] = 20
"""

# pylint: disable=too-many-statements,singleton-comparison, too-many-lines

import itertools as it

import numpy as np
import matplotlib as mpl
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from flamingpy.codes import Stabilizer
from flamingpy.viz.GraphStates import draw_EGraph_matplotlib
from flamingpy.viz.cube_helper import _plot_cubes_at

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
def syndrome_plot(code, ec, index_dict=None, drawing_opts=None):
    """Draw the syndrome plot for a code.

    A comprehensive graphing tool for drawing the error syndrome of code.
    Labelling options are specified with the help of
    drawing_opts, and can include:

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per draw_EGraph
                'label_stabilizers' -> indices of the stabilizers
                'label_boundary'-> indices of the boundary points
                'legend' -> legends for the nodes and stabilizers

    Stabilizers are shown as transparent voxels, green for even parity and
    red for odd. For now, stabilizers on periodic boundaries are not
    drawn in a special way, stabilizers on dual boundaries are unshifted
    from the primal stabilizer location, and incomplete stabilizers
    are still represented as complete cubes.

    Args:
        code (SurfaceCode): the qubit QEC code
        ec (string): the error complex ('primal' or 'dual')
        index_dict (dict): the stabiizer-to-index mapping
        drawing_opts (dict): a dictionary of drawing options, with
            all possibilities described above.

    Returns:
        tuple: figure and axes
    """

    cubes = getattr(code, ec + "_stabilizers")
    # Default drawing options.
    draw_dict = {
        "show_nodes": False,
        "color_nodes": ("state", {"p": None, "GKP": None}),
        "color_edges": "k",
        "label": None,
        "legend": True,
        "show_title": True,
        "show_axes": True,
        "label_stabilizers": True,
        "label_boundary": False,
    }
    # Combine default dictionary with supplied dictionary, duplicates
    # favor supplied dictionary.
    if drawing_opts is None:
        drawing_opts = {}
    drawing_opts = {**draw_dict, **drawing_opts}

    # Shape and font properties from the original graph.
    shape = np.array(code.dims)
    # If show_nodes is True, get the axes object and legend from
    # draw_EGraph (this also plots the graph in the console).
    if drawing_opts["show_nodes"]:
        egraph_args = [
            "color_nodes",
            "color_edges",
            "label",
            "legend",
            "show_axes",
        ]
        egraph_opts = {k: drawing_opts[k] for k in egraph_args}
        fig, ax = draw_EGraph_matplotlib(code.graph, **egraph_opts)
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        plt.xticks(range(0, 2 * shape[0] - 1))
        plt.yticks(range(0, 2 * shape[1] - 1))
        ax.set_zticks(range(0, 2 * shape[2] - 1))
        ax.set_xlabel("x", labelpad=20)
        ax.set_ylabel("z", labelpad=20)
        ax.set_zlabel("y", labelpad=20)
        leg = None
    # Illustrate stabilizers with voxels colored green for even
    # parity and red for odd pariy.
    positions, colors, sizes = [], [], []
    for cube in cubes:

        # Obtain smallest, largest, and middle coordinates for each
        # cube.
        xmin, xmax = np.array(cube.xlims())
        ymin, ymax = np.array(cube.ylims())
        zmin, zmax = np.array(cube.zlims())
        xmid, ymid, zmid = np.array(cube.midpoint())
        # Fill in the color arrays depending on parity.
        if cube.parity:
            color = "#FF000015"
        else:
            color = "#00FF0015"

        # gap defines the distance between adjacent cubes
        gap = 0.15
        min_arr = np.array([xmin, ymin, zmin])
        max_arr = np.array([xmax, ymax, zmax])
        positions.append(min_arr + gap)
        sizes.append(np.abs(max_arr - min_arr) - 2 * gap)
        colors.append(color)

        if drawing_opts["label_stabilizers"] and index_dict:

            if cube in index_dict:
                ax.text(xmid, ymid, zmid, index_dict[cube])

    # draw cubes
    pc = _plot_cubes_at(positions, colors=colors, sizes=sizes)
    ax.add_collection3d(pc)

    # setting plot limits to give some room to the boxes
    ymin = 0 if ec == "primal" else -2
    ax.set_xlim(0, 2 * shape[0])
    ax.set_ylim(ymin, 2 * shape[1])
    ax.set_zlim(0, 2 * shape[2])

    if drawing_opts["label_boundary"] and index_dict:
        bound_points = getattr(code, ec + "_bound_points")
        for point in bound_points:
            ax.scatter(*point, s=5, c="k")
            ax.text(*point, index_dict[point])

    # Define a legend for red/green cubes.
    legend_elements = [
        Patch(facecolor="#00FF0050", label="even parity"),
        Patch(facecolor="#FF000050", label="odd parity"),
    ]
    if drawing_opts["legend"]:
        ax.legend(handles=legend_elements, loc="upper left")
    # Since draw_EGraph legend has been overwritten, re-add
    # it to the plot.
    if leg:
        ax.add_artist(leg)
    if drawing_opts["show_title"]:
        ax.set_title(ec.capitalize() + " syndrome")

    return fig, ax


@mpl.rc_context(plot_params)
def add_recovery_drawing(ax, **kwargs):
    """Plot the recovery."""
    if kwargs.get("show_title"):
        ax.set_title("Syndrome and recovery")
    recovery_set = kwargs.get("recovery_set")
    if recovery_set:
        for point in recovery_set:
            ax.plot(*point, "o", c="k", markersize=6)
    matching = kwargs.get("matching")
    G_match = kwargs.get("matching_graph")
    if matching:
        virtual_points = G_match.virtual_points
        for pair in matching:
            if pair not in it.product(virtual_points, virtual_points):
                xlist, ylist, zlist = [], [], []
                path = G_match.edge_path(pair)
                for node in path:
                    if isinstance(node, Stabilizer):
                        x, y, z = node.midpoint()
                    else:
                        x, y, z = node
                        plt.plot(x, y, z, marker="2", markersize=15, c="k")
                    xlist += [x]
                    ylist += [y]
                    zlist += [z]
                ax.plot(
                    xlist,
                    ylist,
                    zlist,
                    "o-",
                    c=np.random.rand(3),
                    linewidth=plot_params.get("lines.linewidth", None) * 0.9,
                )

    return ax


def draw_decoding(code, ec, dec_objects=None, drawing_opts=None):
    """Draw the stabilizer and matching graphs, and the plot the syndrome."""
    if drawing_opts is None:
        drawing_opts = {}

    # Drawing the stabilizer graph
    G_stabilizer = getattr(code, ec + "_stab_graph")
    G_match = dec_objects.get("matching_graph")
    # An integer label for each node in the stabilizer and matching
    # graphs. This is useful to identify the nodes in the plots.
    if drawing_opts.get("label_stabilizers") or drawing_opts.get("label_boundary"):
        # Node labels for the stabilizer graph (avoid "low" / "high" nodes)
        node_labels = {node: index for index, node in enumerate(list(G_stabilizer.nodes())[2:])}
        # Update node labels to work with the matching graph---needs to be done
        # because virtual boundary nodes are of the form ((x, y, z), i).
        if G_match:
            for virtual_node in set(G_match.graph.nodes()) - set(G_stabilizer.nodes()):
                index = node_labels[virtual_node[0]]
                node_labels[virtual_node] = index
    else:
        node_labels = None
    label_edges = bool(drawing_opts.get("label_edges"))
    show_title = bool(drawing_opts.get("show_title"))
    # title = drawing_opts.get_title()
    fig1, ax1 = code.draw_stabilizer_graph(
        ec,
        title=ec.capitalize() + " stabilizer graph" if show_title else "",
        label_edges=label_edges,
        node_labels=node_labels,
    )

    # Drawing the matching graph
    fig2, ax2 = (None, None)
    if G_match:
        if len(G_match.graph):
            fig2, ax2 = G_match.draw(
                title=ec.capitalize() + " matching graph" if show_title else "",
                label_edges=label_edges,
                node_labels=node_labels,
            )
        else:
            print("\nMatching graph empty!\n")

    # Drawing the syndrome
    fig3, ax3 = syndrome_plot(code, ec, drawing_opts=drawing_opts, index_dict=node_labels)
    if drawing_opts.get("show_recovery"):
        ax3 = add_recovery_drawing(ax3, show_title=drawing_opts.get("show_title"), **dec_objects)

    return (fig1, ax1), (fig2, ax2), (fig3, ax3)
