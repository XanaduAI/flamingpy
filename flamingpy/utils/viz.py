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
"""Helper functions to draw various graphs and generate plots using
Matplotlib."""

# pylint: disable=too-many-statements

import itertools as it
from matplotlib import markers

import numpy as np
import networkx as nx
import matplotlib as mpl
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from flamingpy.codes import Stabilizer
from flamingpy.cv import gkp

plot_params = {
    "font.size": 12,
    "font.family": "sans-serif",
    "axes.labelsize": 13,
    "axes.titlesize": 17,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 13,
    "grid.color": "lightgray",
    "lines.markersize": 5,
    "lines.linewidth": 5,
    "figure.figsize": (8, 6),
}


@mpl.rc_context(plot_params)
def plot_integer_part(xs, ns, alpha, show=True):
    """Plot the integer part of real numbers mod alpha."""
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    fig = plt.figure()
    ax = plt.gca()
    plt.plot(xs, ns, ".")
    plt.title("Integer Part")
    plt.xlabel("$x$")
    plt.xticks(newxticks, newxlabels)
    plt.ylabel(r"$\mathrm{int}(x)$")
    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def plot_fractional_part(xs, fs, alpha, show=True):
    """Plot the fractional part of real numbers mod alpha."""
    plt.title("Fractional Part")
    plt.plot(xs, fs, ".", linewidth=10)
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    newylabels = ["{:.3f}".format(tick) for tick in newyticks[1:-1]]
    newylabels = [gkp.to_pi_string(-alpha / 2)] + newylabels + [gkp.to_pi_string(alpha / 2)]
    plt.xticks(newxticks, newxlabels)
    plt.xlabel("$x$")
    plt.yticks(newyticks, newylabels)
    plt.ylabel(r"$\mathrm{frac}(x)$")
    if show:
        plt.show()


@mpl.rc_context(plot_params)
def plot_GKP_bins(outcomes, bit_values, alpha, show=True):
    """Plot binned real numbers mod alpha."""
    xmin, xmax = alpha * (outcomes[0] // alpha), alpha * (outcomes[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    fig = plt.figure()
    ax = plt.gca()
    plt.plot(outcomes, bit_values, ".")
    plt.title("Binned values")
    plt.xticks(newxticks, newxlabels)
    plt.xlabel("Outcomes")
    plt.yticks([0, 1], [0, 1])
    plt.ylabel("Bit values")
    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def plot_Z_err_cond(hom_val, error, alpha, use_hom_val, show=True):
    """Plot conditional phase probabilities for GKP states."""
    _, frac = gkp.GKP_binner(hom_val, return_fraction=True)
    val = hom_val if use_hom_val else frac
    xmin, xmax = alpha * (hom_val[0] // alpha), alpha * (hom_val[-1] // alpha) + alpha
    print(xmin, xmax, min(val), max(val))
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    fig = plt.figure()
    ax = plt.gca()
    plt.plot(val, error, ".")
    plt.xlabel("Homodyne value")
    plt.ylabel("Error")
    addendum = "Full homodyne value" if use_hom_val else "Central peak"
    plt.title("Conditional phase probabilities: " + addendum)
    plt.xticks(newxticks, newxlabels)
    if show:
        plt.show()

    return fig, ax


@mpl.rc_context(plot_params)
def draw_EGraph(
    egraph,
    color_nodes=False,
    color_edges=False,
    state_colors=None,
    label=None,
    title=False,
    legend=False,
    display_axes=True,
):
    """Draw the graph state represented by the EGraph.

    Args:
        color_nodes (bool or string or dict): Options are:

            True: color the nodes based on the 'color' attribute
                attached to the node. If unavailable, color nodes black.
            'state': color nodes based on the 'state' attribute. Uses
                the color wheel by default, but colors can also be
                specified via a dictionary in the state_colors argument.
            string: color all nodes with the color specified by the stirng
            dict: color nodes based on the 'type' attribute of the node,
                with the dictionary specifying the colours. For example,
                if the type can be primal or dual, the dictionary should
                be of the form:

                    {"primal": primal_color, "dual": dual_color}.

        color_edges (bool):

            True: color the edges based on the 'color' attribute
                attached to the node. If unavailable, color nodes grey.
            string: color all edges with the color specified by the stirng
            dict: color edges based on the 'weight' attribute of the node,
                with the dictionary specifying the colours. For example,
                if the weight can be +1 or -1, the dictionary should
                be of the form:

                    {-1: minus_color, +1: plus_color}.

        label (NoneType or string): plot values next to each node
            associated with the node attribute label. For example,
            to plot bit values, set label to "bit_val". If set to 'index',
            it will plot the integer indices of the nodes. If the attribute
            for some or all of the nodes, a message will print indicating
            for how many nodes the attribute has not been set.
        title (bool): if True, display the title, depending on the label.
            For default labels, the titles are converted from attribute
            name to plane English and capitalized.
        legend (bool): if True and label is set to 'state', display
            the state color legend.
        display_axes (bool): if False, turn off the axes.

    Returns:
        A Matplotib Axes object.
    """
    if state_colors is None:
        state_colors = {}

    # Recommended to be viewed with IPython.
    # Font properties
    dims = egraph.graph.get("dims")
    xmax, ymax, zmax = dims

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if label:
        title_dict = {
            "p_phase": "Phase error probabilities",
            "p_phase_cond": "Conditional phase error probabilities",
            "hom_val_p": "p-homodyne outcomes",
            "hom_val_q": "q-homodyne outcomes",
            "bit_val": "Bit values",
            "weight": "Weights",
            "index": "Indices",
        }
        name = title_dict.get(label) if title_dict.get(label) else label
        n_uncomputed = 0
        if title:
            ax.set_title(name)
            ax.title.set_size(plot_params.get("axes.titlesize"))
        if label == "index":
            indices = egraph.index_generator()

    if color_nodes == "state":
        handles = []
        color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        i = 0
        for state in state_colors.keys():
            color = state_colors.get(state)
            if not color:
                color = color_cycle[i]
            line = mlines.Line2D([], [], color=color, marker=".", label=state)
            handles += [line]
            state_colors[state] = color
            i += 1

    # Plotting points. y and z are swapped in the loops so that
    # z goes into the page; however, the axes labels are correct.
    for point in egraph.nodes:
        x, z, y = point

        # Color nodes based on color_nodes if string, or based on
        # color attribute if True; black otherwise.
        if color_nodes == "state":
            color = state_colors.get(egraph.nodes[point].get("state"))
        elif isinstance(color_nodes, dict):
            color = color_nodes.get(egraph.nodes[point].get("type"))
        elif isinstance(color_nodes, str):
            color = color_nodes
        else:
            color = egraph.nodes[point].get("color") if color_nodes else "k"

        ax.scatter(x, y, z, c=color, s=plot_params.get("lines.markersize", 3))

        if label:
            value = egraph.nodes[point].get(label) if label != "index" else indices[point]
            if value is not None:
                x, z, y = point
                # Raise negative sign above node.
                sign = "^{-}" if value < 0 else " "
                if not isinstance(value, int):
                    value = r"${}{:.2g}$".format(sign, np.abs(value))
                ax.text(
                    x + 0.05,
                    y,
                    z,
                    value,
                    color="MediumBlue",
                    # backgroundcolor="w",
                    zorder=2,
                )
            else:
                n_uncomputed += 1

    if label and n_uncomputed > 0:
        message = "{} at {} node(s) have not yet been computed."
        print(message.format(name.lower(), n_uncomputed))

    # Plotting edges.
    for edge in egraph.edges:

        # Color edges based on color_edges if string, or based on
        # color attribute if True; black otherwise.
        if isinstance(color_edges, str):
            color = color_edges
        elif isinstance(color_edges, dict):
            color = color_edges.get(egraph.edges[edge].get("weight"))
        else:
            color = egraph.edges[edge].get("color") if color_edges else "grey"

        x1, z1, y1 = edge[0]
        x2, z2, y2 = edge[1]
        plt.plot([x1, x2], [y1, y2], [z1, z2], color=color, linewidth=1)

    if color_nodes == "state" and legend:
        ax.legend(handles=handles)

    plt.xticks(range(0, 2 * xmax + 1))
    plt.yticks(range(0, 2 * zmax + 1))
    ax.set_zticks(range(0, 2 * ymax + 1))
    ax.set_xlabel("x", labelpad=15)
    ax.set_ylabel("z", labelpad=15)
    ax.set_zlabel("y", labelpad=15)
    if not display_axes:
        ax.axis("off")
    plt.tight_layout(pad=5)
    plt.draw()
    return fig, ax


@mpl.rc_context(plot_params)
def plot_binary_mat_heat_map(mat, show=True, title=None):
    """Plot the heat map of a binary matrix."""
    fig = plt.figure()
    ax = plt.gca()
    if not isinstance(mat, np.ndarray):
        mat = mat.toarray()
    cmap = mpl.colors.ListedColormap(["C0", "C1"])
    plt.matshow(mat, 0, cmap=cmap)
    if title:
        plt.title(title)
    cbar = plt.colorbar(ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.set_label(
        "value", rotation=270, fontsize=plot_params.get("axes.labelsize", 10) * 1.2, labelpad=20
    )
    if show:
        plt.show()

    axs = [ax, cbar.ax]
    return fig, axs


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


@mpl.rc_context(plot_params)
def draw_dec_graph(graph, label_edges=False, node_labels=None, title=""):
    """Draw a stabilizer or matching graph with a color legend.

    This requires that the graph is implemented with the NetworkX backend.

    Args:
        graph (NxStabilizerGraph or NxMatchingGraph): the graph to draw.
        label_edges (bool, optional): if `True`, label the edges
            of the graph with their weight. Defaults to False.
        node_labels (dict of node to label, optional): if provided, the nodes
            will be identified with the given labels. Else, there will be no
            label for the nodes.
        title (string, optional): add the given title to the plot.
    """
    if not isinstance(graph.graph, nx.Graph):
        raise ValueError("The graph must be implemented with the networkx backend.")
    graph = graph.graph
    layout = nx.circular_layout(graph)

    fig, ax = plt.subplots()
    if title != "":
        ax.set_title(title)
    # NetworkX drawing function for circular embedding of graphs.

    # Color edges based on weight, and draw a colobar.

    # divider = make_axes_locatable(ax)
    ax.axis("off")
    cmap, norm = draw_curved_edges(graph, layout, ax)
    nx.draw_networkx_nodes(graph, pos=layout, node_color="#202020", ax=ax)

    # Draw node labels
    if node_labels is not None:
        node_labels = {node: label for node, label in node_labels.items() if node in graph}
        nx.draw_networkx_labels(
            graph, pos=layout, labels=node_labels, font_color="white", font_weight=100, ax=ax
        )

    if label_edges:
        weight_dict = {edge: "{:.2f}".format(graph.edges[edge]["weight"]) for edge in graph.edges}
        nx.draw_networkx_edge_labels(
            graph,
            layout,
            edge_labels=weight_dict,
            font_size=plot_params.get("font.size", 7) * 0.75,
            clip_on=False,
            alpha=0.7,
            verticalalignment="center_baseline",
            bbox={"alpha": 0},
        )

    cax, kw = mpl.colorbar.make_axes(ax, location="right", fraction=0.15)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, **kw)
    cbar.ax.tick_params(labelsize=plot_params.get("axes.labelsize", 10), rotation=270)
    cbar.set_label(
        "weight", rotation=270, fontsize=plot_params.get("axes.labelsize", 10), labelpad=20
    )

    axs = [ax, cax]
    return fig, axs


def draw_curved_edges(graph, layout, ax, rad=0.5):
    """Draw curved edges using matplotlib directly instead of networkx.

    This requires that the graph is implemented with the NetworkX backend.

    Args:
        graph (NxGraph): The graph to draw.
        layout (dict): A dictionary of positions keyed by node.
        ax (matplotlib.pyplot.Axis): The axis on which arrows should be drawn.
        rad (float, optional): Curvature of the arrows in radians.
    """

    edges = graph.edges
    edge_weights = [edges[edge]["weight"] for edge in edges]

    cmap = mpl.cm.get_cmap("Spectral")
    norm = mpl.colors.Normalize(vmin=np.min(edge_weights), vmax=np.max(edge_weights))
    for edge in graph.edges():
        source, target = edge
        arrowprops = dict(
            arrowstyle="-",
            color=cmap(norm(edges[edge]["weight"])),
            connectionstyle=f"arc3,rad={rad}",
            linestyle="-",
            linewidth=plot_params.get("lines.linewidth", 1) / 2,
            alpha=0.8,
        )
        ax.annotate("", xy=layout[source], xytext=layout[target], arrowprops=arrowprops)
    return cmap, norm


@mpl.rc_context(plot_params)
def syndrome_plot(code, ec, index_dict=None, drawing_opts=None):
    """Draw the syndrome plot for a code.

    A comprehensive graphing tool for drawing the error syndrome of code.
    Labelling options are specified with the help of
    drawing_opts, and can include:

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per draw_EGraph
                'label_cubes' -> indices of the stabilizers
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
        matplotlib.pyplot.axes: the 'axes' object
    """

    cubes = getattr(code, ec + "_stabilizers")
    # Default drawing options.
    draw_dict = {
        "show_nodes": False,
        "color_nodes": "state",
        "color_edges": "k",
        "label": None,
        "legend": True,
        "title": True,
        "state_colors": {"p": None, "GKP": None},
        "display_axes": True,
        "label_cubes": True,
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
            "state_colors",
            "label",
            "legend",
            "display_axes",
        ]
        egraph_opts = {k: drawing_opts[k] for k in egraph_args}
        fig, ax = draw_EGraph(code.graph, **egraph_opts)
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        plt.xticks(range(0, 2 * shape[0] + 1))
        plt.yticks(range(0, 2 * shape[1] + 1))
        ax.set_zticks(range(0, 2 * shape[2] + 1))
        ax.set_xlabel(
            "x",
            labelpad=20,
        )
        ax.set_ylabel(
            "z",
            labelpad=20,
        )
        ax.set_zlabel(
            "y",
            labelpad=20,
        )
        leg = None
    # Illustrate stabilizers with voxels colored green for even
    # parity and red for odd pariy.
    filled = np.zeros(shape, dtype=object)
    for cube in cubes:

        # Obtain smallest, largest, and middle coordinates for each
        # cube. Divided by 2 becaues voxels are 1X1X1.
        xmin, xmax = np.array(cube.xlims(), dtype=int) // 2
        ymin, ymax = np.array(cube.ylims(), dtype=int) // 2
        zmin, zmax = np.array(cube.zlims(), dtype=int) // 2
        xmid, ymid, zmid = np.array(cube.midpoint())
        # Fill in the color arrays depending on parity.
        if cube.parity:
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = "#FF000015"
        else:
            filled[xmin:xmax, ymin:ymax, zmin:zmax] = "#00FF0015"
        if drawing_opts["label_cubes"] and index_dict:
            if cube in index_dict:
                ax.text(
                    xmid,
                    ymid,
                    zmid,
                    index_dict[cube],
                    # fontdict=font_props
                )

    # This portion adapted from a Matplotlib official example to fix
    # an issue with filling in the insides of voxels: the code
    # expands the indices and creates small gaps between the voxels.

    def explode(data):
        size = np.array(data.shape) * 2
        data_e = np.zeros(size - 1, dtype=data.dtype)
        data_e[::2, ::2, ::2] = data
        return data_e

    # upscale the above voxel image, leaving gaps
    filled_e = explode(filled)
    # Shrink the gaps
    x, y, z = np.indices(np.array(filled_e.shape) + 1, dtype=float)
    x[0::2, :, :] += 0.05
    y[:, 0::2, :] += 0.05
    z[:, :, 0::2] += 0.05
    x[1::2, :, :] += 0.95
    y[:, 1::2, :] += 0.95
    z[:, :, 1::2] += 0.95
    ax.voxels(x, y, z, filled_e, facecolors=filled_e)

    if drawing_opts["label_boundary"]:
        bound_points = getattr(code, ec + "_bound_points")
        for point in bound_points:
            ax.scatter(point[0], point[1], point[2], s=70, c="k")
            ax.text(
                point[0],
                point[1],
                point[2],
                index_dict[point],
                # fontdict=font_props
            )

    # Define a legend for red/green cubes.
    legend_elements = [
        Patch(facecolor="#00FF0050", label="even parity"),
        Patch(facecolor="#FF000050", label="odd parity"),
    ]
    if drawing_opts["legend"]:
        ax.legend(
            handles=legend_elements,
            # prop=font_props,
            loc="upper left",
        )
    # Since draw_EGraph legend has been overwritten, re-add
    # it to the plot.
    if leg:
        ax.add_artist(leg)
    if drawing_opts["title"]:
        ax.set_title(ec.capitalize() + " syndrome")

    ax.autoscale()
    return fig, ax


@mpl.rc_context(plot_params)
def draw_matching_on_syndrome_plot(ax, matching, G_match):
    """Plot the matching output by MWPM."""
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
                    plt.plot(x, y, z, marker="2", c="k")
                xlist += [x]
                ylist += [y]
                zlist += [z]
            ax.set_title("Minimum-weight perfect matching")
            ax.plot(xlist, ylist, zlist, "o-", c=np.random.rand(3))
    return ax
