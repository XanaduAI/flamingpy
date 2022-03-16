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
import itertools as it

import numpy as np
import networkx as nx
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from flamingpy.codes import Stabilizer
from flamingpy.cv import gkp


def plot_integer_part(xs, ns, fs, alpha, show=True):
    """Plot the integer part of real numbers mod alpha."""
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    plt.plot(xs, ns, ",")
    plt.title("Integer Part", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    if show:
        plt.show()


def plot_fractional_part(xs, ns, fs, alpha, show=True):
    """Plot the fractional part of real numbers mod alpha."""
    plt.title("Fractional Part", fontsize="medium")
    plt.plot(xs, fs, ",")
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    newylabels = ["{:.3f}".format(tick) for tick in newyticks[1:-1]]
    newylabels = [gkp.to_pi_string(-alpha / 2)] + newylabels + [gkp.to_pi_string(alpha / 2)]
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks(newyticks, newylabels)
    if show:
        plt.show()


def plot_GKP_bins(outcomes, bit_values, alpha, show=True):
    """Plot binned real numbers mod alpha."""
    xmin, xmax = alpha * (outcomes[0] // alpha), alpha * (outcomes[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    plt.plot(outcomes, bit_values, ",")
    plt.title("Binned values", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks([0, 1], [0, 1])
    if show:
        plt.show()


def plot_Z_err_cond(hom_val, error, alpha, use_hom_val, show=True):
    """Plot conditional phase probabilities for GKP states."""
    _, frac = gkp.GKP_binner(hom_val, return_fraction=True)
    val = hom_val if use_hom_val else frac
    xmin, xmax = alpha * (hom_val[0] // alpha), alpha * (hom_val[-1] // alpha) + alpha
    print(xmin, xmax, min(val), max(val))
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [gkp.to_pi_string(tick) for tick in newxticks]
    plt.plot(val, error, ",")
    addendum = "Full homodyne value" if use_hom_val else "Central peak"
    plt.title("Conditional phase probabilities: " + addendum, fontsize="small")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    if show:
        plt.show()


def draw_EGraph(
    egraph,
    color_nodes=False,
    color_edges=False,
    state_colors={},
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
    # Recommended to be viewed with IPython.
    # Font properties
    dims = egraph.graph.get("dims")
    # TODO: automate/improve the following figure size designation
    if dims:
        font_size = 10 * sum(dims) ** (1 / 2)
    else:
        dims = (5, 5, 5)
        font_size = 14
    xmax, ymax, zmax = dims
    # Set plotting options
    plot_params = {
        "font.size": font_size,
        "font.family": "serif",
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "grid.color": "lightgray",
        "lines.markersize": font_size,
    }
    plt.rcParams.update(plot_params)

    fig = plt.figure(figsize=((2 * (sum(dims) + 2), 2 * (sum(dims) + 2))))
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

        ax.scatter(x, y, z, c=color, s=plt.rcParams["lines.markersize"] * 5)

        if label:
            value = egraph.nodes[point].get(label) if label != "index" else indices[point]
            if value is not None:
                x, z, y = point
                # Raise negative sign above node.
                sign = "^{-}" * (-int(np.sign(value)))
                if not isinstance(value, int):
                    value = r"${}{:.2g}$".format(sign, np.abs(value))
                ax.text(
                    x,
                    y,
                    z,
                    value,
                    color="MediumBlue",
                    backgroundcolor="w",
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
        plt.plot([x1, x2], [y1, y2], [z1, z2], color=color)

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
    return ax


def plot_binary_mat_heat_map(mat, show=True):
    """Plot the heat map of a matrix."""
    plt.figure()
    if not isinstance(mat, np.ndarray):
        mat = mat.toarray()
    plt.matshow(mat, 0)
    if show:
        plt.show()


def draw_dec_graph(graph, label_edges=True, node_labels=None, title=""):
    """Draw a stabilizer or matching graph with a color legend.

    This requires that the graph is implemented with the NetworkX backend.

    Args:
        graph (NxStabilizerGraph or NxMatchingGraph): the graph to draw.
        label_edges (bool, optional): if True (the default), label the edges
            of the graph with their weight.
        node_labels (dict of node to label, optional): if provided, the nodes
            will be identified with the given labels. Else, there will be no
            label for the nodes.
        title (string, optional): add the given title to the plot.
    """
    if not isinstance(graph.graph, nx.Graph):
        raise ValueError("The graph must be implemented with the networkx backend.")
    graph = graph.graph
    plt.figure()
    if title != "":
        plt.title(title, family="serif", size=10)
    # NetworkX drawing function for circular embedding of graphs.
    if node_labels is not None:
        node_labels = {node: label for node, label in node_labels.items() if node in graph}
    nx.draw_circular(
        graph,
        edgelist=[],
        with_labels=node_labels is not None,
        labels=node_labels,
        node_color="k",
        font_size=7,
        font_color="w",
        font_family="serif",
    )
    # Color edges based on weight, and draw a colobar.
    weight_list = [graph.edges[edge]["weight"] for edge in graph.edges]
    weight_dict = {edge: "{:.2f}".format(graph.edges[edge]["weight"]) for edge in graph.edges}
    if label_edges:
        nx.draw_networkx_edge_labels(
            graph, nx.circular_layout(graph), edge_labels=weight_dict, font_size=7
        )
    r = nx.draw_networkx_edges(graph, nx.circular_layout(graph), edge_color=weight_list)
    cbar = plt.colorbar(r)
    cbar.ax.tick_params(labelsize=10)


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
    # Font properties
    font_size = 10 * sum(code.dims) ** (1 / 2)
    # Set plotting options
    plot_params = {
        "font.size": font_size,
        "font.family": "serif",
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "grid.color": "lightgray",
        "lines.markersize": font_size,
    }
    plt.rcParams.update(plot_params)

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
    # font_props = state.font_props
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
        ax = draw_EGraph(code.graph, **egraph_opts)
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        # TODO: Initialize axes based on empty ax object from draw_EGraph
        # but prevent from draw_EGraph from plotting.
        fig = plt.figure(figsize=(2 * (np.sum(shape) + 2), 2 * (np.sum(shape) + 2)))
        ax = fig.gca(projection="3d")
        # ax.tick_params(labelsize=font_props['size'])
        plt.xticks(range(0, 2 * shape[0] + 1))
        plt.yticks(range(0, 2 * shape[1] + 1))
        ax.set_zticks(range(0, 2 * shape[2] + 1))
        ax.set_xlabel(
            "x",
            # fontdict=font_props,
            labelpad=20,
        )
        ax.set_ylabel(
            "z",
            # fontdict=font_props,
            labelpad=20,
        )
        ax.set_zlabel(
            "y",
            # fontdict=font_props,
            labelpad=20,
        )
        plt.rcParams["grid.color"] = "lightgray"
        leg = None
    # Illustrate stabilizers with voxels colored green for even
    # parity and red for odd pariy.
    filled = np.zeros(shape, dtype=object)
    for cube in cubes:

        # TODO: Deal appropriately with cubes on periodic and dual
        # boundaries.

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
    return ax


def draw_matching_on_syndrome_plot(ax, matching, G_dec, G_match, label_edges):
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
                    plt.plot(x, y, z, marker="2", ms=50, c="k")
                xlist += [x]
                ylist += [y]
                zlist += [z]
            ax.set_title("Minimum-weight perfect matching", family="serif", size=20)
            ax.plot(xlist, ylist, zlist, "o-", ms=20, linewidth=5, c=np.random.rand(3))
    return ax
