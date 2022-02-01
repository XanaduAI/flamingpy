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
""" 
Helper functions to draw various graphs 
and generate plots using Matplotlib.
"""
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import copy
import networkx as nx
import itertools as it

from matplotlib.patches import Patch
from flamingpy import GKP, RHG


def plot_integer_fractional(xs, ns, fs, alpha):
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(xs, ns, ",")
    plt.title("Integer Part", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.show()

    plt.title("Fractional Part", fontsize="medium")
    plt.plot(xs, fs, ",")
    newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
    newylabels = ["{:.3f}".format(tick) for tick in newyticks[1:-1]]
    newylabels = [GKP.to_pi_string(-alpha / 2)] + newylabels + [GKP.to_pi_string(alpha / 2)]
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks(newyticks, newylabels)
    plt.show()


def plot_GKP_bins(outcomes, bit_values, alpha):
    xmin, xmax = alpha * (outcomes[0] // alpha), alpha * (outcomes[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(outcomes, bit_values, ",")
    plt.title("Binned values", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks([0, 1], [0, 1])
    plt.show()


def plot_Z_err_cond(hom_val, error, alpha, use_hom_val):
    _, frac = GKP.GKP_binner(hom_val, return_fraction=True)
    val = hom_val if use_hom_val else frac
    xmin, xmax = alpha * (hom_val[0] // alpha), alpha * (hom_val[-1] // alpha) + alpha
    print(xmin, xmax, min(val), max(val))
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(val, error, ",")
    addendum = "Full homodyne value" if use_hom_val else "Central peak"
    plt.title("Conditional phase probabilities: " + addendum, fontsize="small")
    plt.xticks(newxticks, newxlabels, fontsize="small")
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
    """Draw the graph.

    Args:
        color_nodes (bool): If True, color nodes based on 'color'
            attributes attached to the node. Black by default.
        color_edges (bool): If True, color edges based on 'color'
            attributes attached to the node. Grey by default.
        label (NoneType): ...
        display_axes (bool): if False, turn off the axes.

    Returns:
        A matplotib Axes object.
    """
    # Recommended to be viewed with IPython.
    # Font properties
    dims = egraph.graph.get("dims")
    if dims:
        # TODO: Store dims as EGraph attributes, rather than a graph
        # attribute?
        font_size = 10 * sum(dims) ** (1 / 2)
    else:
        # TODO: If dims not specified find x, y, z limits of graph,
        # supposing the graph is filled. Alternatively, just change
        # the figure size?
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
            "indices": "Indices",
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
        elif type(color_nodes) == str:
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
                if type(value) != int:
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
        # color ttribute if True; black otherwise.
        if type(color_edges) == str:
            color = color_edges
        else:
            color = egraph.edges[edge].get("color") if color_edges else "k"

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


def plot_binary_mat_heat_map(symplectic):
    plt.figure()
    if type(symplectic) != np.ndarray:
        symplectic = symplectic.toarray()
    plt.matshow(symplectic, 0)
    plt.show()


def draw_code_lattice(
    graph,
    node_color={"primal": "k", "dual": "grey"},
    edge_color={1: "b", -1: "r"},
):
    graph = copy.deepcopy(graph)
    for (_, attr) in graph.nodes.items():
        attr["color"] = node_color[attr["type"]]
    for (_, attr) in graph.edges.items():
        attr["color"] = edge_color[attr["weight"]]
    draw_EGraph(graph, color_nodes=True, color_edges=True)


def draw_dec_graph(G, label_edges=True, title=None):
    """Draw decoding and matching graphs G with a color legend."""
    if title is None:
        try:
            title = G["title"]
        except Exception:
            title = ""
    plt.figure()
    plt.title(title, family="serif", size=10)
    # NetworkX drawing function for circular embedding of graphs.
    nx.draw_circular(
        G,
        edgelist=[],
        with_labels=True,
        node_color="k",
        font_size=7,
        font_color="w",
        font_family="serif",
    )
    # Color edges based on weight, and draw a colobar.
    weight_list = [G.edges[edge]["weight"] for edge in G.edges]
    weight_dict = {edge: "{:.2f}".format(G.edges[edge]["weight"]) for edge in G.edges}
    if label_edges:
        nx.draw_networkx_edge_labels(G, nx.circular_layout(G), edge_labels=weight_dict, font_size=7)
    r = nx.draw_networkx_edges(G, nx.circular_layout(G), edge_color=weight_list)
    plt.colorbar(r)


def syndrome_plot(code, G_dec, index_dict=None, drawing_opts=None):
    """Draw the syndrome plot for a CVGraph G.

    A comprehensive graphing tool for drawing the error syndrome of
    a CVGraph G. Labelling options are specified with the help of
    drawing_opts, and can include:

                'show_nodes' -> the underlying graph displayed
                'label_nodes' -> node labels, as per CVGraph.draw
                'label_cubes' -> indices of the stabilizers
                'label_boundary'-> indices of the boundary points
                'legend' -> legends for the nodes and cubes

    Cubes are shown as transparent voxels, green for even parity and
    red for odd. For now, cubes on periodic boundaries are not shown,
    and stabilizers on dual boundaries occupy are shown to occupy
    the full space of a six-body cube.

    Args:
        G_dec (networkx.Graph): the decoding graph
        G (CVGraph): the encoded CVGraph
        index_dict (dict): the stabiizer-to-index mapping
        drawing_opts (dict): a dictionary of drawing options, with
            all possibilities described above.

    Returns:
        matplotlib.pyplot.axes: the 'axes' object
    """
    # Font properties
    # TODO: Make consistent with EGraph fontprops
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

    bound_inds = G_dec.graph["boundary_points"]
    points = [G_dec.nodes[i]["stabilizer"] for i in bound_inds]

    cubes = code.stabilizers
    # Default drawing options.
    draw_dict = {
        "show_nodes": False,
        "color_nodes": "state",
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
    # CVGraph.sketch (this also plots the graph in the console).
    if drawing_opts["show_nodes"]:
        # TODO: If draw method moved out of CVGraph and into EGraph,
        # the state argument would be unnecessary here.
        egraph_args = [
            "color_nodes",
            "label",
            "legend",
            "title",
            "state_colors",
            "display_axes",
        ]
        egraph_opts = {k: drawing_opts[k] for k in egraph_args}
        ax = draw_EGraph(code.graph, **egraph_opts)
        leg = ax.get_legend()
    # If show_nodes is False, create a new figure with size
    # determined by the dimensions of the lattice.
    else:
        # TODO: Initialize axes based on empty ax object from state.draw()
        # but prevent from state.draw() from plotting.
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
        for point in points:
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
    # Since CVGraph.sketch() legend has been overwritten, readd
    # it to the plot.
    if leg:
        ax.add_artist(leg)
    return ax


def draw_matching_on_syndrome_plot(ax, matching, G_dec, G_match, label_edges):
    virtual_points = G_match.virtual_points
    for pair in matching:
        if pair not in it.product(virtual_points, virtual_points):
            xlist, ylist, zlist = [], [], []
            path = G_match.edge_path(pair)
            for node in path:
                stabe = G_dec.nodes[node]["stabilizer"]
                if isinstance(stabe, RHG.RHGCube):
                    x, y, z = stabe.midpoint()
                else:
                    x, y, z = stabe
                    plt.plot(x, y, z, marker="2", ms=50, c="k")
                xlist += [x]
                ylist += [y]
                zlist += [z]
            ax.set_title("Minimum-weight perfect matching", family="serif", size=20)
            ax.plot(xlist, ylist, zlist, "o-k", ms=20, linewidth=5, c=np.random.rand(3))
    return ax
