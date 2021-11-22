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
""" 
Helper functions to draw various graphs 
and generate plots using Matplotlib.
"""
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import ft_stack.GKP as GKP
import copy


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
    newylabels = (
        [GKP.to_pi_string(-alpha / 2)] + newylabels + [GKP.to_pi_string(alpha / 2)]
    )
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


def plt_Z_err_cond(hom_val, error, alpha, use_hom_val):
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
            value = (
                egraph.nodes[point].get(label) if label != "index" else indices[point]
            )
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


def draw_CVGraph(self, **kwargs):
    """Draw the underlying graph with state colors.

    Run draw_EGraph with state information. State colors optionally
    supplied using state_colors argument; otherwise they are
    determined by the default color cycle.
    """
    default_args = {
        "color_nodes": "state",
        "legend": True,
        "title": True,
        "state_colors": {state: None for state in self._states},
    }
    kwargs = {**default_args, **kwargs}
    draw_EGraph(self.egraph, **kwargs)


def plot_SCZ_mat_heat_map(symplectic):
    print("The symplectic CZ matrix (dark spots 0, bright spots 1):")
    plt.figure()
    if type(symplectic) != np.ndarray:
        symplectic = symplectic.toarray()
    plt.matshow(symplectic, 0)
    plt.show()


def draw_RHG_graph(
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
