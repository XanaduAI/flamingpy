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
"""Helper functions to draw and manipulate graphs that are implemented with the
NetworkX backend."""

# pylint: disable=too-many-statements,singleton-comparison,too-many-lines

import numpy as np
import networkx as nx
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
def draw_dec_graph(graph, label_edges=False, node_labels=None, title=""):
    """Draw a stabilizer or matching graph with a color legend.

    This requires that the graph is implemented with the NetworkX backend.

    Plots are configured via the ``plot_params`` dictionary. These parameters
    are associated with Matplolib's rc settings and are modified within the
    plotting functions using the ``rc_context`` context manager. This approach
    avoids having to modify the global Matplotlib ``rc_params``.

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

    # Remove 'low' and 'high' nodes, which are not important for visualization.
    if graph.__class__.__name__ == "NxStabilizerGraph":
        graph.graph.remove_nodes_from({"low", "high"})

    graph = graph.graph
    layout = nx.circular_layout(graph)

    fig, ax = plt.subplots()
    if title != "":
        ax.set_title(title)

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
