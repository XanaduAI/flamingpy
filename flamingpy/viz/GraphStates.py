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
"""Helper functions to draw EGraphs and Graph States represented by them using
various backends."""

# pylint: disable=too-many-statements,singleton-comparison,too-many-lines

import numpy as np
import plotly.graph_objects as go
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from flamingpy.viz.EGraph_basics import _get_title, _get_node_info, _get_node_color, _get_edge_color

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
def draw_EGraph_matplotlib(
    egraph,
    color_nodes=False,
    color_edges=False,
    label=None,
    title=None,
    legend=False,
    show_axes=True,
    **kwargs,
):
    """Draw the graph state represented by the EGraph.

    Plots are configured via the ``plot_params`` dictionary. These parameters
    are associated with Matplolib's rc settings and are modified within the
    plotting functions using the ``rc_context`` context manager. This approach
    avoids having to modify the global Matplotlib ``rc_params``.

    Args:
        See draw_EGraph for general keyword arguments, see keyword arguments below for
            matplotlib-specific arguments.

    Keyword args:
        dimensions (tuple): Dimensions of the region that should be plotted.
            Should be of the form:

                ``([xmin, xmax], [ymin, ymax], [zmin, zmax])``

            If None, sets the dimensions to the smallest rectangular space
            containing all the nodes.

    Returns:
        tuple: Matplotib Figure and Axes.
    """

    dimensions = kwargs.get("dimensions", None)

    if dimensions is None:
        mins = map(min, zip(*egraph.nodes))
        maxs = map(max, zip(*egraph.nodes))
        mins = map(lambda x: int(np.floor(x)), mins)
        maxs = map(lambda x: int(np.ceil(x)), maxs)
        dimensions = zip(mins, maxs)

    xlim, ylim, zlim = [list(lim) for lim in dimensions]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    title = _get_title(title, label)
    if title:
        ax.set_title(title)
        ax.title.set_size(plot_params.get("axes.titlesize"))

    # plot graph
    ax = _plot_EGraph_nodes(ax, egraph, color_nodes, label, title, legend)
    ax = _plot_EGraph_edges(ax, egraph, color_edges)

    # plot generalities
    plt.xticks(range(xlim[0], xlim[1] + 1))
    plt.yticks(range(zlim[0], zlim[1] + 1))
    ax.set_zticks(range(ylim[0], ylim[1] + 1))

    for lim in [xlim, ylim, zlim]:
        if lim[0] == lim[1]:
            lim[0] -= 1
            lim[1] += 1

    plt.xlim(xlim)
    plt.ylim(zlim)
    ax.set_zlim(ylim)

    ax.set_xlabel("x", labelpad=15)
    ax.set_ylabel("z", labelpad=15)
    ax.set_zlabel("y", labelpad=15)
    if not show_axes:
        ax.axis("off")
    plt.tight_layout(pad=5)
    plt.draw()

    return fig, ax


def draw_EGraph_plotly(
    egraph,
    color_nodes=False,
    color_edges=False,
    label="all",
    title=None,
    legend=False,
    show_axes=True,
    **kwargs,
):
    """Draw the graph state represented by the EGraph with plotly.

    NOTE: Plotly gives you a wide range of options for how and where to
    display your figures. Refer to the `Plotly documentation
    <https://plotly.com/python/renderers/>`_ for more information.

    Args:
        See draw_EGraph for general keyword arguments, see keyword arguments below for
            matplotlib-specific arguments.

    Keyword args:
        showbackground (bool): if True, shows the background of the graph. Default is False.
        showspikes (bool): if True, shows coordinate-lines when hovering over a node. Default is
            False.
        showgrid (bool): if True, shows the grid of the graph. Default is True.
        width (int): width of the graph. Default is 750.
        height (int): height of the graph. Default is 750.

    Returns:
        figure: Plotly Express figure.
    """
    # Layout and axis
    axis = dict(
        showbackground=kwargs.get("showbackground", False),
        showline=True,
        showspikes=kwargs.get("showspikes", False),
        zeroline=False,
        showgrid=kwargs.get("showgrid", True),
        showticklabels=show_axes,
        tickmode="linear",
        tick0=0,
        dtick=1,
    )

    layout = go.Layout(
        width=kwargs.get("width", 750),
        height=kwargs.get("height", 750),
        showlegend=legend,
        hovermode="closest",
        scene=dict(
            xaxis={**dict(title="x" * show_axes), **axis},
            yaxis={**dict(title="y" * show_axes), **axis},
            zaxis={**dict(title="z" * show_axes), **axis},
        ),
        title=_get_title(title, label),
    )

    # Figure object
    fig = go.Figure(layout=layout)

    # Nodes
    nodes = np.array(egraph.nodes)
    x_nodes, y_nodes, z_nodes = nodes[:, 0], nodes[:, 1], nodes[:, 2]

    nodeColors = [_get_node_color(egraph, node, color_nodes) for node in egraph.nodes]
    nodeInfo = [_get_node_info(egraph, node, information=label) for node in egraph.nodes]

    fig.add_traces(
        go.Scatter3d(
            name="nodes",
            x=x_nodes,
            y=y_nodes,
            z=z_nodes,
            mode="markers",
            marker=dict(
                symbol="circle",
                size=5,
                color=nodeColors,
                line=dict(color="black", width=0.5),
            ),
            hovertext=nodeInfo,
            hoverinfo="text",
        )
    )

    # Edges
    x_edges, y_edges, z_edges = [], [], []

    for edge in egraph.edges:
        x0, y0, z0 = edge[0]
        x1, y1, z1 = edge[1]
        x_edges.extend([x0, x1, None])
        y_edges.extend([y0, y1, None])
        z_edges.extend([z0, z1, None])

    edgeColors = [_get_edge_color(egraph, edge, color_edges) for edge in egraph.edges]

    fig.add_traces(
        go.Scatter3d(
            name="edges",
            x=x_edges,
            y=y_edges,
            z=z_edges,
            line=dict(color="black", width=1),
            mode="lines",
            marker=dict(color=edgeColors),
            hoverinfo="none",
        )
    )

    return fig


def draw_EGraph(
    egraph,
    backend="matplotlib",
    **kwargs,
):
    """Draw an EGraph using either matplotlib or plotly as backend.

    Args:
        egraph (EGraph): The EGraph to draw.
        backend (str): The backend to use, either "matplotlib" or "plotly".

    Keyword args:
        color_nodes (bool or string or dict): Options are:
            True: color the nodes based on the 'color' attribute
            attached to the node. If unavailable, color nodes black.
            string: color all nodes with the color specified by the string
            tuple[str, dict]: color nodes based on attribute and defined colour
            string by providing a tuple with [attribute, color_dictionary],
            for example:

                ``["state", {"GKP": "b", "p": "r"}]``

            will look at the "state" attribute of the node, and colour
            according to the dictionary.

        color_edges (bool or string or dict):
            True: color the edges based on the 'color' attribute
            attached to the node. If unavailable, color nodes grey.

            string: color all edges with the color specified by the stirng

            tuple: color edges based on attribute and defined colour
            string by providing a tuple with [attribute, color_dictionary],
            for example: if the edge attribute "weight" can be +1 or -1,
            the tuple should be of the form:

                ``("weight", {-1: minus_color, +1: plus_color})``

        label (NoneType, string or iterable): plot values next to each node
            associated with the node attribute label. For example,
            to plot bit values, set label to "bit_val". If set to 'index',
            it will plot the integer indices of the nodes. If the attribute
            for some or all of the nodes, a message will print indicating
            for how many nodes the attribute has not been set.
        title (bool): if True, display the title, depending on the label.
            For default labels, the titles are converted from attribute
            name to plane English and capitalized.
        legend (bool): if True and color_nodes argument is a tuple(str, dict),
            display the a color legend with node attributes.
        show_axes (bool): if False, turn off the axes.

    See draw_EGraph_matplotlib or draw_EGraph_plotly for backend specific keyword arguments.
    """
    if backend == "matplotlib":
        return draw_EGraph_matplotlib(egraph, **kwargs)
    if backend == "plotly":
        return draw_EGraph_plotly(egraph, **kwargs)
    raise ValueError(f"Unknown backend: {backend}")


def _plot_EGraph_edges(ax, egraph, color_edges):
    """Draw the edges of the graph state represented by the EGraph.

    Args:
        ax (matplotlib.axes.Axes): the axes to draw the lines in
        color_edges (bool or string or dict):

            True: color the edges based on the 'color' attribute
                attached to the node. If unavailable, color nodes grey.
            string: color all edges with the color specified by the stirng
            dict: color edges based on attribute and defined colour
                string by providing a tuple with [attribute, color_dictionary],
                for example: if the edge attribute "weight" can be +1 or -1,
                the tuple should be of the form:
                ``("weight", {-1: minus_color, +1: plus_color})``.

    Returns:
        A Matplotib Axes object.
    """
    # Plotting edges.
    for edge in egraph.edges:
        color = _get_edge_color(egraph, edge, color_edges)

        x1, z1, y1 = edge[0]
        x2, z2, y2 = edge[1]
        ax.plot([x1, x2], [y1, y2], [z1, z2], color=color, linewidth=0.5)

    return ax


def _plot_EGraph_nodes(ax, egraph, color_nodes, label, name, legend):
    """Draw the nodes of the graph state represented by the EGraph.

    Args:
        ax (matplotlib.axes.Axes): the axes to draw the points in
        color_nodes (bool or string or dict): Options are:

            True: color the nodes based on the 'color' attribute
                attached to the node. If unavailable, color nodes black.
            string: color all nodes with the color specified by the string
            tuple[str, dict]: color nodes based on attribute and defined colour
                string by providing a tuple with (attribute, color_dictionary),
                for example: ``("state", {"GKP": "b", "p": "r"})``
                will look at the "state" attribute of the node, and colour
                according to the dictionary.

        label (NoneType or string): plot values next to each node
            associated with the node attribute label. For example,
            to plot bit values, set label to "bit_val". If set to 'index',
            it will plot the integer indices of the nodes.
        name (bool): attribute name to display as title.
        legend (bool): if True and color_nodes argument is a tuple(str, dict),
            display the a color legend with node attributes.

    Returns:
        A Matplotib Axes object.
    """

    if label:
        n_uncomputed = 0
        if label == "index":
            indices = egraph.index_generator()

    # Plotting points. y and z are swapped in the loops so that
    # z goes into the page; however, the axes labels are correct.
    for node in egraph.nodes:
        x, z, y = node
        color = _get_node_color(egraph, node, color_nodes)
        ax.scatter(x, y, z, c=color)

        if label:
            value = egraph.nodes[node].get(label) if label != "index" else indices[node]
            if value is not None:
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
                    zorder=2,
                )
            else:
                n_uncomputed += 1

    if label and n_uncomputed > 0:
        message = "{} at {} node(s) have not yet been computed."
        print(message.format(name.lower(), n_uncomputed))

    # Plotting nodes legend
    if isinstance(color_nodes, tuple) and legend:
        # these two lines are just a handy way to create a legend for
        # the node colors and attributes by generating handles of empty lines
        # with the label and color of the node property
        handles = [
            mlines.Line2D([], [], marker="o", linewidth=0, color=color, label=node_property)
            for node_property, color in color_nodes[1].items()
        ]
        ax.legend(handles=handles)

    return ax
