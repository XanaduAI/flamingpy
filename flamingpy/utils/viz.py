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
"""Helper functions to draw various graphs and generate plots using Matplotlib.

Plots are configured via the ``plot_params`` dictionary. These parameters
are associated with Matplolib's rc settings and are modified within the
plotting functions using the ``rc_context`` context manager. This approach
avoids having to modify the global Matplotlib ``rc_params``.

To modify the plot parameters use, for example,

.. code-block:: python

    from flamingpy.utils.viz import plot_params as fp_plot_params
    fp_plot_params["font.size"] = 20
"""

# pylint: disable=too-many-statements,singleton-comparison, too-many-lines

import itertools as it
import math

import numpy as np
import networkx as nx
import plotly.graph_objects as go
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.ticker import Formatter
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from flamingpy.codes import Stabilizer
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


def _get_title(title=None, label="index"):
    """Get the title for the EGraph plot, depending on the label if title is a
    boolean.

    Args:
        title (string, boolean or NoneType): variable to determine the returned title. If ``title``
            is a string, it will simply return the string. Else, if ``title is None``, it will
            return None. If the ``title`` is a boolean set to ``True``, it will return a title based
            on ``label``. In all other cases, the function will return None (i.e. there will be no
            title on the figure).
        label (string, list or tuple): Only relevant if ``title == True``. In that case, there are
            three options
            - if the label is set to p_phase, p_phase_cond, hom_val_p, hom_val_q, bit_val, weight
                or index, the title will be the label converted to a plane English word.
            - if the label is another string, the title will simply be that string.
            - if the label is a list or tuple of strings, the title will be the list or tuple
                unpacked separated by a comma.
    """
    # Return nothing is not title
    if not title:
        return None

    # Unpack list or tuple...
    if isinstance(label, (tuple, list)):
        if len(label) > 1:
            return ", ".join(label)
        # or convert to a single string
        label = label[0]

    # Return title directly...
    if isinstance(title, str):
        return title
    # ... or base it on label value
    if isinstance(label, str):
        if isinstance(title, bool):
            title_dict = {
                "p_phase": "Phase error probabilities",
                "p_phase_cond": "Conditional phase error probabilities",
                "hom_val_p": "p-homodyne outcomes",
                "hom_val_q": "q-homodyne outcomes",
                "bit_val": "Bit values",
                "weight": "Weights",
                "index": "Indices",
            }
        return title_dict.get(label, label)
    return None


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


def _get_node_color(egraph, node, color_nodes):
    """Color nodes based on ``color_nodes`` arg:

    - if `color_nodes` is a string use the string as color,
    - using the attribute and color dict if `color_nodes` is a tuple(str,dict),
    - or based on color attribute (when available) if `color_nodes` is bool and True;
    - black otherwise.
    """
    default_color = "black"
    if isinstance(color_nodes, str):
        color = color_nodes
    elif isinstance(color_nodes, tuple):
        node_attribute, color_dict = color_nodes
        if not (isinstance(node_attribute, str) and isinstance(color_dict, dict)):
            raise ValueError(
                "Inappropiate value for `color_nodes` argument:"
                "Check that it complies with the type `tuple(str, dict)`,"
                "where the string corresponds to a valid node attribute,"
                "the dictionary keys to valid attribute values and"
                "dictionary values to valid matplotlib color strings."
            )
        node_property = egraph.nodes[node].get(node_attribute)
        color = color_dict.get(node_property, default_color)
    elif color_nodes == True:
        color = egraph.nodes[node].get("color", default_color)
    else:
        color = default_color
    return color


def _get_node_info(egraph, node, information="coordinates"):
    """Information to be displayed when hovering over a node based on
    ``information``

    Arguments:
        egraph (EGraph): the EGraph with the node of interest.
        node (tuple): the node to get the information from.
        information (str, iterable or NoneType): the information to be displayed:
            - if ``information`` is a string, the value of the node attribute ``information``,
            - if ``information`` is an iterable, a list of the values of the node attributes in
                ``information``,
            - if ``information`` is None, return None (nothing will be displayed).
            - if ``information`` contains ``"index"``, include the index from
                ``egraph.to_indices[node]``.
            - if ``information`` is "all" it will display the coordinates, index and all the
                information avaible in the node.
    """
    # information dictionary
    info_dict = egraph.nodes[node].copy()

    # list all available information
    if information == "all":
        info_list = list(info_dict.keys())
        info_list.sort()
        information = ["index"] + info_list

    # add index if desired
    if "index" in information:
        info_dict["index"] = egraph.to_indices[node]

    # returning relevant information
    if information == "coordinates":
        return str(node)
    if information is None:
        return None
    if isinstance(information, str):
        node_property = info_dict.get(information)
        return f"{information}: {node_property}"
    if isinstance(information, (tuple, list)):
        node_info = str(node)
        for key in information:
            node_property = info_dict.get(key, None)
            if node_property is not None:
                node_info += "<br />" + f"{key}: {node_property}"
        return node_info
    raise ValueError(
        "Inappropiate value for `information` argument:"
        "Check that it complies with the type `str`,"
        "`tuple` or `list`, or has value `None`."
    )


def _get_edge_color(egraph, edge, color_edges):
    """Return the color of an EGraph edge.

    Args:
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
        A color (string)
    """
    # Color edges based on `color_edges` choices (see docstring)
    if isinstance(color_edges, str):
        color = color_edges
    elif isinstance(color_edges, tuple):
        edge_attribute, color_dict = color_edges
        if not (isinstance(edge_attribute, str) and isinstance(color_dict, dict)):
            raise ValueError(
                "Inappropiate value for `color_edges` argument:"
                "Check that it complies with the type `tuple(str, dict)`,"
                "where the string corresponds to a valid edge attribute,"
                "the dictionary keys to valid attribute values and"
                "dictionary values to valid matplotlib color strings."
            )
        edge_property = egraph.edges[edge].get(edge_attribute)
        color = color_dict.get(edge_property)
    elif color_edges == True:
        color = egraph.edges[edge].get("color") or "grey"
    else:
        color = "grey"

    return color


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


def _plot_cubes_at(positions, sizes=None, colors=None, **kwargs):
    """Plot cubes with their origin located at ``positions``.

    Note cubes are located by displacing them from the origin. In that sense,
    the location is defined by the corner matching the origin of the coordinate
    system before displacement.

    Args:
        positions (Iterable): An interable of dimension ``(N,3)`` containing the
            position of the corner of the cube.
        sizes (Iterable): An interable of dimension ``(N,3)`` containing the size of
            the cube in the coordinate directions.
        colors (Iterable): An iterable of size ``N`` containing the colors of the cube.
            This can be any of the option allowed by matplolib.

    Keyword args:
        **kwargs: all other parameters are forwarded to
            ```Poly3DColletion``
            <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.art3d.Poly3DCollection.html>`_.

    Returs:
        Poly3DCollection: A collection of 3D polygons defining the cubes.
    """

    g = [_cuboid_data(p, size=s) for p, s in zip(positions, sizes)]
    return Poly3DCollection(np.concatenate(g), facecolors=np.repeat(colors, 6, axis=0), **kwargs)


def _cuboid_data(origin, size=(1, 1, 1)):
    """Return an array with the corners of a cube of size 1."""

    X = np.array(
        [
            [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
            [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
            [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
            [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
            [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
            [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]],
        ]
    ).astype(float)
    # scale the sides of the cube
    for i in range(3):
        X[:, :, i] *= size[i]
    # displace the cube origin
    X += np.array(origin)

    return X


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


def to_pi_string(x, tex: bool = True, d=2):
    """Convert x, a multiple of sqrt(pi)/2, to a pretty string.

    If x is not a multiple of sqrt(pi)/2, return the unmodified string
    of x with `d` integers after the decimal. If tex is True, add LaTeX
    $ signs.
    """
    remainder = math.remainder(x, np.sqrt(np.pi) / 2)
    if np.isclose(remainder, 0):
        integer = round(x / (np.sqrt(np.pi) / 2))
        pref = int(integer * ((1 - integer % 2) / 2 + integer % 2))
        x_str = (not bool(round(x))) * "0" + bool(round(x)) * (
            bool(tex) * "$"
            + (not bool(1 + pref)) * "-"
            + bool(1 - abs(pref)) * str(pref)
            + r"\sqrt{\pi}"
            + (integer % 2) * "/2"
            + bool(tex) * "$"
        )
        return x_str
    return f"{x:.{d}f}"


class PiFormatter(Formatter):
    """Formatter for axis-ticks containing multiples of sqrt(pi)/2."""

    def __init__(self, tex: bool = True, d: int = 2):
        """Initialize the formatter.

        Args:
            tex: Whether to use LaTeX formatting (i.e. adding $ around the string).
            d: Number of decimals to use.
        """
        self.tex = tex
        self.d = d

    def __call__(self, x, pos=None):
        return to_pi_string(x, tex=self.tex, d=self.d)
