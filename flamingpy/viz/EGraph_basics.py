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
"""A series of basic functions to extract the properties of a given EGraph."""

# pylint: disable=too-many-statements,singleton-comparison,too-many-lines


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
