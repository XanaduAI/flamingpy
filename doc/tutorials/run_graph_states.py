"""

.. _graph-states-tutorial:

Graph States
============

"""


######################################################################
# *Author: Ilan Tzitrin*
#


######################################################################
# A graph state is a special kind of entangled state. Certain types of
# graph states, called *cluster states*, are important resources for
# universal fault-tolerant quantum computation. In this tutorial we will
# tell you a bit about graph states, and show you how to define and
# visualize them using FlamingPy.
#


######################################################################
# Definitions
# -----------
#
# There are a few ways of writing down a qubit graph state. We'll split
# these up into *graphical*, *operational*, and *through the stabilizer*.
# The graphical definition is easiest to understand: you draw some
# (undirected) mathematical graph, i.e. a collection of nodes and edges,
# like so:
#
# .. image:: /_static/bell_cluster.svg
#    :width: 200
#    :align: center
#    :alt: Two-qubit cluster states (bell pair)
#
# In the above graph (call it :math:`B`), the nodes represent qubit
# states, and the edges indicate entanglement. This can be made more
# concrete using the operational definition: associate the nodes with
# :math:`\vert + \rangle` states (i.e. superpositions of computational
# basis states:
# :math:`\frac{1}{\sqrt{2}} \left(\vert 0 \rangle + \vert 1 \rangle\right)`),
# and associate the edges with :math:`CZ` gates:
# :math:`\vert 0 \rangle \langle 0 \vert \otimes I + \vert 1 \rangle \langle 1 \vert \otimes Z`.
# This means the graph state :math:`\vert B \rangle` corresponding to the
# graph above is:
#
# .. math::
#
#    \vert B \rangle = CZ (\vert + \rangle_1 \vert + \rangle_2) = \frac{1}{\sqrt{2}} \left( \vert 0 \rangle_1 \vert + \rangle_2 + \vert 1 \rangle_1 \vert -\rangle_2 \right),
#
# where
# :math:`\vert - \rangle = \frac{1}{\sqrt{2}} \left(\vert 0 \rangle - \vert 1 \rangle\right)`.
# We can write down a circuit diagram for this:
#
# .. image:: /_static/bell_circuit.svg
#    :width: 200
#    :align: center
#    :alt: Circuit to produce Bell cluster state
#
# You may recognize :math:`\vert B \rangle` as a type of *Bell* or *EPR
# pair*: a maximally entangled state of two qubits. To create a
# three-qubit entangled state (equivalent to a *GHZ state*), we can follow
# the same process:
#
# .. math::
#
#    \vert \text{GHZ} \rangle =  CZ_{23}CZ_{12} (\vert + \rangle _1 \vert + \rangle _2 \vert + \rangle _3) = \frac{1}{\sqrt{2}}(\vert 0 \rangle _1 \vert + \rangle _2 \vert + \rangle _3 + (\vert 1 \rangle _1 \vert - \rangle _2 \vert - \rangle _3)).
#
# The corresponding graph and circuit are:
#
# .. image:: /_static/GHZ_cluster.svg
#    :width: 200
#    :align: center
#    :alt: Three qubit cluster state (GHZ state)
#
#
#
# .. image:: /_static/GHZ_circuit.svg
#    :width: 200
#    :align: center
#    :alt: Circuit to produce GHZ cluster
#

######################################################################
# This is the general way entanglement is ''generated'' in a graph state:
# the :math:`CZ` gate (edge) entangles each pair of qubits (nodes).
#
# Stabilizer definition
# ^^^^^^^^^^^^^^^^^^^^^
#
# It is possible to completely understand graph states using the
# definitions above. However, it is very convenient---especially for error
# correction---to talk about the *stabilizer* of the graph state. The idea
# is that we can define graph states not through individual qubit states,
# but through operators.
#
# In general, given a state of :math:`n` qubits, you can uniquely
# determine the state by :math:`n` distinct operators called *stabilizer
# elements* (or, colloquially, just *stabilizers*). When applied to the
# state, these operators leave the state unchanged. The stabilizers for
# the above states are:
#
# .. math::
#
#     S_B = \{X_1 Z_2, X_2 Z_1\} \qquad \text{and} \qquad S_\text{GHZ} = \{ X_1 Z_2 Z_3, Z_1 X_2 Z_3, Z_1 Z_2 X_3 \}.
#
# To see this, you can use the fact that the :math:`Z` operator flips
# phases:
#
# .. math::
#
#    Z \vert 0 \rangle = \vert 0 \rangle, \quad Z\vert 1 \rangle = - \vert 1 \rangle, \quad Z\vert + \rangle = \vert - \rangle, \quad Z\vert - \rangle= \vert + \rangle,
#
# while the :math:`X` operator flips bits:
#
# .. math::
#
#    X \vert 0 \rangle = \vert 1 \rangle, \quad X\vert 1 \rangle= \vert 0 \rangle, \quad X\vert + \rangle = \vert + \rangle, \quad X\vert - \rangle = -\vert - \rangle.
#
# If you apply the stabilizer elements to the above states, you'll see
# they don't change. In general, if you're given a graph state with
# :math:`n` qubits, you can come up with :math:`n` stabilizers. Each one
# will apply :math:`X` to a given node, and :math:`Z` to all the
# neighbours.
#


######################################################################
# Creating and visualizing graph states with FlamingPy
# ----------------------------------------------------
#
# In FlamingPy, graph states are represented through the ``EGraph``
# (enhanced graph) class. Let us import this class, along with the
# visualization module:
#

from flamingpy.codes.graphs import EGraph
import flamingpy.utils.viz as viz


######################################################################
# An ``EGraph`` is a type of (inherits from) a ``Graph`` object from the
# ``networkx`` package, but it builds on ``networkx`` graphs with its own
# functionality. EGraphs (like NetworkX graphs) have dictionaries of nodes
# and edges. To properly define a graph state, the ``EGraph`` class
# assumes that the nodes are specified by three-tuples :math:`(x, y, z)`
# corresponding to coordinates in three dimensions.
#
# We can construct a GHZ state using FlamingPy. To do so, we have to place its nodes in 3D
# space. There are infinite choices of coordinates available to us, but
# let us place the points at corners of the unit cube:
#

GHZ_edge_1 = {(0, 0, 0), (0, 0, 1)}
GHZ_edge_2 = {(0, 0, 1), (1, 0, 1)}


######################################################################
# We can give an ``EGraph`` its edges right away, but let us instead first
# initialize an empty graph:
#

GHZ_state = EGraph(dims=(1, 1, 1))
GHZ_state.add_edges_from([GHZ_edge_1, GHZ_edge_2])


######################################################################
# Notice that there is currently nothing in the graph dictionaries except
# for the coordinates of the points:
#

print("Node attributes: ", *GHZ_state.nodes.data())
print("Edge attributes ", *GHZ_state.edges.data())


######################################################################
# Eventually, these dictionaries can be populated attributes useful for
# error correction and visualization. Now, we can create a plot of the state.
# We will first specify some drawing options, and then use the ``draw``
# method of the ``EGraph``. This is as easy as:
#

drawing_opts = {
    "color_nodes": "MidnightBlue",
    "color_edges": "slategray",
}
GHZ_state.draw(**drawing_opts)


######################################################################
# We can also extract some information about the graph state, including
# the *adjacency matrix* :math:`A` of the underlying graph. The indices
# (rows and columns) of this matrix correspond to the nodes of the graph.
# The entry :math:`A_{ij}` is 0 if there is no edge connecting the nodes,
# and 1 (or another number, for weighted graphs) otherwise. We can
# generate the adjacency matrix and then plot its heat map:
#

GHZ_state.adj_generator(sparse=False)
adj = GHZ_state.adj_mat

viz.plot_params["figure.figsize"] = (5.4, 4)
viz.plot_mat_heat_map(adj)
