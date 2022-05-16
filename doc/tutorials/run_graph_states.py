"""
Graph States
============


"""


######################################################################
# A graph state is a special kind of entangled state. Certain types of
# graph states, called cluster states, are important resources for
# universal fault-tolerant quantum computation. In this tutorial we’ll
# tell you a bit about graph states, and show you how to define and
# visualize them using FlamingPy.
#
######################################################################
# Definitions
# -----------
#
# There are a few ways of writing down a qubit graph state. We’ll split
# these up into *graphical*, *operational*, and *through the stabilizer*.
# The graphical definition is easiest to understand: you draw some
# (undirected) mathematical graph, i.e. a collection of nodes and edges,
# like so:
#
# In the above, the nodes represent qubit states, and the edges indicate
# entanglement. This can be made more concrete using the operational
# definition: associate the nodes with :math:`\vert + \rangle` states
# (i.e. superpositions of computational basis states:
# :math:`\frac{1}{\sqrt{2}} \left(\vert 0 \rangle + \vert 1 \rangle\right)`),
# and associated the edges with :math:`CZ` gates:
# :math:`\vert 0 \rangle \langle 0 \vert \otimes I + \vert 1 \rangle \langle 1 \vert \otimes Z`.
# This means the graph state :math:`\vert B \rangle` corresponding to the
# graph above is
#
# .. math::
#   \begin{equation} \vert B \rangle = CZ (\vert + \rangle \vert + \rangle) = \frac{1}{\sqrt{2}} \left( \vert 0 \rangle \vert + \rangle + \vert 1 \rangle \vert -\rangle \right), \end{equation}
#
# where
# :math:`\vert - \rangle = \frac{1}{\sqrt{2}} \left(\vert 0 \rangle - \vert 1 \rangle\right)`.
# We can write down a circuit diagram for this:
#
# You may recognize :math:`\vert B \rangle` as a Bell or EPR pair: a
# maximally entangled state of two qubits. (To be precise, it is a Bell
# pair with a :math:`\pi / 2` or Hadamard rotation on the second qubit).
# To create a three-qubit entangled state(a rotated GHZ) state, we can
# follow the same process:
#
# .. math::
#   \begin{equation}
#   \vert \text{GHZ} \rangle =  CZ_{23}CZ_{21} (\vert + \rangle _1 \vert + \rangle _2 \vert + \rangle _3) = \frac{1}{\sqrt{2}}(\vert 0 \rangle _1 \vert + \rangle _2 \vert + \rangle _3 + (\vert 1 \rangle _1 \vert - \rangle _2 \vert - \rangle _3)).
#   \end{equation}
#

######################################################################
# This is the general way entanglemenet is “generated” in a graph state:
# the :math:`CZ` gate (edge) entangles each pair of qubits (nodes).
#
# Stabilizer definition
# ^^^^^^^^^^^^^^^^^^^^^
#
# It is possible to completely understand graph states using the
# definitions above. However, it is very convenient—especially for error
# correction—to talk about the *stabilizer* of the graph state. The idea
# is that we can define graph states not through individual qubit states,
# but through operators. In general, given a state of :math:`n` qubits,
# you can uniquely detemrine the state by :math:`n` distinct operators
# called *stabilizer elements* (or, colloquially, just “stabilizers”).
# When applied to the state, these operators leave the state unchanged.
# The stabilizers for the above states are:
#
# .. math::
#   \begin{equation} S_B = \{X_1 Z_2, X_2, Z_1\} \qquad \text{and} \qquad S_\text{GHZ} = \{ X_1 Z_2 Z_3, Z_1 X_2 Z_3, Z_1 Z_2 X_3 \}. \end{equation}
#
# To see this, you can use the fact that the :math:`Z` operator changes
# phases:
#
# .. math::
#   \begin{align}
#   Z \vert 0 \rangle = \vert 0 \rangle, \quad Z\vert 1 \rangle = - \vert 1 \rangle, \quad Z\vert + \rangle = \vert - \rangle, \quad Z\vert - \rangle= \vert + \rangle,
#   \end{align}
#
# while the :math:`X` operator flips bits:
#
# .. math::
#   \begin{align}
#   X \vert 0 \rangle = \vert 1 \rangle, \quad X\vert 1 \rangle= \vert 0 \rangle, \quad X\vert + \rangle = \vert + \rangle, \quad X\vert - \rangle = -\vert - \rangle.
#   \end{align}
#
# If you apply the stabilizer elements to the above states, you’ll see
# they don’t change. In general, if you’re given a graph state with
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
# ``NetworkX`` package, but it builds on NetworkX graphs with its own
# functionality. EGraphs (like NetworkX graphs) have dictionaries of nodes
# and edges. To properly define a graph state, the ``EGraph`` class
# assumes that the nodes are specified by three-tuples :math:`(x, y, z)`
# corresponding to coordinates in three dimensions.
#
# Let us construct a Bell pair. To do so, we have to place its nodes in 3D
# space. There are infinitely choices of coordinates available to us, but
# let us place the points at the corners of the unit qube:
#

bell_edge = {(0, 0, 0), (1, 1, 1)}

######################################################################
# We can give an EGraph its edges right away, but let us instead first
# initialize an empty graph:
#
bell_state = EGraph(dims=(1,1,1))
bell_state.add_edge(*bell_edge)

######################################################################
# Notice that there is currently nothing in the graph dictionaries except
# for the coordinates of the points:
#
print("Node attributes: ", bell_state.nodes.data())
print("Edge attributes ", bell_state.edges.data())

######################################################################
# Eventually, these dictionaries can be populated attributes useful for
# error correction and visualization.
#

######################################################################
# Let us create a plot the state. We can first specify some drawing
# options: This is as easy as:
#
drawing_opts = {
    "color_nodes": "MidnightBlue",
    "color_edges": "Firebrick"
}
bell_state.draw(**drawing_opts)

######################################################################
# We can also extract some information about the graph state, including
# the *adjacency matrix* :math:`A` of the underlying graph. This indices
# (rows and columns) of this matrix correspond to the nodes of the graph.
# The entry :math:`A_{ij}` is 0 if there is no edge connecting the nodes,
# and 1 (or another number, for weighted graphs) otherwise. We can
# generate the adjacency matrix and then plot its heat map:
#
bell_state.adj_generator(sparse=False)
adj = bell_state.adj_mat
viz.plot_mat_heat_map(adj)