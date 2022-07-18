""".. _getting-started-tutorial:

Getting Started
===============================
"""

##############################################################################
# There is a vast literature available to understand the theoretical concepts behind FlamingPy. For
# a self-contained description, see
# `Xanadu's blueprint <https://quantum-journal.org/papers/q-2021-02-04-392/>`_ for a fault-tolerant
# photonic quantum computer. To see a sample of what FlamingPy can do, let us first import a few
# important objects:
#

from flamingpy.codes import SurfaceCode
from flamingpy.decoders import correct
from flamingpy.noise import CVLayer

##############################################################################
# Next, let us instantiate an RHG lattice -- the measurement-based version of the surface code:
#

RHG = SurfaceCode(3)

##############################################################################
#
# The integer denotes the code distance. By default, the boundaries are set to "open". Next, let us
# associate the nodes in the RHG lattice with CV states:
#

CVRHG = CVLayer(RHG, delta=0.08, p_swap=0.25)
CVRHG.apply_noise()

##############################################################################
#
# This had the effect of labelling half the lattice (on average) with GKP states and the other half
# with p-squeezed states. Then, a Gaussian random noise model was applied with a squeezing parameter
# of 0.1 to the states in the lattice. Finally, a syndrome measurement (sequence of homodyne
# measurements) was be conducted on the lattice, with the outcomes translated to bit values.
#
# At this point, we are ready to perform error correction on the code. We can choose the decoder,
# but also if (and how) we want to visualize it. Let's do so to get a better understanding of what
# is going on.
#

# Drawing options
node_colors = ("state", {"GKP": "gold", "p": "blue"})

dw = {
    "show_nodes": True,
    "color_nodes": node_colors,
    "show_recovery": True,
    "label_stabilizers": False,
    "label_boundary": False,
    "label_edges": False,
    "label": None,
    "legend": True,
    "show_title": True,
    "show_axes": True,
}

c = correct(RHG, decoder="MWPM", draw=True, drawing_opts=dw)

##############################################################################
#
# Finally, let's print a message identifying success or failure.
#

outcome = "succeeded." * bool(c) + "failed." * (1 - bool(c))
message = "Error correction {}".format(outcome)
print(message)
