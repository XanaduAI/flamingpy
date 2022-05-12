r"""
.. _tutorial_macro_reduce:
Macronodes
=================================

*Author: FlamingPy dev team. Last updated: 12 May 2022.*
"""


######################################################################
# First, we import the relevant modules and functions

#

from flamingpy.codes import SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import BS_network, reduce_macro_and_simulate
from flamingpy.decoders.decoder import correct

import matplotlib.pyplot as plt


######################################################################
# Next, we set the number of trials and the code and noise parameters. You

# are encouraged to play around with these!

#

# Code parameters
d = 3
boundaries = "open"
ec = "primal"

# Noise parameters
delta = 0.1
p_swap = 0.4

# The reduced lattice.
RHG_code = SurfaceCode(d, ec=ec, boundaries=boundaries)
RHG_reduced = RHG_code.graph


######################################################################
# We can generate indices for the nodes (qubits) with ``index_generator``

# as below. Note that I used ``;`` to suppress the output (this can be

# pretty large if there are a lot of qubits).

#

RHG_reduced.index_generator()


######################################################################
# Let's draw the lattice to see what it looks like.

#

RHG_reduced.draw()
plt.show()


######################################################################
# From the reduced RHG lattice, we can generated a so-called macronized

# version of it, see `this

# paper <https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.2.040353>`__

# for more details. Essentially, we replace every node with 4 nodes. This

# is particularly useful, because it allows for a physical implementation

# with **static** optics!

#

# The lattice with macronodes
pad_bool = boundaries != "periodic"
RHG_macro = RHG_code.graph.macronize(pad_boundary=pad_bool)
RHG_macro.index_generator()
RHG_macro.adj_generator(sparse=True)

RHG_macro.draw()
plt.show()


######################################################################
# As you can see, we replaced every node with a group of four nodes,

# sweet!

#

# The empty CV layer, uninitiated with any error model.
CVRHG_reduced = CVLayer(RHG_code)

# Define the 4X4 beamsplitter network for a given macronode.
# star at index 0, planets at indices 1-3.
bs_network = BS_network(4)


######################################################################
# Let's see how well our encoding works. To do so, we run a simulation

# ``n_trials``, and see how many times we corrected it truthfully.

#

# Number of trials
n_trials = 100

# Simulation
successes = 0
for trial in range(n_trials):
    # The empty CV state, uninitiated with any error model.
    reduce_macro_and_simulate(RHG_macro, RHG_reduced, CVRHG_reduced, bs_network, p_swap, delta)

    weight_options = {
        "method": "blueprint",
        "prob_precomputed": True,
    }

    decoder = {"outer": "MWPM"}

    c = correct(code=RHG_code, decoder=decoder, weight_options=weight_options)

    successes += int(c)

# Results
error = (n_trials - successes) / n_trials
print("Trials: ", n_trials)
print("Successes: ", successes)
print("Error rate: ", error)
