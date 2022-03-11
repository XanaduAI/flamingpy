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
"""Example for building and visualizing RHG lattices and surface codes."""
import matplotlib.pyplot as plt
import numpy as np

from flamingpy.codes.surface_code import SurfaceCode, alternating_polarity
from flamingpy.utils import viz
from flamingpy.cv.ops import CVLayer

show = __name__ == "__main__"


# Instantiate an RHG latice of a certian distance, with certain
# boundaries and polarity..
d = 2
boundaries = "periodic"
polarity = None
# polarity = alternating_polarity

RHG = SurfaceCode(d, boundaries=boundaries, polarity=polarity)
RHG_fig = RHG.draw()
RHG_lattice = RHG.graph

# Check edges between boundaries for periodic boundary conditions.
if boundaries == "periodic":
    all_boundaries = []
    for plane in ("x", "y", "z"):
        for i in (0, 2 * d - 1):
            all_boundaries += RHG.graph.slice_coords(plane, i)
    RHG_subgraph = RHG_lattice.subgraph(all_boundaries)
    RHG_subgraph.draw()

# Plot coordinates of stabilizers
syndrome = RHG.stabilizers
# Change [0:1] in the following line to other consecutive indices,
# corresponding to another stabilizer you'd like to plot.
for cube in syndrome[0:1]:
    color = np.random.rand(3)
    for point in cube.egraph:
        x, z, y = point
        RHG_fig.scatter(x, z, y, color=color, s=200)
if show:
    plt.show()
else:
    plt.close()

# Check CV noise model and sampling
delta = 0.001
p_swap = 0
CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
model = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
CVRHG.apply_noise(model)
CVRHG.measure_hom("p", RHG.syndrome_inds)
outcomes = [CVRHG.hom_outcomes()[i] for i in RHG.syndrome_inds]
plt.figure(figsize=(16, 9))
plt.hist(outcomes, bins=100)
CVRHG.draw(label="hom_val_p", legend=True, title=True)

if show:
    plt.show()
else:
    plt.close()
