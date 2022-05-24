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


def surface_code(d, boundaries, err, polarity, show=False):
    """Example for building and visualizing RHG lattices and surface codes."""

    # Instantiate a surface code.
    RHG_code = SurfaceCode(d, ec=err, boundaries=boundaries, polarity=polarity)
    RHG_lattice = RHG_code.graph
    RHG_fig = RHG_code.draw()

    # Check edges between boundaries for periodic boundary conditions.
    if boundaries == "periodic":
        all_boundaries = []
        for plane in ("x", "y", "z"):
            for i in (0, 2 * d - 1):
                all_boundaries += RHG_lattice.slice_coords(plane, i)
        RHG_subgraph = RHG_lattice.subgraph(all_boundaries)
        RHG_subgraph.draw()

    # Plot the stabilizers
    for ec in RHG_code.ec:
        # Stabilizers are available in the attributes primal_stabilizers and/or
        # dual_stabilizers, depending on the error complex.
        stabilizers = getattr(RHG_code, ec + "_stabilizers")
        # Change [0:1] in the following line to other consecutive indices,
        # corresponding to another stabilizer you'd like to plot.
        for stabilizer in stabilizers[0:1]:
            color = np.random.rand(3)
            for point in stabilizer.egraph:
                x, z, y = point
                RHG_fig.scatter(x, z, y, color=color, s=200)

    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":

    params = {
        # Code distance (an integer)
        "d": 2,
        # Boundaries ("open" or "periodic")
        "boundaries": "open",
        # Error complex ("primal" or "dual")
        "err": "primal",
        # Polarity (edge weight pattern in graph state -- all unit weights by default)
        "polarity": None
        # polarity = alternating_polarity'
    }

    surface_code(**params, show=True)
