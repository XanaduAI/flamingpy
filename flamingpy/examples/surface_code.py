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

# pylint: disable=unused-import

import matplotlib.pyplot as plt
import numpy as np

from flamingpy.codes.surface_code import (
    SurfaceCode,
    alternating_polarity,
)


def illustrate_surface_code(d, boundaries, err, polarity, stabilizer_inds=None, show=False):
    """Example for building and visualizing RHG lattices and surface codes.

    See more details about the following arguments in the SurfaceCode class.

    Args:
        d (int): the code distance
        boundaries (str): the code boundaries ("open" or "periodic")
        err (str): the error complex ("primal" or "dual")
        polarity (func): the polarity (edge weight) function
        show (bool): if True, show the plot
        stabilizer_indices (list): indices of the stabilizers to plot (all by
            default).
    """

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
        # Stabilizers are available in the attributes primal_stabilizers and/or dual_stabilizers,
        # depending on the error complex.

        stabilizers = getattr(RHG_code, ec + "_stabilizers")

        if stabilizer_inds is None:
            stabilizer_inds = range(len(stabilizers))
        for ind in stabilizer_inds:
            stabilizer = stabilizers[ind]
            color = np.random.rand(3)
            for point in stabilizer.egraph:
                x, z, y = point
                RHG_fig.scatter(x, z, y, color=color, s=40)

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
        "polarity": None,
        # polarity = alternating_polarity,
        # indices of stabilizer nodes to scatter
        # "stabilizer_inds": [0, 3],
        "show": True,
    }

    illustrate_surface_code(**params)
