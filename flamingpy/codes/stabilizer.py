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
"""The abstract stabilizer class."""
import numpy as np


class Stabilizer:
    """A class for representing a stabilizer element of a code.

    Arguments:
        G (EGraph): the subgraph of the EGraph of the lattice induced by
            the stabilizer vertices.

    Attributes:
        egraph (EGraph): the corresponding EGraph
        physical (NoneType or list): the list of physical coordinates
            associated with the complete stabilizer element, mainly for plotting
            purposes. Specified during syndrome identification in the code
            class. For example, the list will contain six points even for
            incomplete surface code stabilizers, and will have only the
            local points in case of periodic boundaries.
    """

    def __init__(self, G):
        """Initialize the Stabilizer with its associated egraph."""
        self.egraph = G
        self.physical = None

    @property
    def parity(self):
        """Compute total parity of the cube.

        Uses the bit_values found in the 'bit_val' attribute of the nodes of
        the EGraph of the code.

        Returns:
            int: the total parity
        """
        bit_vals = [self.egraph.nodes[node]["bit_val"] for node in self.egraph]
        par = int(np.sum(bit_vals) % 2)
        return par

    def coords(self):
        """Return the coordinates of the syndrome vertices."""
        return list(self.egraph.nodes)

    def xlims(self):
        """Return the smallest and largest x vals of the coordinates."""
        xs = [tup[0] for tup in self.physical]
        xmin, xmax = min(xs), max(xs)
        return xmin, xmax

    def ylims(self):
        """Return the smallest and largest y vals of the coordinates."""
        ys = [tup[1] for tup in self.physical]
        ymin, ymax = min(ys), max(ys)
        return ymin, ymax

    def zlims(self):
        """Return the smallest and largest z vals of the coordinates."""
        zs = [tup[2] for tup in self.physical]
        zmin, zmax = min(zs), max(zs)
        return zmin, zmax

    def midpoint(self):
        """Return the midpoint of the cube."""
        return (
            np.average(self.xlims()),
            np.average(self.ylims()),
            np.average(self.zlims()),
        )
