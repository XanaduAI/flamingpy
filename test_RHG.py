# Copyright 2020 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
""""Unit tests for RHG classes and methods in RHG.py."""

from RHG import RHG_graph, RHGCode, RHGCube

d = 3


class TestRHGGraph:
    """Test the RHG_graph function."""

    def test_periodic_boundaries(self):
        RHG_lattice = RHG_graph(d, "periodic")
        assert len(RHG_lattice) == 6 * (d ** 3)

        all_boundaries = []
        for plane in ("x", "y", "z"):
            for i in (0,  2 * d - 1):
                all_boundaries += RHG_lattice.slice_coords(plane, i)
            # For macronode lattice, change to (-0.1, 2 * d - 0.9)
        for point in all_boundaries:
            assert len(RHG_lattice[point]) == 4


class TestRHGCode:
    """"Test the RHGCode class."""
    pass


class TestRHGCube:
    """Test the RHGCube class."""

    def test_RHG_cube(self):
        pass
