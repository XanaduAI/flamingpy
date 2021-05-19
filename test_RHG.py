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
import itertools as it

d = 3


class TestRHGGraph:
    """Test the RHG_graph function."""

    def test_boundary_combinations(self):
        for boundaries in it.product(['primal', 'dual', 'periodic'], repeat=3):
            assert RHG_graph(d, boundaries)

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

    def test_macronodes(self):
        for boundaries in {"finite", "periodic"}:
            macronode_lattice = RHG_graph(d, boundaries, macronodes=True)
            RHG_macronodes = macronode_lattice.macro
            RHG_reduced = RHG_graph(d, boundaries)
            assert not set(RHG_macronodes) - set(RHG_reduced)
            if boundaries == "periodic":
                assert len(macronode_lattice) == 4 * len(RHG_reduced)


class TestRHGCode:
    """"Test the RHGCode class."""

    def test_stabilizers(self):
        for boundaries in ("finite", "periodic"):
            RHG_code = RHGCode(d, error_complex="primal", boundaries=boundaries)
            cubes = RHG_code.stabilizers
            assert len(cubes) == d ** 3
            for cube in cubes:
                assert len(cube.physical) == 6
            if boundaries == "periodic":
                assert len(cube.egraph) == 6

class TestRHGCube:
    """Test the RHGCube class."""

    def test_RHG_cube(self):
        pass
