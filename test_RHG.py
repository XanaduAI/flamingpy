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
from numpy.random import default_rng as rng
from networkx import fast_gnp_random_graph
import itertools as it
import pytest


@pytest.mark.parametrize("d", range(2, 5))
class TestRHGGraph:
    """Test the RHG_graph function."""

    def test_boundary_combinations(self, d):
        # Check that different boundary conditions produce a nonempty
        # lattice.
        for boundaries in it.product(["primal", "dual", "periodic"], repeat=3):
            assert len(RHG_graph(d, boundaries))
        assert len(RHG_graph(d, "finite"))

    def test_periodic_boundaries(self, d):
        RHG_lattice = RHG_graph(d, "periodic")
        assert len(RHG_lattice) == 6 * (d ** 3)

        all_boundaries = []
        for plane in ("x", "y", "z"):
            for i in (0, 2 * d - 1):
                all_boundaries += RHG_lattice.slice_coords(plane, i)
            # For macronode lattice, change to (-0.1, 2 * d - 0.9)
        # Check that every boundary coordinate has four neighbours.
        for point in all_boundaries:
            assert len(RHG_lattice[point]) == 4

    def test_polarity(self, d):
        pass

    def test_macronodes(self, d):
        for boundaries in {"finite", "periodic"}:
            macronode_lattice = RHG_graph(d, boundaries, macronodes=True)
            RHG_macronodes = macronode_lattice.macro
            RHG_reduced = RHG_graph(d, boundaries)
            # Check that the reduced lattice has the same node set
            # as the nodes stored in the macro attribute of the
            # macronode lattice.
            assert not set(RHG_macronodes) - set(RHG_reduced)
            # Check that the macronode lattice in the case of periodic
            # boundaries is exactly four times larger than the reduced
            # lattice.
            if boundaries == "periodic":
                assert len(macronode_lattice) == 4 * len(RHG_reduced)


code_params = it.product(range(2, 5), ["primal", "dual", "finite", "periodic"])


@pytest.fixture(scope="module", params=code_params)
def RHG_code(request):
    distance, boundaries = request.param
    return RHGCode(distance, error_complex="primal", boundaries=boundaries)


class TestRHGCode:
    """"Test the RHGCode class."""

    def test_init(self, RHG_code):
        # Check proper initialization
        distance = RHG_code.distance
        assert RHG_code.dims == (distance, distance, distance)
        assert RHG_code.complex
        assert RHG_code.boundaries
        assert RHG_code.syndrome_inds

    def test_stabilizers(self, RHG_code):
        cubes = RHG_code.stabilizers
        # Check that there are d^3 stabilizers
        assert len(cubes) == RHG_code.distance ** 3
        # Check that each stabilizer has 6 corresponing physical
        # vertices, no matter if it's a 3, 4, 5, or 6 body stabilizer.
        for cube in cubes:
            assert len(cube.physical) == 6
        # For all primal or periodic boundaries, check that there are
        # only 6-body stabilizers.
        if RHG_code.boundaries in [["primal"] * 3, ["periodic"] * 3]:
            assert len(cube.egraph) == 6
        # TODO: "finite" boundary test

    def test_syndrome(self, RHG_code):
        syndrome_coords = RHG_code.syndrome_coords
        syndrome_from_stabes = [set(cube.egraph) for cube in RHG_code.stabilizers]
        flattened_syndrome = [a for b in syndrome_from_stabes for a in b]
        # Check that syndrome coordinates in the code class match the
        # coordinates in the stabilizers.
        assert not set(syndrome_coords) - set(flattened_syndrome)

    def test_boundary(self, RHG_code):
        if "primal" in RHG_code.boundaries:
            assert RHG_code.boundary_coords
        else:
            assert not RHG_code.boundary_coords


class TestRHGCube:
    """Test the RHGCube class."""

    def test_parity(self):
        # Check that parity is properly computed.
        n = rng().integers(10)
        random_graph = fast_gnp_random_graph(n, 0.5)
        bit_list = []
        for node in random_graph:
            random_bit = rng().integers(2)
            bit_list += [random_bit]
            random_graph.nodes[node]["bit_val"] = random_bit
        random_cube = RHGCube(random_graph)
        assert random_cube.parity == sum(bit_list) % 2
