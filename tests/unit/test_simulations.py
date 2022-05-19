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
"""Unit tests for Monte Carlo simulations for estimating FT thresholds."""

# pylint: disable=no-self-use,protected-access,too-few-public-methods

import itertools as it

import re
import pytest

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.cv.macro_reduce import splitter_symp
from flamingpy.simulations import ec_monte_carlo, run_ec_simulation

code_params = it.product([2, 3, 4], ["primal", "dual"], ["open", "periodic"])


@pytest.fixture(scope="module", params=code_params)
def code(request):
    """A SurfaceCode object for use in this module."""
    distance, ec, boundaries = request.param
    surface_code = SurfaceCode(distance, ec, boundaries, alternating_polarity)
    surface_code.graph.index_generator()
    return surface_code


class TestBlueprint:
    """A class with members to test Monte Carlo simulations for FT threshold
    estimations for Xanadu's blueprint architecture."""

    def test_all_GKP_high_squeezing(self, code):
        """Tests Monte Carlo simulations for FT threshold estimation of a
        system with zero swap-out probability and high squeezing."""
        p_swap = 0
        delta = 0.001
        trials = 10
        errors_py = ec_monte_carlo(code, trials, delta, p_swap, passive_objects=None)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors_py == 0


class TestPassive:
    """A class with members to test Monte Carlo simulations for FT threshold
    estimations for Xanadu's passive architecture."""

    def test_all_GKP_high_squeezing(self, code):
        """Tests Monte Carlo simulations for FT threshold estimation of a
        system with zero swap-out probability and high squeezing."""
        p_swap = 0
        delta = 0.001
        trials = 10

        pad_bool = code.bound_str != "periodic"
        # The lattice with macronodes.
        RHG_macro = code.graph.macronize(pad_boundary=pad_bool)
        RHG_macro.index_generator()
        RHG_macro.adj_generator(sparse=True)

        # The empty CV state, uninitiated with any error model.
        CVRHG_reduced = CVLayer(code)

        # Define the 4X4 beamsplitter network for a given macronode.
        # star at index 0, planets at indices 1-3.
        bs_network = splitter_symp()
        noise_model = {"noise": "grn", "delta": delta, "p_swap": p_swap}
        passive_objects = [RHG_macro, code.graph, CVLayer, noise_model, bs_network]
        errors_py = ec_monte_carlo(code, trials, delta, p_swap, passive_objects=passive_objects)
        # Check that there are no errors in all-GKP high-squeezing limit.
        assert errors_py == 0


@pytest.mark.parametrize("passive", [True, False])
@pytest.mark.parametrize("empty_file", [True, False])
@pytest.mark.parametrize("sim", [run_ec_simulation])
def test_simulations_output_file(tmpdir, passive, empty_file, sim):
    """Check the content of the simulation benchmark output file."""

    expected_header = (
        "distance,passive,ec,boundaries,delta,p_swap,decoder,errors_py,"
        + "trials,current_time,decoding_time,simulation_time"
    )
    dummy_content = "2,True,primal,open,0.04,0.5,2,10,12:34:56,10,20"
    if "benchmark" in sim.__name__:
        expected_header += ",cpp_to_py_speedup"
        dummy_content += ",1"

    f = tmpdir.join("sims_results.csv")
    if not empty_file:
        f.write_text(
            f"{expected_header}\n" + "\n",
            encoding="UTF-8",
        )

    # simulation params
    params = {
        "distance": 2,
        "ec": "primal",
        "boundaries": "open",
        "delta": 0.04,
        "p_swap": 0.5,
        "trials": 100,
        "passive": True,
        "decoder": "MWPM",
    }

    sim(**params, fname=f)

    file_lines = f.readlines()
    # file is created with header and result lines
    assert len(file_lines) > 0

    # contains the expected header
    assert file_lines[0] == expected_header + "\n"

    # contents has the expected number of columns
    assert len(re.split(",", file_lines[-1])) == len(re.split(",", expected_header))
