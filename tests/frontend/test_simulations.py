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

# pylint: disable=protected-access

from datetime import datetime
import logging
import warnings

import itertools as it
import re

from numpy.random import default_rng as rng
import pytest

try:
    from mpi4py import MPI
except ImportError:  # pragma: no cover
    warnings.warn("Failed to import mpi4py libraries.", ImportWarning)

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.noise import CVLayer, CVMacroLayer, IidNoise
from flamingpy.simulations import ec_monte_carlo, run_ec_simulation

now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)

code_args = it.product(
    [rng(int_time).integers(2, 5), rng(int_time).integers(2, 5, 3)],
    ["primal", "dual"],
    ["open", "toric", "periodic"],
)

if "MPI" in locals():
    world_comm = MPI.COMM_WORLD
    mpi_size = world_comm.Get_size()
    mpi_rank = world_comm.Get_rank()
else:
    world_comm = None
    mpi_size = 1
    mpi_rank = 0


@pytest.fixture(scope="module", params=code_args)
def code(request):
    """A SurfaceCode object for use in this module."""
    distance, ec, boundaries = request.param
    surface_code = SurfaceCode(distance, ec, boundaries, alternating_polarity)
    return surface_code


class TestBlueprint:
    """A class with members to test Monte Carlo simulations for FT threshold
    estimations for Xanadu's blueprint architecture."""

    def test_all_GKP_high_squeezing(self, code):
        """Tests Monte Carlo simulations for FT threshold estimation of a
        system with zero swap-out probability and high squeezing."""
        trials = 10
        noise_args = {"delta": 0.001, "p_swap": 0}
        decoder_args = {}
        decoder_args["weight_opts"] = {
            "method": "blueprint",
            "integer": False,
            "multiplier": 1,
            "delta": noise_args.get("delta"),
        }

        noise_instance = CVLayer(code, **noise_args)

        errors = ec_monte_carlo(
            trials,
            code,
            noise_instance,
            "MWPM",
            decoder_args,
            world_comm=world_comm,
            mpi_rank=mpi_rank,
            mpi_size=mpi_size,
        )
        # Check that there are no errors in all-GKP high-squeezing limit for selected decoder:
        assert errors == 0


class TestPassive:
    """A class with members to test Monte Carlo simulations for FT threshold
    estimations for Xanadu's passive architecture."""

    def test_all_GKP_high_squeezing(self, code):
        """Tests Monte Carlo simulations for FT threshold estimation of a
        system with zero swap-out probability and high squeezing."""
        trials = 10
        noise_args = {"delta": 0.001, "p_swap": 0}
        noise_instance = CVMacroLayer(code, **noise_args)
        decoder_args = {"weight_opts": None}

        errors = ec_monte_carlo(
            trials,
            code,
            noise_instance,
            "UF",
            decoder_args,
            world_comm=world_comm,
            mpi_rank=mpi_rank,
            mpi_size=mpi_size,
        )
        # Check that there are no errors in all-GKP high-squeezing limit for selected decoder:
        assert errors == 0


@pytest.mark.parametrize("empty_file", sorted([True, False]))
@pytest.mark.parametrize("noise", sorted(["CVLayer", "CVMacroLayer", "IidNoise"]))
@pytest.mark.parametrize("decoder", sorted(["MWPM", "UF"]))
def test_simulations_output_file(tmpdir, empty_file, noise, decoder):
    """Check the content of the simulation benchmark output file."""
    # The Monte Carlo simulation params
    trials = 100

    code = SurfaceCode
    code_args = {
        "distance": 3,
        "ec": "primal",
        "boundaries": "open",
    }

    if noise in ["CVLayer", "CVMacroLayer"]:
        noise_args = {"delta": 0.09, "p_swap": 0.25}
        expected_header = (
            "code,distance,ec,boundaries,noise,delta,p_swap,decoder,weight_opts,"
            + "errors,trials,current_time,simulation_time,mpi_size"
        )
    else:
        noise_args = {"error_probability": 0.1}
        expected_header = (
            "code,distance,ec,boundaries,noise,error_probability,decoder,weight_opts,"
            + "errors,trials,current_time,simulation_time,mpi_size"
        )
    noise_dict = {"CVLayer": CVLayer, "CVMacroLayer": CVMacroLayer, "IidNoise": IidNoise}
    noise = noise_dict[noise]

    decoder_args = {"weight_opts": {}}

    f = tmpdir.join("sims_results.csv")
    if not empty_file:
        f.write_text(
            f"{expected_header}\n" + "\n",
            encoding="UTF-8",
        )

    args = [trials, code, code_args, noise, noise_args, decoder, decoder_args]
    run_ec_simulation(*args, fname=f)

    file_lines = f.readlines()
    # file is created with header and result lines
    assert len(file_lines) > 0

    # contains the expected header
    assert file_lines[0] == expected_header + "\n"

    # contents has the expected number of columns
    assert len(re.split(",", file_lines[-1])) == len(re.split(",", expected_header))
