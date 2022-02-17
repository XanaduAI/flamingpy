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
""""Unit tests for threshold_plots.py using project's existing processed data."""
import numpy as np
import pandas as pd

from flamingpy.utils.threshold_plots import find_threshold, plot_results


# Some common input options for tests
options = {
    "input_file": "src/flamingpy/test_data/test_processed.csv",
    "save_file": "test_results.pdf",
    "p_swap": 0,
    "unit": "dB",
}
results = pd.read_csv(options["input_file"])


def test_find_thresholds():
    """Tests for find_threshold function."""
    delta_expected = 10.051
    a0_expected = 0.149
    d_p_t = find_threshold(results, options["p_swap"], options["unit"], options["save_file"], False)
    # Check if the expected threshold was returned
    assert round(d_p_t[0], 3) == delta_expected
    # Check if the expected failure probability (the offset) was returned
    assert round(d_p_t[1], 3) == a0_expected


def test_plot_results():
    """Tests for plot_results function."""
    d_p_t = find_threshold(results, options["p_swap"], options["unit"], options["save_file"], False)
    p_swap = np.array([0, 0.36, 0.71])
    sqz = np.array([10.1, 12.4, 60])
    swap_tol_data = zip(sqz, p_swap)
    plot = plot_results(
        results,
        p_swap=options["p_swap"],
        distances=[7, 9, 11],
        unit=options["unit"],
        threshold=d_p_t,
        file_name=options["save_file"],
        show=False,
        swap_tol_plot=swap_tol_data,
        inset=True,
        break_axis=True,
    )
    # Check if plot_results has returned nonempty plot
    assert plot.lines
