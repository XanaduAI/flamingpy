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
"""
Unit tests for flamingpy top _init_.py functions.
"""

import contextlib
import io
import re

import flamingpy as fp


def test_about():
    """Tests if the about string prints correctly."""
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        ff.about()
    out = f.getvalue().strip()

    assert "Platform info:" in out
    assert "Installation path:" in out
    assert "Python version:" in out
    assert "FlamingPy version:" in out
    assert "Numpy version:" in out
    assert "Scipy version:" in out
    assert "NetworkX version:" in out
    assert "RetworkX version:" in out
    assert "Matplotlib version:" in out
    assert "lemonpy shared object:" in out
    assert "cpp_mc_loop shared object:" in out
