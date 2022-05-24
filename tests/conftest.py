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
"""Default parameters, environment variables, fixtures, and shared objects
the unit tests and manual checks."""
import os
import pytest
import flamingpy

# defaults
TOL = 1e-5


@pytest.fixture(scope="session")
def tol():
    """Numerical tolerance for equality tests."""
    return float(os.environ.get("TOL", TOL))
