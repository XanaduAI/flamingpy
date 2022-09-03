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
"""Helper functions to plot and extract properties of the cubes .
"""

# pylint: disable=too-many-statements,singleton-comparison, too-many-lines

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def _plot_cubes_at(positions, sizes=None, colors=None, **kwargs):
    """Plot cubes with their origin located at ``positions``.

    Note cubes are located by displacing them from the origin. In that sense,
    the location is defined by the corner matching the origin of the coordinate
    system before displacement.

    Args:
        positions (Iterable): An interable of dimension ``(N,3)`` containing the
            position of the corner of the cube.
        sizes (Iterable): An interable of dimension ``(N,3)`` containing the size of
            the cube in the coordinate directions.
        colors (Iterable): An iterable of size ``N`` containing the colors of the cube.
            This can be any of the option allowed by matplolib.

    Keyword args:
        **kwargs: all other parameters are forwarded to
            ```Poly3DColletion``
            <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.art3d.Poly3DCollection.html>`_.

    Returs:
        Poly3DCollection: A collection of 3D polygons defining the cubes.
    """

    g = [_cuboid_data(p, size=s) for p, s in zip(positions, sizes)]
    return Poly3DCollection(np.concatenate(g), facecolors=np.repeat(colors, 6, axis=0), **kwargs)


def _cuboid_data(origin, size=(1, 1, 1)):
    """Return an array with the corners of a cube of size 1."""

    X = np.array(
        [
            [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
            [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
            [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
            [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
            [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
            [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]],
        ]
    ).astype(float)
    # scale the sides of the cube
    for i in range(3):
        X[:, :, i] *= size[i]
    # displace the cube origin
    X += np.array(origin)

    return X
