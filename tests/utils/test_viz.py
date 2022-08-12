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
""""Unit tests for functions in the viz module."""

# pylint: disable=no-self-use

import math

import numpy as np
from numpy.random import default_rng as rng
import pytest
import matplotlib
import matplotlib.pyplot as plt
import plotly

from flamingpy.utils import viz
from flamingpy.codes.graphs import EGraph
from flamingpy.codes import SurfaceCode


def test_to_pi_string():
    """Test for the convenience function to_pi_string."""
    # Test +- sqrt(pi) and sqrt(pi)/2.
    assert viz.to_pi_string(np.sqrt(np.pi)) == "$\\sqrt{\\pi}$"
    assert viz.to_pi_string(-np.sqrt(np.pi)) == "$-\\sqrt{\\pi}$"
    assert viz.to_pi_string(np.sqrt(np.pi) / 2) == "$\\sqrt{\\pi}/2$"
    assert viz.to_pi_string(-np.sqrt(np.pi) / 2) == "$-\\sqrt{\\pi}/2$"

    # Test random odd integer multiples of sqrt(pi)/2, eccept 1 and -1
    odd_int = (2 * rng().integers(2, 25) - 1) * (-1) ** rng().integers(2)
    assert viz.to_pi_string(odd_int * np.sqrt(np.pi) / 2) == "${}\\sqrt{{\\pi}}/2$".format(odd_int)

    #  Test random even multiples of sqrt(pi).
    even_int = odd_int + 1
    assert viz.to_pi_string(even_int * np.sqrt(np.pi)) == "${}\\sqrt{{\\pi}}$".format(even_int)

    # Check everything else converted into a str.
    rand_numb = rng().random()
    rand_d = rng().integers(2, 25)
    if not np.isclose(math.remainder(rand_numb, np.sqrt(np.pi) / 2), 0):
        assert viz.to_pi_string(rand_numb, d=rand_d) == "{:.{}f}".format(rand_numb, rand_d)

    # Test for tex=False
    assert viz.to_pi_string(-np.sqrt(np.pi) / 2, tex=False) == "-\\sqrt{\\pi}/2"


class TestDrawEGraphMaplotlib:
    """Tests for visualizing EGraphs."""

    def test_draw_egraph_bell(self):
        """Test for the draw method of EGraph of Bell state."""
        # Bell state EGraph
        edge = [(0, 0, 0), (0, 0, 1)]
        bell_state = EGraph()
        bell_state.add_edge(*edge, color="MidnightBlue")

        # Test for drawing the EGraph
        _, a = viz.draw_EGraph_matplotlib(bell_state)
        plt.close()

        assert len(a.get_xticks()) == 1
        assert a.get_xlim() == (-1, 1)

    def test_wrapper_draw_egraph(self):
        """Tests the returned object of EGraph.draw of EGraph with one node."""
        E = EGraph()
        E.add_node((0, 0, 0))
        f, a = E.draw()
        assert issubclass(type(f), matplotlib.figure.Figure)
        assert issubclass(type(a), matplotlib.axes.Axes)

    @pytest.mark.parametrize("d", (2, 3))
    def test_draw_egraph_rhg(self, d):
        """Test for the draw method of EGraph of RHG lattice."""
        # Bell state EGraph
        RHG = SurfaceCode(d).graph

        # Test for drawing the EGraph
        _, a = viz.draw_EGraph_matplotlib(RHG)
        plt.close()

        n_ticks = 2 * d - 1

        ticks = (a.get_xticks(), a.get_yticks(), a.get_zticks())
        assert [len(tick) for tick in ticks] == [n_ticks] * 3

        actual_lims = (a.get_xlim(), a.get_ylim(), a.get_zlim())
        assert actual_lims == ((0, n_ticks - 1), (1, n_ticks), (1, n_ticks))


class TestDrawEGraphPlotly:
    """Tests for visualizing EGraphs with Plotly."""

    @pytest.mark.parametrize("d", (2, 3))
    def test_draw_egraph_rhg_plotly(self, d):
        """Test for the draw method of EGraph of RHG lattice."""
        # Bell state EGraph
        SC = SurfaceCode(d)
        RHG = SC.graph

        # Test for drawing the EGraph
        fig1 = viz.draw_EGraph_plotly(RHG)
        assert isinstance(fig1, plotly.graph_objs.Figure)

        # Test for drawing the SurfaceCode
        fig2 = SC.draw(backend="plotly")
        assert isinstance(fig2, plotly.graph_objs.Figure)
