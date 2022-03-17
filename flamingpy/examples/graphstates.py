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
"""Example for building and visualizing EGraphs with CV noise."""
import matplotlib.pyplot as plt

from flamingpy.utils import viz
from flamingpy.codes.graphs import EGraph
from flamingpy.cv.ops import CVLayer

show = __name__ == "__main__"

# Bell state EGraph
edge = [(0, 0, 0), (1, 1, 1)]
dims = (1, 1, 1)
bell_state = EGraph(dims=dims)
bell_state.add_edge(*edge, color="MidnightBlue")
# Plot the bell state
bell_state.draw(color_edges="MidnightBlue", color_nodes="magenta", label="index")

if show:
    plt.show()
else:
    plt.close()

bell_state.adj_generator(sparse=False)
print("Adjacency matrix: \n", bell_state.adj_mat, "\n")

CVbell = CVLayer(bell_state, p_swap=0.5)
# Noise model for CVLayer
model = {"noise": "grn", "delta": 1, "sampling_order": "initial"}
CVbell.apply_noise(model)
CVbell.measure_hom("p", [0])
CVbell.measure_hom("q", [1])

# Some parameters to make the graph of CVbell.
CV_graph_params = {"legend": True, "title": True}

CVbell.draw(label="hom_val_p", **CV_graph_params)
CVbell.draw(label="hom_val_q", **CV_graph_params)

if show:
    plt.show()
else:
    plt.close()

print("\nNodes :", bell_state.nodes.data())
print("Edges :", bell_state.edges.data(), "\n")
print("p indices: ", CVbell.p_inds)
print("GKP indices: ", CVbell.GKP_inds)
print("\nSymplectic CZ matrix: ", CVbell.SCZ(), "\n")

viz.plot_binary_mat_heat_map(CVbell.SCZ(), show)

if not show:
    plt.close()
