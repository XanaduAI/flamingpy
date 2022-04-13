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
"""Example of instantiating, applying noise, decoding, recovering, and
visualizing this procedure for the measurement-based surface code."""

# pylint: disable=no-member

import matplotlib.pyplot as plt

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders import decoder as dec
from flamingpy.noise import IidNoise
from flamingpy.utils import viz

show = __name__ == "__main__"

# QEC code parameters
distance = 3
# Boundaries ("open" or "periodic")
boundaries = "open"
# Error complex ("primal", "dual", or "both")
ec = "primal"

# Code and code lattice (cluster state)
RHG_code = SurfaceCode(
    distance=distance,
    ec=ec,
    boundaries=boundaries,
    polarity=alternating_polarity,
    backend="retworkx",
)
RHG_lattice = RHG_code.graph

# Noise model: set to "dv" for iid Z errors; "cv" for Gaussian Random Noise
# over a GKP/sqeezed state architecture
noise = "cv"

if noise == "cv":
    # CV (inner) code / state preparation
    p_swap = 0.05  # probability of having squeezed states (the rest are GKPs)
    CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
    # Noise model
    delta = 0.1  # GKP squeezing parameter
    cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
    # Apply noise, measure syndrome, translate to bit values
    CVRHG.apply_noise(cv_noise)
    CVRHG.measure_hom("p", RHG_code.all_syndrome_inds)
    dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
    # Decoding options
    weight_options = {
        "method": "blueprint",
        "integer": True,
        "multiplier": 100,
        "delta": delta,
    }
    decoder = {"inner": "basic", "outer": "MWPM"}

if noise == "dv":
    # i.i.d Pauli Z errors with probability p_Z
    p_Z = 0.02
    IidNoise(RHG_code, p_Z).apply_noise()
    weight_options = {"method": "unit"}
    decoder = {"outer": "MWPM"}

# Drawing options
node_colors = "state" if noise == "cv" else False
dw = {
    "show_nodes": True,
    "color_nodes": node_colors,
    "label": None,
    "legend": True,
    "title": True,
    "display_axes": True,
    "label_edges": True,
    "label_cubes": False,
    "label_boundary": False,
}

# Manual decoding (to plot intermediate results).
dec.assign_weights(RHG_code, **weight_options)
for ec in RHG_code.ec:
    G_match = dec.build_match_graph(RHG_code, ec)
    matching = G_match.min_weight_perfect_matching()
    G_stabilizer = getattr(RHG_code, ec + "_stab_graph")

    # An integer label for each nodes in the stabilizer and matching graphs.
    # This is useful to identify the nodes in the plots.
    node_labels = {node: index for index, node in enumerate(G_stabilizer.nodes())}

    # The draw_dec_graph function requires the networkx backend. Most backends implement
    # the to_nx() method to perform the conversion if needed.
    RHG_code.draw_stabilizer_graph(
        ec, title=ec.capitalize() + " stabilizer graph", node_labels=node_labels
    )

    ax = viz.syndrome_plot(RHG_code, ec, drawing_opts=dw, index_dict=node_labels)
    viz.draw_matching_on_syndrome_plot(ax, matching, G_match)
    if len(G_match.graph):
        G_match.draw(title=ec.capitalize() + " matching graph", node_labels=node_labels)
    else:
        print("\nMatching graph empty!\n")

    if show:
        plt.show()
    else:
        plt.close()

# Automatic decoding
c = dec.correct(
    code=RHG_code,
    decoder=decoder,
    weight_options=weight_options,
    sanity_check=True,
)
print(f"Success: {c}")
