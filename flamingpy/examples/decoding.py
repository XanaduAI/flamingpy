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
from flamingpy.utils import viz

show = __name__ == "__main__"

# DV (outer) code
distance = 3
# Boundaries ("open" or "periodic")
boundaries = "open"
# Error complex ("primal", "dual", or "both")
ec = "both"

RHG_code = SurfaceCode(
    distance=distance, ec=ec, boundaries=boundaries, polarity=alternating_polarity
)
RHG_lattice = RHG_code.graph

# CV (inner) code/state
p_swap = 0.2
CVRHG = CVLayer(RHG_lattice, p_swap=p_swap)
# Noise model
delta = 0.1
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Apply noise
CVRHG.apply_noise(cv_noise)
# Measure syndrome
CVRHG.measure_hom("p", RHG_code.all_syndrome_inds)

# Decoding options
decoder = {"inner": "basic", "outer": "MWPM"}
weight_options = {
    "method": "blueprint",
    "integer": True,
    "multiplier": 100,
    "delta": delta,
}

# Drawing options
dw = {
    "show_nodes": True,
    "color_nodes": "state",
    "label": None,
    "legend": True,
    "title": True,
    "display_axes": True,
    "label_edges": True,
    "label_cubes": False,
    "label_boundary": True,
}

# Manual decoding to plot intermediate results.
dec.assign_weights(RHG_code, **weight_options)
dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
for ec in RHG_code.ec:
    G_match = dec.build_match_graph(RHG_code, ec)
    matching = G_match.min_weight_perfect_matching()
    G_stabilizer = getattr(RHG_code, ec + "_stab_graph")

    # An integer label for each nodes in the stabilizer and matching graphs.
    # This is useful to identify the nodes in the plots.
    node_labels = {node: index for index, node in enumerate(G_stabilizer.graph)}

    # The draw_dec_graph function requires the networkx backend. Most backends implement
    # the to_nx() method to perform the conversion if needed.
    G_stabilizer.draw(title=ec.capitalize() + " stabilizer graph", node_labels=node_labels)
    ax = viz.syndrome_plot(RHG_code, ec, drawing_opts=dw, index_dict=node_labels)
    # viz.draw_matching_on_syndrome_plot(ax, matching, G_stabilizer, G_match, dw.get("label_edges"))
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
