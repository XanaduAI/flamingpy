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
import matplotlib.pyplot as plt

from flamingpy.codes import alternating_polarity, SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders import decoder as dec
from flamingpy.noise import IidNoise


def decode_surface_code(distance, boundaries, ec, noise, decoder="MWPM", draw=True, show=False):
    """Example of instantiating, applying noise, decoding, recovering, and
    visualizing this procedure for the measurement-based surface code."""

    # Code and code lattice (cluster state)
    RHG_code = SurfaceCode(
        distance=distance,
        ec=ec,
        boundaries=boundaries,
        polarity=alternating_polarity,
    )

    # Noise model: set to "dv" for iid Z errors; "cv" for Gaussian Random Noise
    # over a GKP/sqeezed state architecture
    if noise == "cv":
        # CV (inner) code / state preparation
        p_swap = 0.05  # probability of having squeezed states (the rest are GKPs)
        CVRHG = CVLayer(RHG_code, p_swap=p_swap)
        # Noise model
        delta = 0.1  # GKP squeezing parameter
        cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}
        # Apply noise, measure syndrome, translate to bit values
        CVRHG.apply_noise(cv_noise)
        CVRHG.measure_hom("p", RHG_code.all_syndrome_inds)
        dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
        # Decoding options
        if decoder == "MWPM":
            weight_options = {
                "method": "blueprint",
                "integer": True,
                "multiplier": 100,
                "delta": delta,
            }
        else:
            weight_options = None
        decoder = {"inner": "basic", "outer": decoder}

    if noise == "dv":
        # i.i.d Pauli Z errors with probability p_Z
        p_Z = 0.02
        IidNoise(RHG_code, p_Z).apply_noise()
        weight_options = {"method": "uniform"}
        decoder = {"outer": decoder}

    # Drawing options
    node_colors = "state" if noise == "cv" else False
    dw = {
        "show_nodes": True,
        "color_nodes": node_colors,
        "show_recovery": True,
        "label_stabilizers": False,
        "label_boundary": False,
        "label_edges": False,
        "label": None,
        "legend": True,
        "show_title": True,
        "show_axes": True,
    }

    # Decode and plot
    c = dec.correct(
        code=RHG_code,
        decoder=decoder,
        weight_options=weight_options,
        sanity_check=True,
        draw=draw,
        drawing_opts=dw,
    )
    print(f"Success: {c}")

    if show:
        plt.show()
    else:
        plt.close()

    return c


if __name__ == "__main__":

    params = {
        # QEC code parameters
        "distance": 3,
        # Boundaries ("open" or "periodic")
        "boundaries": "open",
        # Error complex ("primal" or "dual")
        "ec": "primal",
        # Noise model: set to "dv" for iid Z errors; "cv" for Gaussian Random Noise
        # over a GKP/sqeezed state architecture
        "noise": "cv",
        # Decoder: set to "MWPM" for minimum-weight perfect matching, or
        # "UF" for Union-Find
        "decoder": "MWPM",
    }

    c = decode_surface_code(**params, show=True)
