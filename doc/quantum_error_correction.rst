.. default-role:: math

.. _quantum-error-correction:

Quantum Error Correction
========================

.. raw:: html

    <style>
        .breadcrumb {
            display: none;
        }
        h1 {
            text-align: center;
            margin-bottom: 15px;
        }
        p.lead.grey-text {
            margin-bottom: 30px;
        }
        .footer-relations {
            border-top: 0px;
        }
    </style>


    <script type="text/javascript">
        $(function(){
            // Change active navbar element to "Quantum Error Correction".
            $(".nav-item.active").removeClass("active");
            $(".nav-item a:contains('Quantum Error Correction')").parent().addClass("active");
        });
    </script>


Introduction
------------
Quantum systems spontaneously interact with their environment, manifesting as errors in a quantum computation. Quantum error correction, and the theory of fault-tolerant quantum computation in general, has been developed to ensure reliable logical operations in the presence of noise.

The key idea behind quantum error correction is the use of redundancy to represent or *encode* quantum information. The simplest example of an encoding is a *repetition code*, wherein a qubit is associated with a two-dimensional space spanned by `|000\rangle` and `|111\rangle`, i.e., `|\psi\rangle = \alpha|000\rangle + \beta|111\rangle`. For simplicity let us consider the effect of a generic single-qubit Pauli error `P` on `|\psi\rangle`. `P` is incapable of transforming one encoded state to another. Moreover, a simple majority vote can infer the state of the logical qubit given `P|\psi\rangle`. This is quantum error correction in its simplest form.

In general, a `[[n,k]]` quantum error-correcting code `Q` is identified by a `2^{k}` dimensional subspace of a `2^{n}` dimensional Hilbert space. It is often convenient to restrict to a subset of codes that admit a succinct description — as a common eigenspace of a set of mutually commuting Pauli operators — known as *stabilizer* codes. In general, a `[[n,k]]` stabilizer code `Q` is a common eigenspace with eigenvalue `+1`, of a set of `n-k` independent commuting Pauli matrices, which are called stabilizer generators. For our earlier example of a repetition code, it is easy to verify that the corresponding stabilizer generators are `Z_1 Z_2` and `Z_2 Z_3`.

Stabilizer codes are, by far, the most widely used class of codes. Popular examples include the |5q_code|, |Steane_code|, |Shor_code| and |topological codes|. An in-depth survey of stabilizer codes for quantum error correction can be found in Daniel Gottesman's |Gottesman_thesis| and the two popular texts |Nielsen_Chuang| by Michael Nielsen and Isaac Chuang, and |Lidar_Brun| by Daniel Lidar and Todd Brun.

In the stabilizer formalism, quantum error correction comprises two steps. The first is *syndrome extraction*: detecting the presence of an error using quantum measurements. The second, *decoding*, refers to inferring the error from the observed measurement outcomes.

For Pauli errors, inferring the error is tantamount to *recovery*, since errors square to the identity. Ideally, encoded states `|\psi\rangle` should satisfy `S_i |\psi\rangle = |\psi\rangle`. Any deviation from this rule indicates the presence of errors on the encoded state. In other words, a Pauli error `P` on the encoded state will result in `S_i (P |\psi\rangle) = s_i (P |\psi\rangle)`, whenever `S_i P = (-1)^{s_i} P S_i`. The binary sequence `s_1 … s_{n-k}` is called the *error syndrome*. It provides important information on the unknown Pauli error `P`.

Decoding involves inferring the most likely error given its error syndrome. This is a classical algorithm that performs a statistical inference given the physical error model. Popular strategies include *minimum weight* decoding and *maximum likelihood* decoding. For surface codes, the former is implemented by |perfect_matching|. Other types of decoding include the |Union_Find| and the |RG_decoder| (see also |q_Self_Correction|).


.. |5q_code| raw:: html

   <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.5811" target="_blank">5-qubit</a>

.. |Steane_code| raw:: html

   <a href="https://royalsocietypublishing.org/doi/10.1098/rspa.1996.0136" target="_blank">Steane</a>

.. |Shor_code| raw:: html

   <a href="https://journals.aps.org/pra/abstract/10.1103/PhysRevA.52.R2493" target="_blank">Shor’s code</a>

.. |topological_codes| raw:: html

   <a href="https://arxiv.org/abs/1311.0277" target="_blank">topological codes</a>

.. |Gottesman_thesis| raw:: html

   <a href="https://arxiv.org/abs/quant-ph/9705052" target="_blank">Ph.D thesis</a>

.. |Nielsen_Chuang| raw:: html

   <a href="https://www.cambridge.org/highereducation/books/quantum-computation-and-quantum-information/01E10196D0A682A6AEFFEA52D53BE9AE#overview" target="_blank">Quantum Computation and Quantum Information</a>

.. |Lidar_Brun| raw:: html

   <a href="https://www.cambridge.org/us/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-error-correction" target="_blank">Quantum Error Correction</a>

.. |perfect_matching| raw:: html

   <a href="https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.032324" target="_blank">computing a perfect matching</a>

.. |Union_Find| raw:: html

   <a href="https://quantum-journal.org/papers/q-2021-12-02-595/" target="_blank">Union-Find decoder</a>

.. |RG_decoder| raw:: html

   <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.050504" target="_blank">renormalization group decoder</a>

.. |q_Self_Correction| raw:: html

   <a href="https://doi.org/10.1103/PhysRevLett.111.200501" target="_blank">this reference</a>
