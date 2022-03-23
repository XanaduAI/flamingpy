![flamingpy_logo_light](https://user-images.githubusercontent.com/25132802/159598111-fcf6b75a-26a0-4d24-b267-d9d7597bdf39.svg#gh-light-mode-only)
![flamingpy_logo_dark](https://user-images.githubusercontent.com/25132802/159598097-6a16733c-a954-49ba-a29c-ce469ae19fcc.svg#gh-dark-mode-only)

<p align="center">
  <!-- Tests (GitHub actions) -->
  <a href="https://github.com/XanaduAI/flamingpy/actions/workflows/build_tests.yaml">
    <img src="https://github.com/XanaduAI/flamingpy/actions/workflows/build_tests.yaml/badge.svg" />
  </a>
  <!-- ReadTheDocs -->
  <a href="https://flamingpy.readthedocs.io">
    <img src="https://img.shields.io/readthedocs/flamingpy.svg" />
  </a>
  <!-- CodeFactor -->
  <a href="https://www.codefactor.io/repository/github/xanaduai/flamingpy">
    <img src="https://img.shields.io/codefactor/grade/github/XanaduAI/flamingpy/main" />
  </a>
  <!-- CodeCov -->
  <a href="https://codecov.io/gh/XanaduAI/flamingpy">
    <img src="https://codecov.io/gh/XanaduAI/flamingpy/branch/main/graph/badge.svg?token=3FUq4JZL7X" />
  </a>
  <!-- PyPI (Python Version) -->
  <a href="https://pypi.org/project/flamingpy">
    <img src="https://img.shields.io/pypi/pyversions/flamingpy.svg" />
  </a>
  <!-- PyPI -->
  <a href="https://pypi.org/project/flamingpy">
    <img src="https://img.shields.io/pypi/v/flamingpy.svg" />
  </a>
  <!-- License -->
  <a href="https://www.apache.org/licenses/LICENSE-2.0">
    <img src="https://img.shields.io/pypi/l/flamingpy.svg?logo=apache" />    
  </a>
</p>

<p align="center">
 <a href="https://flamingpy.readthedocs.io/en/latest/">FlamingPy</a> is a cross-platform Python library with a variety of backends for efficient simulations of error correction in fault-tolerant quantum computers.
</p>

## Features

<img src="doc/_static/RHG_matching.svg" width="330px" align="right">

* Simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. 
* Supports encoding qubits into GKP states (more precisely, combinations of GKP and squeezed states). 
* Is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features. 
* Provides a host of visualization tools for ease of verifying correctness.
  
## Download and installation 

FlamingPy requires **Python 3.8 or above**. The recommended method to download and install FlamingPy, as well as all dependencies and precompiled C++ binaries, is through `pip` and our [PyPI package](https://pypi.org/project/flamingpy). In your choice of CLI (with a Python environment activated) run the following single command:

```bash
python -m pip install flamingpy
``` 

#### Installation from Source (advanced users)

If you are a developer and wish to manipulate and test FlamingPy source code, you can install the project from Source. First, clone FlamingPy through the Code tab above. Then, create and activate a new virtual environment (if you prefer using an existing environment, you may need to uninstall existing FlamingPy builds). If you use **Conda**, for example, you may run the following:

```bash
conda create -n flamingpy python=3.8
conda activate flamingpy
```

Finally, change to the directory where FlamingPy was cloned and run:

```bash
python -m pip install -r dev_requirements.txt
python setup.py develop # only installs Python libraries
python setup.py build_cython --inplace # [OPTIONAL] compiles Cython-based backends
python setup.py build_cmake --inplace # [OPTIONAL] compiles CMake-based backends
``` 

The purpose of the commands is as follows:
- The first command installs dependencies for building the project and testing purposes, and can be skipped if already satisfied. 
- The second command (develop) installs FlamingPy Python libraries without the compiling the optional backends. 
- The next optional commands compile various FlamingPy backends as required (given you have appropriate compilers pre-installed). 

If you encountered a CMake error, you may need to (re-)install it through `conda install cmake` or other means before re-attempting the above. Furthermore, you may wish to try `conda install git`. For more detailed instructions and recommendations, including how to configure your environments, compilers and resolve errors, see our Frequently Encountered Errors page in the [documentation](https://flamingpy.readthedocs.io/en/stable/help/frequently_encountered_errors.html).

## Getting started and basic usage

There is a vast literature available to understand the theoretical concepts behind FlamingPy. For a self-contained description, see Xanadu's [blueprint](https://quantum-journal.org/papers/q-2021-02-04-392/) for a fault-tolerant photonic quantum computer. You can also visit the documentation, which will be updated with more resources over time.

To see a sample of what FlamingPy can do, let us first import a few important objects:

```
from flamingpy.codes import SurfaceCode
from flamingpy.cv.ops import CVLayer
from flamingpy.decoders import correct
```

Next, let us instantiate an RHG lattice -- the measurement-based version of the surface code:

```
RHG = SurfaceCode(3)
```

The integer denotes the code distance. By default, the boundaries are set to "open". Next, let us associate the nodes in the RHG lattice with CV states:

```
CVRHG = CVLayer(RHG.graph, p_swap=0.5)
```

Now, half the lattice (on average) will be labelled a GKP state, and the other half a p-squeezed state. Next, we can apply a noise model to the states:

```
grn_model = {"noise": "grn", "delta": 0.1}
CVRHG.apply_noise(grn_model)
```

This results in Gaussian random noise model with a squeezing parameter of 0.1 to the GKP states in the lattice. We can now conduct a homodyne measurement on the lattice to measure the syndrome:

```
CVRHG.measure_hom("p", RHG.all_syndrome_inds)
```

At this point, we are ready to perform error correction on the lattice. First, we can specify some options for the decoder:

```
decoder = {"inner": "basic", "outer": "MWPM"}
```

This corresponds to a basic GKP binning function for the inner decoder, and minimum-weight perfect matching (MWPM) for the outer decoder. Lastly, we can detect and correct for errors, and print a message identifying success or failure:

```
c = correct(code=RHG, decoder=decoder)
outcome = "succeeded." * bool(c) + "failed." * (1 - bool(c))
message = "Error correction {}".format(outcome)
print(message)
```

<!-- ## Performance Demos
  -->

## Contribution

We welcome new contributions -- simply fork the FlamingPy repository and make a pull request (PR) containing your contribution. All contributors to FlamingPy will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, functionalities, etc.) may be listed on the arXiv preprints for FlamingPy. See our [release notes and changelogs](https://github.com/XanaduAI/flamingpy/releases) for more details.



## Support

If you are having issues, please let us know by posting the issue on our GitHub issue tracker.

- **Source Code:** https://github.com/XanaduAI/flamingpy
- **Issue Tracker:** https://github.com/XanaduAI/flamingpy/issues

You can also start a general discussion and connect with our community members in our [Discussions Page](https://github.com/XanaduAI/flamingpy/discussions). 

## Attribution for authors

FlamingPy is the work of [many contributors](https://github.com/XanaduAI/ft-stack/graphs/contributors). 

If you are doing research using FlamingPy, please cite our paper below:

> Ilan Tzitrin, Takaya Matsuura, Rafael N. Alexander, Guillaume Dauphinais, J. Eli Bourassa, Krishna K. Sabapathy, Nicolas C. Menicucci, and Ish Dhand,
> Fault-Tolerant Quantum Computation with Static Linear Optics, PRX Quantum, Vol. 2, No. 4, 2021, 
> [DOI:10.1103/prxquantum.2.040353](http://dx.doi.org/10.1103/PRXQuantum.2.040353) 

## License

FlamingPy is **free** and **open source**, and released under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
