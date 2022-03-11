<p align="center">
  <a href="https://github.com/XanaduAI/ft-stack">
    <img width=40% src="https://s10.gifyu.com/images/imagee72e454028813818.png">
  </a>
</p>

<p align="center">
  <!-- Tests (GitHub actions) -->
  <a href="https://github.com/XanaduAI/ft-stack/actions/workflows/build_tests.yaml">
    <img src="https://img.shields.io/github/workflow/status/XanaduAI/strawberryfields/Tests?label=build%20%26%20tests&style=flat-square" />
  </a>
  <!-- CodeFactor -->
  <a href="https://www.codefactor.io/repository/github/xanaduai/strawberryfields">
    <img src="https://img.shields.io/codefactor/grade/github/XanaduAI/strawberryfields/master?style=flat-square" />
  </a>
  <!-- PyPI (Python Version) -->
  <a href="https://pypi.org/project/flamingpy">
    <img src="https://img.shields.io/pypi/pyversions/flamingpy.svg?style=flat-square" />
  </a>
  <!-- PyPI -->
  <a href="https://pypi.org/project/flamingpy">
    <img src="https://img.shields.io/pypi/v/flamingpy.svg?style=flat-square" />
  </a>
  <!-- CodeCov -->
  <a href="https://codecov.io/gh/XanaduAI/strawberryfields">
    <img src="https://img.shields.io/codecov/c/github/xanaduai/strawberryfields/master.svg?style=popout-square" />
  </a>
  <!-- License -->
  <a href="https://www.apache.org/licenses/LICENSE-2.0">
    <img src="https://img.shields.io/pypi/l/flamingpy.svg?logo=apache&style=flat-square" />    
  </a>
</p>

<p align="center">
 <a href="https://github.com/XanaduAI/ft-stack">FlamingPy</a> is a cross-platform Python library with several backends for efficient simulations of error correction in fault-tolerant quantum computers.
</p>

## Features

<img src="https://s10.gifyu.com/images/ftstack_featured.jpg" width="330px" align="right">

* Simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. 
* Supports encoding qubits into GKP states (more precisely, combinations of GKP and squeezed states). 
* Is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features. 
* Provides a host of visualization tools for ease of verifying correctness.
  
## Download and installation 

FlamingPy requires **Python 3.8 or 3.9**. The recommended method to download and install FlamingPy, as well as all dependencies and precompiled C++ binaries, is through `pip` and our [PyPI package](https://pypi.org/project/flamingpy). In your choice of CLI (with a Python environment activated) run the following single command:

```bash
python -m pip install -i https://test.pypi.org/simple/ flamingpy # TODO: TestPyPI cannot properly install dependencies. Please run `python -m pip install matplotlib networkx retworkx numpy pandas scipy thewalrus --upgrade` beforehand manually. Remove this comment when we moved to PyPI.
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

If you encountered a CMake error, you may need to (re-)install it through `conda install cmake` or other means before re-attempting the above. Furthermore, you may wish to try `conda install git`. For more detailed instructions and recommendations, including how to configure your environments, compilers and resolve errors, see our Frequently Encountered Errors page in the documentation [coming soon].

## Getting started and basic usage

> Coming soon ...

## Performance Demos

> Coming soon ...

## Contribution

We welcome new contributions -- simply fork the FlamingPy repository and make a pull request (PR) containing your contribution. All contributers to FlamingPy will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, functionalities, etc.) may be listed on the arXiv preprints for the FlamingPy. See our
changelog for more details.

## Support

- **Source Code:** https://github.com/XanaduAI/ft-stack
- **Issue Tracker:** https://github.com/XanaduAI/ft-stack/issues

If you are having issues, please let us know by posting the issue on our GitHub issue tracker.

## Attribution for authors

FlamingPy is the work of [many contributors](https://github.com/XanaduAI/ft-stack/graphs/contributors). 

If you are doing research using FlamingPy, please cite our paper below:

> Ilan Tzitrin, Takaya Matsuura, Rafael N. Alexander, Guillaume Dauphinais, J. Eli Bourassa, Krishna K. Sabapathy, Nicolas C. Menicucci, and Ish Dhand,
> Fault-Tolerant Quantum Computation with Static Linear Optics, PRX Quantum, Vol. 2, No. 4, 2021, 
> [DOI:10.1103/prxquantum.2.040353](http://dx.doi.org/10.1103/PRXQuantum.2.040353) 

## License

FlamingPy is **free** and **open source**, and released under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
