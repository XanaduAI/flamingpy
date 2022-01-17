<h1 align="center">FT-Stack</h1>

<p align="center">
  <!-- Tests (GitHub actions) -->
  <a href="https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml">
    <img src="https://img.shields.io/github/workflow/status/XanaduAI/strawberryfields/Tests/master?logo=github&style=flat-square" />
  </a>
  <!-- PyPI -->
  <a href="https://pypi.org/project/ft-stack">
    <img src="https://img.shields.io/pypi/v/ft-stack.svg?style=popout-square" />
  </a>
  <!-- PyPI - Python Version -->
  <a href="https://pypi.org/project/ft-stack">
    <img src="https://img.shields.io/pypi/pyversions/ft-stack.svg?style=popout-square" />
  <!-- License -->
  <a href="https://www.apache.org/licenses/LICENSE-2.0">
    <img src="https://img.shields.io/pypi/l/ft-stack.svg?logo=apache&style=flat-square" />    
  </a>
</p>

<p align="center">
  FT-Stack is a hybrid Python/C++ library for efficient threshold estimations of combined DV-CV quantum codes and simulations of fault-tolerant quantum computers.
</p>

## Features
* Simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. 
* Supports encoding qubits into GKP states (more precisely, combinations of GKP and squeezed states). 
* Is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features. 
* Provides a host of visualization tools for ease of verifying correctness.

## Download and installation 

FT-Stack requires **Python 3.8 or higher**. The recommended method to download and install FT-Stack, as well as all dependencies and precompiled C++ binaries, is through `pip` and our [PyPI package](https://pypi.org/project/ft-stack). In your choice of CLI (with a Python environment activated) run the following single line:

```bash
pip install ft-stack
``` 

### Installation from Source (advanced users)

If you are a developer and wish to manipulate and test FT-Stack source code, you can install from Source and appropriately configure C++ compilers and Python environments.

#### Short option

Having cloned `ft-stack`, installed **conda**, opened Anaconda Prompt, and activated a new environment, you may simply navigate to your local `ft-stack` directory and run the lines:

```bash
pip install -e .
pip install -r dev-requirements.txt
``` 

These steps should work as long as there is already a C++ compiler installed on your machine. Otherwise, we recommend the steps in the following section.

#### Thorough (recommended) option

First, clone `ft-stack` through the Code tab above. A clean installation and configuration of compilers and Python environments can be achieved through CLIs. We recommend users employ **BASH** here (Windows users can rely on [WSL 2](https://docs.microsoft.com/en-gb/windows/wsl/install)). 

Next, install C/C++ compilers. Debian users can run `sudo apt update`, `sudo apt upgrade`, and `sudo apt-get install build-essential`, while Windows users should rely on Visual Studio C/C++ or similar. 

Now, set up your interpreters and customized virtual environments. We recommend using the light version of **Conda**, i.e. miniconda, or **virtualenv** for configuring Python environments. In particular, miniconda can be installed in BASH following this [guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Run the following to set up and switch to a new virtual env called, e.g., `ftstack`:
```bash
conda init
conda create -n ftstack python=3.9
conda activate ftstack
```
Lastly, change to the directory where FT-Stack was cloned and run:
```bash
pip install -e .
pip install -r dev-requirements.txt
``` 
The first command installs and compiles the dependencies of FT-Stack while the second installs extra dependencies for testing purposes.    

## Getting started ad basic usage

> Coming soon ...

## Performance Demos

> Coming soon ...

## Contributing to FT-Stack

We welcome all constructive contributions — simply fork an FT-Stack repository and then make a pull request (PR) containing your contributions. All contributors to FT-Stack will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, functionalities, etc.) may be listed on the arXiv preprints pertaining to the FT-Stack.

## Support

- **Source Code:** https://github.com/XanaduAI/ft-stack
- **Issue Tracker:** https://github.com/XanaduAI/ft-stack/issues

If you are experiencing any type of issue or have found bugs, please let us know by posting them our GitHub issue tracker. While we welcome and are committed to respond to all reports, please note FT-Stack is distributed with no guarantee. 

## Attribution for authors

FT-Stack is the work of [many contributors](https://github.com/XanaduAI/ft-stack/graphs/contributors). If you are using FT-Stack for research purposes, please cite the reference below:

```bash
@article{tzitrin2021,
   title={Fault-Tolerant Quantum Computation with Static Linear Optics},
   volume={2},
   ISSN={2691-3399},
   url={http://dx.doi.org/10.1103/PRXQuantum.2.040353},
   DOI={10.1103/prxquantum.2.040353},
   number={4},
   journal={PRX Quantum},
   publisher={American Physical Society (APS)},
   author={Tzitrin, Ilan and Matsuura, Takaya and Alexander, Rafael N. and Dauphinais, Guillaume and Bourassa, J. Eli and Sabapathy, Krishna K. and Menicucci, Nicolas C. and Dhand, Ish},
   year={2021},
   month={Dec}
}

```

## License

FT-Stack is **free** and **open source**, and released under the Apache License, Version 2.0.
