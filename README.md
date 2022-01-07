# FT-Stack
![tests](https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml/badge.svg?branch=main)

## Simulations of fault-tolerant quantum computation
The FT-Stack simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. Among these are the GKP code (in particular, combinations of GKP and squeezed states) concatenated with the surface code. The package is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features, and provides a host of visualization tools for ease of verifying correctness.

## Installation 

### Method #1: Quick installation (pre-compiled binaries)

Our recommended method to download and install FT-Stack is through our PyPi package. In your choice of CLI with a Python environment activated, run the following single-line. Pip will automatically install FT-Stack and fetch all pre-compiled binaries it requires:

> Coming soon ...

### Method #2: Installation from Source through Conda packages  

If you are a developer and need to manipulate and test FT-Stack, you will need to install from Source and appropriately configure C++ compilers and Python environments. 

If you use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html) to manage your virtual environments, installation from Source can be performed easily through a few commands. In your CLI, change to FT-Stack directory, activate your **Conda environment** and then run:
```bash
pip install -e .
pip install -r dev-requirements.txt
``` 
The first command installs and compiles the dependencies of FT-Stack while the second installs extra dependencies for testing purposes.

### Method #3: Full installation from Source using BASH (advanced users only)

1. Clone `ft-stack` through the Code tab above.

2. Clean installation and configuration of compilers and Python environments are achievable through CLIs. While Windows users can rely on Visual Studio C/C++ and/or MinGW compilers and prompts for such purposes, we recommend all Windows/MacOS/Linux users employ **BASH** for concreteness and to avoid some known path-setting and compilation issues. BASH is natively available for non-Windows users -- we recommend Windows 10 users use [WSL 2](https://docs.microsoft.com/en-gb/windows/wsl/install).

3. If you have already set up C/C++ compilers, you can skip the following step. Ubuntu users (including WSL 2 clients) can run `sudo apt update && sudo apt upgrade` in BASH if not already done so. Now simply run `sudo apt-get install build-essential` to get all necessary compilation tools. Non-Debian users need to get C/C++ compilers through their own package installers.

4. If you have already set up Python interpreters and customized virtual environments, you can skip the following step. We recommend using the light version of **Conda**, i.e. miniconda, or **virtualenv** for configuring Python environments. In particular, miniconda can be installed in BASH following this [guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Run the following to set up and switch to a new virtual env called `ft-stack`:
    ```bash
    conda init
    conda create -n ft-stack python=3.9
    activate ft-stack
    ```
5. Next, change to the directory where `ft-stack` was cloned. As before, running `pip install -e .` will install and compile all required libraries and you need to run `pip install -r dev-requirements.txt` to get additional dependencies. 

6. You can now switch to your Python IDE, set `ft-stack` as your interpreter, and enjoy using FT-Stack. For Windows users, popular Python IDEs such as PyCharm and Visual Studio Code support importing Python interpreters through WSL environments natively. If you would like to use an IDE without WSL environments support, such as Spyder 5.2.1, you can do so using a GUI X-server for WSL such as [vcxsrv](https://sourceforge.net/projects/vcxsrv/) -- such IDEs are available to be installed on your Linux subsystem using `conda`.

## Usage

> Coming soon ...

## Performance and benchmarks

> Coming soon ...

## Contribution and support

We welcome all constructive contributions — simply fork an FT-Stack repository, and then make a pull request (PR) containing your contributions. All contributors to FT-Stack will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, functionalities, etc.) may be listed on the existing and upcoming FT-Stack arXiv papers.

While the FT-Stack is distributed with no guarantee, we welcome all issues and bug reports. If you are experiencing any type of issue or have found bugs, please let us know by posting them on our [GitHub issue tracker](https://github.com/XanaduAI/ft-stack/issues).

## Attribution

If you are using FT-Stack for research purposes, please cite the reference below:

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

> Coming soon ...
