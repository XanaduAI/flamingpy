# FT-Stack
![tests](https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml/badge.svg?branch=main)

## Simulations of fault-tolerant quantum computation
The FT-Stack simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. Among these are the GKP code (in particular, combinations of GKP and squeezed states) concatenated with the surface code. The package is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features, and provides a host of visualization tools for ease of verifying correctness.

## Installation 

### Quick installation (pre-compiled binaries)

> Coming soon ...

### Installation from Source (advanced users)  

To manipulate and run FT-Stack you will need appropriately configured C/C++ compilers and Python virtual environments. Follow these recommended steps for a full installation: 

1. Clone `ft-stack` through the Code tab above.

2. Clean installation and configuration of compilers and Python envs are achievable through command line interfaces. While Windows users can rely on Visual Studio C/C++ and/or MinGW compilers and prompts for such purposes, we recommend all Windows/MacOS/Linux users employ **BASH** for concreteness and to avoid some known path-setting and compilation issues. BASH is natively available for non-Windows users -- they can now open a terminal and skip the rest of this installation step.

    We recommend Windows 10 users use WSL 2. If you have never used WSL, first make sure that Windows is up to date via Windows Update. Then, open **Turn Windows feature on or off** from Start and mark **Virtual Machine Platform** and **Windows Subsystem for Linux**. Click OK and await instructions to restart your machine. Now, from the **Microsoft Store**, find and install the **Ubuntu** app (you may choose other Linux distros if you wish). Open a **Windows PowerShell** as an administrator and issue `wsl --set-version Ubuntu 2` (replacing Ubuntu with distro you have previously selected). You can now open and enjoy a native BASH through WSL 2.

3. If you have already set up C/C++ compilers, you can skip the following step. Ubuntu users (including WSL 2 clients) can run `sudo apt update && sudo apt upgrade` in BASH if not already done so you not done yet in your BASH. Now simply run the following to get all necessary compilation tools:
    ```bash
    sudo apt-get install build-essential
    ```
    Other users need to get C/C++ compilers through their own package installers.

4. If you have already set up Python interpreters and customized virtual envs, you can skip the following step. We strongly recommend using **miniconda** for configuring envs, which can be installed in BASH following this [guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Run the following to set up and switch to a new virtual env called `ft-stack`:
    ```bash
    conda init
    conda create -n ft-stack python=3.9
    activate ft-stack
    ```
5. Next, change to the directory where `ft-stack` was cloned. The following line will install and compile all required libraries:
   ```bash
   pip install -e .
   ``` 
   If you are a developer, you also need to run `pip install -r dev-requirements.txt` to get additional dependencies. 

6. You can now switch to your Python IDE, set `ft-stack` as your interpreter, and enjoy using FT-Stack. For Windows users, popular Python IDEs such as PyCharm and Visual Studio Code support importing Python interpreters through WSL virtual envs natively. If you would like to use an IDE without WSL environments support, such as Spyder 5.2.1, you can do so using a GUI X-server for WSL such as [vcxsrv](https://sourceforge.net/projects/vcxsrv/) -- such IDEs are available to be installed on your Linux subsystem using `conda`. In some cases, you may need to name your environment `venv` or `env` and/or add its name to `.gitignore`.

## Usage

> Coming soon ...

## Performance and benchmarks

> Coming soon ...

## Contribution and support

We welcome all constructive contributions — simply fork an FT-Stack repository, and then make a pull request (PR) containing your contributions. All contributors to FT-Stack will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, functionalities, etc.) may be listed on the existing and upcoming FT-Stack arXiv papers.

While the FT-Stack is distributed with no guarantee, we welcome all issue and bug reports. If you are experiencing any type of issue or have found bugs, please let us know by posting them on our [GitHub issue tracker](https://github.com/XanaduAI/ft-stack/issues).

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
