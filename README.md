# FT-Stack
![tests](https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml/badge.svg?branch=main)

## Threshold estimations for concatenated quantum codes
FT-Stack simulates all the error-correction steps (encoding, decoding, recovery, success check) for concatenated measurement-based quantum-error-correcting codes in fault-tolerant bosonic systems. Among these are CV-DV codes, whose inner encoding can be GKP states uniquely, or hybrid variations consisting of both Gaussian and GKP states. The package is conveniently modularized, allowing the user to easily swap encodings, decoders, and other features, and provides a host of visual tools for ease of verifying correctness.

## Installation 

### Quick installation (pre-compiled binaries)

> Coming soon ...

### Installation from source (advanced users)  

To manipulate and run FT-Stack you will need appropriately configured C/C++ compilers and Python virtual environments. Follow these recommended steps for a full installation: 

1. Clone `ft-stack` through the Code tab above.

2. Clean installation and configuration of compilers and Python envs are achievable through terminal command lines. While Windows users can rely on Visual Studio C/C++ and/or MinGW compilers and prompts for such purposes, we recommend all Windows/MacOS/Linux users use BASH to avoid known path setting and compilation issues. BASH is natively available for non-Windows users -- they can now **open a shell** and skip the rest of this installation step. 

    We recommend Windows 10 users use WSL 2. If you never used WSL, first, make sure Windows is updated. Open **Turn Windows feature on or off** from Start and mark **Virtual Machine Platform** and **Windows Subsystem for Linux**. Click OK and await instructions to restart your machine. Now from the Windows Store find and install **Ubuntu** app (you may choose other Linux distros if you wish). Open a Windows PowerShell as an administrator and issue `wsl --set-version Ubuntu 2` (replacing Ubuntu with distro you have previously selected). You can now open and enjoy a native BASH through WSL 2.        

3. If you have already set up C/C++ compilers, you can skip this step. Ubuntu users (including WSL 2 clients) run `sudo apt update && sudo apt upgrade` if you not done yet in your BASH. Now simply run the following to get all necessary compilation components:
    ```bash
    sudo apt-get install build-essential
    ```
    Other users need to get C/C++ compilers through their own package installers.

4. If you have already set up Python interpreters and customized virtual envs, you can skip this step. We strongly recommend using **miniconda** for configuring envs, which can be installed in BASH following this [guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Run the following to set up and switch to a new virtual env called `ft-stack`:
    ```bash
    conda init
    conda create -n ft-stack python=3.9
    activate ft-stack
    ```
5. Next change to the directory `ft-stack` was cloned. The following line will install and compile all required libraries: 
   ```bash
   pip install -e .
   ``` 
   If you are a developer, you also need to run `pip install -r dev-requirements.txt` to get additional dependencies. 

6. You can now switch to your Python IDE, set `ft-stack` as the interpreter, and enjoy using FT-Stack. For Windows users, popular Python IDEs such as PyCharm and Visual Studio Code support importing Python interpreters through WSL virtual envs natively. If you would like to use an IDE without WSL environments support, such as Spyder 5.2.1, you can do so easily using a GUI X-server for WSL such as [vcxsrv](https://sourceforge.net/projects/vcxsrv/) -- such IDEs are available to be installed on your Linux subsystem using `conda`. Note also, in some cases, you need to name your environment `venv` or `env` and/or add its name to `.gitignore`.

## Usage

> Coming soon ...

## Performance and benchmarks

> Coming soon ...

## Contribution and support

We welcome all constructive contributions — simply fork the FT-Stack repository, and then make a pull request (PR) containing your contribution. All contributors to FT-Stack will be listed as authors on the releases. Users who contribute significantly to the code (new plugins, new functionality, etc.) may be listed on the existing and upcoming FT-Stack arXiv paper.

While FT-Stack is distributed with NO GUARANTEE, we welcome all issue reports. If you are having any type of issues or have found bugs, please let us know by posting them on our [GitHub issue tracker](https://github.com/XanaduAI/ft-stack/issues).

## Attribution

If you are using FT-Stack for research purposes, please cite the reference below:

```bash
@misc{tzitrin2021faulttolerant,
      title={Fault-tolerant quantum computation with static linear optics}, 
      author={Ilan Tzitrin and Takaya Matsuura and Rafael N. Alexander and Guillaume Dauphinais and J. Eli Bourassa and Krishna K. Sabapathy and Nicolas C. Menicucci and Ish Dhand},
      year={2021},
      eprint={2104.03241},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

## License

> Coming soon ...
