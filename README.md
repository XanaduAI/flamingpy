# FT-Stack
![Tests](https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml/badge.svg?branch=main)

## Threshold estimations for concatenated quantum codes
FT-Stack simulates all the error-correction steps (encoding, decoding, recovery, success check) for concatenated measurement-based quantum-error-correcting codes. Among these are CV-DV codes, whose inner encoding can be GKP states uniquely, or hybrid variations consisting of both Gaussian states and GKP states. The package is conveniently modularized, allowing the user to easily swap encodings, decoders, and other features, and provides a host of visual tools for ease of verifying correctness.

## Installation  
To run FT-Stack you will need C/C++ compilers and Python virtual env, which have been appropriately set up. Follow these recommended steps for a full installation: 

> If you have already set up a virtual env and C/C++ compilers, you can skip Steps 1-4 below.

1. Clone `ft-stack` through the Code tab above.
2. If never done before, you will need to download and install C/C++ compilers. Windows users should install minimally Visual Studio C/C++ compilers and toolset -- 2019 version and above are recommended. Go to this [link](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022) to set up a minimal Visual Studio 2022 installation -- downloading IDEs and other tools are optional. Non-Windows users can install GNU compilers.
3. If never done before, you will need to set up your Python interpreter and customized virtual envs. We strongly recommend using miniconda for this purpose, which can be installed through this [link](https://docs.conda.io/en/latest/miniconda.html).
4. Open a shell terminal. For Windows users, we recommend **x64 Native Tools prompt**, which can be searched and accessed through Start. Non-Windows users can use BASH. We discourage the use of Anaconda Prompt due to the known issues in setting up system paths. Run the following to set up and switch to a new virtual env called `ft-stack`:
    ```bash
    conda init
    conda create -n ft-stack python=3.9
    activate ft-stack
    ```
5. Next change to the directory `ft-stack` was cloned. The following single line will pip install and compile all required libraries.   
    ```bash
    pip install -e .
    ```
    If you are a developer, you also need to run `pip install -r dev-requirements.txt` to get additional dependencies. 

You can now switch to your Python IDE, set `ft-stack` as the interpreter, and enjoy using FT-Stack. Note in some cases, you need to name your environment `venv` or `env` and/or add its name to `.gitignore`.
