# FT-Stack

![Tests](https://github.com/XanaduAI/ft-stack/actions/workflows/tests.yaml/badge.svg?branch=main)

## Threshold estimations for concatenated quantum codes
Simulate all the error-correction steps--encoding, decoding, recovery, success check--for concatenated measurement-based quantum-error-correcting codes. Among these are CV-DV codes, whose inner encoding can be GKP states uniquely, or hybrid variations consisting of both Gaussian states and GKP states. The package is conveniently modularized, allowing the user to easily swap encodings, decoders, and other features, and provides a host of visual tools for ease of verifying correctness.

## Development setup

In a virtual environment, run 
```bash
pip install .
pip install -r dev-requirements.txt
```
The first command installs the dependencies of `ft_stack`
while the second installs extra dependencies for testing.
Notice that you sometime need to add the `-e` flag to `pip install .` 
when you are adding new modules.

Depending on your Operating System,
you may need to install an extra backend for `matplotlib`.
For example, using Linux with the Plasma desktop, I needed to install Qt bindings.
```bash
pip install PyQt5
```
You should name your environment `venv` or `env` or add its name to `.gitignore`.
