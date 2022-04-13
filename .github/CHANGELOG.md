## Release 0.6.0a3 (development release)

### New features since the last release

### Bug fixes

### Housekeeping changes

* A large number of linting corrections were made to improve the overall pylint report. These were mostly minor but essential modifications including restructuring, re-coding, optimization, updating `.pylintrc`, adding `.coveragerc`. The code quality score is improved to `A` for the released version. Check ["Files changed"](https://github.com/XanaduAI/flamingpy/pull/11/files) for details.
* `TODO` comments have been removed exclusively from files with low code quality grades. The Dev team has created tickets to be actioned for all removed TODO comments on separate (private) FlamingPy boards.
* `CONTRIBUTING.md`, `code_of_conduct.md`, and `CHNAGLELOG.md` were imported and modified from the StrawberryFields project. Dev team plans to extend these with customized details in future PRs.

### Documentation changes

### Contributors

This release contains contributions from (in alphabetical order) ... See full commit details ...

---

## Release 0.4.9a1 (current release)

<h3>Housekeeping changes since the last release</h3>

* Minor updates to `README.md` following the start of the public project. This included:
   * Relative paths cannot be used in README.md logos and were replaced with Github-hosted links.
* C++ imports are now placed within `try` blocks to avoid interrupting non-compiled installations, such as the ones currently used by readthedocs.
* Code coverage (quality) score was improved to a more acceptable `B-` level.

<h3>Bug fixes</h3>

* Fixed a bug in [pull_request_template.md](https://github.com/XanaduAI/flamingpy/pull/2/commits/e30f2cb65daffece08b193ffc4b8fe7a8d90b90e). The template was not loading properly due to a whitespace problem.
* Fixed a bug in [`simulations.py`](flamingpy/simulations.py) and related examples. See [here](https://github.com/XanaduAI/flamingpy/commit/771b0e66e5471c3696ac2e779a2df1cc49e5e684) for details.

<h3>Documentation changes</h3>

* Making Documentation more usable and consistent with other Xanadu projects:
   * API details and inheritance diagrams should be now correctly displayed.
   * "Edit on Github" links were fixed
   * The general style and section structures made more consistent with the company requirements and other packages such as StrawberryFields.
   * Fixed the documentation links in `README.md`
   * Minor updates to `doc/conf.py`, `doc/dev_requirements.txt`, and `doc/Makefile`. 

<h3>Contributors</h3>

This release contains contributions from (in alphabetical order) @nariman87 and Ilan Tzitrin. See full commit details https://github.com/XanaduAI/flamingpy/compare/v0.4.6a1...v0.4.9a1.

## Release 0.4.6a1

<h3>New features since the last private release</h3>

* The first Cython function for Monte Carlo sampling, mostly to provide cythonization samples and testbeds, has been added. See [`cpp_mc_loop.pyx`](flamingpy/cpp/cpp_mc_loop.pyx) and [`simulations.py`](flamingpy/benchmarks/simulations.py) for detailes.

<h3>Housekeeping changes</h3>

* More options for Installation from Source: 
`setup.py` was updated to provide a no-compilation option for only installing (purely) Python libraries and separate options to compile `cmake` and `cython`-based codes. See the new [README.md](https://github.com/XanaduAI/ft-stack/blob/mc-cpp/README.md) for details.

<h3>Contributors</h3>

This release contains contributions from (in alphabetical order) @nariman87 and Ilan Tzitrin.
