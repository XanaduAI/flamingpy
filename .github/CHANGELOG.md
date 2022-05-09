## Release 0.6.1a3 (development release)

### New features since the last release

### Bug fixes

### Improvements

* Test were added to improve the overall test coverage. These included changes to
  `.coveragerc` as well as the refactoring of some examples to allow for proper
  imports from testing modules. Code coverage is now above 95% and
  the fail treshold is bumped accordingly. [(#14)](https://github.com/XanaduAI/flamingpy/pull/14)
* `CVLayer` has been modified to allow for instantiation with a code object
  in addition to an `EGraph`. This makes more semantic sense (applying a noise model
  to a code) and makes it easier for the user. [(#25)](https://github.com/XanaduAI/flamingpy/pull/25)
 * The sometimes failing `test_hybridize` in `test_graphstates.py` has been fixed. [(#25)](https://github.com/XanaduAI/flamingpy/pull/25)
 * PR template has been changed to inform user about 95% + codecov requirement. [(#25)](https://github.com/XanaduAI/flamingpy/pull/25)
 * Introduced `codecov.yml` to customize codecov automated tests. For this version, we have added a `threshold: %0.01` to avoid undesired failures due to just removing a few lines, etc. [(#25)](https://github.com/XanaduAI/flamingpy/pull/25)
* Several changes were made to improve the visualization of MWPM decoding for debugging and understanding purposes. [(#23)](https://github.com/XanaduAI/flamingpy/pull/23)
  * A function was added to the `viz` module and new options added to the `correct` function in the decoder module to be able to simply plot all decoding objects (stabilizer graph, matching graph, matching, syndrome plot) in sync with the actual error correction trial. 
  * Appearances and presence of node labels (specifically the virtual nodes of the matching graph) were fixed. 
  * The `label_cubes` argument was renamed to the more accurate `label_stabilizers`.
  * The argument `show_matching` was added to the drawing options be able to turn the matching plot on or off.
  * One can now plot a non-NetworkX matching graph (by an automatic conversion to a NetworkX graph).
  * The above changes allowed for a significant simplification to the decoding example.

### Documentation changes

### Contributors

This release contains contributions from (in alphabetical order):

See full commit details ...

---


## Release 0.6.1a3 (current release)

### New features since the last release

* Fixed drawing of stabilizer graph for zero syndrome: [(#9)](https://github.com/XanaduAI/flamingpy/pull/9)(backward incompatible)
  * Previously, the drawing function for a stabilizer graph relied on a non-documented feature. That is, it was assumed that when building the matching graph, all edges of a Networkx-based stabilizer graph were assigned a weight. This, however, was not a fair assumption for many reasons.
  * As a solution, we have added a new method to the `SurfaceCode` class to draw the primal or dual stabilizer graph, which makes sure that each edge has a weight. Now, using that method, the user does not have to rely on unfair assumptions.
  * Furthermore, we added a quick check to not add any edges to the matching graph when the syndrome is trivial. In this case, the cost of decoding should be almost zero.
* Pauli noise: have added a new noise model sampling i.i.d Z error for each qubit. [(#8)](https://github.com/XanaduAI/flamingpy/pull/8)(backward incompatible)

### Improvements

* A large number of linting corrections were made to improve the overall pylint report. These were mostly minor but essential modifications including restructuring, re-coding, optimization, updating `.pylintrc`, adding `.coveragerc`. The code quality score is improved to `A` for the released version. Check ["Files changed"](https://github.com/XanaduAI/flamingpy/pull/11/files) for details. [(#11)](https://github.com/XanaduAI/flamingpy/pull/11)

* `TODO` comments have been removed exclusively from files with low code quality grades. The Dev team has created tickets to be actioned for all removed TODO comments on separate (private) FlamingPy boards. [(#11)](https://github.com/XanaduAI/flamingpy/pull/11)

* `CONTRIBUTING.md`, `code_of_conduct.md`, and `CHANGLELOG.md` were imported and modified from the Strawberry Fields project. Dev team plans to extend these with customized details in future PRs. [(#11)](https://github.com/XanaduAI/flamingpy/pull/11)

### Documentation changes

* The new Xanadu Sphinx theme has been applied. Currently, most Xanadu OSS projects include their own version of the Xanadu Sphinx theme; however, the Xanadu Sphinx Theme repository is now publicly available and is the preferred way to access the Xanadu CSS theme and Sphinx directives. [(#17)](https://github.com/XanaduAI/flamingpy/pull/17)
  * Deleted the doc/xanadu_theme directory
  * Updated doc/requirements.txt and doc/conf.py to reference and use the (centralized) Xanadu Sphinx Theme.
  * Replaced the Quantum Error Correction, Install, and FAQ static HTML files with reST ones.

* Updated old FT-Stack links in docs header to correct FlamingPy pages. [(#7)](https://github.com/XanaduAI/flamingpy/pull/7)

### Contributors

This release contains contributions from (in alphabetical order):

[Mikhail Andrenkov](https://github.com/Mandrenkov), [Sebastián Duque Mesa](https://github.com/sduquemesa), Nariman Saadatmand, [Maxime Tremblay](https://github.com/maxtremblay), [Ilan Tzitrin](https://github.com/ilan-tz)

See full commit details [here](https://github.com/XanaduAI/flamingpy/compare/v0.4.9a1...v0.6.1a3).


## Release 0.4.9a1

### Improvements since the last release

* Relative paths cannot be used in README.md logos and were replaced with GitHub-hosted links. [(#5)](https://github.com/XanaduAI/flamingpy/pull/5)

* C++ imports are now placed within `try` blocks to avoid interrupting non-compiled installations, such as the ones currently used by readthedocs. [(#5)](https://github.com/XanaduAI/flamingpy/pull/5)

* Code coverage (quality) score was improved to a more acceptable `B-` level. [(#5)](https://github.com/XanaduAI/flamingpy/pull/5)

### Bug fixes

* Fixed a bug in [`pull_request_template.md`](https://github.com/XanaduAI/flamingpy/pull/2/commits/e30f2cb65daffece08b193ffc4b8fe7a8d90b90e). The template was not loading properly due to a whitespace problem. [(#2)](https://github.com/XanaduAI/flamingpy/pull/2)

* Fixed a bug in [`simulations.py`](flamingpy/simulations.py) and related examples. See [here](https://github.com/XanaduAI/flamingpy/commit/771b0e66e5471c3696ac2e779a2df1cc49e5e684) for commit details. [(#6)](https://github.com/XanaduAI/flamingpy/pull/6)

### Documentation changes

* Making Documentation more usable and consistent with other Xanadu projects [(#5)](https://github.com/XanaduAI/flamingpy/pull/5):
   * API details and inheritance diagrams should be now correctly displayed.
   * "Edit on Github" links were fixed
   * The general style and section structures made more consistent with the company requirements and other packages such as StrawberryFields.
   * Fixed the documentation links in `README.md`
   * Minor updates to `doc/conf.py`, `doc/dev_requirements.txt`, and `doc/Makefile`.

### Contributors

This release contains contributions from (in alphabetical order):

Nariman Saadatmand, [Ilan Tzitrin](https://github.com/ilan-tz)

See full commit details [here](https://github.com/XanaduAI/flamingpy/compare/v0.4.6a1...v0.4.9a1).


## Release 0.4.6a1

### New features since the last private release

* This is the initial public release started from the private template and our sister project [FT-Stack](https://github.com/XanaduAI/ft-stack).
* The first Cython function for Monte Carlo sampling, mostly to provide cythonization samples and testbeds, has been added. See [`cpp_mc_loop.pyx`](flamingpy/cpp/cpp_mc_loop.pyx) and [`simulations.py`](flamingpy/benchmarks/simulations.py) for details. (backward incompatible)

### Improvements

* More options for Installation from Source:
`setup.py` was updated to provide a no-compilation option for only installing (purely) Python libraries and separate options to compile `cmake` and `cython`-based codes. See the new [README.md](https://github.com/XanaduAI/ft-stack/blob/mc-cpp/README.md) for details.

### Contributors

This release contains contributions from (in alphabetical order):

Nariman Saadatmand, [Ilan Tzitrin](https://github.com/ilan-tz)
