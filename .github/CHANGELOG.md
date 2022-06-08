## Release 0.8.2a5 (development release)

### New features since the last release
* Add functions to create different graph states (star and complete graphs, ring graphs, linear clusters, and Bell pairs) in a new module, `utils.graph_states`. [#68](https://github.com/XanaduAI/flamingpy/pull/68). (backward compatible)

### Bug fixes
* Small fix in `viz.draw_EGraph` that raised an error whenever a graph state with non-integer coordinates was plotted. [#68](https://github.com/XanaduAI/flamingpy/pull/68)

### Improvements

* Added tests for `EGraph` plots. [#60](https://github.com/XanaduAI/flamingpy/pull/60)
* Added `fig, ax` returns for the draw methods in `utils/viz.py` and some additional tests. [#55](https://github.com/XanaduAI/flamingpy/pull/55)
* Unit tests have been re-grouped in individual sub-dirs inside `tests/` based on error correction and software layers. This helps manage and target each test unit. [#70](https://github.com/XanaduAI/flamingpy/pull/70)
* `build_tests.yaml` workflow now supports executing unit tests in parallel using `pytest-xdist` package. GitHub runners have at least 2 processors, which helps speed up the pytest blocks by ~1.5 times in practice. [#70](https://github.com/XanaduAI/flamingpy/pull/70)
* Pylint is pinned to stable version `pylint==2.14.0` and added to `dev_requirements.txt`. [#76](https://github.com/XanaduAI/flamingpy/pull/76)
 * pylint `no-self-use` tags are removed as this check has been removed from pylint (see [here](https://github.com/PyCQA/pylint/issues/5502)).
* Added `.gitattributes` to the repository, so git automatically handles consistent `eol`'s for all commits and contributors across different operating systems. [#78](https://github.com/XanaduAI/flamingpy/pull/78)
* Increased the scope of `docformatter` to all `.py` files in the repository. [#79](https://github.com/XanaduAI/flamingpy/pull/79)
* Increased the scope of `black` formatter to include documentation files. [#79](https://github.com/XanaduAI/flamingpy/pull/79)
* Added automatically generated `.svg` files to gitignore. [#84](https://github.com/XanaduAI/flamingpy/pull/84)

### Documentation changes

* Mention the new graph state functions from `flamingpy.utils.graph_states` in the `run_graph_states.py` tutorial. [#68](https://github.com/XanaduAI/flamingpy/pull/68)
* Typo fix and minor changes for README file. [#80](https://github.com/XanaduAI/flamingpy/pull/80)
* non-Xanadu links now open in a new tab, while HTML references are listed scientific-style at the end of a file. [#82](https://github.com/XanaduAI/flamingpy/pull/82)
* Changed the math rendering Sphinx to MathJax (before equations were rendered as png). [#84](https://github.com/XanaduAI/flamingpy/pull/84)


### Contributors

This release contains contributions from (in alphabetical order):

[Joost Bus](https://github.com/soosub), [Sebastián Duque Mesa](https://github.com/sduquemesa), [Luis Mantilla](https://github.com/BestQuark), Nariman Saadatmand, [WingCode](https://github.com/WingCode) 

See full commit details ...

---

## Release 0.8.2a5 (current release)

### New features since the last release

* Node and edge coloring can now be done based on any attribute and personalized colors can be defined via a dictionary: [#32](https://github.com/XanaduAI/flamingpy/pull/32) (backward incompatible)
* The `EGraph` plot legend is not limited to the "state" attribute of the node but to any attribute. [#32](https://github.com/XanaduAI/flamingpy/pull/32) (backward incompatible)
* The `dims` attribute of `EGraph` has been removed. Its function is replaced by the `dimensions` parameter that is passed to the `draw_EGraph` method. This method does not require the `EGraph` to have a `dims` attribute defined anymore. [#42](https://github.com/XanaduAI/flamingpy/pull/42) (backward incompatible)
* Our frontend simulator script, [`simulations.py`](flamingpy/simulations.py), now supports simple and highly-scalable MPI jobs through `mpi4py` libraries in a non-intrusive manner. The users who do **not** have or want MPI, can run `simulations.py` single-threaded as per usual without facing any errors. MPI users can speed up Monte Carlo samplings in EC steps virtually up to as many processors they can throw at it. The script support jobs both on local machines and large-scale clusters. [#47](https://github.com/XanaduAI/flamingpy/pull/47) (backward compatible)
  * MPI users on their local machines can simply run the following for a 4-processor job:
  `mpirun -np 4 python flamingpy/simulations.py`

### Bug fixes
* Fixed the class inheretance diagram displayed in `fp.codes`. [#41](https://github.com/XanaduAI/flamingpy/pull/41)

### Improvements

* The `draw_EGraph` function is refactored. [#32](https://github.com/XanaduAI/flamingpy/pull/32)
  * This reduces the function complexity; ensures nodes, edges, and general plot attributes are handled in different places; and allows for better code maintenance and readability.
  * `display_axes` is changed to `show_axes` for consistency.
* `xlim` in `viz.plot_Z_err_cond` is adjusted to the relevant domain when plotting the central peak. [#33](https://github.com/XanaduAI/flamingpy/pull/33)
* Added `fig, ax` returns for the draw methods in `utils/viz.py`. [#33](https://github.com/XanaduAI/flamingpy/pull/33)
* Both upper and lower axes limits can now be specified for `EGraph` plots. [#42](https://github.com/XanaduAI/flamingpy/pull/42)
* Improvements to the decoding example. [#44](https://github.com/XanaduAI/flamingpy/pull/44)
  * Rename function and add dosctring.
  * Decrease the size of markers for plotting stabilizer nodes.
  * Improve the way to scatter stabilizers via specifying indices.
* Improved codefactor score for several key files. [#51](https://github.com/XanaduAI/flamingpy/pull/51)
* Pandas is removed from the package requirements. [#63](https://github.com/XanaduAI/flamingpy/pull/63)
* `mpi4py` is **not** a development requirement for Linux users. [#64](https://github.com/XanaduAI/flamingpy/pull/64)
* CI test check that code executes properly with and without MPI. [#64](https://github.com/XanaduAI/flamingpy/pull/64)

### Documentation changes

* A pipeline for adding tutorials to the docs was introduced. [#24](https://github.com/XanaduAI/flamingpy/pull/24)
  * To add a tutorial, use the ``gallery-item`` directive from the ``xanadu-sphinx-theme``. For the new document to be compiled its filename should start with `run_`. Thumbnails will be created out of the first figure generated and stored in `tutorials/_out/images/thumb/` with the same name of the tutorial prepended with `sphx_glr_`.
* Brief tutorials about graph states and error correction were added. [#24](https://github.com/XanaduAI/flamingpy/pull/24)
* An introduction to quantum error correction was added. [#24](https://github.com/XanaduAI/flamingpy/pull/24)
* Added UML class and package diagrams for `fp` page. [#41](https://github.com/XanaduAI/flamingpy/pull/41)
* Improved class inheritance diagram for `fp.codes`, `fp.cv`, and `fp.decoders`. [#41](https://github.com/XanaduAI/flamingpy/pull/41)
* Added `libopenmpi-dev` package to the apt list of `.readthedoc.yml` to allow documentation successful builds after adding recent `mpi4py` requirements. [#59](https://github.com/XanaduAI/flamingpy/pull/59)
* Adds a section to `guide_for_devs.rst` explaining how to install and use MPI along with FlamingPy. [#64](https://github.com/XanaduAI/flamingpy/pull/64)

### Contributors

This release contains contributions from (in alphabetical order):

[Joost Bus](https://github.com/soosub), [Sebastián Duque Mesa](https://github.com/sduquemesa), [Luis Mantilla](https://github.com/BestQuark), Nariman Saadatmand, [Ilan Tzitrin](https://github.com/ilan-tz), [Trevor Vincent](https://github.com/trevor-vincent)

See full commit details [here](https://github.com/XanaduAI/flamingpy/compare/v0.7.0a4...v0.8.2a5).


## Release 0.7.0a4

### New features since the last release

* The voxel plotting function has been refactored to allow for easy location in space as well as resizing (the latter being important for stabilizers at boundaries that are represented by incomplete cubes). These changes are reflected in two new functions in the `viz` module: _plot_cube_ and _cuboid_data_. [#20](https://github.com/XanaduAI/flamingpy/pull/20)(backward incompatible)
* `Union-Find` --- a fast new decoder based on [arXiv:1709.06218](https://arxiv.org/abs/1709.06218) and [arXiv:1703.01517](https://arxiv.org/abs/1703.01517) --- has been implemented. Now the user may change between the existing minimum-weight perfect matching decoder ("MWPM" setting) and Union-Find ("UF" setting). We have also temporarily disabled the "both" `ec` option in `SurfaceCode` while we investigate a bug, and make some further minor changes related to the Union-Find decoder. [#37](https://github.com/XanaduAI/flamingpy/pull/37)(backward incompatible)

### Bug fixes

* Voxel plots of dual stabilizers used to be drawn incorrectly since only integer locations and cube sizes were allowed. Furthermore, no cube could be placed in a coordinate less than zero. These have been fixed. [#20](https://github.com/XanaduAI/flamingpy/pull/20)
* The occasionally failing `test_hybridize` in `test_graphstates.py` has been fixed. [#25](https://github.com/XanaduAI/flamingpy/pull/25)
* The `lemon` backend used for MWPM decoding was performing worse compared to the other matching backends. The problem was that missing edges in the graph were associated with 0-entries in the adjacency matrix, leading to them always having the minimal weight and making them indistinguishable from edges with an actual weight of 0. The missing edges are now assigned a very large weight. [#28](https://github.com/XanaduAI/flamingpy/pull/28)

### Improvements

* Tests were added to improve the overall test coverage. These included changes to
  `.coveragerc` as well as the refactoring of some examples to allow for proper
  imports from testing modules. Code coverage is now above 95% and
  the overall fail threshold was bumped accordingly. [#14](https://github.com/XanaduAI/flamingpy/pull/14)
* The visuals produced by FlamingPy have been improved and made more consistent. [#20](https://github.com/XanaduAI/flamingpy/pull/20)

  * The figure, marker, line, label and title size, font family, and colormaps were modified.
  When drawing, FlamingPy no longer changes the global matplotlib's `rcParams`,
  but uses `rc_context` together with the plot parameters defined within the `viz` module.
  To customize such parameters, simply use the following and every new plot produced by FlamingPy will use them accordingly.
  ```python
  from flamingpy.utils.viz import plot_params as fp_plot_params
  fp_plot_params["font.size"] = 20
  ```

  * Most functions in the visualization module now return the figure and axes for further processing.
  * The offered method to draw voxels is much clearer and has an easier-to-use API.
  * Graphs of decoding objects (stabilizer and matching graphs) are prettier and easier
  to parse, thanks partially to a new function, `draw_curved_edges`.
  * `draw_adj` and `draw_SCZ` wrapper methods were added to `EGraph` and `CVLayer`, respectively.
* Several changes were made to improve the visualization of MWPM decoding for debugging and understanding purposes. [#23](https://github.com/XanaduAI/flamingpy/pull/23)
  * A function (`draw_decoding`) was added to the `viz` module and new options were added to the `correct` function in the decoder module to be able to simply plot all decoding objects (stabilizer graph, matching graph, matching, syndrome plot) in sync with the actual error correction trial.
  * The appearance and presence of node labels (specifically the virtual nodes of the matching graph) were fixed.
  * The `label_cubes` argument was renamed to the more accurate `label_stabilizers`.
  * The argument `show_matching` was added to the drawing options to be able to turn the matching plot on or off.
  * One can now plot a non-NetworkX matching graph (by automatic conversion to a NetworkX graph).
  * The above changes allowed for a significant simplification of the decoding example.

* The PR template has been changed to inform the user about the 95%+ codecov requirement. [#25](https://github.com/XanaduAI/flamingpy/pull/25)
* `CVLayer` has been modified to allow for instantiation with a code object
  in addition to an `EGraph`. This makes more semantic sense (applying a noise model
  to a code), making it conceptually easier for the user and avoiding noise layers having to reference the internal mechanics of codes. [#25](https://github.com/XanaduAI/flamingpy/pull/25)
* `codecov.yml` was introduced to customize codecov automated tests. For this version, we have added a `threshold: 0.5%` to avoid undesired delta failures due to just removing a few lines, etc. [#25](https://github.com/XanaduAI/flamingpy/pull/25)
* The Walrus has been re-added as a dependency and its functions are used instead of a verbatim
  copy of the code. [#27](https://github.com/XanaduAI/flamingpy/pull/27)
* Since `retworkx` and `lemon` are the fastest backends and `retworkx` follows the same convention
  as `networkx`, the default backend for stabilizer graphs and MWPM has been changed to `retworkx`. [#28](https://github.com/XanaduAI/flamingpy/pull/28)
* Some more tests were added to `test_matching.py` to compare the output of different matching backends. [#28](https://github.com/XanaduAI/flamingpy/pull/28)
* The _display_axes_ option has been changed to show_axes and title to show_title for consistency. The show_title option is now respected. [#37](https://github.com/XanaduAI/flamingpy/pull/37)
* Decoders have become more organized and compartmentalized. [#37](https://github.com/XanaduAI/flamingpy/pull/37)
  * They are located in a directory with their name, with separate modules for decoding objects and algorithms. The latter -- `algos.py` -- contains
  a cumulative decoding function combining all the steps. This function is imported by `decoder.py`, which is now a more general module.
  * The `draw_decoding` function in `viz` can now accommodate plotting generic decoding procedures: a stabilizer graph, a syndrome plot, and the recovery.


### Documentation changes

* The documentation now mentions that `retworkx` is the default backend. [#28](https://github.com/XanaduAI/flamingpy/pull/28)

### Contributors

This release contains contributions from (in alphabetical order):

[Mikhail Andrenkov](https://github.com/Mandrenkov), [Sebastián Duque Mesa](https://github.com/sduquemesa), [Theodor Isacsson](https://github.com/thisac), [Josh Izaac](https://github.com/josh146), [Priya Nadkarni](https://github.com/PriNad), Nariman Saadatmand, [Maxime Tremblay](https://github.com/maxtremblay), [Ilan Tzitrin](https://github.com/ilan-tz)

See full commit details [here](https://github.com/XanaduAI/flamingpy/compare/v0.6.1a3...v0.7.0a4).


## Release 0.6.1a3

### New features since the last release

* Pauli noise: have added a new noise model sampling i.i.d Z error for each qubit. [#8](https://github.com/XanaduAI/flamingpy/pull/8)(backward incompatible)
* Fixed drawing of stabilizer graph for zero syndrome: [#9](https://github.com/XanaduAI/flamingpy/pull/9)(backward incompatible)
  * Previously, the drawing function for a stabilizer graph relied on a non-documented feature. That is, it was assumed that when building the matching graph, all edges of a Networkx-based stabilizer graph were assigned a weight. This, however, was not a fair assumption for many reasons.
  * As a solution, we have added a new method to the `SurfaceCode` class to draw the primal or dual stabilizer graph, which makes sure that each edge has a weight. Now, using that method, the user does not have to rely on unfair assumptions.
  * Furthermore, we added a quick check to not add any edges to the matching graph when the syndrome is trivial. In this case, the cost of decoding should be almost zero.

### Improvements

* A large number of linting corrections were made to improve the overall pylint report. These were mostly minor but essential modifications including restructuring, re-coding, optimization, updating `.pylintrc`, adding `.coveragerc`. The code quality score is improved to `A` for the released version. Check ["Files changed"](https://github.com/XanaduAI/flamingpy/pull/11/files) for details. [#11](https://github.com/XanaduAI/flamingpy/pull/11)

* `TODO` comments have been removed exclusively from files with low code quality grades. The Dev team has created tickets to be actioned for all removed TODO comments on separate (private) FlamingPy boards. [#11](https://github.com/XanaduAI/flamingpy/pull/11)

* `CONTRIBUTING.md`, `code_of_conduct.md`, and `CHANGLELOG.md` were imported and modified from the Strawberry Fields project. Dev team plans to extend these with customized details in future PRs. [#11](https://github.com/XanaduAI/flamingpy/pull/11)

### Documentation changes

* Updated old FT-Stack links in docs header to correct FlamingPy pages. [#7](https://github.com/XanaduAI/flamingpy/pull/7)
* The new Xanadu Sphinx theme has been applied. Currently, most Xanadu OSS projects include their own version of the Xanadu Sphinx theme; however, the Xanadu Sphinx Theme repository is now publicly available and is the preferred way to access the Xanadu CSS theme and Sphinx directives. [#17](https://github.com/XanaduAI/flamingpy/pull/17)
  * Deleted the `doc/xanadu_theme` directory
  * Updated `doc/requirements.txt` and `doc/conf.py` to reference and use the (centralized) Xanadu Sphinx Theme.
  * Replaced the Quantum Error Correction, Install, and FAQ static HTML files with reST ones.

### Contributors

This release contains contributions from (in alphabetical order):

[Mikhail Andrenkov](https://github.com/Mandrenkov), [Sebastián Duque Mesa](https://github.com/sduquemesa), Nariman Saadatmand, [Maxime Tremblay](https://github.com/maxtremblay), [Ilan Tzitrin](https://github.com/ilan-tz)

See full commit details [here](https://github.com/XanaduAI/flamingpy/compare/v0.4.9a1...v0.6.1a3).


## Release 0.4.9a1

### Improvements since the last release

* Relative paths cannot be used in README.md logos and were replaced with GitHub-hosted links. [#5](https://github.com/XanaduAI/flamingpy/pull/5)

* C++ imports are now placed within `try` blocks to avoid interrupting non-compiled installations, such as the ones currently used by readthedocs. [#5](https://github.com/XanaduAI/flamingpy/pull/5)

* Code coverage (quality) score was improved to a more acceptable `B-` level. [#5](https://github.com/XanaduAI/flamingpy/pull/5)

### Bug fixes

* Fixed a bug in [`pull_request_template.md`](https://github.com/XanaduAI/flamingpy/pull/2/commits/e30f2cb65daffece08b193ffc4b8fe7a8d90b90e). The template was not loading properly due to a whitespace problem. [#2](https://github.com/XanaduAI/flamingpy/pull/2)

* Fixed a bug in [`simulations.py`](flamingpy/simulations.py) and related examples. See [here](https://github.com/XanaduAI/flamingpy/commit/771b0e66e5471c3696ac2e779a2df1cc49e5e684) for commit details. [#6](https://github.com/XanaduAI/flamingpy/pull/6)

### Documentation changes

* Making Documentation more usable and consistent with other Xanadu projects [#5](https://github.com/XanaduAI/flamingpy/pull/5):
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
