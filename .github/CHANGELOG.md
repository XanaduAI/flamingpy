# Release 0.5.1a1 (development release)

### New features since last release

### Housekeeping chores

### Bug fixes

### Documentation changes

### Contributors

This release contains contributions from (in alphabetical order) ...

---

# Release 0.4.9a1 (current release)

<h3>New features since last release</h3>

* `Device.layout` and `Device.gate_parameters` may now return `None`. This can happen
  when a remote simulator device is used.
  [(#661)](https://github.com/XanaduAI/strawberryfields/pull/661)

* A new interferometer decomposition method is implemented following the proposal of the paper
  [_Simple factorization of unitary transformations_](https://doi.org/10.1103/PhysRevA.97.022328).
  [(#665)](https://github.com/XanaduAI/strawberryfields/pull/665)

  ```python
  import numpy as np

  import strawberryfields as sf
  from strawberryfields import ops

  U = np.array([[-0.39302099+0.28732291j,  0.83734522+0.24866248j],
                [ 0.00769051+0.87345344j, -0.3847068 +0.29836325j]])

  prog = sf.Program(2)

  with prog.context as q:
      ops.Interferometer(U, mesh="sun_compact") | q
  ```

* A `Device.certificate` method is added which returns the hardware device certificate.
  [(#679)](https://github.com/XanaduAI/strawberryfields/pull/679)

  ```pycon
  >>> import strawberryfields as sf
  >>> eng = sf.RemoteEngine("X8")
  >>> print(eng.device.certificate)
  {'target': 'X8_01' ... }
  ```

* Setting `shots=None` in the engine or program run options will not execute any measurements
  applied on the circuit.
  [(#682)](https://github.com/XanaduAI/strawberryfields/pull/682)

  ```python
  import strawberryfields as sf
  from strawberryfields import ops

  prog = sf.Program(1)
  eng = sf.Engine("gaussian")

  with prog.context as q:
      ops.Sgate(0.5) | q[0]
      ops.MeasureFock() | q

  results = eng.run(prog, shots=None)

  # samples will output an empty list []
  print(results.samples)

  # the resulting Gaussian state is still accessible
  # via its vector of means and covariance matrix
  print(results.state.means())
  print(results.state.cov())
  ```

* There's a `program_equivalence` function in `strawberryfields/program_utils.py` which checks
  Strawberry Fields programs for equivalence.
  [(#686)](https://github.com/XanaduAI/strawberryfields/pull/686)

* An equality operator is implemented for `strawberryfields.Program`, checking that the exact same
  gates and respective parameters, are applied in order.
  [(#686)](https://github.com/XanaduAI/strawberryfields/pull/686)

  ```python
  import strawberryfields as sf
  from strawberryfields import ops

  prog_1 = sf.Program(1)
  prog_2 = sf.Program(1)

  with prog.context as q:
      ops.Sgate(0.42) | q[0]
      ops.MeasureFock() | q

  with prog.context as q:
      ops.Sgate(0.42) | q[0]
      ops.MeasureFock() | q

  assert prog_1 == prog_2
  ```

* A `Program.equivalence` convenience method is added which calls the `program_equivalence`
  utility function.
  [(#686)](https://github.com/XanaduAI/strawberryfields/pull/686)

  ```python
  prog_1 = sf.Program(1)
  prog_2 = sf.Program(1)

  with prog.context as q:
      ops.Sgate(1.1) | q[0]
      ops.MeasureFock() | q

  with prog.context as q:
      ops.Sgate(0.42) | q[0]
      ops.MeasureFock() | q

  assert prog_1.equivalence(prog_2, compare_params=False)
  ```

* A `Device.validate_target` static method is added which checks that the target in the layout is the same as
  the target field in the specification. This check is also performed at `Device` initialization.
  [(#687)](https://github.com/XanaduAI/strawberryfields/pull/687)

* Tests are run in random order and the seed for NumPy's and Python's random number generators are
  set by `pytest-randomly`.
  [(#692)](https://github.com/XanaduAI/strawberryfields/pull/692)

* Adds support for Python 3.10.
  [#695](https://github.com/XanaduAI/strawberryfields/pull/695)

<h3>Breaking Changes</h3>

* `DeviceSpec` is renamed to `Device`, which now also contains more than only the device specification.
  [(#679)](https://github.com/XanaduAI/strawberryfields/pull/679)

  ```pycon
  >>> import strawberryfields as sf
  >>> eng = sf.RemoteEngine("X8")
  >>> isinstance(eng.device, sf.Device)
  True
  >>> print(eng.device.target)
  X8_01
  ```

<h3>Bug fixes</h3>

* It's now possible to show graphs using the plot apps layer when not run in notebooks.
  [(#669)](https://github.com/XanaduAI/strawberryfields/pull/669)

* `program.compile` now raises an error if the device specification contains gate parameters but no
  circuit layout. Without a layout, the gate parameters cannot be validated against the device
  specification.
  [(#661)](https://github.com/XanaduAI/strawberryfields/pull/661)

* The teleportation tutorial `examples/teleportation.py` now uses the correct value (now `phi = 0`
  instead of `phi = np.pi / 2`) for the phase shift of the beamsplitters.
  [(#674)](https://github.com/XanaduAI/strawberryfields/pull/674)

* `Program.compile()` returns a deep copy of the program attributes, except for the circuit and
  the register references.
  [(#688)](https://github.com/XanaduAI/strawberryfields/pull/688)

<h3>Contributors</h3>

This release contains contributions from (in alphabetical order) Sebastian Duque, Theodor Isacsson, Jon Schlipf, and Hossein Seifoory.
