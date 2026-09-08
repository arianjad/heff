# Package architecture

`heff` uses a term-matrix design. Basis enumeration and angular-momentum
algebra depend on the model; parameter scans do not rebuild those matrices.

```text
TOML or bundled model
        |
        v
model_io.py -> model.py -> backend adapter
                              |       |
                              v       v
                         basis kets  term registry
                              \       /
                               assemble.py
                                   |
                         TermMatrices + ParamSet
                                   |
                    hamiltonian() / engine.sweep()
                                   |
                    eigenvalues, eigenvectors, manifest
```

## Public composition layer

- `heff.model_io` parses schema version 1, freezes definitions, applies
  isotopologue overrides, checks convention records, and converts parameters
  to canonical units.
- `heff.model` exposes `list_models()`, `load_model()`, `MoleculeModel`, and
  `Problem`. A `MoleculeModel` holds one selected manifold and isotopologue. A
  `Problem` holds the backend, basis specification, parameters, and active term
  names.
- `heff.backend_registry` maps a TOML backend ID to its basis builder, term
  registry, units, runtime knobs, and conventions. Bundled adapters load only
  when a model is selected.

The high-level import stays light: importing `heff` does not initialize model
backends, SymPy, SciPy, or plotting libraries. Public model symbols are loaded
on first access.

## Physics backends

- `heff.backends.case_c` adapts a zero- or one-spin Hund's-case-(c) model to
  `StateSpec`; its matrix elements live in `heff.elements_c`.
- `heff.backends.case_c2` adapts the ordered Th-plus-F two-spin chain; its
  matrix elements live in `heff.elements_c2`.
- `heff.backends.amide` builds the restricted C2v metal-amide model. Basis and
  exchange bookkeeping live in `heff.amide_basis`; the inherited core kernels
  live in `heff.elements_amide`; the backend recouples them through the outer
  metal spin and defines the metal contact and Zeeman terms.

Each backend names its canonical units, accepted basis shape, runtime knobs,
and selected term registry. The TOML loader does not infer a basis or operator
from a molecule's name.

## Shared numerical core

- `heff.spec` defines the case-c ket layouts, basis specifications, basis
  enumeration, and signed-`mF` blocking.
- `heff.terms` defines a `Term`: its parameter symbols, selection rules,
  reality and Hermiticity metadata, citation, and unit-coefficient matrix
  element.
- `heff.assemble` evaluates each selected term once to make immutable
  `TermMatrices`. It combines those matrices with canonical parameter values
  and runtime knobs to produce one Hamiltonian or a batch.
- `heff.engine` diagonalizes batches and returns a `SweepResult` containing
  eigenvalues, eigenvectors, kets, active terms, ordering and gauge choices,
  reference index, and the conventions manifest.
- `heff.track`, `heff.observe`, and `heff.spectra` provide state assignment,
  observables, and one-photon spectroscopy on the assembled representation.
  `heff.twophoton` supplies the separate rank-K closure operator.

The matrix cache boundary is the practical center of the package:

```python
import numpy as np
import heff

model = heff.load_model("thf_plus")
problem = model.problem(J_max=2)

kets = problem.kets                    # enumerate once
matrices = problem.term_matrices        # angular algebra once
H = problem.hamiltonian(E_z=0, B_z=0)  # new coefficients only
E = np.linspace(0.0, 100.0, 51)
B = np.zeros_like(E)
scan = problem.sweep(E_z=E, B_z=B)     # batched coefficients and eigh
```

## Data ownership and extension points

Model files own molecular constants, term selection, coupling records, basis
limits, and convention choices. Backends own the translation from those
records to an existing basis and operator family. Matrix-element modules own
the physics formulas. The engine owns coefficient assembly, diagonalization,
ordering, and result metadata.

To add an isotopologue that uses an existing backend, add or override TOML data
and verify its parameter provenance and basis convergence. To add a new model
class, implement its basis and matrix elements, register a backend adapter, and
compare signed matrix elements against an independent limit or construction.
Do not encode missing physics as a silent zero: omit or refuse an unsupported
term, or label a deliberate switch-off or sensitivity value in the model data.

See [Models and TOML files](models.md) for the executable schema and current
backend limits. The canonical ThF+ equations and convention decisions remain
in [the effective-Hamiltonian note](thf-plus-x3delta1-effective-hamiltonian.md).
