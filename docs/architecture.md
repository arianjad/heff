# Package architecture

For a fixed basis and set of interactions, a field sweep changes the
coefficients of the Hamiltonian. We calculate each interaction's matrix once,
then reuse it at every field. The package follows that separation:

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

## Loading a model

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

Importing `heff` loads neither model backends nor SymPy, SciPy, or plotting
libraries. Model symbols load when first accessed.

## Physics backends

- `heff.backends.case_c` adapts a zero- or one-spin Hund's-case-(c) model to
  `StateSpec`; its matrix elements live in `heff.elements_c`.
- `heff.backends.case_c2` adapts the ordered Th-plus-F two-spin chain; its
  matrix elements live in `heff.elements_c2`.
- `heff.backends.amide` builds the restricted C2v metal-amide model. Basis and
  exchange bookkeeping live in `heff.amide_basis`; the inherited core kernels
  live in `heff.elements_amide`; the backend recouples them through the outer
  metal spin and defines the metal contact and Zeeman terms.

Each backend specifies its units, allowed basis, field controls, and
interactions. The TOML file must select a backend that implements the intended
physics.

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

This example shows which work is reused during a scan:

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

## Where to make a change

Edit model files to change constants, active terms, spin coupling, basis
limits, or conventions. Backends translate those records into a basis and
operators; the matrix-element modules contain the formulas. The engine
assembles and diagonalizes the Hamiltonian, orders the states, and records the
calculation's settings.

To add an isotopologue that uses an existing backend, add or override TOML data
and verify its parameter provenance and basis convergence. To add a new model
class, implement its basis and matrix elements, register a backend adapter, and
compare signed matrix elements against an independent limit or construction.
If an interaction is unsupported, omit it or reject the input. If you set a
coefficient to zero for a comparison, record that choice in the parameter
metadata so it cannot be mistaken for a measured value.

See [Models and TOML files](models.md) for the executable schema and current
backend limits. The canonical ThF+ equations and convention decisions remain
in [the effective-Hamiltonian note](thf-plus-x3delta1-effective-hamiltonian.md).
