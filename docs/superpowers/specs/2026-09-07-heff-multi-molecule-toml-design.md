# Design — TOML-configured, multi-molecule effective Hamiltonians

Status: approved architecture direction on 2026-09-07; awaiting review of this
written specification before implementation planning.

## 1. Goal

Turn `heff` into a distributable package that can load, build, diagonalize, and
sweep effective Hamiltonians for several molecular structures through one public
API. The first supported families are:

- ThF⁺ `X ³Δ₁`, using the existing case-(c) implementation;
- RaF, beginning with bosonic `X ²Σ⁺(v=0)` and then `A ²Π`;
- YbOH, including the `X(010)` linear-polyatomic manifold;
- SrNH₂, using the C₂ᵥ asymmetric-top model.

Each molecule has one authoritative TOML file. The file contains its basis
specifications, enabled effective-Hamiltonian terms, parameter values, units,
uncertainties, provenance, conventions, observables, and optional sweep
defaults. A user can copy a bundled file, change values or term selections, and
load the result without editing package code.

The package retains the existing term-matrix architecture:

`H = Σₖ cₖ Mₖ`.

The TOML file selects `Mₖ` and supplies `cₖ`. Registered Python backends define
the basis, symmetry constraints, selection rules, and matrix elements.

## 2. Evidence from the current repositories

The design follows the live code rather than imposing one coupling case on all
four systems.

- `heff` already separates `StateSpec`, `ParamSet`, registered terms, term
  matrices, assembly, sweeps, tracking, and observables. Its current public path
  remains low-level, and its enumerator accepts only case (c), one electronic
  state, and at most two coupled nuclear spins.
- Molecule-Structure stores RaF and YbOH constants in nested dictionaries in
  `Source Code/molecule_parameters.py`. It dispatches `bBJ`, `bBS`, and `aBJ`
  bases through repeated dictionaries in `molecule_library_class.py`.
- C2V-Molecules stores SrNH₂-family constants in production-script dictionaries
  such as `scripts/scan2k.py`. Its basis and fixed Hamiltonian sum live together
  in `atm_core/physics.py`.
- The three target additions span different physical backends: diatomic Hund's
  cases, a linear-polyatomic bending manifold, and an exchange-symmetrized
  asymmetric top. They share the numerical spine but not one matrix-element
  implementation.

The code audit found no external molecule configuration format in either source
repository. This design replaces repeated dictionaries and dispatch tables with
one schema and registered backend identifiers.

## 3. Design principles

### 3.1 TOML contains data, not executable formulas

TOML may specify:

- quantum-number ranges, spins, isotopologues, frames, and block choices;
- a registered backend ID and registered term IDs;
- parameter values and metadata;
- named convention and transformation IDs;
- observables and example sweep defaults.

TOML may not contain:

- matrix-element expressions or Wigner-symbol calls;
- general arithmetic strings such as `"p + 2*q"`;
- conditional expressions such as `"if isotope == 173"`;
- coupling-tree equations or arbitrary state-filter expressions;
- import paths or executable callbacks.

This boundary keeps the configuration portable while preserving one audited
implementation of each physical formula.

### 3.2 Coupling cases remain separate backends

The shared package does not pretend that case (c), `bBJ`, `bBS`, a
linear-polyatomic vibronic basis, and a C₂ᵥ asymmetric top have the same kets or
symmetry rules. Each backend owns exactly three things:

1. its structured ket dtype;
2. its basis enumerator and structural invariants;
3. its registered matrix-element functions.

Every backend uses the same parameter records, term metadata, term-matrix
assembler, sweep engine, tracking code, observable evaluation, and result
manifest.

### 3.3 Configuration errors fail before matrix construction

The loader rejects unknown keys, units, backends, terms, conventions, and
parameters. It also rejects a term that its selected backend does not support.
The resulting error identifies the file and full TOML key path.

Physics validity remains distinct from schema validity. A file can be
well-formed while containing an estimated or placeholder constant; those
statuses remain visible in the model description and result manifest.

## 4. Package layout

```text
heff/
  models/
    thf_plus.toml
    raf.toml
    yboh.toml
    srnh2.toml
  backends/
    case_c.py
    case_a.py
    case_b_betaj.py
    case_b_betas.py
    linear_polyatomic.py
    asym_top_c2v.py
  model.py             # parsed model records and high-level user object
  model_io.py          # TOML loading, overlays, validation, normalization
  backend_registry.py  # stable backend IDs
  terms.py             # shared term metadata and registry machinery
  assemble.py
  engine.py
  observe.py
  ...
```

The exact backend count may shrink during implementation when two backends can
share an enumerator without hiding different conventions. Backend IDs remain
stable public configuration values after the first release.

Bundled TOML files ship as package data and load through `importlib.resources`,
so installed wheels do not depend on the source checkout layout.

## 5. TOML schema

### 5.1 Top-level identity

Each file contains:

```toml
schema_version = 1
model_id = "raf"
display_name = "RaF"
default_isotopologue = "226Ra19F"
default_manifold = "X2Sigma_v0"
```

`model_id` is unique across the installed catalogue. `schema_version` selects
the parser contract; the loader refuses unsupported versions rather than
guessing how to reinterpret them.

### 5.2 Conventions

The file declares decision-relevant conventions explicitly:

```toml
[conventions]
formalism = "R2"
phase = "brown_carrington_hirota"
body_axis = "ligand_to_metal"
parity = "signed_primitive"
zeeman_energy = "E_minus_g_muB_B_mF"
internal_energy_unit = "MHz"
electric_field_unit = "V/cm"
magnetic_field_unit = "G"
```

Backends declare which convention fields they require and which named
transformations they support. For example, a parameter set labelled `N2` can
invoke the registered `n2_to_r2` transformation; a raw formula in TOML cannot.

### 5.3 Isotopologues and spins

```toml
[isotopologues."226Ra19F"]

[[isotopologues."226Ra19F".spins]]
label = "19F"
I = 0.5
couple_to = "J"

[isotopologues."225Ra19F"]

[[isotopologues."225Ra19F".spins]]
label = "225Ra"
I = 0.5
couple_to = "J"

[[isotopologues."225Ra19F".spins]]
label = "19F"
I = 0.5
couple_to = "F1"
```

The backend interprets `couple_to` only from its allowed coupling labels. The
loader rejects an unsupported chain before enumerating kets.

Equivalent-nucleus exchange symmetry is a named backend option, not a TOML
expression. The SrNH₂ file can select a registered rule such as
`exchange_symmetry = "two_equivalent_protons"` and provide the nuclear spins.

### 5.4 Manifolds and basis

One molecule file may contain several electronic or vibronic manifolds. A
molecule does not imply one basis.

```toml
[manifolds.X2Sigma_v0]
label = "X 2Sigma+ (v=0)"
backend = "case_b_betaJ"
terms = [
  "rotation",
  "centrifugal",
  "spin_rotation",
  "hyperfine_F_contact",
  "hyperfine_F_dipolar",
  "stark_z",
  "zeeman_electron_z",
  "zeeman_nuclear_z",
]
observables = ["parity", "dipole_z", "g_factor"]

[manifolds.X2Sigma_v0.basis]
N_min = 0
N_max = 4
M = "blocks"
frame = "lab"
parity_basis = "signed_primitive"
```

Each backend publishes a typed basis schema. The generic TOML loader validates
common fields, then delegates backend-specific fields to that schema. Backend
schemas reject unused fields so a misspelled or irrelevant truncation cannot be
silently ignored.

### 5.5 Effective-Hamiltonian parameters

The TOML file contains the full effective-Hamiltonian parameter record:

```toml
[manifolds.X2Sigma_v0.parameters.B]
value = 5755.56
unit = "MHz"
status = "measured"
source = "source citation"

[manifolds.X2Sigma_v0.parameters.gamma]
value = 175.38
unit = "MHz"
status = "measured"
source = "source citation"

[manifolds.X2Sigma_v0.parameters.d_mf]
value = 3.91
unit = "D"
status = "measured"
source = "source citation"
convention = "center_of_mass"
```

Allowed metadata fields match and extend the existing `Param` record:
`value`, `unit`, `uncertainty`, `status`, `source`, `isotopologue`,
`convention`, and `note`.

Every bundled nonzero physical parameter must specify `unit`, `status`, and
`source`. Zero-valued switches may use `held-fixed` with a note explaining the
choice. Unknown or unsupported units fail during loading.

### 5.6 Composition and override precedence

The schema uses four fixed layers:

1. file-wide shared parameters;
2. manifold parameters;
3. isotopologue-specific manifold parameters;
4. explicit runtime overrides.

Later layers replace whole parameter records; they do not recursively merge
individual metadata fields. A runtime bare float follows the existing `ParamSet`
policy: it receives the previous unit, loses source and uncertainty, and becomes
`status = "unspecified"` with a warning.

The resolved model records all four sources and writes the normalized effective
configuration into its manifest. The loader does not implement arbitrary
inheritance or user-defined merge rules.

### 5.7 Coefficients and derived parameters

Registered terms own their coefficient construction. Most terms use one
parameter or a product of named knobs, as `heff` does now. Named transformations
handle established parameter conventions such as `N² → R²`.

The TOML file stores the sourced parameters in their published convention. The
resolved manifest records both the raw values and transformed canonical values.
No TOML string evaluates arithmetic.

## 6. Python interfaces

### 6.1 Loading and discovery

```python
import heff

heff.list_models()
model = heff.load_model("RaF", manifold="X2Sigma_v0", isotope="226Ra19F")
model.describe()
```

`load_model()` accepts a bundled model ID or a path to a user TOML file.
Identifiers are case-insensitive for lookup but normalize to the canonical
`model_id` in artifacts.

### 6.2 Building and running

```python
problem = model.problem(N_max=5)
H = problem.hamiltonian(E_z=20.0, B_z=0.01)
result = problem.sweep(E_z=E, B_z=B)
```

The high-level object performs the existing low-level sequence:

1. enumerate and validate kets;
2. choose or validate blocking;
3. build parameter-free term matrices;
4. resolve canonical coefficients;
5. assemble and diagonalize the Hamiltonian;
6. evaluate requested observables and stamp the result manifest.

The current low-level functions remain available. The model API composes them;
it does not replace the numerical spine.

### 6.3 Inspectability

Before a costly build, users can inspect:

```python
model.describe()              # basis, terms, parameters, statuses, conventions
model.validate()              # schema and compatibility checks
problem.term_table()          # term ID, coefficient inputs, rules, citation
problem.parameter_table()     # raw and canonical values with provenance
```

This makes a TOML edit reviewable before it produces a spectrum.

## 7. Registries and backend contracts

### 7.1 Backend registry

Each backend registration provides:

- a stable public ID;
- its basis schema and enumerator;
- ket invariants;
- supported term IDs and observable IDs;
- required convention fields;
- conserved quantum numbers and supported blocking modes.

The registry stores callables in Python. TOML can reference only registered IDs.

### 7.2 Term registry

Each term retains the present metadata and adds explicit backend compatibility:

- stable term ID;
- parameter or knob dependencies;
- declared selection rules;
- Hermitian and real/complex properties;
- source citation;
- matrix-element callable.

The model selects a subset of compatible terms. An omitted term is absent from
the model; a selected term with a zero coefficient is present but inactive.
Artifacts preserve that distinction.

### 7.3 Blocking

The backend declares the quantum numbers conserved by each term. The model may
use a blocked basis only when every active term preserves that block label.
Activating a transverse field or another mixing term either selects the full
basis or raises a structural error before assembly.

## 8. Data flow and artifacts

```text
TOML file
  → strict parse and schema validation
  → backend, terms, conventions, and parameter resolution
  → structured basis and block map
  → parameter-free term matrices
  → scalar or broadcast coefficients
  → Hamiltonian batch and eigensolver
  → observables, labels, and provenance-stamped result
```

Each result manifest adds:

- model ID, manifold ID, and isotopologue ID;
- model schema version and package version;
- source TOML path or bundled resource ID;
- a hash of the normalized effective configuration;
- backend and term IDs;
- raw and canonical parameter records;
- convention and transformation IDs;
- basis truncation, blocking, ordering, gauge, and assignment policy.

The hash records provenance. It never acts as a scientific acceptance gate.

## 9. Errors and warnings

The loader raises on structural ambiguity:

- unsupported schema version;
- duplicate model or manifold IDs;
- unknown TOML key, backend, term, observable, unit, or convention;
- missing required parameter;
- incompatible backend and term;
- invalid spin-coupling chain or basis range;
- illegal blocking for an active term;
- convention mismatch between a parameter and its model.

The loader warns but continues for declared epistemic limitations:

- placeholder, stale, estimate, or unspecified parameter status;
- runtime replacement of a sourced parameter by a bare float;
- a requested observable that is mathematically available but lacks a
  species-specific literature benchmark.

Warnings appear in `describe()` and the result manifest rather than only on
stderr.

## 10. Validation strategy

### 10.1 Schema and loader gates

Tests cover every allowed field and one reachable failure for each structural
error. Every bundled TOML file must parse, validate, and resolve without unknown
or ignored keys.

### 10.2 Hamiltonian-agnostic gates

The existing resum, derivative, Hermiticity, selection-rule, Wigner, tracking,
blocking, and import-cost tests remain. New backends must pass these gates before
a species-specific spectrum test can support a claim.

### 10.3 Backend physics gates

Each backend requires symmetry identities, closed-form limits, or independent
cross-basis checks appropriate to that coupling case. The package does not infer
that a passing ThF⁺ case-(c) test validates RaF `bBJ` or SrNH₂ asymmetric-top
matrix elements.

### 10.4 Migration comparisons

During each port, compare the new implementation against the current producing
repository for:

- ordered basis labels;
- each individual term matrix;
- the assembled zero-field Hamiltonian;
- a small longitudinal-field spectrum;
- selected observable matrices.

These differential comparisons validate the transplant. Permanent scientific
tests should prefer identities and primary-source comparisons over pinned
spectral snapshots.

The first non-ThF⁺ acceptance slice is bosonic `²²⁶Ra¹⁹F X ²Σ⁺(v=0)` in a
`bBJ` basis with longitudinal `E_z` and `B_z`. RaF `A ²Π` waits until the stale
`A_SO`/origin inputs in Molecule-Structure have an explicit parameter ruling.

### 10.5 Distribution gates

The release check builds both a wheel and source distribution, installs the
wheel into an isolated environment, loads all bundled TOML resources outside
the checkout, runs a documented RaF or ThF⁺ example, and verifies package
metadata. Publishing to TestPyPI or PyPI remains a separate external action.

## 11. Delivery sequence

1. **Model foundation.** Add the TOML schema, loader, registries, model object,
   resource packaging, and negative validation tests. Move the existing ThF⁺
   composition and parameters into `thf_plus.toml` without changing its physics.
2. **RaF vertical slice.** Port the `bBJ` basis and the minimal `X ²Σ⁺(v=0)`
   terms. Validate against Molecule-Structure and backend identities.
3. **YbOH.** Add the required linear-polyatomic and `bBS` coverage one manifold
   at a time, starting with `X(010)` only after the RaF boundary is stable.
4. **SrNH₂.** Add the C₂ᵥ asymmetric-top backend, equivalent-proton exchange
   rule, and the term subset used by the current production model. Reassess
   dense term-matrix memory before accepting its target truncation.
5. **Shareable release.** Complete README examples, metadata, license decision,
   wheel/sdist checks, and clean-environment installation. Stop before any push,
   package-index upload, or publication without explicit authorization.

## 12. Scope boundaries

This program does not initially add:

- arbitrary formula evaluation in configuration files;
- a general spherical-tensor compiler;
- fitting, dynamics, optical Bloch equations, or GPU support;
- every state in the two source repositories;
- automatic conversion of uncertain or stale constants into trusted values;
- public publishing or a license choice without a separate user decision.

The implementation migrates one verified vertical slice at a time. It preserves
the source repositories unchanged and records unresolved physics as parameter
status or a stopped validation claim.

## 13. Accepted decisions

- Use authoritative TOML files, not Python molecule modules.
- Store both basis specifications and effective-Hamiltonian parameters in TOML.
- Keep TOML data-only; keep formulas and symmetry logic in Python.
- Retain separate coupling-case backends behind one model API.
- Begin non-ThF⁺ validation with bosonic RaF `X ²Σ⁺(v=0)`.
- Preserve the current low-level `heff` API while adding a high-level model API.

## 14. Spec self-review

- Placeholder scan: no unresolved `TBD` or `TODO` fields remain.
- Consistency: the schema, public API, data flow, errors, tests, and delivery
  sequence use the same data-only TOML boundary.
- Scope: the first plan can implement the model foundation and RaF vertical
  slice without implementing YbOH or SrNH₂ simultaneously.
- Ambiguity: parameter layering, backend ownership, term selection, and
  external publishing boundaries are explicit.
