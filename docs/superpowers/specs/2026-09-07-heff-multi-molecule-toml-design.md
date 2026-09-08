# Design — TOML-configured, multi-molecule effective Hamiltonians

Status: approved architecture direction on 2026-09-07. The foundation plan
implements the existing-ThF+ configuration slice only; later backends and basis
transformations remain follow-on work. See the maintained
[foundation acceptance report](../handoffs/2026-09-07-toml-foundation-acceptance.md).

## 1. Goal

Turn `heff` into a distributable package that can load, build, diagonalize, and
sweep effective Hamiltonians for several molecular structures through one public
API. The first supported families are:

- ThF⁺ `X ³Δ₁`, using the existing case-(c) implementation;
- RaF, beginning with bosonic `X ²Σ⁺(v=0)` and then `A ²Π`;
- YbOH, including the `X(010)` linear-polyatomic manifold;
- SrNH₂, using the C₂ᵥ asymmetric-top model.

Each molecule has one authoritative TOML file. Its executable minimum is the
backend, basis specification, and effective-Hamiltonian parameter values (plus
the small amount of identity needed to select a manifold). Terms may be listed
explicitly or supplied by a named backend default. Units, uncertainties,
sources, notes, observables, and sweep defaults are optional. A user can copy a
bundled file, change values or term selections, and load the result without
editing package code.

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

### 3.2 Coupling cases have native backends and shared transforms

The shared package does not pretend that case (c), `bBJ`, `bBS`, a
linear-polyatomic vibronic basis, and a C₂ᵥ asymmetric top have the same native
kets or symmetry rules. A backend is the natural place to enumerate a basis and
implement its matrix elements; it is not an isolated physical model. Each
backend owns exactly three things:

1. its structured ket dtype;
2. its basis enumerator and structural invariants;
3. its registered matrix-element functions.

Every backend uses the same parameter records, term metadata, term-matrix
assembler, sweep engine, tracking code, observable evaluation, and result
manifest.

Registered basis transformations connect compatible representations. In
particular, coupled Hund's-case bases can transform to one another and to a
shared uncoupled/product basis. This allows an operator to be implemented once
in a convenient representation, then transformed, rather than requiring a
separate formula silo for every coupling case.

### 3.3 Validation is narrow and execution-focused

The loader rejects only conditions that prevent an unambiguous calculation:
invalid TOML, an unknown backend, an invalid basis field or range, a missing
parameter required by an active term, an unsupported explicit unit, or an
incompatible term/backend combination. The resulting error identifies the file
and full TOML key path.

Descriptive metadata is not an acceptance gate. Unknown optional metadata and
unused parameters produce an inspection warning rather than preventing a run.
Sources may be ordinary TOML comments. If machine-readable provenance is useful,
`source`, `uncertainty`, `status`, `convention`, and `note` are optional fields.

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
  basis.py             # basis descriptors and transformation paths
  transitions.py       # cross-manifold operators in a shared representation
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

Each backend and term declares documented canonical conventions. A TOML file
may override them explicitly when importing constants in another registered
convention:

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

Omitting this table selects the backend defaults. Backends declare which named
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
common fields, then delegates backend-specific fields to that schema. Basis
fields are execution-critical, so an unknown basis key is an error rather than
a silently ignored truncation.

The basis descriptor includes the coupling tree, quantum-number ordering,
phase convention, truncation, and spectator degrees of freedom. These fields
identify a transformation unambiguously. A runtime basis override may request a
compatible registered representation, including an uncoupled/product basis.

### 5.5 Effective-Hamiltonian parameters

The common form is deliberately small. Bare numbers use the canonical unit
declared by the consuming term or backend:

```toml
[manifolds.X2Sigma_v0.parameters]
# B and gamma are in the backend's canonical frequency unit (MHz).
# Literature source may be recorded here as a human-readable comment.
B = 5755.56
gamma = 175.38

# Rich records remain available when an explicit conversion or metadata helps.
d_mf = { value = 3.91, unit = "D", convention = "center_of_mass" }
```

An optional rich record may contain `value`, `unit`, `uncertainty`, `status`,
`source`, `isotopologue`, `convention`, and `note`. Only `value` is required in
that form. An explicit unit is checked and converted; no unit means the
documented canonical unit.

Neither bundled nor user files require source, status, or uncertainty metadata.
TOML comments are suitable for citations intended only for readers, with the
important limitation that comments are not retained in the resolved runtime
manifest. Machine-readable metadata is preserved when supplied but never acts
as a schema or scientific-acceptance gate.

### 5.6 Composition and override precedence

The schema uses four fixed layers:

1. file-wide shared parameters;
2. manifold parameters;
3. isotopologue-specific manifold parameters;
4. explicit runtime overrides.

Later layers replace parameter values. A bare numeric override keeps the
parameter's canonical unit. A rich override replaces the rich record; the
loader does not recursively merge optional metadata.

The resolved model records the effective values and their originating layer in
its manifest. The loader does not implement arbitrary inheritance or
user-defined merge rules.

### 5.7 Coefficients and derived parameters

Registered terms own their coefficient construction. Most terms use one
parameter or a product of named knobs, as `heff` does now. Named transformations
handle established parameter conventions such as `N² → R²`.

When a file explicitly supplies a noncanonical unit or convention, the resolved
manifest records both the supplied and transformed canonical values. Otherwise,
the bare number is already canonical. No TOML string evaluates arithmetic.

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
problem = model.problem(N_max=5, basis="case_b_betaJ")
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

### 6.3 Basis changes and transitions

A basis transformation has a declared direction. For
`U_target_from_source`, state coefficients and square operators transform as

```text
c_target = U_target_from_source c_source
O_target = U_target_from_source O_source U_target_from_source†.
```

The transformation registry may compose exact paths through a shared
uncoupled/product basis. The API exposes the operation rather than making users
manually align ket orderings:

```python
case_a = model.problem(basis="case_a")
case_b = case_a.to_basis("case_b_betaJ")

amplitude = heff.transition_matrix_element(
    initial=initial_state,
    final=final_state,
    operator="E1",
    common_basis="uncoupled",
)
```

For a transition between different manifolds, the operator can be rectangular.
"Common basis" means a compatible product-label convention and coupling tree
for the degrees of freedom connected by the operator; the initial and final
electronic/vibronic labels need not describe the same Hilbert space. The engine
transforms both eigenvectors and the transition operator into that declared
representation before evaluating the matrix element. It raises a clear error
when no compatible transformation path exists.

An exact basis change requires the source and target truncations to span the
same physical subspace. A rectangular map caused by truncation is labelled a
projection, not a basis change, and cannot support a basis-invariance claim.

Brown and Carrington explicitly give the case-(b)-to-case-(a) expansion in
Eq. (6.149), PDF p. 262 (printed p. 230), and note that the basis choice is in
principle irrelevant when the bases span the same space. That relation is the
first registered cross-case transformation and fixes its phase convention.

### 6.4 Inspectability

Before a costly build, users can inspect:

```python
model.describe()              # basis, terms, values, and optional metadata
model.validate()              # executable-structure compatibility checks
problem.term_table()          # term ID, coefficient inputs, rules, citation
problem.parameter_table()     # effective values and any supplied metadata
problem.basis_paths()         # exact transforms and labelled projections
```

This makes a TOML edit reviewable before it produces a spectrum.

## 7. Registries and backend contracts

### 7.1 Backend registry

Each backend registration provides:

- a stable public ID;
- its basis schema and enumerator;
- ket invariants;
- supported term IDs and observable IDs;
- canonical units and conventions;
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

### 7.3 Basis-transform registry

Each transform registration provides source and target basis descriptors, the
phase and ordering convention, applicability conditions, and a matrix builder.
Exact transforms must satisfy the appropriate identity and composition checks.
The registry distinguishes square unitary transformations from rectangular
projections introduced by truncation.

Terms need not be independently implemented in every coupling case when an
exact registered transformation is available. A second native implementation
is useful as an independent physics check, but is not part of the minimum
backend contract.

### 7.4 Blocking

The backend declares the quantum numbers conserved by each term. The model may
use a blocked basis only when every active term preserves that block label.
Activating a transverse field or another mixing term either selects the full
basis or raises a structural error before assembly.

## 8. Data flow and artifacts

```text
TOML file
  → TOML parse and minimal executable validation
  → backend, terms, basis, and parameter resolution
  → structured basis and block map
  → parameter-free term matrices
  → scalar or broadcast coefficients
  → Hamiltonian batch and eigensolver
  → optional basis transforms and transition operators
  → observables, labels, and configuration-stamped result
```

Each result manifest adds:

- model ID, manifold ID, and isotopologue ID;
- model schema version and package version;
- source TOML path or bundled resource ID;
- a hash of the normalized effective configuration;
- backend and term IDs;
- effective parameter values and any supplied rich records;
- convention and basis-transformation IDs;
- basis truncation, blocking, ordering, gauge, and assignment policy.

The path and hash identify the input configuration used for a run. They do not
require or imply a literature source for every parameter, and never act as a
scientific acceptance gate.

## 9. Errors and warnings

The loader raises on structural ambiguity:

- unsupported schema version;
- duplicate model or manifold IDs;
- unknown backend or explicitly selected term;
- unknown execution-critical basis key;
- unsupported explicit unit or convention;
- missing required parameter;
- incompatible backend and term;
- invalid spin-coupling chain or basis range;
- illegal blocking for an active term;
- incompatible or absent requested basis-transformation path.

The loader warns but continues for non-blocking information:

- unknown optional metadata or a parameter unused by the selected terms;
- a supplied placeholder, stale, or estimate status;
- a requested observable that is mathematically available but lacks a
  species-specific literature benchmark.

Warnings appear in `describe()` and the result manifest rather than only on
stderr.

## 10. Validation strategy

### 10.1 Schema and loader gates

Tests cover the minimal bare-number form, optional rich parameter records, and
one reachable failure for each execution-blocking error. Every bundled TOML
file must parse and resolve all active-term dependencies. Missing sources,
statuses, uncertainties, notes, or convention tables are explicitly tested as
valid configurations.

### 10.2 Hamiltonian-agnostic gates

The existing resum, derivative, Hermiticity, selection-rule, Wigner, tracking,
blocking, and import-cost tests remain. New backends must pass these gates before
a species-specific spectrum test can support a claim.

### 10.3 Backend physics gates

Each backend requires symmetry identities, closed-form limits, or independent
cross-basis checks appropriate to that coupling case. The package does not infer
that a passing ThF⁺ case-(c) test validates RaF `bBJ` or SrNH₂ asymmetric-top
matrix elements.

For compatible case-(a), case-(b), and uncoupled bases spanning the same
truncated physical subspace, permanent tests check:

- transformation unitarity and round-trip ket ordering;
- `H_target = U H_source U†` term by term and after assembly;
- equal spectra up to degeneracy-aware ordering;
- basis-invariant transition strengths after both states and the operator are
  expressed in the same basis.

These tests must fail or downgrade to projection tests when either truncation is
not closed under the transformation.

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

1. **Model foundation.** Add the lean TOML loader, registries, model object,
   resource packaging, and focused validation tests. Move the existing ThF⁺
   composition and parameters into `thf_plus.toml` without changing its physics.
2. **RaF and basis interoperability.** Port the `bBJ` basis and the minimal
   `X ²Σ⁺(v=0)` terms, then register the case-(a)/case-(b)/uncoupled transform
   needed for basis-invariance and transition tests. Validate against
   Molecule-Structure and backend identities.
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
- treating optional source or status metadata as a scientific gate;
- public publishing or a license choice without a separate user decision.

The implementation migrates one verified vertical slice at a time. It preserves
the source repositories unchanged. Optional metadata can record unresolved
physics, but scientific claims remain limited by their actual validation rather
than by TOML completeness.

## 13. Accepted decisions

- Use authoritative TOML files, not Python molecule modules.
- Store both basis specifications and effective-Hamiltonian parameters in TOML.
- Keep TOML data-only; keep formulas and symmetry logic in Python.
- Make bare canonical parameter values the default; keep units and provenance
  metadata optional.
- Treat source comments and metadata as documentation, never strict load gates.
- Retain native coupling-case backends behind one model API, connected by
  registered exact basis transformations and a shared uncoupled representation.
- Evaluate cross-basis transitions only after states and operators are expressed
  in a compatible common basis.
- Begin non-ThF⁺ validation with bosonic RaF `X ²Σ⁺(v=0)`.
- Preserve the current low-level `heff` API while adding a high-level model API.

## 14. Spec self-review

- Placeholder scan: no unresolved `TBD` or `TODO` fields remain.
- Consistency: the schema, public API, data flow, errors, tests, and delivery
  sequence use the same lean, data-only TOML boundary.
- Scope: the first plan can implement the model foundation and RaF vertical
  slice without implementing YbOH or SrNH₂ simultaneously.
- Ambiguity: parameter layering, backend ownership, basis-transform semantics,
  common-basis transition evaluation, term selection, and external publishing
  boundaries are explicit.
