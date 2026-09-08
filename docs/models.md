# Models and TOML files

Before calculating a spectrum, we need a basis, a set of interactions, and
their coefficients. A TOML file records these choices, the isotopologue, and
the conventions. Its backend supplies the basis and matrix elements.
`load_model()` reads the file; `model.problem()` prepares a calculation. You
can change the basis cutoff for that calculation without editing the model.

## Supported model classes

| Backend | Basis | Nuclear-spin chain | Present scope |
|---|---|---|---|
| `case_c` | coupled `(J, Omega, F, mF)` | zero or one spin, coupled to `J` | One axial electronic state with Omega = +/-1 terms used by 232Th19F+ X 3Delta1. Static collinear `E_z` and `B_z`; the stored `frame="rotating"` labels the JILA basis but does not add a rotating-frame Hamiltonian. |
| `case_c2` | coupled `(J, Omega, F1, F, mF)` | exactly two spins, inner spin coupled to `J`, outer spin coupled to `F1` | The native odd-Th ThF+ model. The first spin is thorium and the second is fluorine; this adapter is not a generic two-spin molecule backend. |
| `amide_c2v` | coupled `(N, K, J, F_N, I_T, F_core, F, mF)` | nitrogen, an equivalent proton pair, and optional outer metal spin | Restricted X 2A1 metal-amide model in a Hund's-case-(b) coupled basis. It implements the listed core terms plus isotropic metal contact and metal Zeeman terms. |

The amide backend preserves the exchange filter for two equivalent spin-1/2
protons. For `vibronic_sign = 1`, even K retains the proton singlet `I_T=0`
and odd K retains the triplet `I_T=1`. It omits the difference-spin operator
that can mix those sectors at different K. It also omits metal electric
quadrupole, anisotropic metal hyperfine, shielding, nuclear spin-rotation, and
nuclear spin-spin terms. The example coefficients are synthetic and do not
describe a real CaNH2, SrNH2, or RaNH2 isotopologue.

## Units and parameter records

Every assembled Hamiltonian is in MHz. The loader accepts these input units:

| Quantity | Accepted unit strings | Canonical value used by matrices |
|---|---|---|
| Energy or frequency | `Hz`, `kHz`, `MHz`, `GHz`, `cm-1` | MHz |
| Body-fixed dipole | `D`, `MHz/(V/cm)` | MHz/(V/cm) |
| Effective electric field | `GV/cm`, `MHz/(e cm)` | MHz/(e cm) |
| Dimensionless coefficient | empty string `""` | unchanged |
| Closure polarizability | `(MHz/(V/cm))^2/MHz` | unchanged |

A bare number uses the backend's canonical unit:

```toml
B0 = 7274.3325
```

To record the unit and source, use a parameter record:

```toml
d_mf = { value = 3.37, unit = "D", status = "measured", convention = "center_of_mass", source = "Ng 2022 Table I" }
```

The recognized fields are `value`, `unit`, `uncertainty`, `status`, `source`,
`isotopologue`, `convention`, and `note`. Parameter status is descriptive; it
does not disable a term. Unknown metadata is retained with a warning. An
unknown unit fails when the problem resolves parameters.

The molecular dipole of an ion depends on the coordinate origin. We therefore
record that origin explicitly: the `d_mf` parameter's `convention` must match
`conventions.dipole_origin`. The loader rejects an untagged or mismatched
dipole before constructing matrices.

## Complete case-c example

Save this file as `minimal_case_c.toml`:

```toml
schema_version = 1
model_id = "minimal_case_c"
default_manifold = "X"
default_isotopologue = "test"

[[isotopologues.test.spins]]
label = "19F"
I = 0.5
couple_to = "J"

[manifolds.X]
backend = "case_c"
terms = ["rotation", "omega_doubling", "stark_z"]

[manifolds.X.electronic]
label = "X3Delta1"
Omega = 1.0
S = 1.0
Lambda = 2.0
T0 = 0.0

[manifolds.X.basis]
J_min = 1
J_max = 2
M = "blocks"
frame = "rotating"

[manifolds.X.parameters]
B0 = { value = 7274.3325, unit = "MHz" }
omega_ef = { value = 5.29, unit = "MHz" }
d_mf = { value = 3.37, unit = "D", convention = "center_of_mass" }
```

Load and run it from the repository root:

```python
import numpy as np
from heff import load_model

model = load_model("minimal_case_c.toml")
model.validate()                 # structure only; matrices are still lazy
problem = model.problem(J_max=1) # immutable runtime basis override

H0 = problem.hamiltonian(E_z=0.0)
result = problem.sweep(
    E_z=np.array([0.0, 20.0, 40.0]),
)
print(H0.shape)          # (12, 12)
print(result.evals.shape) # (3, 12)
```

Terms may refer to runtime knobs such as `E_z` and `B_z`; those values do not
belong in `[manifolds.<id>.parameters]`. Every other symbol required by a
selected term must be present. Misspelled terms, basis keys, convention keys,
and missing parameters fail with their TOML path.

## Bundled ThF+ example

Use `heff.list_models()` to list the bundled models. The file
[thf_plus.toml](../heff/models/thf_plus.toml) contains one shared X 3Delta1
manifold and three isotopologues:

```python
from heff import load_model

even = load_model("thf_plus")
odd_229 = load_model("thf_plus", isotope="229Th19F")
odd_227 = load_model("thf_plus", isotope="227Th19F")

even_problem = even.problem(J_max=2)       # case_c
odd_problem = odd_229.problem(J_max=1)     # case_c2
odd_229.validate()

print(even.backend.id, len(even_problem.kets))
print(odd_229.backend.id, len(odd_problem.kets))
# case_c 32
# case_c2 72
```

Isotopologue overrides merge parameters and conventions by name, and replace
`backend` or `terms` when supplied. The base manifold remains immutable.
Identifiers for bundled models, manifolds, and isotopologues are matched
case-insensitively.

The `229Th19F` quadrupole constants are unvalidated sensitivity placeholders.
The `227Th19F` bundled default is the +39.821 GHz Schmidt stress test; it is
distinct from the +1.790176 GHz Minkov theory case used in the field-plot
dataset. Read [the current limitations](open-questions.md) before interpreting
either odd-isotope spectrum.

## Complete amide example

For the amide backend, start with the synthetic model in
[amide_synthetic.toml](../examples/models/amide_synthetic.toml):

```python
import numpy as np
from heff import load_model

model = load_model("examples/models/amide_synthetic.toml")
model.validate()
problem = model.problem()

H = problem.hamiltonian(E_z=0.0, B_z=0.0)
scan = problem.sweep(
    E_z=np.array([0.0, 1.0]),
    B_z=np.zeros(2),
)
print(H.shape)            # (264, 264)
print(scan.evals.shape)   # (2, 264)
```

The order of the spin records specifies the coupling chain:

```toml
[[isotopologues.test.spins]]
label = "nitrogen"
I = 1.0
couple_to = "J"

[[isotopologues.test.spins]]
label = "protons"
I = 0.5
equivalent_count = 2
couple_to = "F_N"

[[isotopologues.test.spins]]
label = "metal"
I = 0.5
couple_to = "F_core"
```

The first two records are required. The metal record is optional; omitting it
sets `I_M=0`. Reordering the records changes or invalidates the coupling chain.
For this backend `frame` must be `"lab"`, `M` must be `"all"` or `"blocks"`,
and `K` contains nonnegative magnitudes. A nonzero K magnitude generates both
signs when the basis is enumerated.

## Inspect before calculating

Use `model.describe()` to inspect the basis, active interactions, and parameter
sources for the selected model and isotopologue. `model.validate()` checks the
input structure. The basis and term matrices are built when you first access
`problem.kets` or `problem.term_matrices`. For a worked calculation, see
[Getting started](getting-started.md).
