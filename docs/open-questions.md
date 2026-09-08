# Scientific limitations and conventions

Use this page to identify which inputs or approximations limit your
calculation. It also records sign and labeling conventions, with `OPEN-*`
identifiers shared by the code and the
[ThF+ Hamiltonian reference](thf-plus-x3delta1-effective-hamiltonian.md).

## ThF+ X 3Delta1 model

| ID | Current status | Consequence |
|---|---|---|
| `OPEN-1` | The observed 3Delta fine-structure intervals cannot be represented by one spin-orbit constant. The shipped model therefore keeps only Omega = +/-1 and fitted effective constants. | A model that includes 3Delta2 and 3Delta3 needs separate electronic origins or a larger coupled-state model. |
| `OPEN-2` | **Resolved.** The package uses Brown's e/f rule: the upper Omega-doublet component has parity `(-1)^J` and is e. | Applying the thesis's `S=1/2` specialization at `S=1` would invert every e/f label. |
| `OPEN-3` | **Resolved for the implementation.** The package uses `+G_par mu_B (J.n)(n.B)`, which reproduces Ng's printed g-factor relation and measured magnitude. Ng thesis Eq. C.6 prints the opposite sign. | The source-level discrepancy remains worth author confirmation, but the implemented branch is fixed and tested. |
| `OPEN-4` | The experiment measured `|g_F|`; theory selects the negative branch used by default. | Signed observables inherit a theory-supported, experimentally unresolved sign. |
| `OPEN-5` | Only `A_par = 2a - b_F - 2c/3` is known for 19F. The separate microscopic constants are unavailable. | The axial Omega = +/-1 model is supported; perpendicular or Delta-Omega = +/-1 hyperfine terms are not. |
| `OPEN-6` | The shipped `c_I(19F) = 20 kHz` is a CsF-scaled sensitivity estimate with unknown sign and unquantified transfer error. | The fitted `A_par` may absorb roughly 40--120 kHz of an omitted or mismodelled spin-rotation contribution. |
| `OPEN-7` | The rotational g-factor `g_r` is omitted and not known for ThF+. | Predicted J = 2--4 g-factors carry an unquantified error of order 1 percent. |
| `OPEN-8` | Hyperfine-dependent Omega doubling `e_Delta` is omitted. The estimate in the Hamiltonian note is only an order-of-magnitude 1--10 kHz scale. | The model cannot predict an F-dependent Omega-doublet splitting at that scale. |
| `OPEN-9` | No ThF+ nuclear-spin-dependent parity-violation coefficient `W_a` or `W_P` was found in the sources reviewed. | NSD-PV requires a new electronic-structure input and operator choice. |
| `OPEN-10` | Parity-dependent Zeeman terms are omitted. | The model does not generate the reported zero-field differential g-factor; use it for the included axial Zeeman physics only. |
| `OPEN-11` | **Resolved as a convention.** The default axis is JILA's `n_hat = F_to_Th`; the alternative `Th_to_F` convention is explicit. | Axis reversal flips signed Omega, dipole, and effective-field quantities together. |
| `OPEN-12` | **Resolved as an output convention.** The default is `Delta g = g^u - g^l`; the thesis convention `delta g = Delta g/2` remains selectable. | Always retain the convention stamp when comparing results. |
| `OPEN-13` | Gresh's printed `k''` sign depends on upper-state branch bookkeeping. | The package takes the Omega-doublet magnitude and ordering from the microwave result rather than inferring them from that sign. |
| `OPEN-14` | The package adopts `E_eff = 35.0 GV/cm` with a 7 percent scale; published calculations give 35.2 and 37.3 GV/cm. | This is an adopted ab-initio input, not a direct measurement. |
| `OPEN-15` | Petrov's printed `D = -0.133 a.u.` is inconsistent with the measured 3.37 D magnitude and is likely missing a factor of ten. | The package uses Ng's measured, center-of-mass-origin dipole and does not use the printed Petrov number. |

## Odd-thorium isotopologues and two-photon operator

### `OPEN-16`: sign of A_parallel(229Th)

The bundled `229Th19F` model uses the negative branch,
`A_par_Th = -1510 MHz`, following Skripnikov and Titov. The positive Denis
branch remains a real disagreement between electronic-structure calculations;
an axis reversal does not resolve it. The older `thf_v2()` API can select the
positive branch explicitly, while the bundled TOML fixes the negative branch.
See section 9.4.3 of the Hamiltonian note.

### `OPEN-17` and `OPEN-18`: thorium quadrupole inputs

The signed normalization bridge from Petrov's `eQq0/eQq2`
definitions to the Brown-and-Carrington matrix elements is not established.
The bundled `229Th19F` values `eQq0_Th = -2600 MHz` and
`eQq2_Th = +300 MHz` are `placeholder` sensitivity points with no numerical
uncertainty. They are not a convention-validated pair or a prediction. The
only implemented `eqq2_norm` is `bc_9p52_q2`; selecting
`petrov2018_eq23` raises rather than guessing. A matched ThF+ electric-field-
gradient calculation is the missing input. See the
[quadrupole audit](superpowers/reports/2026-09-08-thf-quadrupole-estimate-audit.md).

### `OPEN-19`: thorium spin rotation

No ThF+ value for `c_I_Th` is available. The package holds it
at zero for both odd isotopologues. Zero switches off an unknown operator; it
is not a physical estimate. The earlier 1 kHz--1 MHz bracket is a scale study,
not an uncertainty interval.

### `OPEN-20`: 227Th spin and magnetic moment

The bundled `227Th19F` model uses the
`A_par_Th = +39.821 GHz`, `g_N_Th = -3.826` Schmidt single-particle values as
explicit `placeholder` stress-test inputs. The tentative `I^pi = (1/2+)`
assignment does not imply a pure spherical `s_1/2` neutron state.

Minkov et al. (2024) predict
`mu(227Th) = -0.0860 mu_N`, hence `g_N = -0.1720` and
`A_par_Th = +1.790176 GHz` with the same molecular electronic factor. The
field-plot dataset uses this value, but the bundled package model does not.
The calculation quotes no calibrated uncertainty and omits Coriolis and
collective mixing corrections. The ground-state moment and low-lying spin
assignments remain experimental and nuclear-model limitations. See the
[nuclear estimate audit](superpowers/reports/2026-09-08-thf-nuclear-estimate-audit.md).

### `OPEN-21`: rank-one two-photon channel

Exact closure over a complete
opposite-parity intermediate space makes the antisymmetric `K = 1` tensor
zero. The registered closure operator therefore contains only `K = 0, 2`.
`K = 1` can reappear for energy-resolved denominators or a restricted
intermediate manifold; those are different operators and are not implemented.

### `OPEN-22`: rotational truncation

The bundled model uses `J_max = 4`. Convergence must be
checked for the selected odd-isotope parameters because thorium hyperfine can
mix adjacent J manifolds strongly. The separate field-plot calculation used
`J <= 8` and sampled `J <= 7` versus `J <= 8`; that result does not prove the
bundled Schmidt stress-test model converged at `J_max = 4`.

### `OPEN-23`: closure validity for the experiment

The rank-K closure form requires detuning large compared
with the intermediate rotational structure, approximately 7 GHz for the
relevant ThF+ states. The cited JILA detunings are 0.16--1.5 GHz. The package
therefore supplies the closure operator's tensor structure and placeholder
polarizabilities, not a quantitative model of those experiments. A resolved
intermediate-state sum requires a selected ladder, energies, dipoles, and
consistent phase and polarization conventions.
