# heff

Effective Hamiltonians for molecules, term-matrix first. v1 covers Hund's case (c)
for ²³²Th¹⁹F⁺ X ³Δ₁: 96 states (J = 1–4, Ω = ±1, I(¹⁹F) = ½), nine Hamiltonian
terms plus an opt-in PT-odd pair.

## Quickstart

```bash
conda create -n heff python=3.12 numpy scipy sympy pytest jupyter ipykernel nbclient matplotlib -y
conda run -n heff pip install -e .
conda run -n heff python -m pytest tests/ -q
conda run -n heff python notebooks/build_tutorial.py
conda run -n heff jupyter execute --inplace notebooks/ThF_plus_X3Delta1_Tutorial.ipynb
```

The last two lines regenerate the tutorial notebook from its generator script and
run it on a fresh kernel. The committed notebook already carries its outputs, so
reading it needs neither step.

## What it is

`H = Σ_k c_k M_k`. The term matrices `M_k` are built once per basis block and
carry no parameters; a parameter change is a weighted resum (one BLAS call), and
a whole sweep is one `np.tensordot` followed by a batched `eigh`. Because
`∂H/∂knob` is itself a sum of catalogue matrices, derivatives are exact — g-factors
and induced dipoles come from Hellmann–Feynman, not finite differences. Adding a
term is one decorated function plus one parameter entry; its selection rules are
declared as data and its primary source is attached as a citation string.

Every result carries the convention block it was computed under, the list of
terms that were actually active, and — for a sweep — which ordering, gauge and
assignment policy produced its labels.

## Modules

| module | what it does |
| --- | --- |
| `spec.py` | `StateSpec`, the case-(c) ket enumerator, basis invariants, `Blocking` by signed m_F |
| `conventions.py` | the convention block, parity phase and parity operator, e/f labels |
| `params.py` | `Param`/`ParamSet` with unit, uncertainty, source, status and provenance table; `thf_v1()` |
| `formalism.py` | R²↔N² conversion, lifted verbatim from Molecule-Structure |
| `terms.py` | `Ctx`, `Rules`, the `@term` registry, and the selection-rule gate |
| `elements_c.py` | the case-(c) matrix elements: rotation, centrifugal, Ω-doubling, hyperfine (ΔJ = 0 and ±1), `c_I`, Stark, Zeeman, nuclear Zeeman, PT-odd |
| `wigner.py` | 3j/6j/9j behind an `lru_cache`; sympy is imported inside the kernels, not at module scope |
| `assemble.py` | `build_term_matrices`, `hamiltonian`, `hamiltonian_batch`, `sweep_coefficients`, `vertex` |
| `engine.py` | chunked batched `eigh`, `sweep`, and the `SweepResult` artifact |
| `track.py` | eigenvector assignment and sign-gauge kernels, lifted from C2V-Molecules |
| `observe.py` | exact derivatives, g-factors, induced dipoles, expectation values, named-pair differentials |
| `spectra.py` | E1 dipole matrices, line strengths (sum over polarisation, then square), line labelling |

`import heff` pulls numpy only: sympy is lazy inside `wigner`, scipy inside
`track`, and matplotlib is never imported by the package — plotting lives in the
notebook. A test gates this.

## Conventions

Defaults on `Conventions`, stamped onto every result. Each names a real fork in
the literature; the alternative is available and documented in
`docs/thf-plus-x3delta1-effective-hamiltonian.md` §7.

| field | default | meaning |
| --- | --- | --- |
| `n_hat` | `F_to_Th` | JILA orientation: Ω = +1 for Λ = +2, `d_mf = +3.37 D` |
| `ef_rule` | `brown1975` | e ⇔ parity `+(−1)^J` (integer J), `+(−1)^{J−½}` (half-integer) |
| `zeeman_sign` | `plus_Gpar` | `+G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B` |
| `zeeman_energy` | `E_minus_g_muB_B_mF` | `E = −g μ_B B m_F` |
| `dg_def` | `Delta` | report `Δg = g^u − g^ℓ` (`delta` gives half that) |
| `edm_factor` | `ng` | `−(d_e E_eff + W_TP k_TP) Ω/|Ω|`, no Leanhardt ½ |
| `dipole_origin` | `center_of_mass` | `d_mf` is origin-dependent for an ion; the Th-nucleus origin differs by 0.72 D |
| `formalism` | `R2` | rotational operator convention |

Internal unit is MHz throughout; fields are `E_z` in V/cm and `B_z` in G.

## Validation

Gated by the test suite — `conda run -n heff python -m pytest tests/ -q` passes in
full, with the two tier-D literature comparisons skipped unless
`HEFF_RUN_LITERATURE=1` (the opt-in tier described below).

- **Kernel gates (A1–A8)**, none of which depend on a particular Hamiltonian:
  the resum identity and the `tensordot` fast path on random matrices, the exact
  derivative kernel, tracking and gauge invariants, every registered term against
  its own declared selection rules, the cached Wigner backend against uncached
  sympy on a random half-integer sample, basis invariants, and import cost and
  dependency hygiene.
- **Symmetry and structure (B1, B7, B7b, B8)**: the Wigner–Eckart structure of
  the dipole, its sum rule and tensor identity, and the cancellation that only a
  sum-then-square line strength can produce.
- **ThF⁺ closed forms** ([HAM] V-checks, run by default): the Ω-doublet splitting
  law and the parity ordering of the upper component (V4), the g-factor relation
  `g_F = −G∥ γ_F + g_N (μ_N/μ_B) κ_F` and the sign of the Zeeman operator (V5),
  the Zeeman being even in Ω (V6), parity commuting at zero field and not with the
  Stark term (V8), and the enumerator invariants (V11).
- **Comparisons to published numbers** are **opt-in**, because they are
  convention-dependent:
  `HEFF_RUN_LITERATURE=1 conda run -n heff python -m pytest tests/ -q`.
  Each names the convention block it assumes.

Every gate's docstring names the failure mode it uniquely catches, and both PASS
and FAIL are reachable. There are no snapshot, hash or pinned-spectrum tests, and
no constant or tolerance was tuned to make a gate pass. One physics question
surfaced during implementation — the parity ordering of the Ω doublet, OQ-A in
`docs/open-questions.md` — and was confirmed by Arian on 2026-09-05; OQ-A is
closed and the code is unchanged. Both halves of it are hard gates either way.

## Documentation

- Physics: `docs/thf-plus-x3delta1-effective-hamiltonian.md` — every matrix
  element in `heff/elements_c.py` cites it or the primary source it cites.
- Design: `docs/superpowers/specs/2026-09-05-heff-design.md`.
- Open items: `docs/open-questions.md`, and §7 of the physics document.
- Tutorial: `notebooks/ThF_plus_X3Delta1_Tutorial.ipynb`, generated by
  `notebooks/build_tutorial.py`.
