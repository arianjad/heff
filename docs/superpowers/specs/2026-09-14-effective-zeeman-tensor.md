# Effective Zeeman G tensor for ThF⁺ X ³Δ₁ — implementation brief

Status: not started, 2026-09-14. Requested by Arian.
Closes `OPEN-10` of `docs/open-questions.md` ("Parity-dependent Zeeman terms are omitted").

## Read this first

**`docs/lit/lookup-effective-zeeman-tensor.md` is the single source of truth for the
physics, the constants and their provenance.** It was written 2026-09-14 against the
primary sources (arXiv:2503.02840, 1404.4024, 1704.06631, pulled and verified digit by
digit). Do not re-derive its numbers, do not re-open its citations, and do not substitute
your own values for the ones it lists. If you believe one of its numbers is wrong, stop and
say so rather than silently using a different one.

Everything below is scope and mechanics.

## Goal

Replace `heff`'s axial-only Zeeman term with the full three-parameter body-frame tensor

```
H_Z = mu_B B . G . J  with  G_{e/f} = diag(G_perp +- G_Delta, G_perp -+ G_Delta, G_par)
```

so that the model reproduces the measured ThF⁺ Zeeman structure including the
parity-dependent doublet difference Delta g, which it currently cannot generate at all.

## What is already true in the repo (verified 2026-09-14, do not re-check)

- **`elements_c2.py:40` `axial_geometry(bra, ket, ctx, *, k, q=None, p)` is a general
  rank-k, body-q, lab-p molecule-frame tensor geometry** with the full two-spectator
  (I_Th, I_F) recoupling, cited to B&C 5.172/5.174/5.186. It already raises on |q| > k. The
  ΔΩ = ±2 recoupling you need is `k=2, q=±2, p=0` — **it does not need to be derived.**
- `elements_c2.py:284` `zeeman_Gpar` is exactly `sign * mu_B * Om * axial_geometry(k=1, q=0,
  p=0)`. That is the template: a body-frame scalar factor times the geometry.
- `elements_c.py:132` `dipole_geometry(bra, ket, I, p)` is the one-spin twin, but it
  **hardcodes rank 1 and q = 0** (`w3j(J2, 1, J1, -Om, 0, Om)`). It needs the same (k, q)
  generalization; `axial_geometry`'s docstring states it reduces to `dipole_geometry` at
  I_Th = 0, k = 1, q = 0, so the two must stay consistent.
- `elements_c.py:55` `omega_doubling` already declares `rules=Rules(dOm=(-2.0, 2.0), ...)`,
  so the assembler's sparsity mask supports ΔΩ = ±2.
- There are **two backends**: `elements_c.py` (`backend = "case_c"`, one spin, ²³²ThF⁺) and
  `elements_c2.py` (two spins, ²²⁹ and ²²⁷). Every new term needs both, and every new
  parameter needs all three isotopologue blocks in `heff/models/thf_plus.toml`.
- Valid `status` values (`params.py:41`): `measured`, `ab-initio`, `derived`, `estimate`,
  `held-fixed`, `stale`, `unspecified`, `placeholder`. Use `derived`.
- Env: `conda run -n structure python ...`, `heff` installed editable.

## The decomposition

`B·G·J` with G body-diagonal splits into three covariant pieces (derivation and its
verification against Petrov's Eqs. 12–13 are in the lookup doc §1.1 and §5):

```
B.G.J = (G_par - G_perp)(B.n)(J.n) + G_perp (B.J) + G_Delta[(B.ex)(J.ex) - (B.ey)(J.ey)]
```

| term | operator | work |
|---|---|---|
| A | (J·n̂)(n̂·B) | **exists** (`zeeman_Gpar`, both backends). Only its coefficient changes: `G_par` → `G_par - G_perp` |
| B | B·J | new. Rank-1 lab operator on total J; ΔJ = 0, ΔΩ = 0, ΔF = 0,±1 |
| C | body q = ±2 | new. Parity-dependent; `axial_geometry(k=2, q=±2, p=0)` supplies the geometry |

Diagonal elements are Ω²m/[J(J+1)], m, and ±(doublet), giving
`g(J) = -G_perp - (G_par - G_perp)/[J(J+1)] +- G_Delta` — the closed form gate 1 checks.

## Parameters

Replace the single fitted `G_par` with the triple from the lookup doc §4. **In all three
isotopologue blocks** of `heff/models/thf_plus.toml` (the ²²⁹/²²⁷ blocks carry transferred
copies with `status = "estimate"` and a transfer note — preserve that pattern), and in the
matching `params.py` defaults:

| symbol | value | status | source string |
|---|---|---|---|
| `G_par` | 0.0467961 | derived | Petrov & Skripnikov arXiv:2503.02840 Eqs. (3)–(15), second-order reduction — `docs/lit/lookup-effective-zeeman-tensor.md` §2 |
| `G_perp` | 1.02769e-3 | derived | same |
| `G_Delta` | -2.05580e-4 | derived | same; **sign is convention-dependent** — see below |

The current `G_par = 0.04756` is the axial-only fit that forces |g_{F=3/2}| = 0.0149
exactly. It is **not** independent of Petrov's 0.047; both are fits to the same datum. Its
note and `tests/test_observe.py:162` must be updated to say so.

`G_Delta`'s sign depends on the e/f and transverse-axis phase convention (lookup doc §1).
`heff`'s `ef_rule` fork inverts every e/f label at S = 1 (`conventions.py:38`). Add a
`zeeman_tensor_sign` entry to `conventions.py` recording which fork the stored sign belongs
to, rather than asserting a physical sign. **Gate on |Delta g|, never on its sign.**

## Acceptance gates

Each must have both PASS and FAIL reachable — write the test so a deliberate error in the
thing it guards actually breaks it, and say in the docstring which error that is.

1. **Closed form.** `observe.g_factors` at E = 0, B → 0 reproduces
   `-G_perp - (G_par - G_perp)/[J(J+1)] +- G_Delta` for J = 1…5 to round-off.
   Fails on any coefficient-mapping error in the decomposition above.
2. **J = 1 prediction.** |g(J=1, F=3/2)| = 0.014987, inside Ng 2022's 0.0149(3).
   This is now a *prediction*, not a fit. Fails on a wrong term-B projection or a flipped
   nuclear-Zeeman sign.
3. **Doublet difference.** 2|G_Delta| projected to F = 3/2 gives 2.741e-4. Accept the band
   [2.3, 3.6]e-4, whose ends are Petrov's full seven-state calculation and the
   Leanhardt ω_ef/(2B_e) scaling (lookup doc §4). Fails if term C's magnitude is wrong.
4. **Parity.** Under `conventions.parity_operator`, term C is parity-ODD; terms A and B are
   parity-EVEN. Fails on a wrong ΔΩ = ±2 phase.
5. **Regression.** With `G_perp = G_Delta = 0` and `G_par = 0.04756`, every existing number
   reproduces bit-for-bit. The whole current suite is the gate.
6. **Backend parity.** ²³²ThF⁺ run through `elements_c` and through `elements_c2` at
   I_Th = 0 must agree to round-off, for all three new/changed terms.

## Stop conditions — do not guess past these

- **The reduced matrix element for terms B and C.** `axial_geometry` gives the *geometry*;
  what multiplies it is a J-space reduced matrix element, and the body components of **J**
  have anomalous commutation. The sourced operator form is Leanhardt 2011 Eq. (64)
  (`docs/digest-literature-thf-plus.md` §3.2):
  `H_ZeemanDoub = +(1/2) g'_rS mu_B (B+ J+ S+^2 + B- J- S-^2)`.
  If you cannot get from that to a coefficient on `axial_geometry(k=2, q=±2, p=0)` **with a
  citation**, stop and report. Do not fit the coefficient to reproduce gate 3 — that makes
  gate 3 tautological and the operator wrong off-resonance (nonzero E, ΔJ ≠ 0, hyperfine
  mixing), which is precisely where `heff` is used.
- **Anything requiring a new convention fork.** Record it in `conventions.py` and report;
  do not pick one silently.
- `C2V-Molecules/atm_core/physics.py:244` carries the same tensor machinery for an
  asymmetric top and is worth reading for shape (`g_l[k][q]` decomposition, Sears 1984
  eq. 21). **Do not port it**: its tensor contracts B with **S**, a spin that decouples from
  rotation, which is why it has a clean 9j. Ours contracts B with **J**. Substituting S → J
  in that expression is wrong.

## Out of scope

- Third-order / Zeeman-distortion terms. Blocked on ab initio input Petrov does not publish
  (excited-to-excited matrix elements) — lookup doc §5.3.
- Building the explicit seven-state benchmark. It would supersede this model rather than
  test it; separate task.
- Re-deriving or re-sourcing any constant. See "Read this first".

## Working rules

- **Checkpoint at every green milestone.** Commit when a gate first passes — do not batch
  the work into one commit at the end. An interruption should cost at most one milestone.
  Suggested milestones: (i) `axial_geometry`/`dipole_geometry` rank-k generalization with
  gate 6 passing; (ii) term B + gates 1, 2; (iii) term C + gates 3, 4; (iv) parameters and
  TOML across all three isotopologues + gate 5.
- Update this spec's Status line and add a Progress section as you go.
- After the refactor, update the markdown and docstrings that describe the Zeeman term to
  say what the code does *now* — including `docs/open-questions.md` `OPEN-10`,
  `docs/thf-plus-x3delta1-effective-hamiltonian.md`, and the `conventions.py`
  `zeeman_sign` docstring.
- Physics-parameter choices beyond the table above (field ranges, J_max for convergence
  checks, tolerances) go to Arian for a yes before running, not chosen unilaterally.
