# Open questions raised during implementation

## OQ-A — the Ω-doubling sign and the parity of the upper doublet component

Raised: Task 5 of `docs/superpowers/plans/2026-09-05-heff-v1-thf-tutorial.md`.
**Status: CLOSED 2026-09-05 — confirmed by Arian: the upper component has parity (−1)^J at every J (e above f uniformly). Code unchanged.**

### The physical invariant (convention-free)

The **upper** Ω-doublet component has parity `(−1)^J` at every J — equivalently,
**e lies `ω_ef J(J+1)/2` above f, uniformly in J** under Brown's 1975 rule.

Two primary sources, read independently:

- Ng 2022 Fig. 2, rendered and read as an image ([HAM] §2.3). Caption: positive
  (negative) parity levels are black (grey). At **J = 1** the grey line is above
  the black in both F = 3/2 and F = 1/2; at **J = 2** the black is above the grey
  in both F = 5/2 and F = 3/2.
- Gresh 2016's `k″ < 0` for the X ³Δ₁ lower state, which is the same ordering
  ([HAM] §2.3, third overturn note).

### Why heff's element is `−ω_ef J(J+1)/4` and not Ng's `(−1)^J ω_ef J(J+1)/4`

heff's parity operator (`heff.conventions.parity_phase`, thesis Eq. A.15 = B&C
Eq. 6.234) acts as

```
E* |J, Ω⟩ = (−1)^{J−S−ℓ+s} |J, −Ω⟩ = (−1)^{J−1} |J, −Ω⟩     (³Δ: S = 1, ℓ = s = 0)
```

so the combination `(|+1⟩ + σ|−1⟩)/√2` has parity `σ·(−1)^{J−1}`. The eigenvalues
of an off-diagonal `c` in this two-state space are `±|c|`; the **upper** one
(`+|c|`) is the symmetric combination (σ = +1) when `c > 0` and the antisymmetric
one (σ = −1) when `c < 0`.

- With the **J-independent negative** element `c = −ω_ef J(J+1)/4 < 0`: σ = −1 at
  every J, so the upper parity is `−(−1)^{J−1} = (−1)^J` at every J. ✔ invariant.
  (J = 1: phase +1, antisymmetric has parity −1, and it is the upper one, +ω_ef/2.
  J = 2: phase −1, antisymmetric has parity +1, still upper.)
- With **Ng Eq. C.3 transcribed literally**, `c = (−1)^J ω_ef J(J+1)/4`: σ = (−1)^J,
  so the upper parity is `(−1)^J(−1)^{J−1} = −1` at **every** J. ✘ The two
  alternations cancel instead of composing; the figure's alternation is lost.

Ng's `(−1)^J` is the same physics written in a ket-phase convention where
`E*|J, Ω⟩ = |J, −Ω⟩` (no `(−1)^{J−S}`), i.e. his `|J, Ω = −1⟩` carries a
J-dependent phase relative to the ket heff's parity operator acts on. Nothing in
any source read for [HAM] writes that phase down, which is why this was flagged
rather than assumed.

The **splitting magnitude** `ω_ef J(J+1)/2` is identical under either form and is
unaffected by this question.

### Status in code

- `heff/elements_c.py:omega_doubling` returns `−J(J+1)/4` (coefficient `ω_ef`),
  with the transposition recorded in its `cite`.
- Both halves of [HAM] V4 are **hard gates** in
  `tests/test_elements_c_fieldfree.py`: the splitting law
  (`test_V4_omega_doubling_splitting_law`) and the parity ordering
  (`test_V4_upper_doublet_parity_alternates_as_minus_one_to_the_J`). The parity
  gate fails at even J on Ng's literal form and at odd J on a globally flipped
  sign, so both outcomes are reachable.
- No constant, tolerance or parity phase was adjusted.

### The question for Arian

> Confirm that the upper Ω-doublet component of ThF⁺ X ³Δ₁ has parity `(−1)^J` at
> every J (e above f uniformly, Brown 1975). If yes, the code is correct as
> written and this item closes. If the ordering is instead parity `−1` at every J,
> Ng Eq. C.3 is literal in heff's convention and the sign in `omega_doubling`
> flips back — one character, caught either way by the V4 parity gate.

**Answer (Arian, 2026-09-05): yes — (−1)^J on top. Closed; code as written.**

A one-line question to K. B. Ng about the `|J, Ω = −1⟩` ket phase would settle the
convention side of it, and is the same message that settles OPEN-3.

## Erratum — [HAM] §2.8 g_F table

Raised: Task 6 of `docs/superpowers/plans/2026-09-05-heff-v1-thf-tutorial.md`
review, round 1. **Status: erratum in the document; the code is unaffected.**

[HAM] §2.8's printed g_F table is headed `[derived, G_par = 0.04756]`, but its
four rows reproduce only when the closed form
`g_F = −G_par γ_F + g_N (μ_N/μ_B) κ_F` is evaluated at G_par = 0.048 (Ng's
printed, rounded value), not at the header's own 0.04756.

For J = 1, F = 3/2:

- At **G_par = 0.04756**: g_F = −0.0148989 → 20.853 kHz/G. This matches
  [HAM] §2.8's own numerical confirmation line, stated three paragraphs above
  the table, and the measured |g_{F=3/2}| = 0.0149.
- At **G_par = 0.048**: g_F = −0.0150455 → 21.058 kHz/G. This matches the
  table's printed row (−0.015046, i.e. 21.06 kHz/G).

The code uses G_par = 0.04756 (`heff/params.py`, `thf_v1()`), consistent with
the document's own confirmation line, not with the printed table. The
document's table should be regenerated at G_par = 0.04756 (or its header
corrected to read 0.048, whichever Arian intends as the source value); the
code and tests do not pin to the table's numbers either way.

## OPEN-16 through OPEN-23 — v2 isotopologues and two-photon

Raised: `docs/superpowers/specs/2026-09-05-heff-v2-isotopologues-two-photon.md`
§7 (numbering continues [HAM] §7). Each is a place where the code would
otherwise be asserting something no source supports; none blocks any task,
because each is a flag, a default or a figure choice rather than a formula.

### OPEN-16 — the sign of A∥(Th)

Source: [SPEC-v2] §7, [TH] §2.2 (gap G4), `docs/lit/lookup-apar-th-sign-
convention.md`, [HAM] §9.4.3.

Skripnikov & Titov print −4163 (μ/μ_N) MHz, Denis print +1833 MHz from the
same defining equation, agreeing to 2.2 % in magnitude. The sign is what
**orders the F₁ manifold**, so every ²²⁹ThF⁺ level diagram depends on it.
`docs/lit/lookup-apar-th-sign-convention.md` shows the two groups' axis
conventions are opposite but that A∥ is invariant under a consistent
reversal, and that they agree on the analogous HfF⁺ constant — so the
disagreement is a real disagreement between two calculations, not a
convention mismatch, and the audit recommends A∥ < 0. The package therefore
defaults `thf_v2('229', a_par_th_sign='negative')`, trusting Skripnikov &
Titov 2015; `thf_v2('229', a_par_th_sign='positive')` trusts Denis 2015
instead — the keyword records which ab initio calculation you trust, not a
convention.

**Question:** (a) do you accept the audit's recommendation as the shipped
default (`a_par_th_sign='negative'`), given that an ab initio disagreement is
not settled by an audit of conventions? (b) should the notebook draw both
branches side by side anyway, since the F₁ ordering is the most visible
feature of every ²²⁹ThF⁺ figure? [Controller ruling, task 10: the notebook
draws the default `negative` branch only, and names the `a_par_th_sign=`
alternative in one sentence — pending Arian's answer to (b).]

### OPEN-17 — the eQq₂ normalisation bridge

Source: [SPEC-v2] §7, [TH] §3.2, [HAM] §9.4.4.

B&C's (9.52) at q = ±2 versus Petrov 2018 Eq. (23), with its √6 and Y₂₂.
[TH] §3.2 flags this **UNVERIFIED** and says every eQq₂-derived number
inherits the caveat. Two things are not pinned by the printed equations: the
scalar-product pairing in Petrov's Eq. (19), and whether the `√(2π/5)` there
is intended or a typo for `√(4π/5) = C²_q`. The package computes in B&C's own
normalisation (`eqq2_norm='bc_9p52_q2'`, the only implemented value);
`'petrov2018_eq23'` raises `NotImplementedError` naming this item rather than
guessing between the two candidate factors (`−1/√3`, `−1/√6`).

**Question:** if the two normalisations differ by a factor the derivation
cannot pin from the printed equations alone, is asking Petrov the right
move, or do we ship B&C's normalisation and label the HfF⁺-anchored estimate
as order-of-magnitude?

### OPEN-18 — eQq₀ and eQq₂ defaults

Source: [SPEC-v2] §7, [TH] §4.3–4.4 (gap G2).

There is no published ThF⁺ or ThO quadrupole coupling constant and no EFG at
Th, for any isotope or state. `thf_v2('229')` defaults `eQq0_Th` and
`eQq2_Th` to HfF⁺-anchored estimates (`status='estimate'`) because a zero
default would hide a term that [TH] §4.4 argues dominates the Ω-doublet
structure of ²²⁹ThF⁺.

**Question:** estimate-by-default, or zero-by-default with the estimate as
an opt-in knob? (v1 precedent points to estimate-by-default with
`status='estimate'`: that is how `c_I` ships.)

### OPEN-19 — c_I(Th)

Source: [SPEC-v2] §7, [TH] §4.6 (gap G3).

Unconstrained over three decades, ~1 kHz to ~1 MHz; [TH] §4.6 declines to
pick, and at the top of the range it would exceed the entire ¹⁹F hyperfine.
Defaulted to 0 (`status='held-fixed'`) for both ²²⁹Th and ²²⁷Th.

**Question:** does the notebook show a bracket sweep so the reader sees what
is at stake, or is a `note` enough? [This notebook states the bracket in
prose only, per the controller's no-parameter-sweep ruling.]

### OPEN-20 — ²²⁷Th: the spin and the moment

Source: [SPEC-v2] §7, [TH] §1.4, `docs/lit/lookup-227th-nuclear-moment.md`,
[HAM] §9.6.

I = (1/2⁺) is a tentative ENSDF assignment with a 5/2⁺ level 9.3 keV above
it, and no magnetic moment exists in any of the three compilations checked
(Stone INDC(NDS)-0794 (2019), the earlier Stone compilation, the IAEA NDS
live moments database — all three skip A = 227 at Z = 90 entirely). Arian's
ruling (2026-09-05): use the Schmidt single-particle moment
(`A_par_Th = +39821 MHz`, `status='placeholder'`), not `μ(²²⁷Th) = μ(²²⁹Th)`
— the Schmidt value for a tentative s₁/₂ odd-neutron ground state, with the
alternative μ(²²⁷Th) = μ(²²⁹Th) assumption (−7619 MHz, opposite sign, ~5×
smaller) recorded for comparison, and real deformed-actinide moments
typically ~2× smaller than Schmidt.

**Question:** (a) is the Schmidt assumption the one you want, or a different
one (scaled from a neighbouring odd-A actinide)? (b) should the notebook
also draw the I = 5/2 alternative, since it is one keyword argument and it
changes the whole level structure?

### OPEN-21 — does the K = 1 two-photon channel exist? RESOLVED

Source: [SPEC-v2] §7, [HAM] §9.5.3.

**Resolution (Task 1, 2026-09-05): no, not in exact closure.** `K = 1` is the
antisymmetric part of the polarisation dyad; in exact closure (a complete
opposite-parity intermediate manifold) it reduces to `P_X [d_a, d_b] P_X / 2`,
which is identically zero because the Cartesian components of the dipole
operator commute — verified numerically to `2.2 × 10⁻¹⁶` against `||K=0||`,
`||K=2||` of order 0.5 on a complete spherical-harmonic closure model ([HAM]
§9.5.3). `K = 1` reappears only at `O(δ/Δ)` in the resolved sum (a genuinely
different, unregistered operator) or at `O(1)` for a restricted (non-complete)
intermediate manifold — neither is the closure operator this package
implements. **`K ∈ {0, 2}` is registered in `REGISTRY_2G`; there is no
`alpha_K1_*` parameter, and the closure gate `tests/test_twophoton_closure.py`
(V27) is what pins this down.**

### OPEN-22 — J_max for ²²⁹ThF⁺

Source: [SPEC-v2] §7, [TH] §4.2.

The Th ΔJ = ±1 hyperfine is a ~2 GHz off-diagonal element whose second-order
shift is 60–135 MHz and, unlike the ²³²Th case, is not absorbable into `B₀`
because it depends on `F₁`. [TH] infers that J = 1–5 or 1–6 may be needed for
kHz-level J = 4 energies. `J_max` is a `StateSpec` knob and `j_convergence`
reports the convergence at any value; the *default* is a physics parameter
and therefore Arian's.

**Question:** `J_max = 4` for the notebook figures with a convergence table
(this notebook's choice), or `J_max = 6` throughout at 1152 states?

### OPEN-23 — the two-photon validity condition

Source: [SPEC-v2] §7, [2γ] §3.4, [HAM] §9.5.5.

The rank-K form is the `Δ ≫ B_i ≈ 7 GHz` limit (`B_i` the intermediate
state's rotational constant); JILA runs at 0.16–1.5 GHz. The notebook states
this plainly, quoting [HAM] §9.5.5.

**Question:** is that enough, or do you want the resolved-sum path built as a
product (not just a test fixture) against one of the two contested ladders
(Denis 2015 or Petrov 2018), with the choice labelled?
