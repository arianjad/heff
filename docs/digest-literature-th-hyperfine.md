# Literature digest: Th hyperfine structure for ²²⁷ThF⁺ and ²²⁹ThF⁺ (X ³Δ₁)

Assembled 2026-09-05 for the extension of `heff` beyond ²³²Th¹⁹F⁺. Companion to
`docs/thf-plus-x3delta1-effective-hamiltonian.md` ([HAM]), whose §2.11 states the
²³²Th case: I(²³²Th) = 0, so every Th hyperfine and Th quadrupole term is
identically zero. This document collects what is known when I(Th) ≠ 0.

Every claim carries a source. Tags: `[measured]`, `[ab initio]`, `[derived]`
(algebra or arithmetic run here from cited inputs), `[estimate, method]`,
`[inferred]`, `UNVERIFIED`. Units are MHz unless stated. Reference scales from
[HAM]: ¹⁹F A∥ = −20.1(1) MHz, ω_ef = 5.29(5) MHz (J = 1 Ω-doublet splitting),
B₀ = 7.2743325 GHz.

**Headline.** For ²²⁹ThF⁺ the Th magnetic hyperfine constant is ≈ 1.5 GHz — two
orders above the ¹⁹F constant and 21 % of B₀. Two consequences that break the
perturbative logic of [HAM]: the ΔJ = ±1 Th hyperfine term, a 2 kHz nicety for
¹⁹F, becomes a 60–135 MHz second-order shift; and the Th electric quadrupole has
a ΔΩ = ±2 component (`eQq₂`) which, by the ¹⁷⁷HfF⁺ precedent, is expected to
exceed ω_ef by an order of magnitude and to dominate the Ω-doublet mixing.
For ²²⁷Th (I = 1/2, tentative) there is no quadrupole at all and no measured
magnetic moment located in the checked sources. A 2024 nuclear-theory moment now supplies a model-dependent hyperfine alternative; see §1.4 and the current audit.

---

## 1. Nuclear data

### 1.1 Table

| nuclide | Iᵖ | T₁/₂ | μ (μ_N) | Q (e·b) | source |
|---|---|---|---|---|---|
| ²²⁷Th | (1/2⁺) *tentative* | 18.697(7) d | −0.0860 [2024 theory; §1.4], no measured value located | **0 by symmetry** (I = ½) | IAEA Live Chart `ground_states`, ENSDF cut-off 15-Jan-2016, ICTP-2014 Workshop Group; retrieved 2026-09-05 |
| ²²⁹Th | 5/2⁺ | 7880(120) y | 0.46(4) | 3.11(6) | IAEA Live Chart `ground_states`, ENSDF cut-off 1-Jun-2008, Browne & Tuli; retrieved 2026-09-05 |
| ²²⁹Th | 5/2⁺ | — | **0.360(7)** | **3.11(6)** | Safronova, Safronova, Radnaev, Campbell, Kuzmich, PRA 88, 060501(R) (2013) = arXiv:1305.0667, abstract |
| ²²⁹Th | 5/2⁺ | (7917 y quoted) | **0.366(6)** | **3.11(2)** | Porsev, Safronova, Kozlov, arXiv:2107.14723, abstract + Tables 3–4 |
| ²³²Th | 0⁺ | 1.40(1)×10¹⁰ y | 0 | 0 | IAEA Live Chart, ENSDF cut-off 1-Nov-2004 |
| ¹⁹F | 1/2⁺ | stable | 2.628321(4) | 0 (I = ½) | IAEA Live Chart, ENSDF cut-off 31-Oct-1994 |

Live-Chart query used for all four nuclides:
`https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=<A><sym>`
and, for the ²²⁷Th level scheme, `fields=levels&nuclides=227th`.

Half-life caveat: the IAEA/ENSDF evaluated value is 7880(120) y; Flambaum &
Mansour, arXiv:2204.01290 §I.2 quote 7917 y (secondary citation, their Ref. 31).
The two overlap at 1 σ; use 7880(120) y and note that neither value affects any
Hamiltonian term.

### 1.2 The ²²⁹Th magnetic moment: three values, one of them stale

They disagree at the 25 % level, and that propagates directly into every Th
hyperfine constant in §2, because those are tabulated **per unit nuclear moment**.

- **μ = 0.46(4) μ_N** — Gerstenkorn et al., J. Phys. (Paris) 35, 483 (1974), as
  cited by Porsev 2021 and by Flambaum & Mansour 2022 Table 1. Still the value in
  the IAEA/ENSDF ground-state table (cut-off 2008). It is the value that
  Skripnikov & Titov 2015 and Denis et al. 2015 had available.
- **μ = 0.360(7) μ_N** — Safronova et al. 2013 (arXiv:1305.0667), from the Th³⁺
  6d and 5f hyperfine measurements [Campbell et al., PRL 106, 223001 (2011)]
  combined with all-order relativistic atomic theory. Abstract, verbatim: five
  times more accurate and "22% smaller than the best previous value μ=0.46(4)".
- **μ = 0.366(6) μ_N** — Porsev, Safronova, Kozlov 2021 (arXiv:2107.14723),
  CCSDT recomputation of the same Th³⁺ constants. Verbatim: "we find the nuclear
  magnetic moment, μ_I = gI ≈ 0.366(6)". They note the 2013 error bar "did not
  include uncertainty due to the magnetization distribution" (Bohr–Weisskopf).

**Recommendation: μ(²²⁹Th) = 0.366(6) μ_N**, i.e.
**g_N(²²⁹Th) = μ/I = 0.1464(24)** `[derived]`. The 1974 value is superseded;
never use it to rescale a published A∥, and always rescale published A∥ values
that were computed with it (§2.2).

### 1.3 The ²²⁹Th quadrupole moment

Q = 3.11 e·b is stable across all sources; only the uncertainty moves,
3.11(6) → 3.11(2) e·b (Porsev 2021, a factor-three reduction over Safronova
2013; Porsev's recommended Q is the weighted average over four Th³⁺ states).
**Use Q(²²⁹Th) = 3.11(2) e·b.**

### 1.4 ²²⁷Th — the isotope with no magnetic hyperfine input

From the IAEA Live Chart `levels` query (ENSDF cut-off 15-Jan-2016, retrieved
2026-09-05):

- Ground state Jᵖ = **(1/2⁺)** — parenthesised, i.e. a *tentative* ENSDF
  assignment. First two excited states 9.3(3) keV (5/2⁺) and 24.38(3) keV (3/2⁺),
  the pattern of a K = 1/2 band.
- **No magnetic-dipole and no electric-quadrupole moment is tabulated** for the
  ground state or for those excited states (the `magnetic_dipole` and
  `electric_quadrupole` fields are empty).

Consequences and current input status:

1. **²²⁷ThF⁺ has no Th electric quadrupole term.** A rank-2 nuclear operator
   vanishes for I ≤ 1/2 — the 3j symbol (I 2 I; −I 0 I) in the denominator of
   B&C Eq. (9.53) has no allowed triangle. eQq₀ and eQq₂ are identically zero,
   exactly as for ¹⁹F ([HAM] §2.11). Only the Th magnetic dipole and Th
   spin–rotation survive, and the Th hyperfine has the *same algebraic shape* as
   the ¹⁹F structure already coded, with I_Th = 1/2.
2. **A model-dependent hyperfine estimate is available.** [Minkov et al., Phys. Rev. C 110, 034327 (2024), Table IV](https://arxiv.org/abs/2408.11010) predicts μ(²²⁷Th) = −0.0860 μ_N for its octupole-deformed 1/2 ground solution. With I = 1/2 and the existing molecular factor −10408 MHz, this gives g_N = −0.1720 and A∥ = +1790.176 MHz. This is a model prediction without a calibrated uncertainty, not a measurement. See the [current audit](superpowers/reports/2026-09-08-thf-nuclear-estimate-audit.md); the selected Schmidt default is retained separately.

If the tentative 1/2⁺ assignment is wrong (the 9.3 keV 5/2⁺ is only 9 keV away),
the quadrupole conclusion and moment-to-g conversion change. The ENSDF parentheses are a live caveat.

---

## 2. ThF⁺ X ³Δ₁ Th hyperfine constants

### 2.1 What exists, in one table

| quantity | isotope / state | value | status | source |
|---|---|---|---|---|
| A∥(Th) | ²²⁹ThF⁺ X ³Δ₁ | **−4163 (μ_Th/μ_N) MHz** | `[ab initio]` 38e-2c-CCSD(T) + corrections | Skripnikov & Titov, PRA 91, 042504 (2015) = arXiv:1503.01001, Table II row FINAL(ThF⁺); definition their Eq. (10) |
| A∥(Th) | ²²⁹ThF⁺ Ω = 1 | **+1833 MHz** (their μ = 0.45 μ_N, I = 5/2) | `[ab initio]` GASCI ",10" / TZ′ | Denis, Nørby, Jensen, Gomes, Nayak, Knecht, Fleig, NJP 17, 043005 (2015), abstract, §3.2.2, Table 6/7, conclusion; definition their Eq. (2) |
| A∥(Th) | ²²⁹ThO H ³Δ₁ | −2949 (μ_Th/μ_N) MHz | `[ab initio]` | Skripnikov & Titov 2015 Table II row FINAL(ThO), citing their Refs. 28–29 |
| W_M | ²²⁹ThF⁺ X ³Δ₁ | 0.88 ×10³³ Hz/(e·cm²) | `[ab initio]` | Skripnikov & Titov 2015 Table II (already in [HAM] §2.11) |
| eQq₀(Th) | ThF⁺ or ThO, any isotope | **no value found** | — | see §5 gap G2 |
| eQq₂(Th) | ThF⁺ or ThO, any isotope | **no value found** | — | see §5 gap G2 |
| EFG at Th | ThF⁺ or ThO | **no value found** | — | see §5 gap G2 |
| c_I(Th) | ThF⁺, any isotope | **no value found** | — | see §5 gap G3 |
| A∥(Th), eQq | ²²⁷ThF⁺ | A∥ ≈ +1.790 GHz from 2024 nuclear theory and the chosen molecular factor; Q ≡ 0 at I = ½ | — | §1.4 |

Nothing here is measured. **No ThF⁺ Th hyperfine constant of any kind has been
measured**, for any isotope; JILA's spectroscopy (Gresh 2016, Ng 2022) is on
²³²ThF⁺, where I(Th) = 0.

### 2.2 The two ab initio A∥ values agree in magnitude and disagree in sign

Both papers define A∥ by the *same* equation — the axial expectation value of the
electron nuclear-spin–orbit operator divided by IΩ:

- Skripnikov & Titov Eq. (10): `A∥ = (μ_Th /(IΩ)) ⟨Ψ| Σ_i (α_i × r_i / r_i³)_ζ |Ψ⟩`
- Denis et al. Eq. (2): the identical expression (their text: "the parallel
  magnetic hyperfine interaction constant A∥ is defined as the z projection of
  the expectation value of the corresponding perturbative Hamiltonian in Dirac
  theory").

Because 1/I is inside the definition, **A∥ ∝ g_N = μ/I**, and the isotope-free
electronic factor is

```
G_el ≡ A∥ / g_N     Skripnikov: −4163 × (5/2) = −10 408 MHz   [derived]
                    Denis:      +1833 / 0.18  = +10 183 MHz   [derived]
```

`[derived]`; the two agree to **2.2 % in magnitude and differ in sign**.

**Rescaled to the modern moment (μ = 0.366 μ_N, §1.2):**

| source | A∥(²²⁹ThF⁺) |
|---|---|
| Skripnikov & Titov 2015 | **−1524 MHz** `[derived]` |
| Denis et al. 2015 | **+1491 MHz** `[derived]` |
| (Skripnikov at the stale μ = 0.46) | −1915 MHz |

**Use |A∥(²²⁹ThF⁺)| = 1.51(6) GHz** `[derived]`, the mean of the two rescalings
with a spread-based uncertainty; add the authors' own ~7 % ab initio uncertainty
in quadrature if a hard error bar is needed (Skripnikov & Titov quote 7 % for
their ThF⁺ properties, [HAM] §8).

**The sign is UNVERIFIED.** Neither paper states its molecular-axis convention
(which way ζ = n̂ points, Th→F or F→Th) or the sign of Ω, and A∥ changes sign
under n̂ → −n̂ combined with Ω → −Ω only if the convention is applied
inconsistently. Since `heff` fixes `n_hat = F_to_Th` ([HAM] §7) and since the
sign of A∥ is what orders the F₁ manifold (§4), this must be settled before any
²²⁹ThF⁺ level ordering is published. Options: (i) ask the authors; (ii) re-derive
the sign from the ³Δ₁ = σδ configuration and the sign of ⟨1/r³⟩ (the orbital
term `a` dominates for a δ electron, so A∥ should follow the sign of Λ·Ω, which
in the [HAM] convention Λ = +2, Ω = +1 is positive-definite for `a > 0` — this is
an argument, not a verification; do not code it). Recorded as gap G4.

### 2.3 The transferable anchors: ¹⁷⁷/¹⁷⁹HfF⁺ X ³Δ₁

Petrov, Skripnikov, Titov & Flambaum, "Evaluation of CP violation in HfF⁺",
PRA 98, 042502 (2018), received 24 July 2018 = arXiv:1809.06701, is the closest
published treatment of
exactly the problem this digest is about: a ³Δ₁ molecular ion with a heavy
nuclear spin *and* ¹⁹F, including the quadrupole. Local copy:
`docs/lit/petrov2018-CP-violation-HfFplus.txt`.

Their §III constants (their Eqs. 20–23 define them):

| constant | ¹⁷⁷HfF⁺ (I=7/2, g=0.2267, Q=3.365 b) | ¹⁷⁹HfF⁺ (I=9/2, g=−0.1424, Q=3.793 b) | status |
|---|---|---|---|
| A∥(F) | −62.0 MHz | −62.0 MHz | `[measured]`, Cairncross et al. 2017 (their Ref. 1) |
| A∥(Hf) | −1429 MHz | +898 MHz | `[ab initio]`, their Ref. 15 |
| eQq₀ | **−2100 MHz** | **−2400 MHz** | `[ab initio]`, computed in that paper, CCSD(T)/Dyall CVTZ at R = 3.41 a₀ |
| eQq₂ | **+110 MHz** | **+125 MHz** | `[estimate]` from Eq. (24), see below |

Verbatim from the paper: "The ratios for the magnetic dipole and electric
quadrupole hyperfine structure constants correspond to the ratios for the nuclear
g-factors and the quadrupole moments of the ¹⁷⁷Hf and ¹⁷⁹Hf nuclei." — i.e. the
A∥ ∝ g_N, eQq ∝ Q scaling of §2.2 is stated explicitly by the source.

Their eQq₂ is *not* a direct calculation: their Eq. (23) shows eQq₂ has no matrix
element within the nonrelativistic ³Δ term, so it comes entirely from the
spin–orbit admixture of a ³Π state, and they estimate

```
eQq₂ = 483 w Q ⟨1/r³⟩_5d  MHz          (their Eq. 24)
G∥ ≈ 2 − 2.002319 + w                   (their Eq. 25)
```

with ⟨1/r³⟩_5d = 4.86 a.u. (Dirac–Hartree–Fock, Hf⁺) and w = 0.014 fixed from the
measured G∥(HfF⁺) = 0.011768.

**Why eQq₂ matters, in the source's own words** (their §V): "the MQM shift
strongly depends on eQq₂ and decreases as eQq₂ increases. … the electric
quadrupole hyperfine interaction causes the sublevels Ω = +1 and Ω = −1, which
have different signs for T,P-odd shifts, to mix." For HfF⁺ the Ω-doubling is
sub-MHz and eQq₂ is 110 MHz, so the quadrupole, not the Ω-doubling operator,
sets the Ω-mixing at low field. §4.4 carries this over to ²²⁹ThF⁺.

### 2.4 What each source does *not* contain

- Skripnikov & Titov 2015 (`docs/lit/skripnikov2015-…txt`): no quadrupole, no
  EFG, no spin–rotation. Explicitly: "In the present paper we do not consider
  fluorine nuclear spin."
- Denis et al. 2015 (`docs/lit/denis2015-…txt`): no quadrupole, no EFG, no
  spin–rotation, no ¹⁹F spin. Grepped for `gradient|quadrupole` — the only hits
  are the words in reference titles.
- Flambaum & Mansour 2022 (arXiv:2204.01290) lists ThF⁺ among MQM candidates and
  tabulates the ²²⁹Th MQM in terms of nuclear constants (their Table 1, using the
  stale μ = 0.46(4)), but publishes no molecular hyperfine or quadrupole
  constants for ThF⁺.
- Flambaum 2018 (arXiv:1808.03629) names ²²⁹ThF⁺ as a Schiff-moment candidate;
  no hyperfine constants.

---

## 3. Operator forms with two nuclear spins

### 3.1 Coupling scheme

**Use F₁ = J + I_Th, F = F₁ + I_F, basis |((J I_Th) F₁, I_F) F, m_F⟩** with J, Ω
as in [HAM] §1.1.

Justification, two independent ones:
- Magnitude: |A∥(Th)| ≈ 1.5 GHz ≫ |A∥(F)| = 20.1 MHz (§2.2, [HAM] §2.4), so F₁ is
  a near-exact quantum number and F is a small splitting on top of it — a factor
  ~75. `[derived]`
- Precedent: Petrov 2018 labels every HfF⁺ result by F₁ = J + I_Hf and states the
  T,P-odd shifts are "independent of the F and m_F quantum numbers" for a given
  F₁ (their §II). Their computational basis is the *decoupled* one,
  `|Ω⟩ θ^J_{M,Ω} U^Hf_{I₁M₁} U^F_{I₂M₂}` (their Eq. 15), with F₁ used only for
  labelling — a legitimate alternative for `heff` if m_F blocking is preferred
  over F₁ blocking.

For ²²⁷ThF⁺ (I_Th = 1/2) the two spins are comparable in the *nuclear* sense but
not in the *coupling* sense, and the same scheme applies; there F₁ = J ± 1/2.

### 3.2 The Th operators: everything reduces to the one-spin formulas with F → F₁

**Key structural result.** H_hf(Th) and H_Q(Th) are scalars built entirely from
the rotational degrees of freedom and I_Th — i.e. from the *inner* pair of the
coupled scheme. By **B&C Eq. (5.176)** (PDF p. 205 / book p. 173, verbatim:
"If the scalar product is formed from spherical tensor operators which both act
on the same inner part of a coupled scheme … = δ_{j₁₂j₁₂′} δ_{m₁₂m₁₂′}
δ_{j₂j₂′} ⟨j₁|T^k(A₁)·T^k(B₁)|j₁′⟩"), their matrix elements in
|((J I_Th)F₁, I_F)F⟩ are **exactly** the single-spin matrix elements in
|J, Ω, I_Th, F₁⟩, diagonal in F, m_F and I_F and independent of them.

So, with the substitution **I → I_Th, F → F₁** and no other change:

**(a-Th) Axial magnetic hyperfine, ΔJ = 0** — B&C Eq. (9.50), PDF p. 636 /
book p. 604 (the same equation `heff/elements_c.py::hyperfine_A_par` already
cites):

```
⟨η, J, Ω, I_Th, F₁| H_hf^Th |η, J, Ω, I_Th, F₁⟩
   = A∥^Th [F₁(F₁+1) − J(J+1) − I_Th(I_Th+1)] / [2 J(J+1)]
```
with `A∥^Th = {aΛ + (b_F + (2/3)c)Σ} Ω` evaluated for the Th nucleus.
Selection rules: ΔJ = ΔΩ = ΔF₁ = ΔF = Δm_F = 0. Parity-even.

**(b-Th) Axial magnetic hyperfine, ΔJ = ±1** — B&C Eq. (9.51), same pages,
same substitution (cited by `hyperfine_A_par_dJ1`):

```
⟨η, J, Ω, I_Th, F₁| H_hf^Th |η, J−1, Ω, I_Th, F₁⟩
   = −{aΛ + (b_F+(2/3)c)Σ}
     × (J² − Ω²)^{1/2} {(F₁−I_Th+J)(F₁+I_Th+J+1)(J+I_Th−F₁)(F₁−J+I_Th+1)}^{1/2}
       / [2J(4J²−1)^{1/2}]
```
(J is the larger of the two, brace = A∥^Th/Ω, exactly as in the existing code.)
ΔF₁ = 0, ΔF = 0.

**(c-Th) Electric quadrupole** — B&C Eq. (9.52) general q, PDF p. 636 / book
p. 604, and Eq. (9.53) for q = 0, PDF p. 637 / book p. 605:

```
⟨η, Λ; S, Σ; J, Ω, I, F | H_Q | η′, Λ′; S, Σ; J′, Ω′, I, F ⟩
  = −½ eQ Σ_q (−1)^{J′+I+F+J−Ω} {(2J+1)(2J′+1)}^{1/2}
      × {J′ I F; I J 2} (J 2 J′; −Ω q Ω′) (I 2 I; −I 0 I)^{−1}
      × ⟨η, Λ| T²_q(∇E) |η′, Λ′⟩                                        (9.52)

q = 0 piece:
  = (eq₀Q/4) (−1)^{J′+I+F+J−Ω} {(2J+1)(2J′+1)}^{1/2}
      × {J′ I F; I J 2} (J 2 J′; −Ω 0 Ω) (I 2 I; −I 0 I)^{−1}           (9.53)
```
with I → I_Th, F → F₁. **B&C's stated sign convention** (same page, verbatim):
"q₀ is the negative of the electric field gradient". This is a live trap — see
the verification in §3.5.

Selection rules from the 3j (J 2 J′; −Ω q Ω′): **ΔJ = 0, ±1, ±2** and **Δ(−Ω) +
q + Ω′ = 0, i.e. q = Ω − Ω′**. In the Ω = ±1 block of ³Δ₁ this gives two
physically distinct terms:
- **q = 0**, ΔΩ = 0: the ordinary axial quadrupole, constant eq₀Q. Petrov 2018
  Eq. (22), `eQq₀ = 2eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₀(θ₁ᵢ,φ₁ᵢ)/r₁ᵢ³ |³Δ₁⟩`.
- **q = ∓2**, Ω = ±1 ↔ Ω = ∓1: a ΔΩ = 2 term, constant eq₂Q. Petrov 2018
  Eq. (23), `eQq₂ = 2√6 eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₂(θ₁ᵢ,φ₁ᵢ)/r₁ᵢ³ |³Δ₋₁⟩`.
  This is a **new kind of term for `heff`**: parity-even, ΔΩ = ±2, diagonal in
  J, F₁, F, m_F. It sits in the same matrix position as the Ω-doubling operator.

> **UNVERIFIED — normalisation bridge.** B&C define eq₀Q by the q = 0
> specialisation of (9.52) (comparing (9.52) and (9.53) gives
> eq₀Q = −2eQ⟨η,Λ|T²₀(∇E)|η,Λ⟩ `[derived]`). Petrov 2018 defines eQq₀ and eQq₂
> by their Eqs. (22)–(23). I have **not** verified that Petrov's eQq₂ maps onto
> B&C's q = ±2 element with the same 1/4 prefactor as eq₀Q does for q = 0 — the
> √6 and the Y₂₂ vs T²₂ normalisation must be checked before a published eQq₂ is
> substituted into (9.52). Every eQq₂-derived number in §4.4 inherits this
> caveat and is order-of-magnitude only.

**(d-Th) Nuclear spin–rotation** — B&C Eq. (8.7), PDF p. 410 / book p. 378,
`H_nsr = c_I T¹(J)·T¹(I)`, coupled-basis element Eq. (8.20), PDF p. 414 /
book p. 382 (both already cited by `spin_rotation_cI`). Both operators sit inside
F₁, so B&C (5.176) applies again:

```
⟨…F₁| c_I^Th T¹(J)·T¹(I_Th) |…F₁⟩
   = c_I^Th [F₁(F₁+1) − I_Th(I_Th+1) − J(J+1)] / 2
```
ΔJ = ΔF₁ = ΔF = 0.

### 3.3 The ¹⁹F operators: no longer diagonal in F₁

These are the ones that change shape. T¹(I_F) acts on the *outer* spin, the
rotational factor acts inside F₁. Two cited steps:

1. **B&C Eq. (5.173)** (PDF p. 205), scalar product across the two parts of a
   coupled pair, with j₁ = F₁, j₂ = I_F, j₁₂ = F, k = 1:
   ```
   ⟨F₁, I_F, F| T¹(A₁)·T¹(I_F) |F₁′, I_F, F⟩
      = (−1)^{F₁′+F+I_F} δ_{FF′} {I_F F₁′ F; F₁ I_F 1}
        ⟨F₁‖T¹(A₁)‖F₁′⟩ ⟨I_F‖T¹(I_F)‖I_F⟩
   ```
2. **B&C Eq. (5.174)** (PDF p. 205), operator on the first part with I_Th a
   spectator, j₁ = J, j₂ = I_Th, j₁₂ = F₁:
   ```
   ⟨J, I_Th, F₁‖T¹(A₁)‖J′, I_Th, F₁′⟩
      = (−1)^{F₁′+J+1+I_Th} [(2F₁+1)(2F₁′+1)]^{1/2} {J′ F₁′ I_Th; F₁ J 1}
        ⟨J‖T¹(A₁)‖J′⟩
   ```
   with ⟨J‖T¹(J)‖J⟩ = [J(J+1)(2J+1)]^{1/2} — B&C Eq. (5.179), PDF p. 205.

Applied to `H_nsr^F = c_I^F T¹(J)·T¹(I_F)` this is the complete matrix element.
Applied to the ¹⁹F axial hyperfine, the ΔJ = 0 block of B&C (9.50) is
`A∥^F T¹(J)·T¹(I_F)/[J(J+1)]`, so the same two steps give it. Selection rules for
both: **ΔF₁ = 0, ±1**, ΔF = 0, Δm_F = 0, ΔΩ = 0, ΔJ = 0 (for the projected form).

The ΔJ = ±1 ¹⁹F term (B&C 9.51) needs the same treatment with the ΔJ ≠ 0 reduced
element; that is a mechanical extension of step 2 with J′ ≠ J, and I have **not**
evaluated it — mark as derivable-not-derived, and note it is a ~2 kHz term
([HAM] §2.5) so it is the last thing to add.

### 3.4 Reduction check against the existing ¹⁹F code — PASSES

Requirement from the brief: with I_Th = 0 the two-spin forms must reproduce
`heff/elements_c.py::hyperfine_A_par`. I evaluated §3.3 steps 1–2 numerically
(own 3j/6j implementation) and formed the F = F₁+½ minus F = F₁−½ splitting in
units of A∥^F:

| J | I_Th = 0 result | [HAM] §2.4 closed form (2J+1)/[2J(J+1)] |
|---|---|---|
| 1 | +0.75000 | +0.75000 |
| 2 | +0.41667 | +0.41667 |
| 3 | +0.29167 | +0.29167 |
| 4 | +0.22500 | +0.22500 |

`[derived]` — exact agreement at J = 1–4, so the recoupling in §3.3 reduces
correctly. The check can fail: a wrong phase in (5.173) or (5.174), or the wrong
6j column order, changes these numbers or their sign.

A second, independent check of the quadrupole implementation of B&C (9.53):
setting Ω = 0 and J′ = J, (9.53) must reproduce the textbook Casimir function
`[¾C(C+1) − I(I+1)J(J+1)] / [2I(2I−1)(2J−1)(2J+3)]`, C = F(F+1)−I(I+1)−J(J+1).
It does, for (J,I) = (1,1), (2,1), (2,3/2), (3,5/2), all F — **with a uniform
ratio of exactly −1** `[derived]`. That factor of −1 *is* B&C's stated convention
that q₀ is the negative of the field gradient. **Code it explicitly or the sign
of every quadrupole splitting flips.**

Cross-check of the same angular structure from a third source: Skripnikov,
Petrov, Titov & Flambaum, arXiv:1408.5368 ("Manifestations of nuclear
CP-violation in ThO"), Eq. (5), give the MQM shift in ³Δ₁ as
`δ(J,F) = (−1)^{Ω+I+F+1} C(J,F) W_M M` with
`C(J,F) = [(2J+1)/2] (J 2 J; −Ω 0 Ω)/(I 2 I; −I 0 I) {J I F; I J 2}` —
the identical rank-2 case-(c) skeleton as B&C (9.53), differing only by the
overall constant and by an overall sign consistent with B&C's q₀ convention
`[derived]`. Three independent sources, one angular structure.

### 3.5 Summary of new terms `heff` would need

| term | operator source | selection rules | new machinery? |
|---|---|---|---|
| `hyperfine_A_par_Th` | B&C 9.50, I→I_Th, F→F₁ | ΔJ=ΔΩ=ΔF₁=0 | no — existing kernel, new index |
| `hyperfine_A_par_Th_dJ1` | B&C 9.51, I→I_Th, F→F₁ | ΔJ=±1, ΔF₁=0 | no |
| `quadrupole_eQq0_Th` | B&C 9.53 | ΔJ=0,±1,±2; ΔΩ=0; ΔF₁=0 | **yes** — first rank-2 nuclear term |
| `quadrupole_eQq2_Th` | B&C 9.52 at q=±2 | ΔJ=0,±1,±2; **ΔΩ=±2**; ΔF₁=0 | **yes** — first ΔΩ=2 term |
| `spin_rotation_cI_Th` | B&C 8.7/8.20, I→I_Th, F→F₁ | ΔJ=ΔF₁=0 | no |
| `hyperfine_A_par` (¹⁹F) | B&C 9.50 + 5.173 + 5.174 | ΔJ=0, **ΔF₁=0,±1**, ΔF=0 | **yes** — F₁ becomes off-diagonal |
| `spin_rotation_cI` (¹⁹F) | B&C 8.7 + 5.173 + 5.174 | ΔJ=0, ΔF₁=0,±1, ΔF=0 | yes, same kernel |
| Stark, Zeeman, Ω-doubling, PT-odd | unchanged operators | + spectator recoupling through B&C 5.174/5.175 for the extra spin | mechanical |

---

## 4. Scale estimates

All using **A∥(²²⁹Th) = −1524 MHz** (Skripnikov rescaled to μ = 0.366 μ_N; use
+1491 if Denis's sign is right — the magnitudes below are sign-independent),
B₀ = 7.2743325 GHz, Ω = 1, I_Th = 5/2, I_F = 1/2. Compare against ¹⁹F
A∥ = −20.1 MHz and ω_ef = 5.29 MHz.

### 4.1 Th magnetic hyperfine, diagonal (²²⁹ThF⁺)

Manifold pattern in units of A∥^Th, from B&C 9.50 `[derived]`:

| J | F₁ values | coefficients | total spread × A∥ |
|---|---|---|---|
| 1 | 3/2, 5/2, 7/2 | −1.7500, −0.5000, +1.2500 | **4572 MHz** |
| 2 | 1/2 … 9/2 | −1.1667, −0.9167, −0.5000, +0.0833, +0.8333 | **3048 MHz** |
| 3 | 1/2 … 11/2 | −0.8333, −0.7083, −0.5000, −0.2083, +0.1667, +0.6250 | **2223 MHz** |
| 4 | 3/2 … 13/2 | −0.6250, −0.5000, −0.3250, −0.1000, +0.1750, +0.5000 | **1715 MHz** |

Largest adjacent F₁ gap: 2667 MHz (J = 1, F₁ = 5/2 ↔ 7/2). Smallest: 190 MHz
(J = 3 and J = 4, lowest pair). **The Th hyperfine manifold at J = 1 spans 63 %
of B₀ and 16 % of the J = 1→2 interval (4B = 29.1 GHz).** Every ¹⁹F structure
(20 MHz) and the whole Ω-doubling (5–53 MHz) live inside one F₁ level.

For ²²⁷ThF⁺ (I = 1/2) the same formula gives the ¹⁹F-shaped pattern: splitting
`A∥^Th(2J+1)/[2J(J+1)]`, i.e. 0.75, 0.4167, 0.2917, 0.2250 × A∥^Th for
J = 1, 2, 3, 4 — magnitude unknown (gap G1).

### 4.2 Th magnetic hyperfine, ΔJ = ±1 — promoted from kHz to 100 MHz

B&C 9.51 elements and their second-order shift |ME|²/(2BJ_upper) `[derived]`:

| J ↔ J′ | F₁ | matrix element | 2nd-order shift |
|---|---|---|---|
| 1 ↔ 2 | 3/2 | 1352 MHz | 63 MHz |
| 1 ↔ 2 | 5/2 | 1928 MHz | 128 MHz |
| 1 ↔ 2 | 7/2 | 1980 MHz | **135 MHz** |
| 2 ↔ 3 | 7/2 | 2172 MHz | 108 MHz |
| 3 ↔ 4 | 9/2 | 2231 MHz | 86 MHz |

Compare the ²³²ThF⁺ case, where the same term gives 1.6–2.6 kHz ([HAM] §2.5).
**Five orders of magnitude larger.** Three consequences:
- The ΔJ = ±1 hyperfine is not optional and is not perturbative in the sense
  [HAM] §2.5 uses; it must be in the matrix, and the J = 1–4 truncation must be
  re-tested (the J = 4 ↔ 5 element is of the same 2 GHz size, so a J = 1–5 or
  1–6 basis is likely needed for kHz-level J = 4 energies). `[inferred]`
- The 9/5 ratio between the J = 1 and J = 2 hyperfine splittings that JILA uses
  to fit A∥ ([HAM] §2.4) is broken at the several-percent level for Th.
- The "second-order shift is linear in J(J+1), so it just renormalises B"
  argument of [HAM] §1.2 does **not** apply here: the shift depends on F₁, so it
  is not absorbable into B.

### 4.3 Th electric quadrupole, ΔΩ = 0 (²²⁹ThF⁺ only)

**No ThF⁺ eQq₀ exists** (§2.1). Estimate by the HfF⁺ anchor:

```
eQq₀(²²⁹ThF⁺) ≈ eQq₀(¹⁷⁷HfF⁺) × [Q(²²⁹Th)/Q(¹⁷⁷Hf)] × R_el
             = −2100 MHz × (3.11/3.365) × R_el = −1.94 GHz × R_el
```
`[estimate, method: HfF⁺ anchor × Q ratio × electronic-density ratio]`, with
R_el the ratio of the valence-δ EFG at Th to that at Hf. Proxy for R_el: the same
ratio for the magnetic hyperfine electronic factor,
|A∥/g_N|(ThF⁺)/|A∥/g_N|(HfF⁺) = 10 408 / (1429/0.2267) = 10 408/6303 = **1.65**
`[derived]`. The two operators weight r differently (1/r² vs 1/r³), so treat this
as bracketing rather than a derivation:

**eQq₀(²²⁹ThF⁺) ≈ −2 to −3.3 GHz, take −2.6(1.0) GHz** `[estimate]`.

Level pattern in units of eq₀Q, B&C 9.53 diagonal `[derived]`:

| J | coefficients over F₁ | spread | × (−2.6 GHz) |
|---|---|---|---|
| 1 | +0.0700, −0.0800, +0.0250 | 0.1500 | **390 MHz** |
| 2 | −0.1000, −0.0357, +0.0357, +0.0607, −0.0357 | 0.1607 | 418 MHz |
| 3 | −0.1500, −0.0825, +0.0050, +0.0750, +0.0750, −0.0625 | 0.2250 | 585 MHz |
| 4 | −0.1518, −0.0607, +0.0320, +0.0916, +0.0734, −0.0773 | 0.2434 | **633 MHz** |

So the quadrupole is a **few hundred MHz** effect — 15–30 % of the Th magnetic
hyperfine, and it does not fall with J the way the magnetic term does. It also
has non-zero ΔJ = 1 and ΔJ = 2 elements (0.01–0.18 and 0.05–0.14 × eq₀Q
`[derived]`, i.e. 25–470 MHz), reinforcing §4.2's J-truncation warning.

Zero for ²²⁷ThF⁺ (§1.4).

### 4.4 Th electric quadrupole, ΔΩ = ±2 — the term that would out-muscle ω_ef

Geometric factor ⟨J, Ω=−1, F₁| H_Q |J, Ω=+1, F₁⟩ in units of the q = ∓2 constant,
same normalisation as eq₀Q in (9.53) `[derived]`, subject to the UNVERIFIED
normalisation bridge of §3.2:

| J | coefficients over F₁ (ascending F₁) |
|---|---|
| 1 | +0.1715, −0.1960, +0.0612 |
| 2 | +0.2449, +0.0875, −0.0875, −0.1487, +0.0875 |
| 3 | +0.2449, +0.1347, −0.0082, −0.1225, −0.1225, +0.1021 |
| 4 | +0.2187, +0.0875, −0.0461, −0.1320, −0.1058, +0.1113 |

Estimate of the constant, by the HfF⁺ route (Petrov 2018 Eqs. 24–25):
w(ThF⁺) = G∥ + 0.002319 = 0.0476 + 0.0023 = **0.0499**, using Ng's measured
|G∥| = 0.0476 ([HAM] §2.8), against w(HfF⁺) = 0.014. Then
```
eQq₂(²²⁹ThF⁺) ≈ 110 MHz × (0.0499/0.014) × (3.11/3.365) × [⟨1/r³⟩_6d,Th / 4.86 a.u.]
             ≈ 362 MHz × [⟨1/r³⟩_6d,Th / 4.86 a.u.]
```
`[estimate]`. ⟨1/r³⟩ for the Th⁺⁺ 6d orbital is not in any source read (gap G2),
and 6d is more diffuse than Hf⁺ 5d, so the bracket is probably 0.6–1.0. Take
**eQq₂(²²⁹ThF⁺) ~ 200–400 MHz**.

Then the Ω = +1 ↔ Ω = −1 matrix element at J = 1 is
`0.17–0.20 × (200–400 MHz) ≈ 35–80 MHz`, against the Ω-doubling off-diagonal
element `ω_ef J(J+1)/4 = 2.65 MHz` at J = 1 ([HAM] §2.3). **The quadrupole would
be 13–30× the Ω-doubling operator, and would set the parity-doublet structure of
²²⁹ThF⁺ at zero and low field.** This is not a novel claim: Petrov 2018 §V
reports precisely this behaviour in ¹⁷⁷HfF⁺, where eQq₂ = 110 MHz dwarfs a
sub-MHz Ω-doubling and visibly changes the T,P-odd shifts.

Caveats, both real: (i) the normalisation bridge of §3.2 is UNVERIFIED, so the
ratio could be off by a small factor; (ii) at the JILA operating field
E_rot = 60 V/cm the Stark polarisation energy is 50.9 MHz per level ([HAM] §2.7),
comparable to the eQq₂ element, so the two compete rather than one simply
winning. Both make this a **first-class reason to get a real ab initio eQq₂ for
ThF⁺** before designing a ²²⁹ThF⁺ measurement.

### 4.5 ¹⁹F hyperfine, projected onto F₁

Once I_Th ≠ 0 the ¹⁹F doublet splitting inside each F₁ is no longer
A∥^F(2J+1)/[2J(J+1)]. From §3.3 `[derived]`, in units of A∥^F, at J = 1:

| F₁ | splitting E(F=F₁+½) − E(F=F₁−½) |
|---|---|
| 3/2 | **−0.4000** (sign inverted relative to I_Th = 0) |
| 5/2 | +0.1714 |
| 7/2 | +0.5714 |

(I_Th = 0 reference: +0.7500.) At A∥^F = −20.1 MHz these are +8.0, −3.4 and
−11.5 MHz. So **the ¹⁹F doublet ordering flips between F₁ manifolds** — a clean,
falsifiable prediction for a ²²⁹ThF⁺ spectrum and a good acceptance test for the
recoupling code.

Off-diagonal ΔF₁ = ±1 ¹⁹F elements are 0.13–0.37 × A∥^F/[J(J+1)] `[derived]`,
i.e. 1.3–3.8 MHz, against F₁ spacings of 190–2670 MHz (§4.1) — second-order
shifts ≲ 25 kHz. **F₁ is a good quantum number**, as §3.1 assumed.

### 4.6 Nuclear spin–rotation

- ¹⁹F: c_I ≈ 20 kHz with a factor-3 uncertainty ([HAM] §2.6, OPEN-6). Unchanged
  in magnitude; now off-diagonal in F₁ (§3.3).
- Th: **no value anywhere** (gap G3). The only bracketing I will commit to:
  c_I is roughly proportional to g_N × (an electronic rotational-field factor).
  g_N(²²⁹Th)/g_N(¹⁹F) = 0.1464/5.2566 = 0.0279 `[derived]`, so if the electronic
  factors were equal, c_I(Th) ≈ 0.6 kHz. They are not equal — the unpaired
  σδ density sits on Th, and the analogous hyperfine electronic factor is
  2700× larger there (10 408 MHz vs 20.1/5.2566 = 3.82 MHz per unit g_N)
  `[derived]`. Those two factors point in opposite directions by ~10²–10³, so
  the honest statement is **c_I(Th) is unconstrained between ~1 kHz and ~1 MHz**;
  I decline to pick a value. At the top of that range it would be a
  c_I[F₁(F₁+1)−I(I+1)−J(J+1)]/2 pattern of up to ~10 MHz at J = 4 — i.e.
  potentially larger than the whole ¹⁹F hyperfine. This deserves a real
  calculation before it is assumed negligible.

### 4.7 Ordering summary (²²⁹ThF⁺, J = 1)

| term | size | vs [HAM] ²³²ThF⁺ |
|---|---|---|
| B₀ J(J+1) | 14.5 GHz | same |
| **Th magnetic hyperfine spread** | **4.6 GHz** | new |
| **Th quadrupole eQq₀ spread** | **~0.4 GHz** `[estimate]` | new |
| **Th ΔJ=±1 hyperfine, 2nd order** | **60–135 MHz** | 2.6 kHz → 10⁵× |
| **Th quadrupole eQq₂, Ω-mixing** | **~35–80 MHz** `[estimate]` | new |
| Stark at 60 V/cm | 50.9 MHz | same |
| ¹⁹F hyperfine (F₁-projected) | 3–12 MHz | was 15.1 MHz |
| Ω-doubling ω_ef | 5.29 MHz | same |
| Th spin–rotation | 1 kHz – 1 MHz?? | unconstrained |
| ¹⁹F spin–rotation | ~30 kHz | same |
| ¹⁹F ΔJ=±1 hyperfine | ~2.6 kHz | same |

---

## 5. Gaps within the searched sources

**G1. μ(²²⁷Th) and the ²²⁷Th ground-state spin.**
The earlier search missed a published nuclear-theory calculation. [Minkov et al., Phys. Rev. C 110, 034327 (2024), Table IV](https://arxiv.org/abs/2408.11010) predicts μ(²²⁷Th) = −0.0860 μ_N for its octupole-deformed 1/2 ground solution. With I = 1/2 and the existing molecular factor −10408 MHz, this gives g_N = −0.1720 and A∥ = +1790.176 MHz. This is a model prediction without a calibrated uncertainty, not a measurement.
The 2026-09-05 compilation/API queries are preserved in
[the lookup record](lit/lookup-227th-nuclear-moment.md); no experimental moment
was located there. The spin assignment (1/2⁺) remains tentative. Nuclear-model
uncertainty, rather than the existence of any estimate, is the unresolved input.

**G2. eQq₀ and eQq₂ for ThF⁺ (or ThO), and the EFG / ⟨1/r³⟩ at Th.**
No published value, for any Th isotope or state.
Queries run across Skripnikov & Titov 2015 and Denis et al. 2015 for
`gradient|quadrupole|eQq` (only reference-title hits); arXiv API
`abs:"ThF" AND abs:"quadrupole"` (5 hits, all MQM papers, none
with an electric-quadrupole coupling constant); `all:"electric field gradient"
AND all:"ThO"` → 0; `all:"quadrupole coupling constant" AND all:"thorium" AND
all:"molecule"` → 0; DuckDuckGo `"ThF+" "eQq" thorium quadrupole coupling
constant` → "No results found".
The §4.3/§4.4 numbers are HfF⁺-anchored estimates, nothing more. This is a
one-paper-sized hole and the obvious ask of Skripnikov/Petrov or Fleig.

**G3. c_I for Th in ThF⁺ (and, still, for ¹⁹F).**
No value for either nucleus in ThF⁺.
Queries run: `all:"thorium" AND all:"nuclear spin-rotation"` → 0 results; a
full-text search of the reviewed source set for `spin-rotation` returns only the ¹⁹F estimate
already recorded in [HAM] §2.6. The ¹⁹F value remains the CsF-anchored 20 kHz
estimate (OPEN-6). The Th value is unconstrained over three decades (§4.6).

**G4. The sign of A∥(Th).**
Skripnikov & Titov 2015 print −4163, Denis et al. 2015 print +1833, from the
same defining equation, agreeing to 2.2 % in magnitude (§2.2). The later
[sign-convention audit](lit/lookup-apar-th-sign-convention.md) found that both
papers state their axes elsewhere, but a consistent axis reversal leaves A∥
invariant. The audit recommends the negative branch; it does not resolve the
underlying disagreement between the electronic-structure calculations.

**G5. Any measured Th hyperfine constant in any ThF⁺ isotopologue.**
None exists. All JILA spectroscopy (Gresh 2016, Ng 2022, Zhou 2020, Roussy 2023)
is ²³²ThF⁺, I(Th) = 0.
Query run: arXiv API `all:"ThF+"` sorted by date, 60 results scanned — the only
ThF⁺ molecular-physics papers are Barker 2019/2012, Gresh 2015, Ng 2022,
Skripnikov 2015, Petrov 2025, the 2025 REMPAD detection paper, and solid-state
²²⁹ThF₄ clock films. No ²²⁷/²²⁹ThF⁺ spectroscopy.

**G6. A ²²⁹ThF⁺-specific structure/measurement proposal.**
Flambaum 2018 (arXiv:1808.03629) and Flambaum & Mansour 2022 (arXiv:2204.01290)
both name ²²⁹ThF⁺ as a Schiff-moment / MQM candidate, but neither computes a
molecular hyperfine or quadrupole constant for it. The nearest complete treatment
is Petrov 2018 on ¹⁷⁷/¹⁷⁹HfF⁺ (§2.3), which is the template this extension
should follow.
Query run: arXiv API `all:"229ThF"` → 1 result (Flambaum 2018); `au:Petrov AND
all:"ThF"` → 1 result (the 2025 g-factor paper, ²³²ThF⁺ only).

**Coverage limit of the 2026-09-05 search:** Stone's INDC(NDS)-0794 compilation
was not part of this digest's original source set, and the JILA theses (Ng,
Stutz, Gresh, Cairncross, Grau, and Loh) were searched selectively rather than
read in full. The later [²²⁷Th lookup](lit/lookup-227th-nuclear-moment.md)
records the targeted compilation search.
