# `heff`'s two-photon notation vs the Pipi papers

Assembled 2026-09-17 for `heff`. Task: check the effective two-photon (2 × E1)
notation in `heff/twophoton.py` and `docs/thf-plus-x3delta1-effective-hamiltonian.md`
§9.5 against Anastasia Pipi's two arXiv papers and their appendices.

**Addendum 2026-09-18:** the package default flipped to `reading='raman'`
(commit 24ccf28); the caller instruction below is now the default behaviour.

**Reading rules for this document** (same as `docs/digest-literature-two-photon.md`).
Every factual claim about a paper carries an inline cite to an equation and a
PDF page I read at source in this session, in true column order. Anything I
derived is tagged `[derived]`. Anything I could not verify is tagged
**UNVERIFIED** and is *not* used to support a downstream conclusion. Where a
convention is simply absent from both papers the row says **NOT IN PAPER** and
§5 records the query whose zero I am reporting. I have not guessed a single
matrix element, phase, or polarization convention.

**Sources read.**

- **[P24]** A. Pipi, X. Tao, A. Wu, P. Narang, D. R. Leibrandt, *Molecular
  Quantum Control Algorithm Design by Reinforcement Learning*,
  [arXiv:2410.11839v5](https://arxiv.org/abs/2410.11839), 41 PDF pages. Main
  text pp. 1-9; Supplementary Material from p. 10 (Sec. SA p. 10, SB p. 11,
  SC p. 12, SD p. 13, SE figures and tables p. 17 onward). Read in full, with
  pp. 6, 11, 12, 18, 28, 39 read in true reading order.
- **[P26]** A. Pipi, V. Duruisseaux, E. Been, X. Tao, T. L. Patti,
  A. Anandkumar, P. Narang, *Inverse Design of Quantum Control Sequences with
  Fourier Neural Operators*, [arXiv:2608.03702v2](https://arxiv.org/abs/2608.03702),
  19 PDF pages. Main text pp. 1-13; **Appendix A: Hamiltonian** p. 16
  (A 1 Hyperfine Hamiltonian, Eqs. A1-A10; A 2 Block Diagonal form p. 17,
  Eqs. A11-A13); Appendix B p. 17; Appendix C p. 18. Read in full, with
  pp. 2, 4, 16, 17 read in true reading order.
- `heff` working tree at `aca8f32`: `heff/twophoton.py`,
  `docs/thf-plus-x3delta1-effective-hamiltonian.md` §9.5 (**[HAM]**),
  `docs/lit/2026-09-15-twophoton-operator-audit.md` (**[AUDIT]**),
  `docs/digest-literature-two-photon.md` §3 (**[2γ]**), `tests/test_twophoton.py`.

Layout text extracts sit **untracked** in `docs/lit/` as
`pipi2024-arXiv2410.11839-molecular-control.txt` and
`pipi2026-arXiv2608.03702-fno-control.txt`. They are not committed: `.gitignore`
line 17 is `docs/lit/**/*.txt`, and the repo's existing Petrov extracts are
untracked for the same reason (`git ls-files docs/lit/` lists no `.txt`). To
re-run the §5 queries, regenerate them with
`curl -sL https://arxiv.org/pdf/2410.11839v5 -o pipi2024.pdf` (same for
`2608.03702v2`) then `pdftotext -layout`. Both papers are two-column and
`pdftotext -layout` interleaves the
columns, so **every equation quoted below was re-read with
`pdftotext -f N -l N <pdf> -`** before being cited.

---

## 0. The headline, before the detail

1. **Neither paper is a notation source for `heff`.** Both compute the
   two-photon coupling as an *explicit numerical sum over intermediate states*
   of products of dipole matrix elements. Neither paper contains a
   polarizability tensor, a rank decomposition, a Wigner-Eckart or
   Clebsch-Gordan statement, a Condon-Shortley phase statement, a
   polarization-vector component definition, or an absolute field amplitude
   (§5 carries the queries). So the rows of §2 that ask "which phase
   convention" are NOT IN PAPER, not DIFFERS — there is nothing to disagree
   with.
2. **The one equation that maps cleanly onto `heff` is [P26] Eq. (A5) = [P24]
   Eq. (S10)**, and its *first* term matches `heff`'s `reading='raman'`
   structure exactly, including the sign of the intermediate-state detuning.
3. **The two papers state opposite `σ±` ↔ `Δm_F` maps for the same physical
   pair** (π absorbed, σ± emitted): [P26] Eq. (A10) p. 16 says `σ+ → Δm_F = +1`,
   while [P24] Fig. S11 caption p. 39 says `σ−`(emit.) → `Δm = +1`. `heff`
   (Raman reading) agrees with [P24] and contradicts [P26]. Neither paper
   defines `σ±` in components or states the beam geometry, so the conflict is
   **not resolvable from the papers**. **STOP item 1.**
4. **[P26] Eq. (A5) retains a rank-1 (antisymmetric) two-photon part that
   `heff` excludes by construction.** Its two time orderings carry *unequal*
   denominators, and [P24] p. 12 explicitly declines the rotating-wave
   approximation precisely because those denominators are comparable in
   magnitude. The exchange-antisymmetric part of such a sum is nonzero
   `[derived]`, and its geometric weight for a (π, σ±) pair equals that of the
   rank-2 part `heff` does keep (measured in §4). This is the same physics as
   the existing **OPEN-21**, now with an external published example.
   **STOP item 2.**
5. **No absolute-scale cross-calibration is possible.** Every Rabi rate in
   both papers is normalized to a reference transition ([P26] Eqs. (A8)/(A9)
   p. 16; [P24] p. 12), never to a stated field amplitude. `heff`'s
   `alpha_K2_dOm0` and `alpha_K2_dOm2` cannot be pinned from these papers.

---

## 1. Where the two-photon coupling is defined

### 1.1 [P26] arXiv:2608.03702v2 — Appendix A 1, PDF p. 16

This is the defining block. Quoted verbatim from PDF p. 16 (reading order):

> The allowed Raman couplings are restricted by the electric-dipole-induced
> two-photon selection rules
>
> ```
> ΔJ = 0, ±2,      ΔK = 0,      p_J = p_J′ .                          (A4)
> ```
>
> In this molecule we can drive transitions using laser fields that couple the
> different molecular states in the subspace. For each allowed transition
> `|J⟩ → |J′⟩`, the effective two-photon Raman Rabi rate is computed by summing
> over the intermediate excited states `|M⟩`,
>
> ```
>                1      ⎡  A_{J M J′}         B_{J M J′}   ⎤
> Ω_{J,J′}  =  ----  Σ  ⎢ ------------  +  --------------- ⎥ ,          (A5)
>              4ℏ²    M ⎣  ω_{J M} − ω₁     ω_{J M} + ω₂   ⎦
> ```
>
> ```
> A_{J M J′} = ⟨J′|d · E₂|M⟩⟨M|d · E₁|J⟩ ,
> B_{J M J′} = ⟨J′|d · E₁|M⟩⟨M|d · E₂|J⟩ ,                             (A6)
> ```
>
> where `E₁` and `E₂` are the two laser fields with frequencies `ω₁` and `ω₂`,
> `d` is the molecular dipole operator, and
>
> ```
> ω_{J M} = (E_M − E_J)/ℏ                                              (A7)
> ```
>
> is the transition frequency between the initial molecular state `|J⟩` and the
> intermediate state `|M⟩`. **The two terms account for the two possible time
> orderings of the Raman process.**

and, same page, the normalization and the polarization rule:

> The field amplitudes are normalized such that the reference transition
> `{J,K,p,m_F,ξ}_i = (2, 2, −, 3/2, 1)`, `{…}_f = (2, 2, −, 5/2, 1)`   (A8)
> has Rabi rate `Ω_ref/2π = 2.000 kHz`.                                 (A9)
>
> The two fields driving each pulse have fixed polarizations: `E₁` is
> π-polarized, while `E₂` is chosen to be either `σ+`- or `σ−`-polarized. The
> allowed transitions satisfy
>
> ```
>                                    ⎧ +1,  σ+ ,
> Δm_F = m_{F,J′} − m_{F,J}    =     ⎨                                  (A10)
>                                    ⎩ −1,  σ− .
> ```

Supporting definitions elsewhere in [P26]:

- **Eq. (A1)-(A3), p. 16**: state label `𝒥 ≡ {J, K, p, m_F, ξ}`, with `K` "its
  projection onto the molecule-fixed symmetry axis", `p = ±` the inversion
  parity, and `m_F` "the projection of the total angular momentum onto the
  laboratory-frame quantization axis".
- **Eq. (2), p. 4**: `Ĥ_int = Σ_{J→J′} (ℏΩ_{J,J′}/2)[e^{i[η(a+a†)−ωt]}|J′⟩⟨J| + h.c.]`,
  with `ω = ω₁ − ω₂` "the effective two-photon angular frequency" and `η = 0.09`.
- **Eqs. (3), (4), p. 4**: `ω^(+) = ν_f + ω_{J′J}`, `ω^(−) = ν_f − ω_{J′J}`, and
  the *sideband* detuning `Δ_{J′J}(ω, σ) = ω^(±) − ω`. Fig. 1 caption p. 2 gives
  the same as `Δ_{J′J} = ν_f + ω_{J′J} − ω`.
- **p. 17, Appendix A 2**: "Within the defined MHz frequency range and the
  selection rules of two-photon transitions, the Hamiltonian becomes
  block-diagonal", with the RWA invoked to drop THz-scale couplings.

### 1.2 [P24] arXiv:2410.11839v5 — Sec. SC, PDF p. 12

Quoted verbatim:

> **Raman Rabi rates: CaH+** — For a given set of two pump/Stokes pulses (with
> known amplitudes, polarization, frequencies, and duration), the Raman-Rabi
> frequency is given by
>
> ```
>            1      ⎡ ⟨f|d·E₂|M⟩⟨M|d·E₁|i⟩     ⟨f|d·E₁|M⟩⟨M|d·E₂|i⟩ ⎤
> Ω_if  =  ----  Σ  ⎢ --------------------- + --------------------- ⎥   (S10)
>          4ℏ²    M ⎣     ω_iM − ω₁                ω_iM + ω₂        ⎦
> ```
>
> where `E₁, E₂` are the two driving fields with respective frequencies
> `ω₁, ω₂`, `d` is the dipole operator, `ω_iM = (E_M − E_i)/ℏ` is the frequency
> difference of the initial, `|i⟩`, and intermediate `|M⟩`, states. **The
> absorption pulse produces a π-polarized field and the stimulated emission
> pulse produces a `σ+/σ−`-polarized field.** More details on the above
> expression can be found in [1, 5], and pg. 22-23 of [6]. It is worth
> mentioning that **we do not apply the rotating wave approximation in
> Eq. S10, because two-photon Raman transitions can utilize a pump/Stokes laser
> that is far detuned from the intermediate states (i.e. `|ω_iM − ω₁|` is
> comparable to `|ω_iM + ω₂|`).** … The amplitudes of the laser pulses are set
> the same as in the experiment [1] such that the Rabi rate for transition
> `|1, −3/2, −⟩ → |1, −1/2, −⟩` is `2π × 2.087 kHz`.

Eq. (S10) is Eq. (A5) with `(i, f) → (𝒥, 𝒥′)`. Supporting locations in [P24]:

- **Eq. (S4)/(S5), p. 11**: `H_int(t) = Σ (Ω_{𝒥,𝒥′}/2)[e^{i[λ_LD(a†+a)−ωt)]}|𝒥′⟩⟨𝒥| + h.c.]`,
  "where `Ω_{𝒥,𝒥′}` is the Rabi frequency for the two-photon stimulated Raman
  transition `|𝒥⟩ → |𝒥′⟩` and is obtained by adiabatically eliminating the
  intermediate states [1] (see Eq. S10)".
- **Table S4 caption, p. 28**: "The selection rules for the electric
  dipole-induced Raman transition is `ΔJ = 0, ±2`, `ΔK = 0`,
  `parity_i = parity_f`. Two Raman pulses are π-polarized and `σ+/σ−`-polarized,
  respectively, thus only `Δm_F = ±1` transitions are allowed. The table only
  lists the transitions with `Δm_F = +1` since the transitions i→f and f→i share
  the same Rabi rate. The rates are obtained by assuming laser pulse amplitudes
  such that the Rabi rate of the * transition is `2.000 × 2π kHz`." The starred
  row is `(2,2,−,3/2,1) → (2,2,−,5/2,1)`, i.e. the same reference transition and
  the same 2.000 kHz as [P26] Eqs. (A8)/(A9).
- **Fig. S11 caption, p. 39** (CaH+): "The Raman pulses are with **π (abs.) and
  `σ−` (emit.)** polarizations, thus only one direction of the population
  transfer (namely, **those with `Δm = 1`**) can be driven." The four listed
  pulse pairs confirm it: the driven (primed) member of each pair is the one
  with `Δm = +1`, e.g. `1′: |2, 3/2, −⟩ → |2, 5/2, +⟩`.
- **Table S2 caption, p. 18**: `D = π/(λ_{L−D} Ω)`, the pulse duration relation.
- **Eq. (S9), p. 12** (CaH+ Hamiltonian, a ¹Σ-type ion with `|J, m, ξ⟩` labels)
  and **Eq. (S11), p. 12** (H₃O⁺, `Ĥ_inv-rot + Ĥ_Zeeman + Ĥ_s-r`, `B = 0.36 mT`).
  Both papers defer the H₃O⁺ level and coupling computation to "a subsequent
  article" ([P24] ref. 7, [P26] ref. 48), which is where any spherical-tensor
  or phase convention would live.

### 1.3 The molecules, and what that costs the comparison

[P24] treats **CaH⁺** (a ¹Σ-type ion; Eq. (S9) has only `R Ĵ²`, rotational and
nuclear Zeeman, and `c_IJ Î·Ĵ`, so `Ω = 0`) and **H₃O⁺** (a C₃ᵥ symmetric top
with umbrella inversion doubling, labels `J, K, p, m_F, ξ`). [P26] treats
H₃O⁺ only. **Neither is a case (c) `Ω = ±1` molecule**, so `heff`'s registered
channels — case `c2`, `Ω = ±1`, `q = Ω′ − Ω ∈ {0, ±2}` — cannot be pointed at
either paper's states. What *can* be compared is the lab-frame polarization
geometry and the selection-rule content, which is what §2 and §4 do. What
cannot be compared is any molecule-frame reduced element, any `ΔΩ = ±2`
statement, and any absolute `α` scale.

---

## 2. Correspondence table

`heff` symbols are as in `heff/twophoton.py` and [HAM] §9.5. Verdicts:
**CONFIRMED** (same object, same convention), **DIFFERS** (stated how, and
whether convention or physics), **NOT IN PAPER**.

### (a) The effective two-photon operator and the detuning sign

| | |
|---|---|
| `heff` | `T_eff = Σ_i (d·ε₂*)\|i⟩⟨i\|(d·ε₁)/Δ_i`, `Δ_i = E_i − E_g − ℏω₁` ([HAM] §9.5.1, after Cossel Eqs. (6.29)/(6.34)); then closure `Δ_i → Δ` makes the intermediate sum the reduced element of one rank-K operator, B&C (5.142) read backwards. Amplitude `⟨f\|T_eff\|g⟩`, summed over channels **then** squared (`heff.spectra._strengths_from_matrices`). |
| Pipi | `Ω_{J,J′} = (1/4ℏ²) Σ_M [A/(ω_{JM} − ω₁) + B/(ω_{JM} + ω₂)]`, [P26] Eq. (A5) p. 16 = [P24] Eq. (S10) p. 12, with `A`, `B` from Eq. (A6) and `ω_{JM} = (E_M − E_J)/ℏ` from Eq. (A7). |

**Verdict: CONFIRMED on the first term, DIFFERS on the structure.**
`heff`'s `Δ_i = E_i − E_g − ℏω₁` is `ℏ(ω_{JM} − ω₁)` — **same quantity, same
sign**, with `heff`'s `i` = Pipi's `M` and `heff`'s `g` = Pipi's `J`. So the
`A`-term of Eq. (A5) is `heff`'s operator term for term. The difference is that
Pipi keep the **second time ordering** `B/(ω_{JM} + ω₂)`, which `heff` has no
representation for: `heff` closes the intermediate sum with one common `Δ`, and
a common denominator is exactly what makes the `A` and `B` numerators combine
into their exchange-**symmetric** part alone. Consequences in row (h).

Pipi's *other* detuning, `Δ_{J′J}(ω, σ) = ω^(±) − ω` ([P26] Eq. (4) p. 4), is a
**two-photon sideband** detuning used to drive pulse dynamics. `heff` computes
static line positions and strengths and has no analogue. **NOT IN `heff`**, not
a conflict.

### (b) Polarizability tensor components

| | |
|---|---|
| `heff` | `alpha_K0_dOm0` (scalar), `alpha_K2_dOm0` (`T²₀ ~ 2a_zz − a_xx − a_yy`), `alpha_K2_dOm2` (`T²_{±2} ~ a_xx − a_yy`); `alpha_{K,dΩ} = ⟨η′‖α^K_q‖η⟩ = (1/Δ)⟨η′‖T^K_q(d,d)‖η⟩`, `q = Ω′ − Ω` ([AUDIT] §6). |
| Pipi | nothing. |

**Verdict: NOT IN PAPER.** No polarizability, Raman tensor, AC-Stark tensor, or
`α_∥/α_⊥/α^{(K)}_P` component appears anywhere in either paper (§5 query 1).
Every `α` in both papers is the reinforcement-learning / gradient-descent
control-parameter vector `α = (ω, τ, σ)` ([P26] p. 4, Fig. 2 caption). Pipi use
the resolved dipole-sum form, which is the *other* of the two code shapes
[2γ] §4.2 names; there is no scalar to map `alpha_K2_dOm0` or `alpha_K2_dOm2`
onto, in either direction.

### (c) Spherical-component and Clebsch-Gordan / 3j phase convention

| | |
|---|---|
| `heff` | Condon-Shortley, stated and gated: `_spherical` returns `{+1: −(v_x+iv_y)/√2, 0: v_z, −1: (v_x−iv_y)/√2}` (B&C Table 5.2); `leg_weights` gives `c[p] = (−1)^p ε_{−p}` (B&C 5.111); `_cg` is `(−1)^{j₁−j₂+M}√(2J+1) w3j(j₁,j₂,J,m₁,m₂,−M)` (B&C 5.141). All three CONFIRMED against B&C at source in [AUDIT] §2-§3. |
| Pipi | nothing. |

**Verdict: NOT IN PAPER.** Zero hits for Clebsch, Wigner, Eckart, Condon,
Shortley, "reduced matrix", irreducible, or any 3j/6j/9j symbol in either paper
(§5 query 2). Both papers evaluate `⟨J′|d·E|M⟩` numerically in a diagonalized
hyperfine basis, so whatever convention they use lives in the deferred
"subsequent article" and in unpublished code. **Nothing in `heff` is
contradicted, and nothing is corroborated.**

### (d) Polarization-vector definitions and the helicity sign

| | |
|---|---|
| `heff` | `POLARIZATIONS['sigma+'] = (−1, −i, 0)/√2 = ê_{+1} = −(x̂+iŷ)/√2`, Condon-Shortley; `leg_weights(ê_{+1})` gives `c[+1] = +1` and nothing else, i.e. **absorbing a σ⁺ photon raises `m_F` by one** — [HAM] §9.5.1(3) calls this "the sanity anchor for the whole convention", and [AUDIT] §2 re-measured `c[+1] = 0.9999999999999999`. |
| Pipi | labels only: "π-polarized", "`σ+`/`σ−`-polarized". The only quantitative statement is the `Δm_F` assignment. |

**Verdict: components NOT IN PAPER; the `Δm_F` assignment DIFFERS, and the two
papers DIFFER FROM EACH OTHER.**

Both papers use the same physical pair — `E₁` = π = **absorbed**, `E₂` = σ± =
**emitted** ([P24] p. 12 verbatim; [P26] via `ω = ω₁ − ω₂` on p. 4, so `ω₂` is
the Stokes field). Given that:

| source | statement | page |
|---|---|---|
| [P26] Eq. (A10) | `σ+ → Δm_F = +1`, `σ− → Δm_F = −1` | 16 |
| [P24] Fig. S11 caption | `σ−`(emit.) `→ Δm = +1` (so `σ+ → −1`) | 39 |
| `heff`, `reading='raman'` | `(ε₁ = π, ε₂ = σ+) → Δm_F = −1` (measured, §4) | — |

`heff` agrees with [P24] and contradicts [P26]. [P24]'s assignment is the
physically natural one under the standard helicity convention `[derived]`:
emitting a `σ−` photon removes `−ℏ` of projection from the molecule, so the
molecule's `m_F` rises by one. [P26] Eq. (A10) is the opposite, and is
*internally* self-consistent — its Eq. (3) pairs `σ+` with `ω^(+) = ν_f + ω_{J′J}`
and `σ−` with `ω^(−) = ν_f − ω_{J′J}`, and p. 4 states the reflection relation
"for every `σ+` coupling `|J⟩ → |J′⟩`, the corresponding `σ−` coupling drives
`|J′⟩ → |J⟩`" — so nothing inside [P26] breaks. The label is simply attached the
other way round.

**Not resolvable from the papers.** Neither paper defines `σ±` in Cartesian or
spherical components, and neither states the beam propagation direction relative
to the quantization axis (§5 queries 3 and 4) — a σ± label referred to the beam
`k`-vector inverts relative to the quantization axis when the beam is reversed,
which is the obvious candidate resolution and is unstated. **STOP item 1 for
Arian; a one-sentence question to the authors settles it.**

### (e) Which leg carries which beam, and is the second photon conjugated

| | |
|---|---|
| `heff` | slot `a` (left/bra-side operator factor) carries `ε₂`, slot `b` carries `ε₁`, following B&C (5.141)/(5.142)'s ordering ([HAM] §9.5.1(3), [AUDIT] §3). `reading='ladder'` (**the default**) takes `ε₂` unconjugated, both photons absorbed; `reading='raman'` conjugates `ε₂`. |
| Pipi | `A = ⟨J′\|d·E₂\|M⟩⟨M\|d·E₁\|J⟩` ([P26] Eq. (A6) p. 16): `E₁` acts on the initial state, `E₂` on the bra side. No conjugate is written; the fields are real and the RWA is explicitly declined ([P24] p. 12). |

**Verdict: CONFIRMED on the slot order; the conjugation is NOT IN PAPER
explicitly but is forced.** Pipi's `A` numerator puts `E₂` on the bra side and
`E₁` on the ket side, which is `heff`'s slot assignment exactly. Because Pipi
write a real field with no RWA, `d·E₂` carries **both** Fourier components; which
one contributes to the `A` term is fixed by that term's own denominator
`(ω_{JM} − ω₁)` together with the net two-photon frequency `ω = ω₁ − ω₂`, i.e.
`E₂` must supply the **emission** (conjugated) component `ε₂*` `[derived]`. So
Pipi's `A` term is `heff`'s `reading='raman'`, and

> **any comparison of `heff` to either Pipi paper must pass
> `reading='raman'`** to `dyad_weights` / `two_photon_line_strengths`. The
> package default is `'ladder'`, which inverts the σ± ↔ `Δm_F` map (measured in
> §4). [HAM] §9.5.1(3) already requires every polarization-labelled result to
> name its reading; this is an external case where it matters.

### (f) Field-amplitude factors and the absolute scale of `α`

| | |
|---|---|
| `heff` | **none.** `two_photon_line_strengths` "multiplies geometry by `alpha` and by the polarisation dyad, nothing else", so strengths come out in units of `alpha²` (`two_photon_norm = 'bc_5p142_reduced'`, the rank-K reduced element unity per channel). |
| Pipi | `1/(4ℏ²)` in Eq. (A5)/(S10) with `E₁`, `E₂` as field vectors, and `ℏΩ/2` in `Ĥ_int` ([P26] Eq. (2) p. 4; [P24] Eqs. (S4)/(S5) p. 11). |

**Verdict: the prefactors are a convention difference with no conflict; the
absolute amplitudes are NOT IN PAPER.**

No `E₀ cos ωt`, no `E_rms`, no intensity-to-field conversion, and no
`W/cm²` anywhere in either paper (§5 query 4). The amplitudes are *defined
implicitly* by normalizing one transition: `Ω_ref/2π = 2.000 kHz` for H₃O⁺
([P26] Eqs. (A8)/(A9) p. 16; the starred row of [P24] Table S4 p. 28) and
`2π × 2.087 kHz` for CaH⁺ ([P24] p. 12). Since `heff` carries no amplitude
either, the two conventions do not collide — but **no cross-calibration of
`alpha_K2_dOm0` / `alpha_K2_dOm2` against these papers is possible**, in either
direction. A caller wiring `heff` geometry into a Rabi rate in Pipi's convention
must supply `1/(4ℏ²)` and the two field magnitudes, and then halve it for
`Ĥ_int`'s `Ω/2`.

### (g) Lab vs molecule frame, and the `q = Ω′ − Ω` rule

| | |
|---|---|
| `heff` | molecule-frame component `q = Ω′ − Ω`, entering B&C (5.186)'s 3j `(J′ K J; −Ω′ q Ω)`, which vanishes unless `\|q\| ≤ K`; with `k₁ = k₂ = 1` bounding `K ≤ 2`, `ΔΩ = 0` admits `K = 0, 1, 2` and `ΔΩ = ±2` admits `K = 2` only ([HAM] §9.5.2, verified numerically there). Lab projection `Δm_F = P`. |
| Pipi | both frames named, no transformation written: `K` is "its projection onto the molecule-fixed symmetry axis", `m_F` "the projection of the total angular momentum onto the laboratory-frame quantization axis" ([P26] Eq. (A1) and text, p. 16). Selection rule `ΔK = 0` ([P26] Eq. (A4) p. 16; [P24] Table S4 caption p. 28). |

**Verdict: CONFIRMED in structure for the `ΔΩ = 0` analogue; the `ΔΩ = ±2`
channel is NOT IN PAPER.**

`ΔK = 0` **is** the `q₁ + q₂ = 0` statement for a molecule whose two E1 legs both
carry `q = 0` `[derived]`: H₃O⁺ is C₃ᵥ and its ground-vibrational-state dipole
lies along the symmetry axis, so the molecule-frame dipole has no perpendicular
component and each leg is `q = 0`. That is the direct analogue of `heff`'s rule
at `Ω′ = Ω`. Pipi never state this reason, so the `[derived]` tag is mine, not
theirs.

`heff`'s `alpha_K2_dOm2` channel has **no Pipi counterpart**: `|q| = 2` requires
two perpendicular (`q = ±1`) legs through an `Ω = 0` intermediate ([HAM] §9.5.2,
[2γ] §3.3), which the axial-dipole symmetry of both Pipi molecules forbids. This
is a molecule difference, not a notation difference, and needs no `heff` change.

### (h) The rank-1 (antisymmetric) part

| | |
|---|---|
| `heff` | `K = 1` is **not registered** (OPEN-21). `dyad_weights` returns `K = 1` weights, but no operator consumes them. [HAM] §9.5.3 proves `α¹ ∝ P_X ½[d_a, d_b] P_X = 0` in **exact closure over a complete intermediate manifold**, measured at `2.2 × 10⁻¹⁶`, and shows it returns at `O(δ/Δ)` for unequal denominators and at `O(1)` for a restricted manifold. |
| Pipi | no rank language at all; but Eq. (A5)/(S10) is `A/D₁ + B/D₂` with `D₁ = ω_{JM} − ω₁`, `D₂ = ω_{JM} + ω₂`, and `B` is `A` with the two legs exchanged. |

**Verdict: DIFFERS — physics, not convention.**

`[derived]`, elementary algebra: writing `X = A`, `Y = B` (leg-exchanged),

```
A/D₁ + B/D₂  =  ½(1/D₁ + 1/D₂)(X + Y)  +  ½(1/D₁ − 1/D₂)(X − Y)
```

Under B&C (5.141) at `k₁ = k₂ = 1`, leg exchange carries `(−1)^{1+1−K} = (−1)^K`,
so `(X + Y)` is the `K = 0, 2` content and `(X − Y)` is the `K = 1` content —
exactly [HAM] §9.5.3's antisymmetry statement. The rank-1 coefficient is

```
½(1/D₁ − 1/D₂) = ½ (D₂ − D₁)/(D₁ D₂) = ½ (ω₁ + ω₂)/(D₁ D₂)
```

which is **zero only if `D₁ = D₂`**, i.e. only in the common-denominator closure
that `heff` assumes and Pipi do not. Its size relative to the retained rank-0/2
coefficient is

```
(ω₁ + ω₂) / (2 ω_{JM} − ω₁ + ω₂)     [derived]
```

and [P24] p. 12 states outright that `|ω_{JM} − ω₁|` is comparable to
`|ω_{JM} + ω₂|`, i.e. that they are in the regime where neither term dominates.

**The magnitude is UNVERIFIED** and is not used downstream: evaluating that
ratio needs `ω_{JM}` and the laser frequencies, and both papers defer the
H₃O⁺ level structure to the unpublished "subsequent article". What *is* verified
is that the rank-1 part is present with a nonzero coefficient, and that its
**geometric** weight for a (π, σ±) pair equals the rank-2 weight `heff` keeps:
`|w¹_P| = |w²_P| = 0.707107` exactly (measured, §4).

For `heff`'s own ThF⁺ regime this source of `K = 1` is separately negligible
`[derived]`: with `ω₁ ≈ ω₂` optical and `D₁ ~ 0.16-1.5 GHz` ([HAM] §9.5.5),
`D₂ ≈ 2ω_optical`, so the `B` term is down by `~D₁/2ω_optical ~ 10⁻⁵`. **`heff`'s
one-term Raman reading is precisely the RWA that [P24] declines** — appropriate
for a near-resonant Raman, not for Pipi's. This is a *different* mechanism from
the `O(1)` restricted-manifold `K = 1` that [HAM] §9.5.3 and [AUDIT] defect 5
already escalate; the two must not be conflated. **STOP item 2**, folded into
the existing OPEN-21 escalation rather than opened as a new one.

### (i) Selection rules

| rule | `heff` | Pipi | verdict |
|---|---|---|---|
| `ΔJ` | `\|ΔJ\| ≤ K = 2`, so `{0, ±1, ±2}` ([HAM] §9.5.4); `ΔJ = ±1` measured live inside one parity eigenspace at 0.144-0.189 (§4) | `ΔJ = 0, ±2` ([P26] Eq. (A4) p. 16; [P24] Table S4 caption p. 28) | **DIFFERS**, traced to the molecule — see below |
| `Δm_F` reach | `Δm_F = P`, `\|P\| ≤ 2`, never `±3` (gate `test_channel_reach`); `(π, σ±)` gives `\|Δm_F\| = 1` | `Δm_F = ±1` for the `(π, σ±)` pair (Eq. (A10); Table S4 caption) | **CONFIRMED** on the reach |
| `Δm_F` sign | see row (d) | see row (d) | **DIFFERS** (and the two papers disagree) |
| parity | even; conserves the `parity_operator` eigenvalue, `[T, Par] = 0` measured at `0.00e+00` (§4, gate V24) | `p_J = p_J′` ([P26] Eq. (A4)) | **CONFIRMED** |
| body-frame projection | `ΔΩ ∈ {0, ±2}` | `ΔK = 0` | **CONFIRMED** for `ΔΩ = 0`; `±2` **NOT IN PAPER** (row g) |
| `ΔF`, `ΔF₁` | `\|ΔF\| ≤ K`, `\|ΔF₁\| ≤ K` ([HAM] §9.5.4) | not stated | **NOT IN PAPER** |

**On `ΔJ`.** Pipi's `ΔJ = 0, ±2` and their `p_J = p_J′` are one statement, not
two `[derived]`: for both their molecules the total parity of a level is the
inversion (or `Σ`-state) label times `(−1)^J`, so a parity-even operator that
conserves the label forces `ΔJ` even. `heff`'s operator is *also* parity-even
(measured `[T, Par] = 0` exactly) yet still reaches `ΔJ = ±1`, because ThF⁺
X³Δ₁ is an **Ω doublet**: both parities exist at **every** `J`, measured in §4
as equal-dimension parity subspaces (36/36 at `J = 1`, 60/60 at `J = 2`, 84/84 at
`J = 3`). So `ΔJ = ±1` connects two states of the *same* total parity and is
parity-allowed. **This is a molecule difference, and `heff` needs no change**;
`heff`'s `ΔJ ∈ {0, ±1, ±2}` is correct for its molecule and Pipi's
`ΔJ = 0, ±2` is correct for theirs.

One wording caveat, for Arian, not a code issue: [HAM] §9.5.4's parity row reads
"even — at zero field it does not connect e to f". §4 shows the operator does
connect `J` to `J ± 1` within one `parity_operator` eigenvalue. That is only
consistent if "e to f" there means the `parity_operator` eigenvalue rather than
the conventional `(−1)^J`-absorbed spectroscopic `e`/`f` label, under which those
same pairs *are* `e ↔ f`. Flagged once; outside this note's write scope.

---

## 3. Concrete discrepancies, and what `heff` would have to change

| # | row | what | `heff` change | observable? |
|---|---|---|---|---|
| D1 | (e) | `heff`'s default `reading='ladder'` is the two-absorbed-photon reading; Pipi's process is Raman, and only `reading='raman'` matches the `A` term of Eq. (A5). | **None to code.** Callers must pass `reading='raman'` to `heff.twophoton.dyad_weights` / `two_photon_line_strengths` for any Pipi comparison. Optional doc-only proposal: cite [P24] Fig. S11 p. 39 in [HAM] §9.5.1(3) as an external anchor for the Raman-reading assignment. | **Labels only.** The reachable `Δm_F` *set* is identical between readings ([HAM] §9.5.4); only the beam-pair → `Δm_F` mapping inverts. |
| D2 | (d) | The σ± ↔ `Δm_F` sign. `heff` (Raman) matches [P24] Fig. S11 and contradicts [P26] Eq. (A10); the two papers contradict each other. | **None. STOP item 1** — a physics/convention question for the authors, not a fix. `heff`'s own convention is internally anchored and gated (σ⁺ = `ê_{+1}`, `c[+1] = +1`, [AUDIT] §2), so there is nothing to change until the papers are reconciled. | **Yes, operationally.** It decides *which laboratory beam* drives a wanted `Δm_F`. Taking [P26] Eq. (A10) at face value with `heff`'s Raman reading targets the opposite `Δm_F`, i.e. the wrong transition. Line *strengths* are unaffected (the magnitudes are equal). |
| D3 | (h) | `heff` has no representation for Pipi's second time ordering, and therefore drops a rank-1 part whose geometric weight equals the retained rank-2 weight for a (π, σ±) pair. | **None here. STOP item 2.** If `heff` ever needs Pipi's regime it is a new operator with its own parameter, exactly as [HAM] §9.5.3 rules ("either requires a distinct operator and parameter definition"). Same escalation as [AUDIT] defect 5 / OPEN-21. | **Yes, in principle** — a dropped amplitude channel changes line strengths. Magnitude **UNVERIFIED** for Pipi's numbers; negligible (`~10⁻⁵`) for `heff`'s own ThF⁺ detunings `[derived]`. |
| D4 | (b), (f) | No polarizability tensor and no absolute field amplitude in either paper. | **None** — this is a GAP, not a discrepancy. `alpha_K2_dOm0` and `alpha_K2_dOm2` cannot be pinned or cross-checked from these papers. | No. |
| D5 | (g) | `alpha_K2_dOm2` has no Pipi analogue (their molecules' dipoles are axial, so `ΔK = 0` strictly). | **None.** | No. Scope only. |
| D6 | (i) | `ΔJ = ±1` is allowed in `heff` and forbidden in Pipi. | **None** — traced to the Ω-doublet structure of X³Δ₁ versus single-parity-per-`J` levels, and measured in §4. A doc-wording item on [HAM] §9.5.4's "does not connect e to f" is flagged in row (i); outside this note's write scope. | No, once the molecule difference is recognised. |

**Nothing in §2 warrants a code change to `heff/twophoton.py`.** Two items are
STOP items for Arian (D2, D3); one is a caller instruction (D1); three are
scope or wording (D4, D5, D6).

---

## 4. Numerical check

**What is and is not checkable.** Neither paper gives a closed-form two-photon
matrix element. Eq. (A5)/(S10) is a numerical sum over an intermediate manifold
`{|M⟩}` whose members, energies `ω_{JM}`, and dipole elements are **not
published** (both papers defer them), and whose absolute scale is fixed by
normalizing a reference transition. The Table S4 magnitudes therefore cannot be
reproduced or even ratio-checked: the H₃O⁺ eigenvectors are Zeeman/spin-rotation
mixtures from an unpublished Hamiltonian, and `heff`'s registered channels are
case (c2) with `Ω = ±1`, which H₃O⁺ and CaH⁺ are not (§1.3). **No magnitude
check is run, and none is invented.**

What *is* checkable are the two **selection-rule statements** the papers make
explicitly and `heff`'s geometry evaluates independently.

### CHECK 1 — the σ± ↔ `Δm_F` map, [P26] Eq. (A10) p. 16 vs [P24] Fig. S11 p. 39

**Stated before running.** PASS against a given paper = `heff`'s
`dyad_weights(eps1='pi', eps2=σ, reading='raman')` puts its single nonzero
weight at the `Δm_F` that paper assigns to `σ`, for both signs of `σ`. FAIL =
the opposite assignment. **Both are reachable**, and the demonstration is built
into the check: flipping to `reading='ladder'` inverts the assignment, so the
same code returns PASS and FAIL for the same paper depending on the reading.
The two papers' claims are mutually exclusive, so exactly one PASS per reading
is the only possible outcome.

### CHECK 2 — `ΔJ` reach, [P26] Eq. (A4) p. 16

**Stated before running.** PASS = `heff`'s registered `(K = 2, ΔΩ = 0)` channel
has no `|ΔJ| = 1` matrix element inside a single parity eigenspace, matching
`ΔJ = 0, ±2`. FAIL = it does. **Negative control**: the raw (non-parity-adapted)
`|ΔJ| = 1` column must be nonzero, so that a zero in the parity columns is
parity's doing and not basis truncation. Run at `J_max = 3` so `|ΔJ| = 2` pairs
genuinely exist in the basis.

### Output, verbatim

`conda run -n structure python chk.py` (scratch, session tmp):

```
CHECK 1  sigma+- <-> Delta m_F for Pipi's pair (E1 = pi absorbed, E2 = sigma+- emitted)
  heff dyad_weights(eps1=pi, eps2=sigma+-); P = Delta m_F
    reading=raman  eps2=sigma+  -> Delta m_F = [-1]   |w^1_P| = 0.707107  |w^2_P| = 0.707107
    reading=raman  eps2=sigma-  -> Delta m_F = [1]    |w^1_P| = 0.707107  |w^2_P| = 0.707107
    reading=ladder eps2=sigma+  -> Delta m_F = [1]    |w^1_P| = 0.707107  |w^2_P| = 0.707107
    reading=ladder eps2=sigma-  -> Delta m_F = [-1]   |w^1_P| = 0.707107  |w^2_P| = 0.707107
  Pipi 2026 Eq.(A10): sigma+ -> dmF=+1, sigma- -> dmF=-1
  Pipi 2024 Fig.S11 : sigma-(emit) -> dm=+1  (so sigma+ -> dm=-1)
    reading=raman : vs Pipi2026(A10) = FAIL   vs Pipi2024(Fig.S11) = PASS
    reading=ladder: vs Pipi2026(A10) = PASS   vs Pipi2024(Fig.S11) = FAIL
```

`conda run -n structure python chk2.py`:

```
CHECK 2 (extended, J_max=3 so |dJ|=2 pairs exist in the basis)
  basis: 360 kets, J in [1.0, 2.0, 3.0], ctx.S = 1.0; max|Par^2 - I| = 0.0e+00
  parity eigenvalues: [-1.0, 1.0]
    J=1: dim 72, tr(P+) = 36.0, tr(P-) = 36.0
    J=2: dim 120, tr(P+) = 60.0, tr(P-) = 60.0
    J=3: dim 168, tr(P+) = 84.0, tr(P-) = 84.0
  max |T| by |dJ|, channel (K=2, dOm=0), raw basis and inside one parity eigenspace:
  P=-2  raw |dJ|=0,1,2,3: 0.2000 0.3780 0.4140 0.0000 | par+1: 0.1000 0.1890 0.2070 0.0e+00 | par-1: 0.1000 0.1890 0.2070 0.0e+00
  P=-1  raw |dJ|=0,1,2,3: 0.1958 0.3162 0.2954 0.0000 | par+1: 0.0979 0.1581 0.1477 0.0e+00 | par-1: 0.0979 0.1581 0.1477 0.0e+00
  P=+0  raw |dJ|=0,1,2,3: 0.2500 0.2887 0.2791 0.0000 | par+1: 0.1250 0.1443 0.1396 0.0e+00 | par-1: 0.1250 0.1443 0.1396 0.0e+00
  P=+1  raw |dJ|=0,1,2,3: 0.1958 0.3162 0.2954 0.0000 | par+1: 0.0979 0.1581 0.1477 0.0e+00 | par-1: 0.0979 0.1581 0.1477 0.0e+00
  P=+2  raw |dJ|=0,1,2,3: 0.2000 0.3780 0.4140 0.0000 | par+1: 0.1000 0.1890 0.2070 0.0e+00 | par-1: 0.1000 0.1890 0.2070 0.0e+00
  commutator (gate V24 re-run here): max|[T,Par]| over P=-2..2 = 0.00e+00
```

### Reading the results

**CHECK 1.** `heff`'s Raman reading — the one that structurally matches Pipi's
`A` term, row (e) — gives `(π, σ+) → Δm_F = −1` and `(π, σ−) → Δm_F = +1`. That
**PASSES against [P24] Fig. S11** and **FAILS against [P26] Eq. (A10)**. The
ladder reading inverts both, so PASS and FAIL are both reachable from the same
code path and the check is not vacuous. `heff` is not wrong here; the two papers
disagree with each other, and `heff` sides with the one whose assignment follows
from the standard helicity convention (row d). The `|w¹_P| = |w²_P| = 0.707107`
columns are the row-(h) measurement: for this polarization pair the dropped
rank-1 channel carries exactly the same geometric weight as the retained
rank-2 one.

**CHECK 2.** `[T, Par] = 0.00e+00` exactly — parity is conserved — and yet
`|ΔJ| = 1` elements of 0.144-0.189 survive inside a single parity eigenspace, of
the same order as the `|ΔJ| = 0` (0.072-0.125) and `|ΔJ| = 2` (0.140-0.207)
elements. So the result is **FAIL against [P26] Eq. (A4)**, and the negative
control holds (the raw `|ΔJ| = 1` column is 0.289-0.378, so the basis does hold
such pairs). The cause is read off the same output: the parity subspaces have
**equal dimension at every `J`** (36/36, 60/60, 84/84), i.e. both parities exist
at every `J` because X³Δ₁ is an Ω doublet, so `ΔJ = ±1` between same-parity
states is allowed. Pipi's molecules carry one total parity per `(J, label)`, which
forces `ΔJ` even. **`heff` is right for its molecule and Pipi are right for
theirs** — row (i), item D6.

One honesty note on this output: the `|ΔJ| = 3` column is zero **by basis
truncation**, since `J ∈ {1, 2, 3}` admits no `|ΔJ| = 3` pair. It is not evidence
for the `|ΔJ| ≤ K` bound and is not used as such; the `|P| ≤ 2` half of that
bound is separately gated in `tests/test_twophoton.py::test_channel_reach`.

---

## 5. What I could not verify, with the query I ran

Every zero below was checked against the control term `Rabi` in the same two
extracts (27 hits in `pipi2024.txt`, 7 in `pipi2026.txt`), so the query shape
reaches populated content. Command shape:
`grep -c -iF '<term>' pipi2024.txt pipi2026.txt`.

1. **Any polarizability tensor.** `polariz` returns 6 / 25 hits, but
   `grep -o -iE 'polariz[a-z]*' | sort | uniq -c` shows **every** hit is
   `polarization` / `polarizations` / `polarized`; the substring
   `polarizability` appears **zero** times in either paper. Also zero:
   `Raman tensor`, `AC-Stark`, `light-shift tensor`. **No map from
   `alpha_K2_dOm0` / `alpha_K2_dOm2` to a Pipi symbol exists.**
2. **Any angular-momentum-algebra convention.** Zero hits in both papers for
   `Clebsch`, `Wigner`, `Eckart`, `Condon`, `Shortley`, `reduced matrix`,
   `irreducible`, and any 3j/6j/9j symbol. **Whether Pipi use Condon-Shortley is
   UNVERIFIED**; it is in the deferred "subsequent article" ([P24] ref. 7,
   [P26] ref. 48) and in unpublished code.
3. **Any polarization-vector component definition.** Zero hits for `x + iy`,
   `handed`, `helic`, and no statement of which of σ± raises `m`. Combined with
   item 4, this is why the row-(d) conflict is unresolvable.
4. **Beam geometry.** Zero hits for `beam` in both papers; the single
   `quantiz` hit ([P26] p. 16) names the lab quantization axis but not the beam
   direction relative to it. `propagat` hits are all about numerical propagation
   of the dynamics. **The σ± reference frame is UNVERIFIED**, and a beam-referred
   label is the leading candidate explanation for the [P24]/[P26] conflict.
5. **Any explicit field amplitude.** Zero hits for `E_0`, `E0 `, `intensity`,
   `W/cm`. Absolute `α` scale unrecoverable (row f).
6. **Any rank decomposition.** Zero hits for `antisymmet`, `irreducible`,
   `rank` in [P24]; the 6 `rank` hits in [P26] are all FNO tensor-factorization
   hyperparameters (p. 18). Pipi never decompose their two-photon operator by
   rank, so row (h)'s rank-1 statement is `[derived]` from their Eq. (A5), not
   quoted from them.
7. **The magnitude of the rank-1 contribution in Pipi's regime.**
   **UNVERIFIED** — needs `ω_{JM}` and the laser frequencies, which are not in
   either paper. Only the *presence* of the term is established, and only the
   presence is used.
8. **Any reproducible two-photon magnitude.** [P24] Table S4 (pp. 28-34),
   Table S2 (p. 18) and Fig. S11 (p. 39) give numbers, but all are normalized to
   a reference transition and all depend on unpublished eigenvectors and
   intermediate-state data. **Not checkable**, and no check was invented.
9. **`Ω` or `Λ` language.** Zero hits for `Lambda`; one `Hund` hit in [P24]
   (a bibliography entry). Neither paper's molecule has an electronic `Λ` or
   `Ω`, consistent with §1.3.

---

## 6. STOP items for Arian

1. **The σ± ↔ `Δm_F` sign conflict** between [P26] Eq. (A10) p. 16 and [P24]
   Fig. S11 p. 39, for the same (π absorbed, σ± emitted) pair. `heff` matches
   [P24]. Not resolvable from the papers: neither defines σ± in components and
   neither states the beam geometry. A one-sentence question to the authors — *is
   σ± referred to the quantization axis or to the beam propagation direction, and
   which of the two papers' `Δm_F` assignments is the intended one?* — settles it.
   Until then, `heff` changes nothing.
2. **The rank-1 part retained by Pipi's two-time-ordering form** and excluded by
   `heff` (OPEN-21). This is an external published instance of exactly the
   limitation [HAM] §9.5.3 and [AUDIT] defect 5 already flag, in a regime
   ([P24] p. 12) where the authors explicitly refuse the approximation `heff`
   relies on. It does not change `heff`'s ThF⁺ validity — that regime suppresses
   this source by `~10⁻⁵` `[derived]` — but it does mean the closure operator
   cannot be used to reproduce a Pipi-style calculation.
