# The effective Zeeman G tensor for ²³²ThF⁺ X ³Δ₁ — source audit and implementation handoff

Audits the working summary `ThF_Zeeman_effective_hamiltonian_summary.html` (2026-09-11)
against its primary sources, and specifies what it takes to replace `heff`'s axial-only
Zeeman term with the full three-parameter tensor. Resolves `OPEN-10` of
`docs/open-questions.md` ("Parity-dependent Zeeman terms are omitted").

Assembled 2026-09-14. Tags: `[verbatim]` quoted from a source read this session,
`[derived]` algebra run here, `[inferred]` an argument, not a citation.

Local source texts (extracted from ar5iv this session; untracked per `.gitignore`, like
every other `docs/lit` source text — regenerate from the arXiv IDs below):
- `docs/lit/petrov2025-arXiv2503.02840-ThFplus-gfactors.txt`
- `docs/lit/petrov2014-arXiv1404.4024-ThO-zeeman.txt`
- `docs/lit/petrov2017-arXiv1704.06631-HfFplus-zeeman.txt`

**Headline.** Every number in the summary's §3–§4 checks out against Petrov's papers:
the three second-order sums C₀₊, C₀₋, C₂ reproduce to all printed digits from the
published matrix elements, and the (G∥, G⊥, G_Δ) inversion is algebraically exact. The
one flagged contradiction — a 10× discrepancy in Δg — **does not exist**; it was a
misreading of Table I. The real finding is a framing correction: **Petrov's G∥ = 0.047 is
itself fitted to the measured g = 0.0149**, so the 0.58 % residual measures second-order
perturbation-theory truncation, not theory-versus-experiment agreement.

---

## 1. The structural formula is Petrov 2014, verbatim

The summary's §3 presents `g_{e/f}(J) = −G₀/[J(J+1)] + C₀± + C₂(J+2)(J−1)/[2J(J+1)]` as
"the ThO calculation". That notation is the summary's own; Petrov 2014 writes the same
content with explicit state labels. Eqs. (12)–(13), p. 4:

> `[verbatim]` g_e(J) = −G∥/[J(J+1)] + G⊥⁽²⁾Δ⁽²⁾/[T_e(H³Δ₁) − T_e(A³Π₀₊)]
> + G⊥⁽¹⁾Δ⁽¹⁾/[T_e(H³Δ₁) − T_e(Q³Δ₂)] · (J+2)(J−1)/[2J(J+1)]   (12)
>
> g_f(J) = −G∥/[J(J+1)] + G⊥⁽²⁾Δ⁽²⁾/[T_e(H³Δ₁) − T_e(³Π₀₋)]
> + G⊥⁽¹⁾Δ⁽¹⁾/[T_e(H³Δ₁) − T_e(Q³Δ₂)] · (J+2)(J−1)/[2J(J+1)]   (13)

So the summary's `G₀ ≡ G∥`, `C₀± ≡ Σ_{Ω=0±} G⊥ₙΔₙ/(E_X − Eₙ)`, `C₂ ≡ Σ_{Ω=2} G⊥ₙΔₙ/(E_X − Eₙ)`.
The generalisation from ThO (one state per Ω class) to ThF⁺ (two at Ω=0⁺, one at Ω=0⁻,
three at Ω=2) is a sum, nothing more `[inferred]`.

**e/f assignment hazard.** Petrov's Eq. (12) — the *e* level — carries ³Π₀₊; Eq. (13) —
*f* — carries ³Π₀₋. Which physical parity that corresponds to depends on the convention,
and `heff`'s `ef_rule` fork (`brown1975` vs `thesis_S_half`) **inverts every label at
S = 1** (`conventions.py:38`). The sign of G_Δ must be pinned against whichever fork is
active; its magnitude and the doublet difference 2|G_Δ| are invariant.

### 1.1 The inversion to a constant tensor is exact `[derived]`

Using (J+2)(J−1)/[2J(J+1)] = 1/2 − 1/[J(J+1)], Petrov's formula becomes

```
g(J) = −[G₀ + C₂]/[J(J+1)] + C̄₀ + C₂/2 ± (C₀₊−C₀₋)/2
```

Matching term by term against `g(J) = −G⊥ − (G∥−G⊥)/[J(J+1)] ± G_Δ`:

| | |
|---|---|
| G⊥ | = −(C̄₀ + C₂/2) |
| G∥ | = G₀ − C̄₀ + C₂/2 |
| G_Δ | = (C₀₊ − C₀₋)/2 |

All three reproduce the summary's §3 exactly. Through second order the ThO/ThF⁺ g-factor
spectrum **is** the spectrum of a J-independent rank-2 body-frame tensor, and a separate
leading rotational g_R is redundant.

---

## 2. Petrov's ThF⁺ inputs, verified

Electronic term values, arXiv:2503.02840 Eq. (3), p. 2 `[verbatim]`:
³Δ₁ T₀=0; ¹Σ⁺ T₀=314; 1³Δ₂ T₀=1052; 2Ω=2 T_v=5806; 3Ω=2 T_v=6402; ³Π₀₋ T₀=3044;
³Π₀₊ T₀=3395 cm⁻¹. (¹Σ⁺, 1³Δ₂, ³Π₀₊ from Barker 2012; the rest vertical at R = 3.75 a₀.)

Matrix elements, Eqs. (4)–(15), p. 3 `[verbatim]`, with
G⊥ₙ ≡ ⟨…|L̂ᵉ₊ − g_S Ŝᵉ₊|…⟩ and Δₙ ≡ 2B_rot⟨…|Ĵᵉ₊|…⟩, B_rot = 0.243 cm⁻¹:

| Eq | state | G⊥ | Eq | Δ (cm⁻¹) |
|---|---|---|---|---|
| 4 | ¹Σ⁺ | 0.48 | 10 | 0.250 |
| 5 | 1³Δ₂ | −2.41 | 11 | −0.510 |
| 6 | 2Ω=2 | −0.06 | 12 | −0.113 |
| 7 | 3Ω=2 | 0.93 | 13 | 0.427 |
| 8 | ³Π₀₋ | 1.07 | 14 | 0.587 |
| 9 | ³Π₀₊ | 1.01 | 15 | 0.555 |

> `[verbatim]` "We found that matrix element (4) provides the main contribution to the
> Ω-doubling effect. It was adjusted to reproduce the experimental value of 5.29 MHz for
> Ω-doubling measured in Ref. [Ng et al. 2022]. To achieve this, we increased the
> calculated value by 8.8%. Due to the similarity between matrix elements (4) and (10),
> matrix element (10) was also increased by 8.8%."

**The printed 0.48 and 0.250 are the pre-enhancement values** `[derived]`: applying ×1.088
to both reproduces the summary's C₀₊ exactly, while the raw values give −5.473×10⁻⁴, which
does not match. Recomputing all three sums from the table above:

| contribution | value |
|---|---|
| ¹Σ⁺ (Ω=0⁺, ×1.088²) | −4.5239 × 10⁻⁴ |
| ³Π₀₊ (Ω=0⁺) | −1.6511 × 10⁻⁴ |
| **C₀₊** | **−6.17497 × 10⁻⁴** (summary: −6.175 × 10⁻⁴) |
| **C₀₋** = ³Π₀₋ | **−2.06337 × 10⁻⁴** (summary: −2.063 × 10⁻⁴) |
| 1³Δ₂ | −1.16835 × 10⁻³ |
| 2Ω=2 | −1.1678 × 10⁻⁶ |
| 3Ω=2 | −6.2029 × 10⁻⁵ |
| **C₂** | **−1.231543 × 10⁻³** (summary: −1.232 × 10⁻³) |

Every digit the summary prints is confirmed.

### 2.1 G∥ = 0.047 is a fitted parameter, not an ab initio one

> `[verbatim]` (p. 3) "We also used G∥ = 0.047 to reproduce the g = 0.0149 value for the
> g-factor of the F = 3/2 state. Our value is slightly different from G∥ = 0.048 estimated
> in Ref. [Ng et al. 2022] since we account for nonadiabatic interactions between ³Δ₁ and
> other states in the basis set."

This is the framing correction. The independent theory content of the summary's §4 is
C₀±/C₂ — the nonadiabatic corrections. G₀ absorbs whatever remains to land on the measured
0.0149. So comparing the reconstructed tensor's J = 1 prediction to 0.0149 does **not**
test theory against data; it tests second-order perturbation theory against Petrov's full
numerical diagonalization `[inferred]`.

---

## 3. The Δg "contradiction" does not exist

`docs/digest-literature-thf-plus.md` carried two Petrov numbers for Δg that appeared to be
10× apart (line 188: 10⁷Δg₀ = 233.8; line 189: Δg(E=0) = 2.3 × 10⁻⁴). Resolved:

> `[verbatim]` (§IV Results, p. 4) "The calculated difference Δg = gᵘ − gˡ = 2.3 × 10⁻⁴
> between the g-factors of the upper (gᵘ) and lower (gˡ) levels of the Ω doublets **for a
> zero electric field** is in agreement with the experimental value |Δg| = 3(3) × 10⁻⁴
> [Ng et al. 2022] (the sign of Δg was not determined in Ref. [Ng et al. 2022])."

Table I is indexed **by electric field**, one row per E, each with its own local
linearisation Δg = Δg₀ + Δg₁E: at E = 40 V/cm, 10⁷Δg₀ = 349.4 and 10⁷Δg₁ = 6.3; at
E = 60 V/cm, 233.8 and 8.7. **Δg₀ is a per-field coefficient, not an E = 0 intercept**, so
there is no conflict with 2.3 × 10⁻⁴ at zero field. Δg falls steeply off zero field
(≈6 × 10⁻⁵ by 40 V/cm `[derived]` from the E = 40 row) then rises slowly — the same
Stark-suppression shape Petrov 2017 reports for HfF⁺, whose minimum Δg = 3 × 10⁻⁶ sits at
7 V/cm.

Correct `digest-literature-thf-plus.md:188` to say Δg₀/Δg₁ are per-field linearisation
coefficients.

---

## 4. How far the second-order tensor misses

Substituting §2's sums (G₀ = 0.047, I(¹⁹F) = 1/2, g_N μ_N/μ_B = 2.86345 × 10⁻³):

```
G∥ = 0.0467961      G⊥ = 1.02769e-3      G_Δ = -2.05580e-4
G_X^eff ≈ diag(8.221e-4, 1.2333e-3, 4.6796e-2)      G⁽²⁾₀ = 3.737e-2
```

| quantity | 2nd-order tensor | Petrov, full 7-state | measured (Ng 2022) |
|---|---|---|---|
| g(J=1, F=3/2) | −0.0149868 | −0.0149 (fit target) | 0.0149(3), sign undetermined |
| Δg(E=0), F=3/2 | 2.741 × 10⁻⁴ | 2.3 × 10⁻⁴ | \|Δg\| = 3(3) × 10⁻⁴ |

- **J = 1 g factor: 0.58 % high.** Inside the measurement's ±2 % band (0.29σ), but as §2.1
  says, the comparison that matters is against Petrov's 0.0149, and 0.58 % is the
  truncation error of stopping at second order.
- **Δg: 19 % high.** This is the number that actually probes the truncation, and the
  measurement cannot adjudicate it (100 % error bar, 0.09σ).
- **Independent magnitude check** `[derived]`: Leanhardt 2011 Eq. (65) gives
  |g′_rS| ≈ ω_ef/(2B_e) for the parity-dependent term. For ThF⁺,
  5.29 MHz / (2 × 0.243 cm⁻¹ × 29979.2458) = 3.63 × 10⁻⁴, against 2|G_Δ| = 4.11 × 10⁻⁴
  (rovibronic, before F projection) — 13 % apart, same ¹Σ⁺ mixing mechanism.

### 4.1 ThF⁺ data cannot separate G∥ from G⊥

At J = 1 the tensor collapses: `g_J(1) = −(G∥ + G⊥)/2`. Ng's single measurement therefore
constrains **only the sum**, to ±1.9 %. Splitting the two needs J ≥ 2:

```
g_J(1) − g_J(2) = −(G∥ − G⊥)/3 = −0.015256      (64 % of g_J(1))
```

enormous and easily measurable — and not measured. `digest-literature-thf-plus.md:213`
records the search: rotational g-factor for ThF⁺, **not found**; Ng 2022 §II D says only
"it might be of interest to compute" it. ThO is where multi-J g factors exist, and is where
Petrov 2014 validated the mechanism (his Table 2: pure case-(a) 1/[J(J+1)] scaling is
"badly violated", and the Ω = 0,2 mixing restores agreement).

**Consequence for the model:** fitting to measured shifts gives one equation for three
unknowns. It cannot replace the literature input; it can only rescale the sum. Take all
three constants from §2 and let g(J=1, F=3/2) be a prediction.

---

## 5. Implementation plan for `heff`

`B·G·J` with G body-diagonal decomposes into three manifestly covariant pieces `[derived]`:

```
B·G·J = (G∥ − G⊥)(B·n̂)(J·n̂) + G⊥ (B·J) + G_Δ [(B·êx)(J·êx) − (B·êy)(J·êy)]
```

whose diagonal elements are Ω²m/[J(J+1)], m, and ±(the doublet term) — i.e. exactly
`g(J) = −G⊥ − (G∥−G⊥)/[J(J+1)] ± G_Δ`. Three operators, one of which already exists.

| term | operator | status |
|---|---|---|
| A | (J·n̂)(n̂·B) | **exists** — `elements_c.py:182` `zeeman_Gpar`; coefficient `G_par` → `G∥ − G⊥` |
| B | B·J | new; rank-1 lab operator on J, ΔJ=0, ΔΩ=0, ΔF=0,±1 |
| C | body q=±2, ΔΩ=±2 | new; the parity-dependent term |

**Term C's source is Leanhardt 2011 Eq. (64)** (`digest-literature-thf-plus.md` §3.2),
in Brown/Nelis operator form:

```
H_ZeemanDist = −½ g_rS  μ_B (B₊J₋S₊S₋ + B₋J₊S₋S₊)      parity-even
H_ZeemanDoub = +½ g′_rS μ_B (B₊J₊S₊²  + B₋J₋S₋² )      parity-odd  ← term C
```

The ΔΩ=±2 sparsity is already supported: `omega_doubling` (`elements_c.py:55`) declares
`dOm=(-2.0, 2.0)`.

**Molecule-frame tensor formalism, for the shape of the implementation.** `C2V-Molecules`
carries the same physics for an asymmetric top in `atm_core/physics.py:244`
(`spin_aniso_zeeman_hamiltonian_element`, Sears 1984 eq. 21 form): the tensor lives as
`g_l[k][q]` over body rank k ∈ {0,1,2} and component q ∈ [−k,k], and the body-frame 3j
`w3j(N,k,N',-K,q,K')` enforces ΔK = q. `atm_core/transverse.py:69` is the generic-lab-q
twin. Two caveats `[inferred]`:
- Their tensor contracts B with **S** (a spin that decouples from rotation, hence the clean
  9j separating N-space rank k from S-space rank 1). Ours contracts B with **J**, whose
  body components have anomalous commutation. Do not substitute S → J in that expression.
- `physics.py:559` records a bug worth not repeating: the prefactor is μ_B with **no** g_s,
  because the tensor carries the absolute anisotropy. An extra g_s made their anisotropic
  Zeeman 2× too large.

### 5.1 Parameters

Replace the single `derived` `G_par = 0.04756` (currently fitted so |g_{F=3/2}| = 0.0149
exactly — `tests/test_observe.py:162`) with the §4 triple, status `derived-from-literature`,
cited to arXiv:2503.02840 Eqs. (3)–(15). `G_par = 0.04756` and Ng's 0.048 are the same
axial-only fit; they are not independent of Petrov's 0.047.

### 5.2 Acceptance gates

Both PASS and FAIL must be reachable for each `[inferred]`:

1. **Closed form.** `g_factors` at E=0 must reproduce −G⊥ − (G∥−G⊥)/[J(J+1)] ± G_Δ for
   J = 1…5 to round-off. Fails on any coefficient-mapping error in §5.
2. **J = 1 prediction.** |g(J=1,F=3/2)| = 0.014987, inside Ng's 0.0149(3). Fails if the
   term-B projection or the nuclear Zeeman sign is wrong.
3. **Doublet difference.** 2|G_Δ| projected to F=3/2 = 2.741 × 10⁻⁴, band [2.3, 3.6] × 10⁻⁴
   spanning Petrov's full calculation and the ω_ef/2B_e scaling. Gate on |Δg|, not its sign.
4. **Parity.** Term C must be parity-odd under `parity_operator`; terms A and B parity-even.
   Fails on a wrong ΔΩ=±2 phase.
5. **Reduction.** Setting G⊥ = G_Δ = 0 and G∥ = 0.04756 must reproduce today's numbers
   bit-for-bit. The regression gate.

### 5.3 What remains open

- **Third order.** The summary's appendix lists the inputs needed (excited-to-excited
  Coriolis and magnetic matrix elements). Petrov publishes only ³Δ₁-to-X elements, so this
  is blocked on new ab initio input, not on algebra.
- **The benchmark the summary asks for** — exact diagonalization of the published
  seven-state Hamiltonian for J = 1…5 — is buildable in `heff` from §2's table alone
  (add Ω=0±, 2 kets plus the Coriolis and electronic-Zeeman couplings). That is the only
  way to measure the truncation without waiting for a J ≥ 2 measurement. Larger change,
  and it would supersede the three-constant model rather than test it.
- **`D = −0.133 a.u.`** is what arXiv:2503.02840 p. 3 prints `[verbatim]`. It is 10× too
  small for ThF⁺ (−1.33 a.u. ≡ −3.38 D is the physical value); the digest's §7.2 flag on
  the printed digit is confirmed as a typo in the paper, not a transcription error.
