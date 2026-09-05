# Effective Hamiltonian: ²³²Th¹⁹F⁺ X ³Δ₁, v = 0, J = 1–4

Assembled 2026-09-04. Every number below carries a primary-source cite, a `[derived]` tag (algebra or numerics I ran, from cited inputs), or an `[estimate, method]` tag. Values from HfF⁺ are never carried into a ThF⁺ row; where HfF⁺ is used as a scaling anchor it says so.

Units: energies as frequencies (`h·ν`). 1 cm⁻¹ = 29 979.2458 MHz. Nuclear spins: ²³²Th I = 0 (no Th hyperfine, no Th quadrupole — see §2.11); ¹⁹F I = ½.

**Primary sources read at source for this document** (not through the digest): Ng PhD thesis App. B–C (PDF pp. 316–325) and Ch. 2.4/4.1; Ng et al. PRA 105, 022823 (2022) incl. Fig. 2 rendered and read as an image; Leanhardt et al. JMS 270, 1 (2011) Sec. IV A–B, G; Gresh et al. JMS 319, 1 (2016) §3 and Table 1 (rendered); Petrov & Skripnikov arXiv:2503.02840; Petrov, Skripnikov & Titov arXiv:2302.02856; Skripnikov & Titov PRA 91, 042504 (2015); Brown & Carrington (B&C) Eqs. 5.172, 5.174, 5.175, 5.179, 5.186, 6.234, 8.7, 8.20, 9.50, 9.51, 9.70 and Table 8.12. Full list with pages in §8.

**Three results in this document overturn or correct what the sources print.** They are flagged in place and collected in §7: (i) Ng thesis Eq. C.6's G∥ sign contradicts the g_F formula printed three lines below it and the G∥ = 0.048(2) that Ng quotes — `[derived]`, resolved in favour of the g_F formula, verified numerically; (ii) Arian's thesis e/f rule `P = (−1)^{J−S−ℓ}` is the S = ½ specialization and gives the **opposite** e/f label to B&C at S = 1 — a live trap for this molecule; (iii) Gresh's printed sign of k″ tracks the upper state's Ω = 0± label, not the physical ordering of the ThF⁺ Ω-doublet.

---

## 1. Basis and truncation

### 1.1 The recommended working basis

**Ω = ±1 block only, in the JILA product basis `|J, Ω, F, m_F⟩`**, with `F = J + I` (I = ½), Ω = ±1 the signed projection of the total electronic angular momentum on the internuclear axis, and J, F, m_F quantised along the instantaneous **E-field** direction.

This is the basis of Ng thesis App. C.3.1 (p. 323): "{|J, Ω, F, m_F⟩}", 12 states for J = 1, 32 for J = 1&2. It is a Hund's case (c)-like basis: Λ and Σ never appear, the rotational energy is `B J(J+1)` with no −Ω² term (Ng thesis Eq. 2.2a p. 37, Eq. C.7c p. 322), and the electronic structure enters only through the effective constants A∥, ω_ef, d_mf, G∥, E_eff. Equivalent case-(a) labels for bookkeeping: Λ = ±2, S = 1, Σ = ∓1, Ω = Λ+Σ = ±1.

**Two frames, and this is the most common way to get the code silently wrong.** Ng thesis App. C.1 (p. 318): the quantisation frame **rotates with E_rot** and defines F, m_F, J; the molecule frame has its axis along the internuclear vector and defines Ω. Petrov's papers use `M_F = M_I + M` with M the projection on the **lab** ẑ and warn explicitly "MF = MI + M is not equal to mF" (arXiv:2302.02856 p. 3). For a static-field calculation the two coincide; the moment a rotating E_rot is introduced they do not.

### 1.2 Why Ω = ±1 alone is enough for J ≤ 4 — the numbers

The Ω = 1 component of ³Δ can mix with Ω = 2 (i.e. ³Δ₂) through the S-uncoupling term `−B(J₊S₋ + J₋S₊)` in `H_tum = B_e(J − S)²` (Leanhardt Eq. 12, p. 13). The case-(a) matrix element is

```
⟨Λ, S, Σ+1, J, Ω+1| −B(J₊S₋ + J₋S₊) |Λ, S, Σ, J, Ω⟩
      = −B √(J(J+1) − Ω(Ω+1)) · √(S(S+1) − Σ(Σ+1))
```

which for ³Δ₁ → ³Δ₂ (Ω: 1→2, Σ: −1→0, S = 1) is `−B√2 · √(J(J+1) − 2)`.

With B = 0.242647 cm⁻¹ (§3) and ΔE(³Δ₂ − ³Δ₁) = 1052.5(1.0) cm⁻¹ (Gresh 2016 Table 2, from Barker et al.):

| J | matrix element (cm⁻¹) | mixing amplitude c | \|c\|² | 2nd-order shift |
|---|---|---|---|---|
| 1 | 0 (Ω = 2 needs J ≥ 2) | 0 | 0 | 0 |
| 2 | −0.6863 | −6.52 × 10⁻⁴ | 4.25 × 10⁻⁷ | 13.4 MHz |
| 3 | −1.0852 | −1.031 × 10⁻³ | 1.06 × 10⁻⁶ | 33.5 MHz |
| 4 | −1.4559 | −1.383 × 10⁻³ | 1.91 × 10⁻⁶ | 60.4 MHz |

`[derived]`, all four rows. Leanhardt makes the J = 1 statement independently: "Off-diagonal couplings in Ω are zero since J·S preserves the value of J (there is no level with J = 1 and Ω = 2)" (p. 14).

**The 13–60 MHz shifts are not a reason to enlarge the basis, because they are already inside the measured constants.** The second-order shift is `2B²(J(J+1) − 2)/ΔE`, exactly **linear in J(J+1)**, so it renormalises B and the band origin and nothing else:

- `δB = 2B²/ΔE = 1.119 × 10⁻⁴ cm⁻¹ = 3.354 MHz` `[derived]` — absorbed into the fitted B₀ = 7.2743325 GHz.
- The J⁴ residue it induces is `(2B²/ΔE)²/ΔE = 1.19 × 10⁻¹¹ cm⁻¹`, four orders below the measured D₀ = 1.30 × 10⁻⁷ cm⁻¹ `[derived]`. D₀ is instead reproduced by the vibrational Kratzer relation `D ≈ 4B³/ω_e² = 1.332 × 10⁻⁷ cm⁻¹` — agreeing with the measurement to 2 % `[derived]`, which is the cleanest available evidence that no anomalous J⁴ physics is hiding in the ³Δ₁ manifold.

What the Ω = 2 admixture **does** change, and cannot be absorbed: the differential g-factor. Ng's own 32-level Ω = ±1 model gives δg/g = −0.00223 against a measured −0.00255(6), and he attributes the 15 % gap to "coupling between the X ³Δ₁ and the ³Δ₂ states that is not accounted for in the 32-level Hamiltonian" (thesis p. 85). Petrov's ThF⁺ model, which does include 1³Δ₂, lands at Δg = 7.56 × 10⁻⁵ ⇒ δg/g = −0.00254 `[derived from arXiv:2503.02840 Table I]`, essentially on the measurement. **I reproduced Ng's −0.00223 exactly** with a 96-state Ω = ±1 model built from the operators in §2 (see §6, check V7). So: the Ω = ±1 block is right for level energies to well below a kHz; it is wrong at the 15 % level for δg, which matters for eEDM systematics but not for spectroscopy.

**Recommendation: v1 = Ω = ±1 block with effective (fitted) parameters.** Provide the full ³Δ manifold and the Ω = 0 perturbers as an opt-in extension (§1.4), needed only for ab initio Δg or for generating ω_ef rather than fitting it.

### 1.3 Dimension per M_F block

`[derived]`, exhaustive enumeration.

**Ω = ±1 block** (2 Ω values × F = J ± ½):

| basis | total states | \|M_F\| = ½ | 3/2 | 5/2 | 7/2 | 9/2 |
|---|---|---|---|---|---|---|
| J = 1 | 12 | 4 | 2 | – | – | – |
| J = 1–2 | 32 | 8 | 6 | 2 | – | – |
| J = 1–3 | 60 | 12 | 10 | 6 | 2 | – |
| **J = 1–4** | **96** | **16** | **14** | **10** | **6** | **2** |

(Each |M_F| value has two signed blocks, ±M_F; the column gives the dimension of one signed block. 2 × (16+14+10+6+2) = 96 ✓. The J = 1 and J = 1–2 totals match Ng's 12 and 32 exactly, Ng thesis p. 323.)

**Full ³Δ manifold** (Ω = ±1, ±2, ±3, i.e. Σ = −1, 0, +1 with Λ = ±2; |Ω| ≤ J):

| basis | total | \|M_F\| = ½ | 3/2 | 5/2 | 7/2 | 9/2 |
|---|---|---|---|---|---|---|
| J = 1–4 | 244 | 36 | 34 | 28 | 18 | 6 |

Per-M_F blocking is exact for a collinear (E ∥ B ∥ ẑ) static problem; it is broken by the rotating-frame term `ℏω_rot F_x` (Ng Eq. C.4) and by any B⊥, both of which couple Δm_F = ±1.

### 1.4 The three tiers, so the decision is explicit

| tier | basis | dim (J = 1–4) | what it buys | what it costs |
|---|---|---|---|---|
| **T1 (recommended v1)** | Ω = ±1, effective ω_ef, A∥, G∥, d_mf | 96 | every energy in §4 to ≪ 1 kHz; δg to 15 % | ω_ef and G∥ are inputs, not outputs |
| T2 | + ³Δ₂ (Ω = ±2), + ³Δ₃ (Ω = ±3) | 244 | δg to the measurement; the A_SO/λ fine structure | needs A_SO, λ, and the ³Δ₁–³Δ₂ off-diagonal element |
| T3 | + Ω = 0 perturbers (a ¹Σ⁺, ³Π₀₋, ³Π₀₊) as in Petrov 2025 Eq. 3 | ~500+ | ω_ef **generated**, not fitted (see §2.3) | six ab initio electronic matrix elements |

For T2, the ³Δ fine-structure constants are not published for ThF⁺. Fitting `E(³Δ_Ω) = 2AΣ + (2/3)λ(3Σ² − S²)` to T₀(³Δ₂) − T₀(³Δ₁) = 1052.5 cm⁻¹ and T₀(³Δ₃) − T₀(³Δ₂) = 2097.5 cm⁻¹ (Gresh 2016 Table 2) gives **A ≈ 787.5 cm⁻¹, λ ≈ 261 cm⁻¹** `[derived]`. The two intervals are not equal (ratio 2.0, not 1.0), so a single spin-orbit constant cannot describe the ³Δ manifold; λ/A = 0.33 is large and the a ¹Σ⁺ state 314 cm⁻¹ from ³Δ₁ certainly contributes. **OPEN-1**.

---

## 2. Term-by-term Hamiltonian

Convention contract (`digest-thesis-effective-hamiltonian.md` §4, verified against B&C where cited):

- Wigner–Eckart: `⟨η,J,M|T^k_p(A)|η′,J′,M′⟩ = (−1)^{J−M} (J k J′; −M p M′) ⟨η,J‖T^k(A)‖η′,J′⟩` — B&C Eq. (5.172), PDF p. 205 = Arian Eq. 2.6.
- Spectator theorem, operator on the **first** coupled part (here J, with I the spectator) — **B&C Appendix Eq. (5.174), PDF p. 205 / book p. 173** (the main-text form is not used, per the brief):
  `⟨j₁,j₂,j₁₂‖T^{k₁}(A₁)‖j₁′,j₂′,j₁₂′⟩ = δ_{j₂j₂′} (−1)^{j₁₂′+j₁+k₁+j₂} [(2j₁₂+1)(2j₁₂′+1)]^{1/2} {j₁′ j₁₂′ j₂; j₁₂ j₁ k₁} ⟨j₁‖T^{k₁}(A₁)‖j₁′⟩`
  and, for the operator on the **second** part (here I), Eq. (5.175).
- Rotation-matrix reduced element: `⟨J,Ω‖D^{(k)}_{·q}(ω)*‖J′,Ω′⟩ = (−1)^{J−Ω} (J k J′; −Ω q Ω′) [(2J+1)(2J′+1)]^{1/2}` — B&C Eq. (5.186), PDF p. 207.
- `⟨J‖T¹(J)‖J′⟩ = δ_{JJ′}[J(J+1)(2J+1)]^{1/2}` — B&C Eq. (5.179).
- Condon–Shortley throughout; anomalous commutation for J (not for L, S, I).
- Parity: `E* = σ_xz R_y(π)`; case (a) composite phase `(−1)^{J−S−ℓ+s}`, which for ³Δ (S = 1, ℓ = 0, s = 0) is `(−1)^{J−1}`. B&C states the same phase as `(−1)^{J−S+s}` in Eq. (6.234), PDF p. 283.

**Result stated once, because it saves work: Ng's Appendix C matrix elements are already in Arian's convention.** Deriving H_Stark from B&C (5.174) + (5.172) + (5.186) reproduces Ng Eq. C.5 term for term, including every phase, with the two 6j symbols related by column exchange `{J′ F′ I; F J 1} = {J F I; F′ J′ 1}` `[derived]`. No conversion layer is needed.

Throughout, `|J, Ω, F, m_F⟩` with Ω′ ≡ Ω wherever an operator is diagonal in Ω.

---

### 2.1 Rotation and centrifugal distortion

**Operator** (Ng thesis Eq. C.7c–d, p. 322): `H_rot = B_v J(J+1) − D [J(J+1)]²`, with `B_v = B_e − α_e(v + ½)`.

Note there is **no −Ω² term**: this is the case-(c) convention of Ng Ch. 2.4 Eq. 2.2a (p. 37), `E(J) = B J(J+1), J ≥ Ω`, with the Ω-dependent constant absorbed into the band origin. It is why "4B = 29.09733(4) GHz" is a literal statement about the J = 1 → 2 interval (Ng thesis p. 78, Ng 2022 §II B).

**Matrix element**: diagonal, `δ_{JJ′} δ_{ΩΩ′} δ_{FF′} δ_{m_F m_F′}`.

**Selection rules**: ΔJ = ΔΩ = ΔF = Δm_F = 0.

**Size** `[derived from B₀ = 7.2743325 GHz, D₀ = 3.8973 kHz]`:

| J | B J(J+1) | D [J(J+1)]² |
|---|---|---|
| 1 | 14.5487 GHz | 15.6 kHz |
| 2 | 43.6460 GHz | 140.3 kHz |
| 3 | 87.2920 GHz | 561.2 kHz |
| 4 | 145.4866 GHz | 1.559 MHz |

---

### 2.2 Spin–orbit A

**Operator** (Leanhardt Eq. 11, p. 13): `H_SO = A Λ Σ`. Diagonal in Ω. Within the Ω = ±1 block it is a **constant** (Λ Σ = −2 for both Ω = +1 and Ω = −1), and is absorbed into the electronic origin. It never appears in the v1 Hamiltonian.

It matters only for T2/T3 (§1.4), where A ≈ 787.5 cm⁻¹ `[derived, §1.4]` sets the ³Δ₁–³Δ₂ separation and hence the S-uncoupling admixture of §1.2. Leanhardt's justification for dropping it, verbatim (p. 13): the spin-orbit, tumbling, spin-spin and spin-rotation terms "primarily describe an overall shift of the ³Δ₁ J-level, and can be ignored in evaluating energy differences in the states we care about."

---

### 2.3 Ω-doubling

**Operator** (Ng thesis Eq. C.3, p. 319, verbatim):

```
H_Ω-doub = ⊕_{J ≥ Ω}  ((−1)^J / 2) ℏ ω_ef Ω_x^(J),
   Ω_x^(J) = [J(J+1)/2] ( |J, Ω=+1⟩⟨J, Ω=−1| + |J, Ω=−1⟩⟨J, Ω=+1| )
```

So the off-diagonal element between the two Ω components at fixed J, F, m_F is `(−1)^J ω_ef J(J+1)/4`, and the **splitting is `ω_ef J(J+1)/2`**, with ω_ef *defined* as the J = 1 splitting `[derived from the operator; confirmed against Ng 2022 Fig. 2, which draws ω_ef at J = 1 and 3ω_ef at J = 2]`.

**Selection rules**: ΔJ = ΔF = Δm_F = 0, Ω → −Ω. Diagonal in parity (it is the operator that *splits* parity).

**Which component is upper — settled from the primary figure.** I rendered Ng 2022 p. 3 and read Fig. 2 as an image. Caption: "The energy levels with positive (negative) parity are denoted by black (grey) lines." In the zero-field panel: at **J = 1** the grey (negative-parity) line is above the black in both F = 3/2 and F = 1/2; at **J = 2** the black (positive-parity) line is above the grey in both F = 5/2 and F = 3/2. So the parity of the upper component is `(−1)^J`, i.e. the **physical parity ordering alternates with J**, exactly as Ng's `(−1)^J` prefactor encodes. This is also what the case-(a) parity phase `(−1)^{J−S} = (−1)^{J−1}` for S = 1 requires: the symmetric combination `(|Ω=+1⟩+|Ω=−1⟩)/√2` is the positive-parity state only for odd J. `[derived, cross-checked two ways]`

In the **B&C e/f convention** — "For integral J values, levels with parities (−1)^J or (−1)^{J+1} are designated as e or f respectively" (B&C PDF p. 283 / book p. 251) — upper parity `(−1)^J` means **the e component lies `ω_ef J(J+1)/2` above the f component, uniformly in J**. That is the J-independent statement to code against.

> **TRAP (OPEN-2).** Arian's thesis defines e by `P = (−1)^{J−S−ℓ}` (thesis Fig. 4.14 caption). At S = ½ that is `(−1)^{J−1/2}` and agrees with B&C. At **S = 1 it is `(−1)^{J−1} = −(−1)^J`, the opposite of B&C.** Applying the thesis rule verbatim to ³Δ₁ inverts every e/f label in this molecule. Use B&C's J-only rule here.

**Microscopic origin, and it is worth having both pictures.**

*Case-(a) picture* (Leanhardt Eqs. 16–17, p. 13): `H_LD = ½(o_Δ + 3p_Δ + 6q_Δ)(S₊²J₊² + S₋²J₋²)`, a ΔΛ = ∓4, ΔΣ = ±2 effective operator whose electronic ΔΛ = 4 part (`e^{∓4iφ_e}`) is folded into the constant, with `õ_Δ = |o_Δ + 3p_Δ + 6q_Δ| ≈ Σ C_{Π,Σ,Π′} A²B_e² / [(E_Δ−E_Π)(E_Δ−E_Σ)(E_Δ−E_Π′)]` (third order in spin-orbit × rotation) and `ω_ef = 4õ_Δ` (Leanhardt p. 14). The `J₊²` matrix element `⟨Ω=−1|J₊²|Ω=+1⟩ = J(J+1)` gives the J(J+1) law directly `[derived]`.

*Case-(c) picture* (Petrov & Skripnikov 2025): ω_ef is generated at **second order** in the rotational coupling through the Ω = 0 states, `Ω=+1 → Ω=0 → Ω=−1`. The coupling operator is the off-diagonal part of `Ĥ_rot = B J² − 2B(J·Jᵉ)`, i.e. `−B(J₊Jᵉ₋ + J₋Jᵉ₊)`, and Petrov tabulates `Δ^(i) = 2B_rot⟨i|Ĵᵉ₊|³Δ₁⟩` (Eqs. 10, 14, 15). Second-order perturbation theory gives

```
H_eff(Ω=+1, Ω=−1) = −(J(J+1)/4) Σ_i σ_i (Δ^(i))² / T_i,     σ_i = +1 for 0⁺, −1 for 0⁻
```

`[derived]` — the reflection-parity sign σ_i is the standard statement that only the *imbalance* between 0⁺ and 0⁻ perturbers produces Ω-doubling. With Petrov's ThF⁺ numbers (a ¹Σ⁺: Δ = 0.250 cm⁻¹, T = 314; ³Π₀₋: 0.587, 3044; ³Π₀₊: 0.555, 3395):

| J | Σ contributions (cm⁻¹) | predicted splitting | ω_ef J(J+1)/2 |
|---|---|---|---|
| 1 | +9.952e−5 − 5.660e−5 + 4.536e−5 = 8.829e−5 | **5.294 MHz** | 5.290 MHz |
| 2 | 2.649e−4 | 15.881 MHz | 15.870 |
| 3 | 5.297e−4 | 31.762 MHz | 31.740 |
| 4 | 8.829e−4 | 52.937 MHz | 52.900 |

`[derived]`. Two-line caveat: Petrov *fitted* Δ^(1) upward by 8.8 % to reproduce 5.29 MHz, so the J = 1 agreement is partly by construction — but it confirms that a closed-form second-order expression reproduces their full diagonalisation to 0.07 % and reproduces the J(J+1) law exactly. With the unscaled ab initio Δ^(1) = 0.2298 cm⁻¹ the prediction would be 4.37 MHz `[derived]`, i.e. ab initio is ~17 % low.

**Sizes** (ω_ef = 5.29(5) MHz): J = 1 → 5.29 MHz, J = 2 → 15.87, J = 3 → 31.74, J = 4 → 52.90 MHz.

**No centrifugal correction to ω_ef exists for X ³Δ₁**: Gresh 2016 Table 1 lists `k_D″ = 0` for every ³Δ₁ lower state (verified from the rendered table). Fine to J = 4; would need care at high J.

---

### 2.4 ¹⁹F magnetic hyperfine — the axial term A∥

**Operator** (Ng thesis Eq. C.2, p. 319): `H_hf = A∥ (F² − I² − J²)/(2J²)`, obtained by projecting `(I·n̂)(J·n̂)` with the projection theorem.

**Matrix element** (diagonal): `A∥ [F(F+1) − I(I+1) − J(J+1)] / (2 J(J+1))`.

**Definition of A∥ — this is the row the brief asks for.** B&C Eq. (9.50), PDF p. 636 / book p. 604, verbatim:

```
⟨η, J|H_hf|η, J⟩ = {aΛ + (b_F + (2/3)c)Σ} Ω {F(F+1) − J(J+1) − I(I+1)} / [2J(J+1)]
```

so **`A∥ ≡ {aΛ + (b_F + (2/3)c) Σ} Ω`**, which for ³Δ₁ (Λ = +2, Σ = −1, Ω = +1) is **`A∥ = 2a − b_F − (2/3)c`** `[derived from B&C 9.50]`. a is the electron orbital (nuclear-spin–orbit) constant, b_F the Fermi contact, c the axial dipolar. None of a, b_F, c is separately known for ThF⁺ (§7, OPEN-5); only the combination is measured.

**Consequences, all `[derived]` and cross-checked against Ng 2022 Fig. 2:**
- splitting `E(F=J+½) − E(F=J−½) = A∥ (2J+1) / (2J(J+1))`
- J = 1: `(3/4) A∥`; J = 2: `(5/12) A∥`; ratio 9/5 (stated by Ng thesis p. 319 and drawn in Ng 2022 Fig. 2)
- J = 3: `(7/24) A∥`; J = 4: `(9/40) A∥`

**Selection rules**: ΔJ = 0 (in Ng's form), ΔΩ = 0, ΔF = 0, Δm_F = 0. Parity-even.

**Sizes** (A∥ = −20.1(1) MHz, so **F = J−½ lies above F = J+½**):

| J | splitting | E(F = J+½) | E(F = J−½) |
|---|---|---|---|
| 1 | −15.075 MHz | −5.025 MHz | +10.050 MHz |
| 2 | −8.375 | −3.350 | +5.025 |
| 3 | −5.863 | −2.513 | +3.350 |
| 4 | −4.523 | −2.010 | +2.513 |

---

### 2.5 ¹⁹F hyperfine, ΔJ = ±1 — a term Ng drops that is above the kHz floor

B&C Eq. (9.51), PDF p. 636 / book p. 604, verbatim:

```
⟨η, J|H_hf|η, J−1⟩ = −{aΛ + (b_F + (2/3)c)Σ}
      × (J² − Ω²)^{1/2} {(F−I+J)(F+I+J+1)(J+I−F)(F−J+I+1)}^{1/2} / [2J(4J²−1)^{1/2}]
```

i.e. the same A∥ combination, off-diagonal in J. Ng's Eq. C.2 is strictly diagonal in J and omits it.

**Selection rules**: ΔJ = ±1, ΔΩ = 0, ΔF = 0, Δm_F = 0, parity-preserving.

**Sizes** `[derived]`, with A∥ = −20.1 MHz, Ω = 1, and the second-order shift `∓|ME|²/ΔE_rot`:

| pair | common F | matrix element | shift on lower J | shift on upper J |
|---|---|---|---|---|
| J = 1 ↔ 2 | 3/2 | +8.704 MHz | −2.603 kHz | +2.603 kHz |
| J = 2 ↔ 3 | 5/2 | +9.475 | −2.057 kHz | +2.057 kHz |
| J = 3 ↔ 4 | 7/2 | +9.731 | −1.627 kHz | +1.627 kHz |

(F = J−½ levels have no ΔJ = ±1 partner: the formula's factor `(F−J+I+1)` vanishes. Verified numerically.)

So this term shifts J = 1 F = 3/2 by −2.6 kHz relative to F = 1/2, changing the J = 1 hyperfine splitting by 0.017 % — invisible in Ng's 0.5 %-precision fit, but **above the kHz threshold this document is written to**, and it is the single largest term missing from Ng's Appendix C. Include it in v1; it is one extra matrix element from a formula already needed.

---

### 2.6 Nuclear spin–rotation c_I — not measured, estimated

**Operator** (B&C Eq. 8.7, PDF p. 410 / book p. 378): `H_nsr = c_I T¹(J) · T¹(I)`.

**Matrix element in the coupled basis** (B&C Eq. 8.20, PDF p. 414 / book p. 382, verbatim):

```
⟨η; J, I, F, M_F| c_I T¹(J)·T¹(I) |η; J, I, F′, M_F⟩
   = c_I (−1)^{J+F+I} δ_{FF′} {I J F; J I 1} {J(J+1)(2J+1) I(I+1)(2I+1)}^{1/2}
   = c_I [F(F+1) − I(I+1) − J(J+1)] / 2
```

**Selection rules**: ΔJ = 0, ΔF = 0, Δm_F = 0, ΔΩ = 0. Parity-even.

**Why it is not degenerate with A∥.** The two have opposite J-dependence:
- A∥ splitting `= A∥(2J+1)/(2J(J+1))` — falls as ~1/J
- c_I splitting `= c_I(2J+1)/2` — grows as ~J

**Estimate `[estimate, method: empirical ¹⁹F anchor + linear-in-B scaling]`.** No ThF⁺ value exists in any source read. B&C Table 8.12 (PDF p. 513 / book p. 481) gives the ¹⁹F spin–rotation constant in **CsF X¹Σ⁺**: `c₂ = 15.1 kHz at B_v = 0.183782 cm⁻¹` (v = 0; 15.0 at v = 1, 14.7 at v = 2). CsF is the closest available analogue: an essentially ionic heavy-metal fluoride, F⁻-like, with a rotational constant within 35 % of ThF⁺'s. The rotational magnetic field at the nucleus scales with the rotation rate, so scaling linearly in B:

`c_I(ThF⁺) ≈ 15.1 kHz × (0.242646 / 0.183782) ≈ 20 kHz`, **assign 20 kHz with a factor-of-3 uncertainty** (the electronic/paramagnetic part of c_I is the dominant and the least transferable piece; B&C p. 453 notes it "contains both first- and second-order parts which are opposite in sign"). ThF⁺ is also open-shell, which CsF is not.

**Sizes at c_I = 20 kHz** `[derived]`: hyperfine splitting contribution `c_I(2J+1)/2` = 30, 50, 70, 90 kHz for J = 1, 2, 3, 4; largest single-eigenvalue shift over J ≤ 4 is 50 kHz (numerically confirmed).

**Consequence for A∥ that nobody in the literature states.** Ng 2022 fits A∥ and ω_ef *simultaneously* to the six J = 1 ↔ J = 2 microwave lines (Ng 2022 §II B, verbatim: "just six distinct resonant frequencies ... from which we perform a fit to the spectroscopic constants A∥ and ω_ef"), assuming the 9/5 ratio built into the operator. With c_I ≠ 0 present but unmodelled, the J = 1 and J = 2 splittings imply `A∥^fit = A∥^true + 2c_I` and `A∥^true + 6c_I` respectively — inconsistent by `4c_I ≈ 80 kHz`, i.e. **0.4 % of A∥, against a quoted uncertainty of 0.5 %** `[derived]`. So c_I up to ~25 kHz is invisible in their data, and **A∥ = −20.1(1) MHz may carry an unmodelled c_I bias comparable to its own error bar.** OPEN-6.

---

### 2.7 Stark

**Operator** (Ng thesis Eq. C.5, p. 320): `H_Stark = −d_mf n̂ · E`, with d_mf > 0 and n̂ the internuclear axis. `q = 0` because the dipole lies along the axis.

**Matrix element** (Ng Eq. C.5; independently re-derived here from B&C (5.174) + (5.172) + (5.186) and found identical):

```
⟨J′,Ω′,F′,m′_F| H_Stark |J,Ω,F,m_F⟩
  = −d_mf Σ_{p∈{0,±1}} E_p
      × (−1)^{F+J′+1+I} { J  F  I ;  F′ J′ 1 }
      × (−1)^{F′−m′_F} √((2F+1)(2F′+1)) ( F′ 1 F ; −m′_F  p  m_F )
      × (−1)^{J′−Ω′} √((2J+1)(2J′+1)) ( J′ 1 J ; −Ω′  0  Ω )
```

(round = 3j, curly = 6j.) Ng's own derivation pointer is B&C p. 167 (book) = PDF p. 199.

**Selection rules**: ΔΩ = 0 (forced by q = 0); ΔJ = 0, ±1; ΔF = 0, ±1; Δm_F = p (so Δm_F = 0 for E ∥ ẑ). **Parity-odd** — within a J it connects the two parity components of the Ω-doublet, which is the whole polarisation mechanism.

**Closed form for the ΔJ = 0 piece, verified numerically.** For every (J = 1–4, F = J±½, two m_F, Ω = ±1) I checked that the expression above with J′=J, F′=F, Ω′=Ω, p=0 equals exactly

```
−Ω m_F γ_F d_mf E,      γ_F = [J(J+1) + F(F+1) − I(I+1)] / [2 F(F+1) J(J+1)]
```

`[derived]` — which is Leanhardt Eq. 20–21 (p. 14). γ_{F=3/2}(J=1) = 1/3, γ_{F=1/2}(J=1) = 2/3; γ_{5/2}(J=2) = 2/15, γ_{3/2}(J=2) = 1/5; γ_{7/2}(J=3) = 1/14, γ_{5/2}(J=3) = 2/21; γ_{9/2}(J=4) = 2/45, γ_{7/2}(J=4) = 1/18.

**Sign convention for "upper" and "lower".** Leanhardt p. 14, verbatim: "The electric field therefore raises the energy of the states with m_F Ω < 0 (denoted 'upper' states, superscript u), and lowers the energy of states with m_F Ω > 0 ('lower')." Petrov's equivalent geometric statement (arXiv:2302.02856 Fig. 1 caption): with n̂ from metal to F, "For upper (lower) Stark doublet n is parallel (antiparallel) to the external rotating electric field", and the experimental doublet switch is `D̃ = +1` for the **lower** doublet.

**Regimes and sizes** `[derived; d_mf = 3.37(9) D = 1.6965 MHz/(V/cm)]`. For the J = 1, F = 3/2, |m_F| = 3/2 stretched pair the exact two-level Ω-doublet + Stark eigenvalue is `±√((ω_ef/2)² + (γ_F m_F d_mf E)²)`:

| E | d_mf E | linear shift `γ_F m_F d_mf E` | exact | polarisation regime |
|---|---|---|---|---|
| 0.1 V/cm | 0.170 MHz | 0.085 MHz | 2.646 MHz | quadratic (unpolarised) |
| 1 V/cm | 1.697 | 0.848 | 2.778 | crossover |
| 3.1 V/cm | 5.26 | 2.645 | 3.74 | `γ_F m_F d_mf E = ω_ef/2` — the crossover field |
| 24 V/cm | 40.72 | 20.358 | 20.529 | polarised (Ng 2022 Fig. 4 field) |
| 60 V/cm | 101.79 | 50.895 | 50.964 | polarised (JILA E_rot) |

Leading correction to the fully-polarised limit: `+(ω_ef J(J+1)/4)² / (2 γ_F m_F d_mf E)` = 69 kHz at 60 V/cm `[derived]`. Quadratic-Stark coefficient in the weak-field limit: `(m_F γ_F d_mf)²/(ω_ef J(J+1)/2)` = 0.136 MHz/(V/cm)² for J = 1, F = 3/2, |m_F| = 3/2 `[derived]`.

Cross-check against the source: `(3/2) γ_{F=3/2} d_mf E_rot` at 60 V/cm = **50.895 MHz**, against Ng thesis Table B.2's "2π × 50.9 MHz" ✓; `d_mf E_rot` = 101.79 MHz against Table B.2's "2π × 102 MHz" ✓.

**ΔJ = ±1 Stark (J-mixing).** Diagonalising the same Hamiltonian in a J = 1-only basis and in a J ≤ 4 basis at 60 V/cm, the lowest level moves by **53 kHz** `[derived]`; this is not a nuisance, it is the mechanism that generates δg (§2.9) and it is why Ng needs 32 states rather than 12 (his footnote 1, thesis p. 85: "coupling with the J = 2 states through the Stark Hamiltonian modifies the g-factor significantly").

**Electronic polarizability**: `−½αE²` with α ~ 20 a.u. gives ~10⁻⁸ Hz at 60 V/cm `[derived]`. Twelve orders below the floor; drop without comment.

> **Origin dependence — ThF⁺ is an ion, so d_mf is origin-dependent and this is a 20 % trap.** The value 3.37(9) D and the ab initio 3.46 D are referenced to the **centre of nuclear mass**. Skripnikov & Titov 2015 Table II tabulates 2.74 D **with respect to the Th nucleus** (their footnote, p. 8); the difference is `e·r(Th→c.m.) = 0.72 D` and 2.74 + 0.72 = 3.46 D, exactly the value Ng attributes to them `[derived, confirms the conversion]`. The code must assert its origin.

---

### 2.8 Zeeman — and a sign in Ng's Eq. C.6 that does not survive its own consistency check

**Operator, as printed** (Ng thesis Eq. C.6, p. 321):
`H_Zeeman = −G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B`

**The relation Ng prints three lines below it** (same page; identical to Petrov & Skripnikov arXiv:2503.02840 Eqs. 2–3 and Petrov et al. arXiv:1704.06631):

```
g_F = −G∥ [F(F+1) + J(J+1) − 3/4] / [2F(F+1)J(J+1)]
      + g_N (μ_N/μ_B) [F(F+1) − J(J+1) + 3/4] / [2F(F+1)]
    ≡ −G∥ γ_F + g_N (μ_N/μ_B) κ_F
```

with the Zeeman shift defined as `E = −g μ_B B m_F` (Petrov arXiv:2503.02840 Eq. 16; "This definition matches the ones in [Ng 2022] and [Petrov 2017]").

**These two are inconsistent by one sign, and the operator is the one that is wrong.** `[derived]`

- `(J·n̂)(n̂·B) = Ω × B ⟨n̂_z⟩`, and `⟨n̂_z⟩ = Ω m_F γ_F` in this basis (I verified this numerically: the D-matrix structure in Ng's own C.5/C.6 evaluates to exactly `Ω m_F γ_F`, which is also what makes H_Stark reduce to Leanhardt Eq. 20). Hence `(J·n̂)(n̂·B) → Ω² m_F γ_F B = m_F γ_F B`.
- Ng's printed operator therefore gives `E = −G∥ μ_B B m_F γ_F ≡ −g μ_B B m_F` ⟹ `g = +G∥ γ_F`, the **opposite sign** to his own formula.
- Ng's quoted numbers settle it. He writes: "converting |g_{F=3/2}| into G∥, we get −0.042(2) if g_{F=3/2} > 0 and 0.048(2) otherwise" (thesis p. 87). Solving `g = −G∥/3 + g_N(μ_N/μ_B)/3` gives G∥ = −0.0418 for g = +0.0149 and G∥ = +0.0476 for g = −0.0149 `[derived]` — reproducing both of his printed values. The formula, not the operator, is what he used.
- Numerical confirmation: with `H = +G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B` and G∥ = 0.04756, diagonalising the J = 1 block at B = 1 G gives adjacent-m_F spacings of **20.85 kHz**, matching the closed form `|(1/3)(−G∥ + g_N μ_N/μ_B)| μ_B = 20.853 kHz/G` to 4 digits. With Ng's printed sign the model gives 23.52 kHz and |g| = 0.0168 — irreconcilable with the measurement.

**Adopt: `H_Zeeman = +G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B`, with `G∥ = (1/Ω)⟨Ψ| L̂ᵉ_n̂ − g_S Ŝᵉ_n̂ |Ψ⟩`, `g_S = −2.0023`** (Skripnikov & Titov 2015 Eq. 11, p. 3; identically arXiv:2302.02856 Eq. 10). Skripnikov notes G∥ "is close to zero for the ³Δ₁ state (and equal to zero when both the scalar-relativistic approximation is applied and the radiation corrections to the free-electron g-factor are ignored)" — for Λ = 2, Σ = −1 the naive value is `2 − 2.0023 = −0.0023`, twenty times smaller than the measured 0.048 and of the opposite sign. G∥ is a small difference of large numbers; treat it as a fitted parameter, not a derived one. **OPEN-3** records the alternative reading (that Ng's operator is right and his formula's relative sign is wrong), which the measurement excludes but which someone should confirm with the JILA authors.

**Matrix elements.** The G∥ term has exactly the Stark structure of §2.7 with `−d_mf E_p → +G∥ μ_B B_p Ω`; Ng writes both out and notes the similarity himself. The nuclear term uses B&C's spectator theorem (5.175) on the I part:

```
⟨…|T¹_{p}(I)|…⟩ = (−1)^{F′+J+I+1} √(I(I+1)(2I+1)) { I  F′ J ;  F  I  1 }
                    × (−1)^{F′−m′_F} √((2F+1)(2F′+1)) ( F′ 1 F ; −m′_F  p  m_F )
```

**Selection rules**: G∥ term — ΔΩ = 0, ΔJ = 0, ±1, ΔF = 0, ±1, Δm_F = 0 for B ∥ ẑ, **parity-even** (it is quadratic in n̂). Nuclear term — ΔJ = 0, ΔΩ = 0, ΔF = 0, ±1, Δm_F = 0, parity-even.

**Structural consequence worth a test (§6, V6): at leading order the Zeeman shift is even in Ω** (the two factors of Ω multiply out), so `g^u = g^ℓ` exactly in this term. All of δg comes from higher order. This is precisely why the eEDM term (odd in Ω) separates from the Zeeman background in the four-way chop.

**g_F values and sizes** `[derived, G∥ = 0.04756, g_N = 5.25773, μ_N/μ_B = 1/1836.15267, μ_B = 1.3996245 MHz/G]`:

| J | F | γ_F | κ_F | g_F | g_F μ_B (kHz/G) |
|---|---|---|---|---|---|
| 1 | 3/2 | 1/3 | 1/3 | −0.015046 | −21.06 |
| 1 | 1/2 | 2/3 | −1/3 | −0.032954 | −46.12 |
| 2 | 5/2 | 2/15 | 1/5 | −0.005827 | −8.16 |
| 2 | 3/2 | 1/5 | −1/5 | −0.010173 | −14.24 |
| 3 | 7/2 | 1/14 | 1/7 | −0.003020 | −4.23 |
| 3 | 5/2 | 2/21 | −1/7 | −0.004980 | −6.97 |
| 4 | 9/2 | 2/45 | 1/9 | −0.001815 | −2.54 |
| 4 | 7/2 | 1/18 | −1/9 | −0.002985 | −4.18 |

`g_{F=1/2} = 2 g_{F=3/2}` holds exactly in the G∥-only limit (Leanhardt Eq. 24, since γ_{1/2} = 2γ_{3/2}); the nuclear term breaks it — the actual ratio is 2.19, a 10 % deviation `[derived]`. That deviation is the cleanest experimental handle on the nuclear contribution, if the F = 1/2 g-factor is ever measured.

Sanity checks at the JILA operating point: `3 g_{F=3/2} μ_B B_rot` at B_rot = 799 µG = **50.0 Hz**, against Ng Table B.2's "2π × 50 Hz" ✓. Nuclear Zeeman alone: `g_N μ_N = 4.008 kHz/G` — the ¹⁹F nuclear moment is a **kHz-scale effect at 1 G** and must be kept.

**Sign of g_F is not measured.** Ng 2022 §II E and thesis p. 87 both say the spectroscopy is insensitive to it; theory (G∥ > 0 from Skripnikov 0.034, Cheng 0.035, Petrov 0.047) forces g_F < 0, and Petrov's Fig. 2 shows g^{u,ℓ} ≈ −1.49 × 10⁻² directly. Expose the sign as a switch. OPEN-4.

---

### 2.9 The differential g-factor δg / Δg, and its E dependence

**Two definitions in circulation, differing by exactly 2.** Ng thesis App. C.4 (p. 324): `δg = (g^u − g^ℓ)/2`. Petrov & Skripnikov arXiv:2503.02840 Eq. 18: `Δg = g^u − g^ℓ`. Both are self-consistent; a code that mixes them is 100 % wrong on a systematic-error budget. Pin one with a comment.

**Three physical contributions, in the order they dominate:**

1. **E-field-induced (dominant at the operating field).** Stark J-mixing gives the J = 1 upper and lower Stark doublets slightly different admixtures of J = 2, whose γ_F is smaller. Leanhardt Eq. 66 (the full 3j/6j sum, p. 24) reduces for J = 1, F = 3/2 to Eq. 67: `|δg_{F=3/2}/g_{F=3/2}| = 9 d_mf E_rot / (40 B_e)` = 0.00315 at 60 V/cm `[derived]`, against a measured 0.00255(6). The sign is negative in the JILA data and in every model.
2. **Zero-field parity-dependent Zeeman.** Leanhardt Eq. 64 (p. 24), from Brown et al. and Nelis et al., is the ³Δ analogue of the parity-dependent terms (vi) and (vii) of B&C's complete Zeeman Hamiltonian Eq. (9.70) (PDF p. 652 / book p. 620): `H_ZeemanDist = −½ g_rS μ_B (B₊J₋S₊S₋ + B₋J₊S₋S₊)` (parity-independent) and `H_ZeemanDoub = +½ g′_rS μ_B (B₊J₊S₊² + B₋J₋S₋²)` (parity-dependent), with `|g′_rS| ≈ ω_ef/(2B_e)` (Eq. 65) = 3.6 × 10⁻⁴ and a zero-field g-difference "twice the value in Eq. 65" ≈ 7.3 × 10⁻⁴ `[derived]`. Petrov's ab initio zero-field value for ThF⁺ is `Δg(E=0) = 2.3 × 10⁻⁴` (arXiv:2503.02840 p. 4), the same order, 3× smaller. This term is **absent from Ng's Appendix C** and from any Ω = ±1 effective model that lacks it.
3. **Rotating-field contribution** from finite ω_rot/(d_mf E_rot): Leanhardt Eq. 68, `δg/g = √6 γ²_{F=3/2} ω²_rot/(d_mf E_rot E_hf)` — a few 10⁻⁴; only relevant with a rotating field.

**Measured and modelled values, all at E = 60 V/cm:**

| quantity | value | source |
|---|---|---|
| δg/g (measured) | −0.00255(6) | Ng thesis Table 4.1 p. 87, Fig. 4.9b |
| δg/g (Ng 32-level model) | −0.00223 | Ng thesis p. 85 |
| δg/g (this document's 96-state model) | **−0.00223** | `[derived]`, §6 V7 |
| Δg = g^u − g^ℓ (Petrov, ³Δ₂ included) | +7.56 × 10⁻⁵ ⟹ δg/g = −0.00254 | arXiv:2503.02840 Table I |
| \|δg\| (Ng 2022 paper, Table I) | 0.0003(3) | Ng 2022 p. 5 — note the *paper* uses the unhalved convention, so this is Petrov's Δg |

**Petrov's Δg₀, Δg₁ are a local linearisation, not global constants.** `Δg = Δg₀ + Δg₁E` (Eq. 18) with Δg₀ falling from 349.4 × 10⁻⁷ at 40 V/cm to 93.8 × 10⁻⁷ at 150 V/cm, and Δg₁ rising from 6.3 to 10.4 × 10⁻⁷ cm/V (Table I, p. 5). Evaluate the pair at the field you are working at; do **not** read Δg₀ as the E = 0 value (which is 2.3 × 10⁻⁴).

**Size in energy**: at 1 G, `Δg μ_B B = 0.106 kHz` per unit m_F (using Δg = 7.56 × 10⁻⁵); at B_rot = 799 µG, 0.08 Hz.

---

### 2.10 Rotational g-factor g_r

Not measured or computed for ThF⁺. Ng 2022 (p. 5): "It might be of interest to compute the rotational g-factor of ThF⁺, because the rotational g-factor has been shown to contribute to about 6 % of the total g-factor of a similar molecular species" (ThO, Petrov et al. PRA 89, 062505). His G∥ is extracted "neglecting the rotational contribution", so the measured |g_{F=3/2}| already contains g_r.

**Operator**: B&C Eq. (9.70) term (iii), `−g_r μ_B B_Z T¹_{p=0}(J − L − S)`. In the J = 1 g-factor it enters as Leanhardt Eq. 22, `γ_F[((g_L+g_r)Λ + (g_S+g_r)Σ)Ω − g_r J(J+1)]`, whose net g_r piece at J = 1, F = 3/2 is `−g_r/3` `[derived]`.

**Estimate `[estimate, method: rigid-rotor nuclear-charge contribution only]`.** For a diatomic AB about the centre of mass, `g_r^nuc = m_p [Z_A m_B² + Z_B m_A²] / (M m_A m_B)`. With Z(Th) = 90, Z(F) = 9, m = 232.0381 and 18.99840 u:

`g_r^nuc = 0.471 μ_N = 2.56 × 10⁻⁴ μ_B` `[derived]`

⟹ contribution to g_{J=1,F=3/2} of `−8.5 × 10⁻⁵`, i.e. **0.6 % of |g_F| = 0.0149**. Energy: 0.36 kHz at 1 G (m_F = ±3/2 splitting), 0.29 Hz at B_rot.

**Uncertainty is large and one-sided.** The electronic (paramagnetic) contribution has the opposite sign and in light molecules largely cancels the nuclear part; in heavy molecules with low-lying states of the right symmetry it can dominate and exceed it. Treat 2.6 × 10⁻⁴ μ_B as an order-of-magnitude floor, not a value. The ThO comparison (6 %) is against a much smaller total g, so it is consistent with this estimate. **OPEN-7.**

The practical consequence: fitting G∥ to reproduce the measured g at J = 1, F = 3/2 gives the right answer there by construction, but g_r has a different (J, F) dependence, so the **predicted g at J = 2, 3, 4 carries an unquantified ~1 % error.** Say so wherever the code reports a J > 1 g-factor.

---

### 2.11 What is absent, and why

- **²³²Th quadrupole and ²³²Th magnetic hyperfine**: I(²³²Th) = 0. Both identically zero. (Skripnikov & Titov tabulate A∥ = −4163 μ_Th/μ_N MHz and W_M = 0.88 × 10³³ Hz/(e·cm²) — those are **²²⁹Th** quantities for the MQM programme and must never enter a ²³²Th row.)
- **¹⁹F quadrupole**: I = ½, no quadrupole moment. Zero.
- **Nuclear spin–spin**: one spin only. Zero.
- **Electron spin–spin λ and spin–rotation γ_SR** (Leanhardt Eqs. 13–14): diagonal in Ω and constant within the Ω = ±1 block; absorbed into the origin. λ ≈ 261 cm⁻¹ `[derived, §1.4]` matters only for T2.
- **Hyperfine-induced Ω-doublet asymmetry, e_Δ** (Leanhardt Eq. 15, the `½e_Δ(J₊I₊S₊² + J₋I₋S₋²)` term): "a previously unreported term ... expected to be even smaller than the already small Λ-doublet splitting itself, however, and will be ignored" (Leanhardt p. 13). It makes ω_ef depend on F. Scaling `e_Δ/õ_Δ ~ A∥/B₀ = 2.8 × 10⁻³` (replacing one rotational `B J₊` by a hyperfine factor) gives `e_Δ ~ 3.7 kHz` and a hyperfine dependence of the Ω-doublet splitting of **order 1–10 kHz** `[estimate, method: ratio of hyperfine to rotational constant applied to õ_Δ = ω_ef/4]`. That is at the v1 threshold, with an order-of-magnitude uncertainty. **OPEN-8.**
- **Off-diagonal-in-Ω hyperfine** (ΔΩ = ±1, connecting ³Δ₁ to ³Δ₂ and to the Ω = 0 states): suppressed by hyperfine/spin-orbit, "hence a factor of 10⁻⁶" (Leanhardt p. 14). With A∥ = 20 MHz and ΔE = 1052 cm⁻¹ the ratio is 6 × 10⁻⁷ `[derived]`, so amplitude ~10⁻⁶ and energy ~10⁻⁵ Hz. Drop.
- **Vibrational structure**: out of scope for v = 0. Record only that `B₀ = B_e − α_e/2 = 0.24261 cm⁻¹` from Gresh's B_e = 0.24311(7), α_e = 1.00(4) × 10⁻³, against the directly measured `B₀ = 29.09733(4)/4 GHz = 0.2426456 cm⁻¹` — agreement to 1.5 × 10⁻⁴ relative `[derived]`. Use the measured 4B.
- **Second-order / perpendicular Zeeman**: `±(3/4)(g_F μ_B B⊥)²/(γ_F d_mf E_rot)` (Leanhardt Eq. 72). At B⊥ = 1 mG and E_rot = 60 V/cm: (0.0209 kHz)²·0.75/33.93 MHz ≈ 7 × 10⁻⁹ Hz `[derived]`. Negligible unless B⊥ reaches ~1 G.

---

### 2.12 PT-odd terms

**eEDM.** Ng thesis Eq. C.8 (p. 322): `H_eEDM = −d_e E_eff Ω/|Ω|`. In code this is an **added diagonal-in-Ω, odd-in-Ω matrix**: `+d_e E_eff` on Ω = −1, `−d_e E_eff` on Ω = +1. It is parity-odd and T-odd; in the parity basis it is purely off-diagonal within each Ω-doublet, structurally identical to the ΔJ = 0 Stark element.

> **The factor of 2 — pin it.** Leanhardt Eq. 26 (p. 15): `E_EDM = −d_e E_eff Ω/(2|Ω|)`. Ng removes the ½, with footnote 3 on p. 322, verbatim: "There is an extra factor of 1/2 in the cited paper [Leanhardt], which uses a different convention for how E_eff is defined in papers that calculate the value. The factor of 1/2 is removed for consistency with the cited value of E_eff." **Use Ng's form** (no ½) with a calculators' E_eff. The observable relation everyone now uses is `f^BD = 2 d_e E_eff` (arXiv:2302.02856 Eq. 3; arXiv:2503.02840 p. 4).

**Definitions from the calculators** (Skripnikov & Titov 2015, PRA 91, 042504, Eqs. 1–4, p. 2, verbatim):
`W_d = (1/Ω)⟨Ψ|Σ_i H_d(i)/d_e|Ψ⟩` with `H_d = 2d_e[[0,0],[0,σ·E]]`, `E_eff = W_d |Ω|`, and `Ω = ⟨Ψ|J·n|Ψ⟩` **with n directed from Th to F** (Ω = +1 for ³Δ₁).
`W_{T,P} = (1/Ω)⟨Ψ|Σ_i H_{T,P}(i)/k_{T,P}|Ψ⟩` with `H_{T,P} = i(G_F/√2)Z k_{T,P} γ⁰γ⁵ ρ_N(r)`; the tabulated W_{T,P} is in kHz, so the energy shift is `W_{T,P} k_{T,P}`.

**Scalar–pseudoscalar.** Operator form in the code is the same added matrix as the eEDM with `d_e E_eff → W_{T,P} k_{T,P}`; the two are experimentally degenerate in a single species.

**MQM (W_M).** Requires I ≥ 1 on the heavy nucleus. Zero for ²³²Th. Skripnikov's `H^MQM_eff = −W_M M/(2I(2I−1)) S T̂ n` (their Eq. 8) applies to ²²⁹ThF⁺ only.

**W_a / W_P (nuclear-spin-dependent PV):** **not found** for ThF⁺ in any source read. Searched Skripnikov & Titov 2015 (full), Denis 2015, Petrov 2025, Petrov 2023, Ng 2022, Ng thesis App. C. OPEN-9.

**Sizes** `[derived]`:

| quantity | value |
|---|---|
| d_e E_eff at d_e = 10⁻³¹ e·cm, E_eff = 35 GV/cm | 0.846 µHz (Ng Table B.2 prints 0.851, using 35.2) ✓ |
| d_e E_eff at the current bound \|d_e\| < 4.1 × 10⁻³⁰ e·cm (Roussy 2023, HfF⁺) | 34.7 µHz |
| f^BD = 2 d_e E_eff at that bound | 69.4 µHz |
| W_{T,P} k_{T,P} at k_{T,P} = 10⁻⁹ | 50 µHz |

Ten orders of magnitude below the kHz floor. They are in the term list because they are the observable, not because they are large.

---

### 2.13 Rotating-frame term (context, not part of the static Hamiltonian)

`H_rot-frame = ℏ ω_rot F_x` (Ng thesis Eq. C.4, p. 320), with ẑ along E_rot and **x̂ anti-parallel to the (counter-clockwise) rotation vector**. Couples Δm_F = ±1 within the same F, J, Ω, breaks per-M_F blocking, and is what generates Berry's phase and the |m_F| = 3/2 avoided crossings `Δ^u, Δ^ℓ` (third order in ω_rot/d_mf E_rot; Leanhardt Eq. 42: `Δ ≈ 170 ω_ef (ω_rot/d_mf E_rot)³`). Ng thesis Table B.2 lists `Δ^u = 2π × 0.511 Hz`, `Δ^ℓ = 2π × 1.43 Hz` at E_rot = 60 V/cm, ω_rot = 2π × 147.5 kHz — model output, not measurement.

Out of scope for a static-field v1, but note it if the code ever claims to reproduce a JILA Ramsey frequency.

---

## 3. Parameter table

Isotopologue ²³²Th¹⁹F⁺, X ³Δ₁, v = 0 throughout. Conversion factors used: 1 cm⁻¹ = 29 979.2458 MHz; 1 D = 0.5034118 MHz/(V/cm); μ_B = 1.3996245 MHz/G; μ_N = μ_B/1836.15267 = 0.7622593 kHz/G.

| Symbol | Value ± unc. as printed | Units as printed | Converted | Source or estimate method | Status |
|---|---|---|---|---|---|
| **B₀** | 4B = 29.09733(4) | GHz | B₀ = 7.2743325(10) GHz = 0.2426456 cm⁻¹ `[29.09733/4]` | Ng 2022 §II B p. 2; Ng thesis §4.1.1 p. 78 | measured |
| B₀ (independent) | 0.24264(3) | cm⁻¹ (95 % CI) | 7274.1(9) MHz `[×29979.2458]` | Gresh 2016 Table 1, X ³Δ₁ rows (rendered) | measured |
| B_e | 0.24311(7) | cm⁻¹ | 7288.2(21) MHz | Gresh 2016 Table 2 p. 11 | measured |
| α_e | 1.00(4) × 10⁻³ | cm⁻¹ | 29.98(12) MHz | Gresh 2016 Table 2 | measured |
| **D₀** | 1.30(4) × 10⁻⁷ | cm⁻¹ | 3.897(120) kHz `[×29979.2458]` | Gresh 2016 Table 1, X ³Δ₁ (D″ column; other bands 1.27–1.40 × 10⁻⁷) | measured |
| D₀ (cross-check) | 4B_e³/ω_e² = 1.332 × 10⁻⁷ | cm⁻¹ | 3.993 kHz | `[derived]`, Kratzer relation from B_e, ω_e | derived |
| **ω_ef** | 5.29(5) | MHz (as ω_ef/2π) | 5.29(5) MHz; splitting at J is ω_ef J(J+1)/2 | Ng 2022 Table I p. 5; Ng thesis Table 4.1 p. 87 | measured |
| ω_ef (independent) | \|k″\| = 0.869–0.892 × 10⁻⁴ | cm⁻¹ | ω_ef = 2\|k\| = 5.21–5.35 MHz `[×2×29979.2458]` | Gresh 2016 Table 1 (rendered; 10 X ³Δ₁ bands) | measured |
| ω_ef (2nd-order prediction) | 5.294 | MHz | — | `[derived]`, §2.3, from Petrov 2025 Eqs. 10, 14, 15 (partly circular: Δ⁽¹⁾ was fitted to ω_ef) | derived |
| **A∥** | −20.1(1) | MHz (as A∥/2π) | −20.1(1) MHz | Ng 2022 Table I p. 5; Ng thesis Table 4.1 p. 87 | measured |
| A∥ (theory) | −21.5 | MHz | — | Ng 2022 Table I (X2CAMF-CCSD(T), L. Cheng); no uncertainty quoted | ab initio |
| E_hf(J=1) = ¾A∥ | −15.1 | MHz (2π×) | −15.075 MHz `[¾ × 20.1]` | Ng thesis Table B.2 p. 317 | derived from measured |
| **d_mf** | 3.37(9) | D | 1.6965 MHz/(V/cm) `[3.37 × 0.5034118]`; Ng Table B.2 prints 1.70 ✓ | Ng 2022 Table I p. 5; **origin = centre of nuclear mass** | measured |
| d_mf (theory, Skripnikov) | 2.74 (w.r.t. **Th nucleus**) | D | 3.46 D at c.m. `[+0.72 D = 4.80320 D/Å × 0.1499 Å]` | Skripnikov & Titov 2015 Table II p. 8 + footnote | ab initio |
| d_mf (theory, Cheng) | 3.46 | D | — | Ng 2022 Table I | ab initio |
| d_mf (theory, Denis) | 4.03 | D (c.m.) | — | Denis 2015 abstract, Table 10 — 16 % above experiment | ab initio |
| D (Petrov's signed body-frame dipole) | −0.133 **[probable typo for −1.33]** | a.u. | −1.33 a.u. = −3.381 D `[×2.541746]` matches d_mf; −0.133 would be 8× too small | arXiv:2503.02840 p. 3; sign is Petrov's n̂ (Th→F) | ab initio input |
| **\|g_{F=3/2}\|** | 0.0149(3) | — | \|g\|μ_B = 20.85 kHz/G | Ng 2022 Table I p. 5 — **sign not measured** | measured |
| **G∥** | 0.048(2) if g_F < 0; −0.042(2) if g_F > 0 | — | 0.04756 reproduces \|g\| = 0.0149 exactly `[derived]` | Ng 2022 §II E p. 5; Ng thesis p. 87; rotational contribution neglected | derived from measurement |
| G∥ (theory) | 0.034 / 0.035 / 0.047 | — | — | Skripnikov & Titov 2015 Table II / Cheng in Ng 2022 / Petrov 2025 p. 3 (with non-adiabatic mixing) | ab initio |
| g_N(¹⁹F) | 5.25773 | — | g_N μ_N = 4.008 kHz/G; g_N μ_N/μ_B = 2.8634 × 10⁻³ | Petrov et al. arXiv:1704.06631 Eq. 3 discussion | constant |
| δg/g at E_rot = 60 V/cm | −0.00255(6) | — | δg = (g^u−g^ℓ)/2 = +3.80 × 10⁻⁵ | Ng thesis §4.1.3 p. 85, Table 4.1 p. 87 | measured |
| Δg₀, Δg₁ at 60 V/cm | 233.8 × 10⁻⁷, 8.7 × 10⁻⁷ cm/V | — | Δg = 7.558 × 10⁻⁵ `[233.8e-7 + 8.7e-7×60]` | Petrov & Skripnikov arXiv:2503.02840 Table I p. 5 (local linearisation; full table 40–150 V/cm) | ab initio |
| Δg at E = 0 | 2.3 × 10⁻⁴ | — | — | arXiv:2503.02840 p. 4; "accuracy not worse than 10 %" | ab initio |
| **E_eff** | 37.3 (±7 %) | GV/cm | — | Skripnikov & Titov 2015 Table II p. 8, uncertainty p. 7 | ab initio |
| E_eff | 35.2 | GV/cm | — | Denis 2015 abstract | ab initio |
| E_eff (JILA adopts) | ≈ 35 | GV/cm | — | Ng thesis Eq. C.8 discussion p. 322 | ab initio |
| W_{T,P} | 50 (±7 %) | kHz | shift = W_{T,P} k_{T,P} | Skripnikov & Titov 2015 Table II p. 8 | ab initio |
| W_{T,P} | 48.4 | kHz | — | Denis 2015 abstract | ab initio |
| **c_I (¹⁹F spin–rotation)** | ~20 (factor 3) | kHz | splitting c_I(2J+1)/2 = 30–90 kHz for J = 1–4 | `[estimate, method: B&C Table 8.12 PDF p.513 gives c₂(¹⁹F, CsF) = 15.1 kHz at B_v = 0.183782 cm⁻¹; scaled linearly in B to ThF⁺]` | **estimate** |
| **g_r (rotational)** | ~0.47 μ_N = 2.6 × 10⁻⁴ μ_B (floor) | — | contributes −g_r/3 = −8.5 × 10⁻⁵ to g_{J=1,F=3/2} | `[estimate, method: rigid-rotor nuclear-charge term g_r^nuc = m_p(Z_Th m_F² + Z_F m_Th²)/(M m_Th m_F); electronic part omitted]` | **estimate** |
| **e_Δ (hyperfine Ω-doubling)** | ~4 kHz (order of magnitude) | kHz | ω_ef acquires an F dependence of order 1–10 kHz | `[estimate, method: e_Δ/õ_Δ ~ A∥/B₀ = 2.8 × 10⁻³ applied to õ_Δ = ω_ef/4]` | **estimate** |
| A (³Δ spin–orbit) | ≈ 787.5 | cm⁻¹ | — | `[derived]` from T₀(³Δ₂)−T₀(³Δ₁) = 1052.5 and T₀(³Δ₃)−T₀(³Δ₂) = 2097.5 (Gresh 2016 Table 2) with E = 2AΣ + (2/3)λ(3Σ²−S²) | derived, T2 only |
| λ (spin–spin) | ≈ 261 | cm⁻¹ | — | same fit `[derived]` | derived, T2 only |
| T₀(a ¹Σ⁺) | 314.282(7) | cm⁻¹ | 9.4224 THz | Gresh 2016 Table 2 p. 11 | measured |
| T₀(1³Δ₂) | 1052.5(1.0) | cm⁻¹ | — | Gresh 2016 Table 2 (from Barker et al.) | measured |
| T₀(³Δ₃) | 3150(30) | cm⁻¹ | — | Gresh 2016 Table 2 | measured |
| T₀(³Π₀₋), T₀(³Π₀₊) | 3044, 3395 | cm⁻¹ | — | arXiv:2503.02840 Eq. 3 p. 2 (vertical, R = 3.75 a₀) | ab initio |
| Δ⁽¹⁾, Δ⁽³ᵃ⁾, Δ⁽³ᵇ⁾ (Ω = 0 rotational couplings) | 0.250, 0.587, 0.555 | cm⁻¹ | — | arXiv:2503.02840 Eqs. 10, 14, 15; Δ⁽¹⁾ scaled +8.8 % to fit ω_ef | ab initio, T3 only |
| Δ⁽²ᵃ⁾ (³Δ₁–1³Δ₂) | −0.510 | cm⁻¹ | — | arXiv:2503.02840 Eq. 11 | ab initio, T2 only |
| ω_e, ω_eχ_e | 656.96(1), 1.920(3) | cm⁻¹ | — | Gresh 2016 Table 2 p. 11 | measured |
| R_e | 3.75 (theory), 3.74(4) (exp) | a₀ | 1.984 Å | Skripnikov & Titov 2015 Table I p. 7 | ab initio / measured |

**HfF⁺ values encountered and deliberately not transferred**: A∥(¹⁹F) = −62.0 MHz, G∥ = 0.011768, D∥ = −1.53(2) a.u., B_rot = 0.2989 cm⁻¹, E_eff = 22.5(9) GV/cm, ³Δ₁ T_e = 976.930 cm⁻¹, ³Δ₂ T_e = 2149.432, ³Π₀₋ = 10212.623, ³Π₀₊ = 10401.723 cm⁻¹ (arXiv:2302.02856 Eqs. 9–13, Table I; arXiv:1704.06631 Eq. 9). The **only** HfF⁺-derived quantity used anywhere in this document is the current eEDM bound |d_e| < 4.1 × 10⁻³⁰ e·cm, used solely to size the PT-odd row.

---

## 4. Energy hierarchy for J = 1–4, largest to smallest

Two field points: **E = 1 V/cm, B = 1 G** (a generic laboratory scale) and **E_rot = 60 V/cm, B_rot = 799 µG** (the JILA Gen. III point, Ng Table B.2). All entries `[derived]` from §2 and §3.

| # | Effect | J = 1 | J = 2 | J = 3 | J = 4 | notes |
|---|---|---|---|---|---|---|
| 1 | Rotation B J(J+1) | 14.549 GHz | 43.646 GHz | 87.292 GHz | 145.487 GHz | sets J spacing |
| 2 | Stark, polarised, at 60 V/cm (max \|shift\|) | 50.895 MHz | 33.930 MHz | 25.447 MHz | 20.358 MHz | `γ_F m_F d_mf E`; max `m_F γ_F = 1/(J+1)` exactly `[derived]` |
| 3 | Ω-doubling ω_ef J(J+1)/2 | 5.29 MHz | 15.87 MHz | 31.74 MHz | 52.90 MHz | grows with J; overtakes the 60 V/cm Stark shift at J = 3 |
| 4 | Hyperfine A∥ splitting | 15.075 MHz | 8.375 MHz | 5.863 MHz | 4.523 MHz | F = J−½ above F = J+½ |
| 5 | Centrifugal D [J(J+1)]² | 15.6 kHz | 140 kHz | 561 kHz | 1.559 MHz | |
| 6 | Stark at 1 V/cm (stretched, exact 2-level shift) | 132.7 kHz | 20.1 kHz | 5.7 kHz | 2.2 kHz | quadratic/crossover regime; crossover field `ω_ef J(J+1)(J+1)/(4 d_mf)` = 3.12 V/cm at J = 1 |
| 7 | 2nd-order Stark J-mixing at 60 V/cm | 53 kHz | 53 kHz | — | — | measured by diagonalising the J = 1-only vs J ≤ 4 basis `[derived]`; source of δg |
| 8 | Correction to full polarisation at 60 V/cm | 69 kHz | — | — | — | (ω_ef J(J+1)/4)²/(2γ_F m_F d_mf E) |
| 9 | Zeeman g_F μ_B B at 1 G (per unit m_F) | 21.1 kHz (F=3/2) / 46.1 (F=1/2) | 8.2 / 14.2 | 4.2 / 7.0 | 2.5 / 4.2 | |
| 10 | Nuclear spin–rotation c_I ≈ 20 kHz | 30 kHz | 50 kHz | 70 kHz | 90 kHz | **estimate**; grows with J |
| 11 | ¹⁹F nuclear Zeeman g_N μ_N B at 1 G | 4.01 kHz | 4.01 | 4.01 | 4.01 | J-independent |
| 12 | Hyperfine ΔJ = ±1 (B&C 9.51), 2nd order | 2.60 kHz | 2.60 / 2.06 | 2.06 / 1.63 | 1.63 kHz | omitted by Ng App. C |
| 13 | e_Δ hyperfine Ω-doubling | ~1–10 kHz | ~ | ~ | ~ | **estimate**, order of magnitude only |
| 14 | Zero-field Δg (parity-dependent Zeeman) at 1 G | 0.97 kHz (3 m_F) | — | — | — | Δg(0) = 2.3 × 10⁻⁴, Petrov |
| 15 | Δg at 60 V/cm, 1 G | 0.32 kHz (3 m_F) | — | — | — | Δg = 7.56 × 10⁻⁵ |
| 16 | Rotational g_r at 1 G (3 m_F) | 0.36 kHz | — | — | — | **estimate**, floor |
| 17 | Zeeman at B_rot = 799 µG (3 g μ_B B, F = J+½) | 50.5 Hz | 19.5 | 10.1 | 6.1 Hz | JILA point; Ng Table B.2 prints "2π × 50 Hz" ✓ |
| 18 | ¹⁹F nuclear Zeeman at B_rot (g_N μ_N B) | 3.2 Hz | 3.2 | 3.2 | 3.2 Hz | |
| 19 | Δg at 60 V/cm, B_rot (3 m_F) | 0.25 Hz | — | — | — | |
| 20 | Ω = 2 admixture, non-absorbable part | 0 | 8.5 Hz | 21 Hz | 38 Hz | \|c\|² × (A∥ scale, 20 MHz); \|c\|² ≤ 1.9 × 10⁻⁶ (§1.2) |
| 21 | Perpendicular Zeeman, B⊥ = 1 mG at 60 V/cm | 7 × 10⁻⁹ Hz | — | — | — | Leanhardt Eq. 72 |
| 22 | d_e E_eff at the current bound | 34.7 µHz | — | — | — | f^BD = 2 d_e E_eff = 69.4 µHz |
| 23 | W_{T,P} k_{T,P} at k_{T,P} = 10⁻⁹ | 50 µHz | — | — | — | |
| 24 | Ω-doubling avoided crossings Δ^u, Δ^ℓ | 0.511, 1.43 Hz | — | — | — | rotating field only; model output |
| 25 | Electronic polarizability at 60 V/cm | ~10⁻⁸ Hz | — | — | — | drop |

**Where the 1 kHz line falls.** Rows 1–15 all reach or exceed ~1 kHz somewhere in J = 1–4 at the stated field points, and are in scope. Rows 16–19 sit below 1 kHz there, but each is the *same operator* as a row above it (16 ↔ 9, 18 ↔ 11, 19 ↔ 15) evaluated at a weaker field, so they cost nothing extra and must not be dropped by construction — only by field. Rows 20, 21, 24, 25 are genuinely droppable for a static-field v1. Rows 22–23 stay because they are the observable, not because they are large.

Note the ordering is field-dependent and crosses over: at 60 V/cm the Stark shift dominates the Ω-doubling for J ≤ 2 and is dominated by it for J ≥ 3; at 1 V/cm the Ω-doubling dominates the Stark shift at every J. Any code that assumes a fixed ordering of terms will be wrong somewhere in this range.

---

## 5. Recommended v1 term list

**Basis**: Ω = ±1, J = 1–4, F = J ± ½, all m_F. 96 states; block-diagonal in M_F (16/14/10/6/2 per signed block) for collinear static fields.

**Include (8 operators):**

| # | Term | Form | Parameters | Why |
|---|---|---|---|---|
| 1 | Rotation + centrifugal | `B₀J(J+1) − D₀[J(J+1)]²` | B₀, D₀ | GHz–MHz |
| 2 | Ω-doubling | Ng Eq. C.3 | ω_ef | MHz |
| 3 | Axial hyperfine, ΔJ = 0 | Ng Eq. C.2 = B&C 9.50 | A∥ | MHz |
| 4 | Axial hyperfine, ΔJ = ±1 | B&C 9.51 | A∥ (same parameter) | 1.6–2.6 kHz; free once (3) is coded |
| 5 | Stark | Ng Eq. C.5 | d_mf | MHz; also generates δg via J-mixing |
| 6 | Zeeman, G∥ | `+G∥ μ_B (J·n̂)(n̂·B)` — **sign per §2.8, not as Ng prints it** | G∥ | tens of kHz/G |
| 7 | Nuclear Zeeman | `−g_N μ_N I·B` | g_N (fixed) | 4 kHz/G |
| 8 | Nuclear spin–rotation | `c_I T¹(J)·T¹(I)`, B&C 8.7 / 8.20 | c_I (**estimate**) | 30–90 kHz; must be a knob, default 20 kHz, and it biases A∥ |

**Add as an opt-in PT-odd block** (not part of the structural Hamiltonian): `H_PT = −(d_e E_eff + W_{T,P} k_{T,P}) Ω/|Ω|`, following Ng Eq. C.8 (no Leanhardt ½).

**Drop, each with the size that justifies it:**

| Dropped | Size | Justification |
|---|---|---|
| Ω = ±2, ±3 (³Δ₂, ³Δ₃) | mixing \|c\| ≤ 1.4 × 10⁻³; energy effect linear in J(J+1), absorbed into B₀ | §1.2. Non-absorbable residual ≤ 38 Hz at J = 4. Exception: 15 % of δg — make it an opt-in tier |
| Ω = 0 perturbers (a¹Σ⁺, ³Π₀±) | generate ω_ef, which is fitted | §2.3. Needed only to *predict* ω_ef |
| Spin–orbit A, spin–spin λ, spin–rotation γ_SR | constant within Ω = ±1 | Leanhardt p. 13 |
| ²³²Th hyperfine, ²³²Th quadrupole, ¹⁹F quadrupole | identically 0 | I(²³²Th) = 0, I(¹⁹F) = ½ |
| Off-diagonal-in-Ω hyperfine | ~10⁻⁵ Hz | ratio hyperfine/spin-orbit = 6 × 10⁻⁷ |
| e_Δ hyperfine Ω-doubling | ~1–10 kHz **estimate** | **borderline** — dropped in v1, flagged OPEN-8; the only dropped term above the kHz line |
| Parity-dependent Zeeman (g_rS, g′_rS) | Δg(E=0) = 2.3 × 10⁻⁴ ⟹ ≈ 1 kHz at 1 G | **borderline** — dropped for level energies, required for a zero-field Δg. OPEN-10 |
| Rotational g_r | ≤ 0.6 % of g_F, 0.36 kHz at 1 G **estimate** | absorbed into the fitted G∥ at J = 1; introduces ~1 % error at J > 1 |
| Perpendicular / 2nd-order Zeeman | 7 × 10⁻⁹ Hz at B⊥ = 1 mG | Leanhardt Eq. 72 |
| Electronic polarizability | ~10⁻⁸ Hz | §2.7 |
| Rotating-frame `ℏω_rot F_x` | not a static-field term | §2.13; add when Berry phase or Δ^{u,ℓ} is wanted |
| Vibrational structure | 653 cm⁻¹ | v = 0 only |
| Higher-order rotational (H, k_D) | no published value; Gresh has k_D″ = 0 for all ³Δ₁ | fine to J = 4 |

---

## 6. Verification checks the code should implement

Each names the failure mode only it catches. V1–V6 and V8–V11 are Hamiltonian-agnostic or closed-form (both PASS and FAIL reachable, no pinned spectra). V7, V12, V13 are labelled **opt-in** comparisons to published numbers.

| id | check | passes iff | failure mode it uniquely catches |
|---|---|---|---|
| **V1** | Hermiticity and reality: `H = Hᵀ` for real E, B; eigenvalues real | exact | a transposed 3j/6j argument order, or a bra/ket swap in the spectator-theorem phase |
| **V2** | **Stark closed form.** For every (J ≤ 4, F, m_F, Ω), the ΔJ = 0, ΔF = 0, p = 0 Stark element equals `−Ω m_F γ_F d_mf E` with `γ_F = [J(J+1)+F(F+1)−I(I+1)]/[2F(F+1)J(J+1)]` | to machine precision | a wrong 6j column order or a wrong `(−1)^{J′−Ω′}` phase — both leave the *magnitude* right and only the (J,F)-dependence wrong. **I ran this; it passes.** |
| **V3** | **Hyperfine J-scaling.** Splitting = `A∥(2J+1)/(2J(J+1))`; J=1 : J=2 = 9/5 exactly | exact | using `F² − I² − J²` over `2J(J+1)` vs `2J²`, or mis-coupling I to J |
| **V4** | **Ω-doubling J-scaling.** Splitting = `ω_ef J(J+1)/2`; J=2 : J=1 = 3 exactly. And the **upper** component has parity `(−1)^J` (e in the B&C convention) at every J | exact | a missing `(−1)^J` — invisible at J = 1, wrong at every even J. Catches the S = 1 e/f trap (OPEN-2) |
| **V5** | **g-factor closed form.** With B alone, `g(J,F) = −G∥ γ_F + g_N(μ_N/μ_B) κ_F`; specifically `g(1,3/2) = (1/3)(−G∥ + g_N μ_N/μ_B)` and, with g_N set to 0, `g(1,1/2) = 2 g(1,3/2)` | to machine precision | **the Ng Eq. C.6 sign error of §2.8.** With the wrong sign the model gives 23.5 kHz/G instead of 20.85 kHz/G — a 13 % error that no other check sees |
| **V6** | **Zeeman is even in Ω.** At E = 0, B ≠ 0, and with the parity-dependent Zeeman terms off, `g^{Ω=+1} = g^{Ω=−1}` exactly; equivalently δg = 0 | exact | an odd-in-Ω contamination of the Zeeman operator, which would fake an eEDM signal in the four-way chop |
| **V7** | **δg from Stark J-mixing** *(opt-in)*. At E = 60 V/cm the Ω = ±1, J = 1–2 model gives δg/g = −0.0022 | within 5 % | Stark ΔJ = ±1 elements wrong in magnitude or relative phase. **I ran this: −0.00223, identical to Ng thesis p. 85's 32-level value.** Measurement is −0.00255(6); the 15 % gap is the known missing ³Δ₂ |
| **V8** | **Parity block structure.** At E = 0 and with the PT-odd block off, H commutes with `E* = σ_xz R_y(π)` and splits into two blocks of equal dimension | exact | a parity-odd term (Stark, eEDM) leaking into the field-free Hamiltonian |
| **V9** | **Kramers/time-reversal degeneracy.** At B = 0 (any E, no eEDM), `E(F, m_F, Ω) = E(F, −m_F, −Ω)` — Leanhardt Eq. 34, which he derives as exact | exact | a sign error in an m_F-odd term. This is the degeneracy the whole experiment measures against |
| **V10** | **eEDM parity.** Turning on d_e splits each Kramers pair by `2 d_e E_eff` with **opposite** sign for the upper and lower Stark doublets, and the Zeeman splitting is unchanged; `f^BD = 2 d_e E_eff` | exact | the Leanhardt ½ (§2.12) — a factor-2 error in the extracted d_e that nothing else catches |
| **V11** | **Zero-field degeneracy count.** At E = B = 0, each (J, F, parity) level is exactly (2F+1)-fold degenerate; total dimension 96 for J = 1–4; per-|M_F| block dimensions 16/14/10/6/2 | exact | a basis-enumeration bug (double-counting Ω, or admitting F outside \|J−I\| … J+I) |
| **V12** | **Kratzer** *(opt-in)*. `D₀ ≈ 4B_e³/ω_e² = 1.33 × 10⁻⁷ cm⁻¹` against the measured 1.30(4) × 10⁻⁷ | within 5 % | a wrong D₀ unit conversion (cm⁻¹ vs MHz vs kHz) |
| **V13** | **Ω-doubling from perturbers** *(opt-in, tier T3 only)*. `−(J(J+1)/4) Σ_i σ_i (Δ^(i))²/T_i` with σ = ±1 for 0± reproduces ω_ef = 5.29 MHz and the J(J+1) law | within 10 % | dropping the 0⁺/0⁻ reflection-parity sign, which turns a 5.29 MHz answer into 12.1 MHz. Note this check is partly circular (Δ⁽¹⁾ was fitted) |
| **V14** | **Traceless Stark.** `Σ_{m_F} E_Stark(J, F, Ω, m_F) = 0` at fixed J, F, Ω in the linear regime | exact | a missing `(−1)^{F−m_F}` in the Wigner–Eckart 3j |
| **V15** | **c_I / A∥ separability.** With A∥ = 0 and c_I ≠ 0, splitting = `c_I(2J+1)/2` (grows with J); with c_I = 0 and A∥ ≠ 0, splitting ∝ `(2J+1)/(2J(J+1))` (falls with J) | exact | conflating the two hyperfine operators, which are indistinguishable at a single J |

Explicitly **not** a check: any comparison to a stored spectrum for one parameter set. Every entry above is either a symmetry, an exact closed form, or a labelled opt-in comparison, and every one is reachable in both PASS and FAIL.

---

## 7. Open items for Arian

**OPEN-1 — ³Δ fine-structure constants for tier T2.** The observed ³Δ intervals are 1052.5 and 2097.5 cm⁻¹, ratio 2.0, so a single A cannot describe the manifold. Fitting `2AΣ + (2/3)λ(3Σ²−S²)` gives A ≈ 787.5, λ ≈ 261 cm⁻¹ `[derived]`, λ/A = 0.33. Options: (a) use these two derived constants; (b) treat the three Ω components as independent case-(c) origins with the measured T₀'s (no A, no λ); (c) go to T3, where a ¹Σ⁺ at 314 cm⁻¹ is an explicit perturber of ³Δ₁ and the apparent λ is an artefact of squeezing case (c) into case (a). I recommend (b) for T2 and (c) if the fine structure is ever the point.

**OPEN-2 — the e/f convention at S = 1 (must be decided before any parity label is written).** Arian's thesis rule `e ⇔ P = (−1)^{J−S−ℓ}` is the S = ½ specialisation; at S = 1 it gives the opposite label to B&C's J-only rule (PDF p. 283 / book p. 251). For ThF⁺, B&C says the **e** component is the **upper** one of every Ω-doublet; the thesis rule would say **f**. Options: (a) adopt B&C's rule globally and add a note in the thesis-convention layer that item 21 is S = ½-only; (b) keep the thesis rule and relabel every literature comparison. The physical statement is convention-free and verified: the **upper** component has parity `(−1)^J`.

**OPEN-3 — the sign of the G∥ Zeeman operator.** Ng thesis Eq. C.6 prints `H = −G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B`. That is inconsistent with the g_F relation printed on the same page, with the G∥ = 0.048(2) / −0.042(2) pair Ng quotes, and with the measured |g| = 0.0149 — all three of which require `+G∥`. Options: (a) `+G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B` (my recommendation; verified numerically to reproduce 20.85 kHz/G and δg/g = −0.00223); (b) Ng's printed operator with the *relative* sign of the g_N term flipped in the g_F formula — which then fails to reproduce his own G∥ = 0.048. Worth a one-line question to K. B. Ng.

**OPEN-4 — the sign of g_F.** Never measured. Every theory route (G∥ = +0.034, +0.035, +0.047; Petrov Fig. 2 showing g^{u,ℓ} ≈ −1.49 × 10⁻²) implies g_F < 0. Expose it as a switch, default negative, and make any observable that depends on it report which branch it used.

**OPEN-5 — the microscopic ¹⁹F hyperfine constants a, b_F, c.** Only the combination `A∥ = 2a − b_F − (2/3)c` is known. No source read gives them separately for ThF⁺, and no ab initio calculation of them exists. This matters only if you ever want the perpendicular/d-type hyperfine or the ΔΩ = ±1 elements; for the Ω = ±1 block A∥ is sufficient.

**OPEN-6 — c_I, and the A∥ bias it implies.** My estimate is ~20 kHz with a factor-3 uncertainty, from CsF. Two questions: (a) do you want c_I as a free knob (my recommendation) or fixed at 0 to match Ng? (b) Should the parameter table carry the caveat that Ng's `A∥ = −20.1(1) MHz` may absorb `2c_I`–`6c_I` ≈ 40–120 kHz, comparable to its own error bar? An ab initio c_I from L. Cheng or the Skripnikov group would settle it; it is a much easier calculation than G∥.

**OPEN-7 — g_r.** Not measured, not computed, absorbed into the fitted G∥. My rigid-rotor nuclear-charge estimate (0.47 μ_N) is a floor: the electronic part has the opposite sign and in a heavy molecule can dominate. Consequence: predicted g-factors at J = 2, 3, 4 carry an unquantified ~1 % error. Ng 2022 flags this as worth computing.

**OPEN-8 — e_Δ, the hyperfine-dependent Ω-doubling.** My scaling puts it at 1–10 kHz, i.e. straddling the v1 threshold; Leanhardt dismisses it without a number. It would make ω_ef depend on F, which is in principle measurable in the existing J = 1 / J = 2 microwave data. Do you want it in v1 as a zero-defaulted knob?

**OPEN-9 — W_a / W_P (nuclear-spin-dependent parity violation).** No ThF⁺ value found in any source read (Skripnikov & Titov 2015 full text, Denis 2015, Petrov 2025, Petrov 2023, Ng 2022, Ng thesis App. C). If you want NSD-PV in the engine, this parameter has to come from a new calculation.

**OPEN-10 — parity-dependent Zeeman terms.** Leanhardt Eq. 64 / B&C Eq. 9.70 terms (vi)–(vii). They produce Petrov's Δg(E = 0) = 2.3 × 10⁻⁴ (≈ 1 kHz at 1 G) and are absent from every Ω = ±1 effective model including Ng's. Include them in v1 with `|g′_rS| ≈ ω_ef/(2B_e) = 3.6 × 10⁻⁴` as the default, or leave the model unable to produce a zero-field Δg?

**OPEN-11 — n̂ direction (a convention, but it flips signs everywhere).** JILA: n̂ points **from F to Th** ("from the more negative atom to the more positive one", Leanhardt p. 15; "pointing towards thorium", Ng thesis p. 318). Petrov/Skripnikov: n̂ is "directed **from Th to F**" (Skripnikov & Titov 2015 p. 2; arXiv:2503.02840 p. 3). The signs of Ω, the signed d, E_eff, and "upper/lower doublet" all flip with the choice. Pick one and assert it at the top of the code. I have used the **JILA convention** throughout this document (n̂ from F to Th, d_mf = +3.37 D, Ω = +1 for ³Δ₁ with Λ = +2).

**OPEN-12 — δg vs Δg.** Ng thesis: `δg = (g^u − g^ℓ)/2`. Petrov and the Ng 2022 *paper*: `Δg = g^u − g^ℓ`. A factor of 2, and both appear in sources you will compare against. Which does the code report?

**OPEN-13 — Gresh's k″ sign.** Magnitudes agree beautifully across ten bands (0.869–0.892 × 10⁻⁴ cm⁻¹, giving ω_ef = 5.21–5.35 MHz against the microwave 5.29(5)). But the printed **sign** tracks the upper state's label — negative for every Ω = 0⁻ band, positive for every Ω = 0⁺ band (verified by reading the rendered Table 1). Since `k′ = 0` for all these bands, `s′` never enters the fit and the e/f assignment of the lower state depends on branch bookkeeping that flips with the upper state's reflection parity. The Ω = 0⁻ sign (k″ < 0, ⟹ e above f) is the one consistent with Ng 2022 Fig. 2. Confirm with the authors, or just take ω_ef and its sign from the microwave measurement.

**OPEN-14 — E_eff: 35.2 (Denis) vs 37.3 (Skripnikov) GV/cm**, both with ~7 % claimed uncertainty; Ng adopts ≈35, Petrov adopts 37.3. Which does the code default to? (Both supersede Meyer & Bohn's 90 GV/cm, which Skripnikov calls "more than twice overestimated".)

**OPEN-15 — Petrov's printed `D = −0.133 a.u.`** (arXiv:2503.02840 p. 3). −0.133 a.u. = −0.338 D, eight times smaller than the d_mf = 3.37 D they cite from the same reference in the same sentence; −1.33 a.u. = −3.381 D matches exactly, and matches the companion HfF⁺ paper's format (D∥ = −1.53(2) a.u.). Almost certainly a typo. Not load-bearing for us — we use Ng's measured 3.37(9) D — but worth knowing if their Δg table is ever re-derived.


**v2 open items.** OPEN-16 through OPEN-23 are defined in the v2 spec (`docs/superpowers/specs/2026-09-05-heff-v2-isotopologues-two-photon.md` §7) and collected in `docs/open-questions.md`. Task 1's rulings on four of them are written into §9 and are the citable source for the code: **OPEN-16** (A∥(Th) sign) → §9.4.3, resolved, default `a_par_th_sign = 'negative'`; **OPEN-17** (eQq₂ normalisation bridge) → §9.4.4, **not** resolved, escalated with two candidate factors and the one pinned ratio; **OPEN-21** (does K = 1 survive) → §9.5.3, resolved: no, not in exact closure; **OPEN-20** (²²⁷Th spin and moment) → §9.6, still open, with a labelled Schmidt placeholder.

---

## 8. Method log

**Read at source, in full or by named page range:**

| source | what was read |
|---|---|
| `ng-thesis-JILA_p316-325.txt` | App. B Table B.1–B.2 (PDF pp. 316–317) and **all of App. C** (pp. 318–325): Eqs. C.1–C.8, C.3.1 basis, C.4 two-level reduction |
| `ng-thesis-JILA_p35-40.txt` | Ch. 2.4: Hund's-case summary Fig. 2.4, Eq. 2.2a case (c) `E(J) = BJ(J+1)`, Ω = J_a·n̂ definition, Eqs. 2.3–2.5 |
| `ng-thesis-JILA_p76-96.txt` | §4.1.3 (Ramsey g-factor, δg/g global fit Fig. 4.9b), §4.1.4 Table 4.1 and the G∥ sign paragraph (pp. 84–87) |
| `ng2022-…txt` (arXiv:2202.01346) | pp. 1–5 in full: §II A–E, Table I, the "six distinct resonant frequencies" fit statement |
| `ng2022-…_p3.png` | **rendered and read as an image**: Fig. 2 parity colouring and the J = 1 / J = 2 doublet ordering — the load-bearing parity check of §2.3 |
| `leanhardt2011-arXiv1008.2997.txt` | Sec. IV A–B (PDF pp. 13–16): Eqs. 10–35; Sec. IV G–I (pp. 24–25): Eqs. 64–73 |
| `gresh2016-…txt` | §2.2 and §3 (PDF p. 3): Eqs. 1–3 including the `s = +1 for e` convention |
| `gresh2016-…_p10.png` | **rendered and read as an image**: Table 1, all X ³Δ₁ rows — B″, D″, k″, k_D″ and the sign pattern of §OPEN-13 |
| `2025-arXiv2503.02840-…txt` | read in full: Eqs. 1–20, Table I, Fig. 2 description, the D = −0.133 a.u. sentence, the electronic basis Eq. 3 |
| `2023-arXiv2302.02856-…txt` | pp. 2–3: Fig. 1 caption (n̂ ↔ upper/lower), Eqs. 1–10 including `f^BD = 2d_eE_eff` and the `M_F ≠ m_F` warning |
| `skripnikov2015-arXiv1503.01001-…txt` | pp. 1–3 and Table II: Eqs. 1–11 (W_d, H_d, W_{T,P}, H_{T,P}, H_MQM, A∥, **G∥ Eq. 11**), the 7 % uncertainty statement, FINAL(ThF⁺) row |
| `digest-thesis-effective-hamiltonian.md` §4 | Arian's 21-item convention contract (used as the target convention; items 2, 7, 9, 17, 18, 21 load-bearing here) |
| `digest-literature-thf-plus.md` | used as a map only; every sign, definition and number it lists was re-read at source |
| `synthesis-draft.md` §1 | Molecule-Structure's confirmed sign bug in the case (a) Λ-doubling q operator (identically zero: `K1 = K0+2q` vs the 3j needing `P1 = P0−2q`) and the absence of any ΔΛ = ±4 operator — noted, but no repo file was opened |

**Brown & Carrington** (`C:\Users\Arian\Zotero\storage\CKZKCGXY\…`, 1045 pages; PDF page = book page + 32). Read via PyMuPDF (pdf-mcp was down, HTTP 401). Page ranges extracted to `notes/lit/bc-pages/`:

| PDF pp. | book pp. | what was used |
|---|---|---|
| 195–210 | 163–178 | Appendix 5.1: **Eq. (5.174)** spectator theorem (used, per the brief, in place of the main-text form), (5.175), (5.172) Wigner–Eckart, (5.179), (5.184)–(5.186) D-matrix reduced element |
| 283 | 251 | §6.9.4 Eq. (6.234) parity combinations, the `(−1)^{J−S+s}` phase, and the **e/f convention definition** |
| 376–404 | 344–372 | Ch. 4/7 operator glossary: `c_I` = nuclear spin–rotation constant, `g_r`, `g_l`, anisotropic Zeeman |
| 410, 414 | 378, 382 | **Eq. (8.7)** `H_nsr = c_I T¹(J)·T¹(I)` and **Eq. (8.20)** its coupled-basis matrix element |
| 513 | 481 | **Table 8.12**: CsF molecular parameters, `c₂(¹⁹F) = 15.1 kHz` at B_v = 0.183782 cm⁻¹ — the c_I anchor |
| 636 | 604 | **Eqs. (9.50) and (9.51)**: the case-(a) axial hyperfine `{aΛ + (b_F + (2/3)c)Σ}Ω` and its ΔJ = ±1 partner |
| 636–652 | 604–620 | **Eq. (9.70)**: the complete seven-term Zeeman Hamiltonian, including the two parity-dependent terms |
| 686–700 | 654–668 | (o + p + q) Λ-doubling operator forms — ³Σ worked example; **no ³Δ (o_Δ, p_Δ, q_Δ) treatment found in B&C**, which is why §2.3 cites Leanhardt/Brown & Merer for the ΔΛ = 4 case |

**Page renders read as images: 2** (Ng 2022 p. 3; Gresh 2016 p. 10), out of a budget of 30. Both were pre-existing renders in `notes/lit/`; no new renders were needed.

**Numerics I ran** (scripts under the session scratchpad, conda env `claude-code`, Python 3.12, NumPy + `sympy.physics.wigner`; nothing written into `notes/` except this file):

1. Unit conversions and every size in §2 and §4.
2. Verification that Ng Eq. C.5's ΔJ = 0 element equals `−Ω m_F γ_F d_mf E` for all (J ≤ 4, F, m_F, Ω) — check V2, **passes**.
3. B&C (9.51) ΔJ = ±1 hyperfine matrix elements and second-order shifts, J = 1–4.
4. Second-order Ω-doubling from the Ω = 0 perturbers with 0± reflection signs — 5.294 MHz at J = 1, exact J(J+1) law.
5. Exhaustive basis enumeration, total and per-M_F, for both Ω = ±1 and the full ³Δ manifold.
6. **A 96-state Ω = ±1 Hamiltonian (J = 1–4) built from the §2 operators and diagonalised**, used to: confirm Hermiticity; confirm the zero-field level pattern (0, ω_ef, 3A∥/4, 3A∥/4 + ω_ef); settle the Eq. C.6 Zeeman sign (20.85 vs 23.52 kHz/G); measure the hyperfine ΔJ mixing (max 2.603 kHz); measure the c_I effect (50 kHz at c_I = 20 kHz); and extract `g^u`, `g^ℓ` from the m_F = ±3/2 splittings of the two Stark doublets at E = 40, 60, 100 V/cm, giving **δg/g = −0.00223 at 60 V/cm**, identical to Ng's own 32-level value, with convergence already reached at J_max = 2.

**What was not available.** `pdf-mcp` and the `zotero` MCP server both failed to connect (HTTP 401); B&C and all PDFs were read with PyMuPDF instead. Brown & Merer 1979 (the ³Δ effective-Hamiltonian paper Leanhardt cites as [118]) and Nelis et al. [120] are not in `notes/lit/` and were not obtained — so the `o_Δ, p_Δ, q_Δ` and `g_rS, g′_rS` operator definitions in §2.3 and §2.9 rest on Leanhardt's transcription of them, not on the originals. Loh et al., Science 342, 1220 (2013) supplement was likewise not obtained. Petrov, PRA 108, 062804 (2023) (the source of the S₂, S₃ formulas) and Caldwell et al., PRA 108, 012804 (2023) were not obtained; neither affects any parameter used here.


### 8.1 Addendum for §9 (v2 derivations, 2026-09-05)

**Read at source for §9, beyond what §8 lists.** Brown & Carrington via PyMuPDF from the same Zotero copy (`pdf-mcp` and the `zotero` MCP server were down again this session — the same HTTP 401 / ConnectionRefused as before, so no PDF-native search was available):

| PDF pp. | book pp. | what was used for §9 |
|---|---|---|
| 164 | 132 | **Eqs. (4.29)–(4.32)**: `H_Q = −e T²(∇E)·T²(Q)`, the definition of `T²(∇E)` and of `eT²(Q)` — the inputs to the §9.4.4 normalisation bridge |
| 198 | 166 | **Eqs. (5.141), (5.142)**: the coupled tensor `T^K(A₁,B₁)` and its reduced element, i.e. the two-photon operator of §9.5 |
| 205 | 173 | **Eqs. (5.172), (5.173), (5.174), (5.175), (5.176)** — Wigner–Eckart, the coupled scalar product, the two spectator forms, and the same-inner-part theorem |
| 206 | 174 | **Eq. (5.179)** `⟨j‖T¹(j)‖j′⟩` — note [HAM] §2 cites this as PDF p. 205; **the correct PDF page is 206** |
| 207 | 175 | **Eqs. (5.185), (5.186)** — the D-matrix reduced element |
| 410, 414 | 378, 382 | **Eqs. (8.7), (8.19), (8.20)** re-read; (8.19) is the (5.174) skeleton `heff`'s `dipole_geometry` implements, and it confirms the phase and 6j column order term for term |
| 636–637 | 604–605 | **Eqs. (9.50)–(9.53)** re-read in full, including B&C's sentence that "q₀ is the negative of the electric field gradient" |

**Other sources read at source for §9:**

| source | what was read |
|---|---|
| `petrov2018-CP-violation-HfFplus.pdf` | **p. 2 rendered and read as an image** to recover Eq. (19) (the text layer mangles it); pp. 2–3 text for Eqs. (15)–(25), the `Q = 2⟨Q̂²₀⟩` definition, the HfF⁺ constants, and §IV's |5s5dπ| spin–orbit-admixture model |
| `docs/lit/lookup-apar-th-sign-convention.md` (commit `0b3e5fa`) | headline, §1.1 and the conclusion — recorded in §9.4.3, not re-derived |
| `docs/lit/lookup-227th-nuclear-moment.md` | all three compilations' Z = 90 rows and the conclusion — the input to §9.6 |
| `docs/digest-literature-th-hyperfine.md` §§1.4, 2.2, 3, 4 | coupling scheme, the two ab initio A∥ values, `G_el = −10 408 MHz`, and both target tables |
| `docs/digest-literature-two-photon.md` §3 | Cossel Eqs. (6.29)/(6.34)/(6.37), the ΔΩ and Δm_F probes, the detuning caveat |
| arXiv:2201.13247 §1 (fetched 2026-09-05) | the Schmidt-moment statement quoted verbatim in §9.6 |

**Page renders read as images for §9: 1** (Petrov 2018 p. 2, Eq. 19).

**Numerics run for §9** (scratch script `verify_s9.py` under the session scratchpad, conda env `heff`, `heff.wigner` + `heff.elements_c.dipole_geometry` + NumPy; **nothing written into the package**):

1. §9.1 reduction of the two-spectator element to `dipole_geometry` at `I_Th = 0` (1720 elements, max dev 2.220 × 10⁻¹⁶) and at `I_F = 0` (same), plus two 6j-corruption falsification probes.
2. §9.3 Tables 1 and 2, each against two independent routes (the full 6j recoupling and a 6j-free projection-theorem closed form), plus three falsification probes on Table 2.
3. §9.3 ΔJ = ±1 axial hyperfine against B&C (9.51) at `I_Th = 0` (max dev 5.551 × 10⁻¹⁷), and `c_I^F` against the (8.20) closed form.
4. §9.3 independent rebuild of the whole recoupling in the **decoupled** `|J,m_J⟩|I_Th,m₁⟩|I_F,m₂⟩` basis with explicit Clebsch–Gordan — the check that settles the B&C (5.173) phase, which Hermiticity cannot.
5. §9.4 B&C (9.53) against the Casimir function for (J, I) = (1,1), (2,1), (2,3/2), (3,5/2), all F — uniform ratio −1.000000.
6. §9.4.4 reconstruction of Petrov Eq. (24)'s 483 MHz prefactor from Eq. (23) under both readings of `√(2π/5)` — **inconclusive** (−284.8 and −402.8 MHz), reported as such.
7. §9.5 rank decomposition of `Σ_i d_a|i⟩⟨i|d_b` in a spherical-harmonic model space with and without a common denominator — `‖K=1‖ = 2.2 × 10⁻¹⁶` under exact closure, non-zero and `∝ 1/Δ` otherwise.
8. §9.6 arithmetic: `μ_Schmidt`, `g_N`, `A∥(²²⁷Th)`, the `μ(²²⁷) = μ(²²⁹)` alternative, and the resulting J = 1 splitting and ΔJ = ±1 element against `4B`.

**What was not available for §9.** The two questions that would close OPEN-17 need the authors of PRA 98, 042502 (2018); no erratum or later restatement of their Eqs. (19)/(22)/(23) was found. Long, *The Raman Effect* (2002) — the conventional home of the Placzek polarisability derivation named in [2γ] §3.4 — was again not accessible, so §9.5's closure form rests on B&C (5.141)/(5.142) and Cossel's thesis, not on that textbook. Bonin & McIlrath, JOSA B **1**, 52 (1984) remains abstract-level only. Kälber et al., Z. Phys. A **334**, 103 (1989), the one paper that might carry a ²²⁷Th hyperfine measurement, is still paywalled.

---

## 9. Two nuclear spins, and the rank-K two-photon operator (v2 derivations)

Written 2026-09-05 for `heff` v2 (²²⁷ThF⁺, ²²⁹ThF⁺, and the two-photon operator). **This section is the single source Tasks 4, 5 and 7 copy from; they cite §9 and do not re-derive.** Every formula carries its Brown & Carrington equation number with both the PDF page and the book page of the copy at `C:\Users\Arian\Zotero\storage\CKZKCGXY\Brown and Carrington, Rotational Spectroscopy of Diatomic Molecules.pdf` (**PDF page = book page + 32**). Steps that are algebra or numerics run here are tagged `[derived]`; anything a source does not support is tagged **UNVERIFIED** and escalated in §7, never guessed.

**Notation contract for all of §9, stated once because it is the commonest way to get a phase backwards.**

- **A prime means the bra.** `⟨J′, Ω′, F₁′, F′, m′_F | … | J, Ω, F₁, F, m_F⟩`. This is the convention of [HAM] §2, of `heff/elements_c.py`, and of Ng thesis Eq. C.5. **B&C's own convention is the opposite** — in (5.172)–(5.176), (5.186), (9.50)–(9.53) the *primed* labels are the **ket**. Every B&C equation reproduced below is first quoted verbatim in B&C's convention and then rewritten with the primes moved onto the bra; the rewrite is nothing but a relabelling, but skipping it transposes 6j columns and flips phases.
- **Coupling scheme**: `F₁ = J + I_Th`, `F = F₁ + I_F`, basis `|((J I_Th) F₁, I_F) F, m_F⟩` ([TH] §3.1; justified there by |A∥(Th)| ≈ 1.5 GHz ≫ |A∥(F)| = 20.1 MHz and by Petrov 2018's use of F₁ as the HfF⁺ label). Ω is the signed molecule-frame projection, as in §1.1. Setting `I_Th = 0` must return every v1 formula of §2 exactly; that is the master gate of §9.1.
- **Molecule-frame component**: `q = Ω′ − Ω = Ω_bra − Ω_ket`. This is forced by the 3j `(J′ k J; −Ω′ q Ω)`, whose projections must sum to zero. B&C write the same physical statement as `q = Ω − Ω′` in (9.52) **because their primes are on the ket** — the two are the same rule, not two conventions. Cairncross thesis p. 141 ("we can simply set q = Ω′ − Ω") uses primed = bra and agrees with the form written here.
- Wigner 3j in round brackets, 6j in curly brackets, both in the `heff.wigner` argument order: `(j₁ j₂ j₃; m₁ m₂ m₃)` and `{j₁ j₂ j₃; j₄ j₅ j₆}`.

---

### 9.1 The two-spectator axial geometry

**The operator.** Any molecule-frame tensor of rank k with a lab component p — the E1 dipole, the axial hyperfine vector n̂, the rank-2 field gradient, the rank-K two-photon polarisability — enters the lab frame through B&C's space-to-molecule transformation, so its matrix element factorises into (i) a lab Wigner–Eckart 3j on F, (ii) a spectator reduction of I_F out of F, (iii) a spectator reduction of I_Th out of F₁, and (iv) a rotation-matrix reduced element. Four factors, four cited equations, in that order.

**(i) Lab Wigner–Eckart — B&C Eq. (5.172), PDF p. 205 / book p. 173, verbatim:**

```
⟨j, m|T^k_p(A)|j′, m′⟩ = (−1)^{j−m} ( j k j′ ; −m p m′ ) ⟨j‖T^k(A)‖j′⟩
```

With primes on the bra and j = F:

```
m1 = (−1)^{F′−m′_F} ( F′  k  F ; −m′_F  p  m_F )
```

**(ii) I_F spectator — B&C Eq. (5.174), PDF p. 205 / book p. 173, verbatim** (this is the Appendix 5.1 form, the one [HAM] §2 already uses; the main-text form is not used anywhere in this document):

```
⟨j₁, j₂, j₁₂‖T^{k₁}(A₁)‖j₁′, j₂′, j₁₂′⟩
   = δ_{j₂ j₂′} (−1)^{j₁₂′ + j₁ + k₁ + j₂} [(2j₁₂+1)(2j₁₂′+1)]^{1/2}
     × { j₁′  j₁₂′  j₂ ;  j₁₂  j₁  k₁ } ⟨j₁‖T^{k₁}(A₁)‖j₁′⟩
```

Applied with `j₁ = F₁, j₂ = I_F, j₁₂ = F, k₁ = k`, and with the primes moved onto the bra:

```
m2 = (−1)^{F + F₁′ + k + I_F} √((2F′+1)(2F+1)) { F₁  F  I_F ;  F′  F₁′  k }
```

**The 6j column order is literal**: upper row `(F₁, F, I_F)` = (ket's inner angular momentum, ket's total, spectator), lower row `(F′, F₁′, k)` = (bra's total, bra's inner, operator rank). Transposing the two upper entries `F₁ ↔ F` breaks the reduction check below by 0.65 in a quantity of order 1 (§9.1 self-check, falsification row 1).

**(iii) I_Th spectator — B&C Eq. (5.174) again**, now with `j₁ = J, j₂ = I_Th, j₁₂ = F₁, k₁ = k`:

```
m3 = (−1)^{F₁ + J′ + k + I_Th} √((2F₁′+1)(2F₁+1)) { J  F₁  I_Th ;  F₁′  J′  k }
```

Same column order, one level in: upper row `(J, F₁, I_Th)`, lower row `(F₁′, J′, k)`.

**Why (5.174) twice and (5.175) never.** B&C (5.175), PDF p. 205 / book p. 173, is the companion form for an operator acting on the **second** constituent of a coupled pair. Here both reductions have the operator on the **first** constituent — `F₁` inside `F = F₁ + I_F`, then `J` inside `F₁ = J + I_Th` — because the coupling scheme was built inner-first for exactly this reason. (5.175) is needed only for an operator acting on a nuclear spin directly; §9.3 handles that case with the scalar-product form (5.173) instead, so (5.175) is not used anywhere in §9. `[derived]`

**(iv) Rotation-matrix reduced element — B&C Eq. (5.186), PDF p. 207 / book p. 175, verbatim:**

```
⟨J, Ω‖D^{(k)}_{·q}(ω)*‖J′, Ω′⟩ = (−1)^{J−Ω} ( J k J′ ; −Ω q Ω′ ) [(2J+1)(2J′+1)]^{1/2}
```

With primes on the bra:

```
m4 = (−1)^{J′−Ω′} √((2J′+1)(2J+1)) ( J′  k  J ; −Ω′  q  Ω ),      q = Ω′ − Ω
```

**The master element.** `[derived]`

```
⟨J′,Ω′,F₁′,F′,m′_F| T^k_p(axial, molecule-frame component q) |J,Ω,F₁,F,m_F⟩
  = (−1)^{F′−m′_F} ( F′  k  F ; −m′_F  p  m_F )                                   ← (5.172)
  × (−1)^{F + F₁′ + k + I_F} √((2F′+1)(2F+1)) { F₁  F  I_F ;  F′  F₁′  k }        ← (5.174), I_F spectator
  × (−1)^{F₁ + J′ + k + I_Th} √((2F₁′+1)(2F₁+1)) { J  F₁  I_Th ;  F₁′  J′  k }    ← (5.174), I_Th spectator
  × (−1)^{J′−Ω′} √((2J′+1)(2J+1)) ( J′  k  J ; −Ω′  q  Ω ) ⟨η′‖T^k‖η⟩             ← (5.186)

with q = Ω′ − Ω, and Δm_F = p forced by the first 3j.
```

This is the same four-line layout §2.7 uses for the Stark element, with one extra line — line 3 — and with `k` and `q` left general instead of frozen at `k = 1, q = 0`.

**Analytic collapse at I_Th = 0** `[derived]`. Then `F₁ = J`, `F₁′ = J′`, and line 3 contains `{J J 0; J′ J′ k}`. A 6j with a zero in the upper-right slot is `{a a 0; c c f} = (−1)^{a+c+f} / √((2a+1)(2c+1))`, so line 3 becomes

```
(−1)^{J + J′ + k} √((2J′+1)(2J+1)) × (−1)^{J+J′+k} / √((2J+1)(2J′+1)) = 1
```

exactly, and the remaining three lines at `k = 1, q = 0` are `heff.elements_c.dipole_geometry`'s `m1 · m2 · m3` term for term. The same argument at `I_F = 0` (so `F = F₁`) collapses line 2 to 1 and leaves `dipole_geometry` with `I → I_Th`, `F → F₁`. Both limits are therefore analytic identities, not approximations.

**Self-check as run** (scratch script `verify_s9.py`, `conda run -n heff python`, `heff.wigner` + `heff.elements_c.dipole_geometry`; script kept out of the repo). Printed verbatim:

```
S9.1  reduction of the two-spectator element
  (a) I_Th = 0, k=1,q=0 vs dipole_geometry : 1720 elements, max |dev| = 2.220e-16
  (b) I_F  = 0, k=1,q=0 vs dipole_geometry(I->I_Th, F->F1): max |dev| = 2.220e-16
  falsification: transpose the I_F 6j  -> max |dev| = 6.455e-01
  falsification: transpose the I_Th 6j -> max |dev| = 2.220e-16   (blind at I_Th = 0, see (b)/S9.3)
```

Row (a) is the check the brief asks for: all 1720 elements over `J ≤ 4`, `J′ ≤ 4`, `Ω = ±1`, every `F = J ± ½`, every `m_F`, every `p ∈ {0, ±1}`, **max deviation 2.220 × 10⁻¹⁶** — machine epsilon, i.e. exact.

**The check's blind spot, stated because a reviewer must know it.** Transposing the *I_F* 6j's upper pair is caught loudly (0.65). Transposing the *I_Th* 6j's upper pair is **not** caught at I_Th = 0, because there `F₁ = J` makes the corrupted symbol identical to the correct one. Two further checks close that hole and both are reported below: row (b) here (`I_F = 0`, `I_Th = 5/2`, `F₁ ≠ J`, max deviation 2.220 × 10⁻¹⁶), and the §9.3 Table 2 falsification, where the same corruption collapses the whole table to zero. Do not treat row (a) alone as validating line 3.

---

### 9.2 The Th operators are the v1 formulas with F → F₁, exactly

**The structural theorem — B&C Eq. (5.176), PDF p. 205 / book p. 173, verbatim:**

```
(viii) Tensor operators acting on the same inner part of a coupled system

If the scalar product is formed from spherical tensor operators which both act on the
same inner part of a coupled scheme, it is intuitively obvious that

  ⟨j₁, j₂, j₁₂, m₁₂| T^k(A₁) · T^k(B₁) |j₁′, j₂′, j₁₂′, m′₁₂⟩
     = δ_{j₁₂ j₁₂′} δ_{m₁₂ m′₁₂} δ_{j₂ j₂′} ⟨j₁| T^k(A₁) · T^k(B₁) |j₁′⟩.        (5.176)

This result can be proved formally by application of the Wigner–Eckart theorem,
equation (5.172), followed by equation (5.174).
```

**Application.** In `|((J I_Th) F₁, I_F) F, m_F⟩` take `j₁ = F₁` (which already contains J and I_Th), `j₂ = I_F`, `j₁₂ = F`. Every Th interaction — the Th magnetic hyperfine, the Th electric quadrupole, the Th nuclear spin–rotation — is a scalar built from the rotational/electronic degrees of freedom and I_Th alone, i.e. from the **inner** part. B&C (5.176) then says its matrix element is diagonal in F and m_F, **independent of them and of I_F**, and equal to the matrix element evaluated inside `|J, Ω, I_Th, F₁⟩`. `[derived from B&C (5.176)]`

**Consequence, stated as the substitution rule Tasks 4 and 5 implement.** With `I → I_Th`, `F → F₁` and no other change, the following v1 expressions are the **exact** Th matrix elements in the two-spin basis, times `δ_{FF′} δ_{m_F m′_F}`:

| v1 formula | B&C source | becomes |
|---|---|---|
| §2.4 axial hyperfine, ΔJ = 0 | (9.50), PDF p. 636 / book p. 604 | `A∥^Th [F₁(F₁+1) − J(J+1) − I_Th(I_Th+1)] / [2J(J+1)]` |
| §2.5 axial hyperfine, ΔJ = ±1 | (9.51), PDF p. 636 / book p. 604 | same expression with `F → F₁`, `I → I_Th`, J the larger of the two |
| §2.6 nuclear spin–rotation | (8.7) PDF p. 410 / book p. 378, element (8.20) PDF p. 414 / book p. 382 | `c_I^Th [F₁(F₁+1) − I_Th(I_Th+1) − J(J+1)] / 2` |
| §9.4 electric quadrupole | (9.52)/(9.53), PDF pp. 636–637 / book pp. 604–605 | same, `I → I_Th`, `F → F₁` |

Selection rules for all four: `ΔF₁ = 0`, `ΔF = 0`, `Δm_F = 0`, and the v1 ΔJ/ΔΩ rules unchanged.

**What this means for the code, so no new algebra is written for the Th terms.** `heff`'s existing `hyperfine_A_par`, `hyperfine_A_par_dJ1` and `spin_rotation_cI` are already the bare Casimir/(9.51) kernels — none of them reads Ω except `hyperfine_A_par_dJ1`, which divides by the signed Ω because B&C's brace is `A∥/Ω`. The Th versions are **those same functions read through a field adapter** that hands them `(I_Th, F₁)` where the v1 call hands them `(I_F, F)`. That is a one-line indirection, not a new matrix element, and it is the reason §9.3 — not §9.2 — is where the work is.

**Cross-check run** `[derived]`: at `I_Th = 5/2` the diagonal (9.50) coefficients over F₁ at J = 1 are `−1.7500, −0.5000, +1.2500`, reproducing [TH] §4.1 exactly.

---

### 9.3 The ¹⁹F operators, recoupled

These are the terms that genuinely change shape: `T¹(I_F)` acts on the **outer** spin, the rotational factor acts inside F₁, so B&C (5.176) does not apply and F₁ becomes off-diagonal.

**(a) Scalar product across the coupled pair — B&C Eq. (5.173), PDF p. 205 / book p. 173, verbatim:**

```
⟨j₁, j₂, j₁₂, m| T^k(A₁) · T^k(A₂) |j₁′, j₂′, j₁₂′, m′⟩
  = (−1)^{j₁′ + j₁₂ + j₂} δ_{j₁₂, j₁₂′} δ_{m, m′}
    { j₂′  j₁′  j₁₂ ;  j₁  j₂  k } ⟨j₁‖T^k(A₁)‖j₁′⟩ ⟨j₂‖T^k(A₂)‖j₂′⟩
```

Applied with `j₁ = F₁, j₂ = I_F, j₁₂ = F, k = 1`, and with primes moved onto the bra:

```
s1 = (−1)^{F₁ + F + I_F} { I_F  F₁  F ;  F₁′  I_F  1 }
     × ⟨J′,F₁′‖T¹(A₁)‖J,F₁⟩ × √(I_F(I_F+1)(2I_F+1))
```

**The phase carries the KET's F₁, not the bra's.** This matters only for the `ΔF₁ = ±1` elements, where the two readings differ by a factor −1 — and **Hermiticity does not discriminate between them**, because the reduced element `⟨J′,F₁′‖T¹(A₁)‖J,F₁⟩` itself changes sign under bra↔ket exchange when `ΔF₁ = ±1`, so both readings give a symmetric matrix. The discriminating check is the independent decoupled-basis rebuild reported below.

`⟨I_F‖T¹(I_F)‖I_F⟩ = [I_F(I_F+1)(2I_F+1)]^{1/2}` is B&C Eq. (5.179), PDF p. 206 / book p. 174 (note: [HAM] §2 cites this equation as PDF p. 205; the correct PDF page is 206).

**(b) I_Th spectator — B&C Eq. (5.174)** as in §9.1 step (iii), `j₁ = J, j₂ = I_Th, j₁₂ = F₁, k₁ = 1`:

```
⟨J′,F₁′‖T¹(A₁)‖J,F₁⟩ = (−1)^{F₁ + J′ + 1 + I_Th} √((2F₁′+1)(2F₁+1)) { J  F₁  I_Th ;  F₁′  J′  1 }
                        × ⟨J′,Ω′‖T¹(A₁)‖J,Ω⟩
```

with two choices of the innermost reduced element:

- `A₁ = J`: `⟨J′‖T¹(J)‖J⟩ = δ_{JJ′}[J(J+1)(2J+1)]^{1/2}` — **B&C (5.179), PDF p. 206 / book p. 174**.
- `A₁ = n̂`: `⟨J′,Ω′‖T¹(n̂)‖J,Ω⟩ = (−1)^{J′−Ω′}√((2J′+1)(2J+1)) (J′ 1 J; −Ω′ 0 Ω)` — **B&C (5.186), PDF p. 207 / book p. 175**, i.e. exactly the `m4` line of §9.1 at `k = 1`, `q = 0`, with `J′ ≠ J` allowed.

**The three ¹⁹F terms.**

1. **ΔJ = 0 axial hyperfine.** B&C (9.50)'s ΔJ = 0 block is `A∥^F T¹(J)·T¹(I_F)/[J(J+1)]` — the projection-theorem form Ng thesis Eq. C.2 uses. Substituting `A₁ = J` above gives the complete two-spin element. Selection rules: `ΔJ = 0`, `ΔΩ = 0`, **`ΔF₁ = 0, ±1`**, `ΔF = 0`, `Δm_F = 0`.
2. **Nuclear spin–rotation.** `H_nsr^F = c_I^F T¹(J)·T¹(I_F)` — B&C (8.7), PDF p. 410 / book p. 378. Identical skeleton with `A₁ = J` and no `1/[J(J+1)]`. Same selection rules.
3. **ΔJ = ±1 axial hyperfine** — the term [TH] §3.3 leaves "derivable-not-derived". It is the **same scalar product with `A₁ = n̂`**:

   ```
   H_hf^F(axial) = (A∥^F / Ω) T¹(n̂) · T¹(I_F)
   ```

   `[derived]` The normalisation is fixed, not fitted: within a fixed J the projection theorem gives `⟨J,Ω‖T¹(n̂)‖J,Ω⟩ = Ω ⟨J‖T¹(J)‖J⟩ / [J(J+1)]` exactly, because `(J 1 J; −Ω 0 Ω) = (−1)^{J−Ω} Ω / √(J(J+1)(2J+1))`, so `(A∥^F/Ω) T¹(n̂)` and `A∥^F T¹(J)/[J(J+1)]` have **identical** reduced elements on the ΔJ = 0 block and therefore identical matrix elements at every F₁ and F. The two forms differ only off-diagonal in J, which is precisely where (9.51) lives. Selection rules: `ΔJ = 0, ±1`, `ΔΩ = 0`, `ΔF₁ = 0, ±1`, `ΔF = 0`, `Δm_F = 0`.

**Closed form used as an independent cross-check** `[derived]`. Applying the projection theorem twice — `T¹(J)·T¹(I_F)` projected onto F₁, then `F₁·I_F` projected onto F — the `F = F₁+½` minus `F = F₁−½` splitting of term 1, in units of `A∥^F`, is

```
Δ(J, F₁) = [J(J+1) + F₁(F₁+1) − I_Th(I_Th+1)] (2F₁+1) / [4 F₁(F₁+1) J(J+1)]
```

which at `I_Th = 0` (so `F₁ = J`) reduces to `(2J+1)/[2J(J+1)]`, the §2.4 result. This closed form uses **no 6j at all**, so agreement between it and the full recoupling is a genuine test of the phases and column orders, not a restatement.

#### Self-checks as run

Scratch script `verify_s9.py`, `conda run -n heff python`. Printed verbatim:

```
S9.3  the 19F operators, recoupled
  (5.173) phase uses KET F1: max |M - M^T| = 0.000e+00
  (5.173) phase uses BRA F1: max |M - M^T| = 0.000e+00
  dJ=0: (A/[J(J+1)]) T1(J).T1(I_F) vs (A/Om) T1(n).T1(I_F): max |diff| = 2.776e-17

  TABLE 1  I_Th = 0, splitting E(F=J+1/2)-E(F=J-1/2) in units of A_par^F
     J | recoupled | (2J+1)/[2J(J+1)] | projection thm
     1 |  +0.75000  |     +0.75000      |   +0.75000
     2 |  +0.41667  |     +0.41667      |   +0.41667
     3 |  +0.29167  |     +0.29167      |   +0.29167
     4 |  +0.22500  |     +0.22500      |   +0.22500

  TABLE 2  I_Th = 5/2, J = 1, splitting E(F=F1+1/2)-E(F=F1-1/2), units A_par^F
      F1 | recoupled | projection thm | [TH] S4.5
     1.5 |  -0.40000  |    -0.40000     |   -0.4000
     2.5 |  +0.17143  |    +0.17143     |   +0.1714
     3.5 |  +0.57143  |    +0.57143     |   +0.5714

  falsification of TABLE 2 (must NOT reproduce -0.4000/+0.1714/+0.5714):
    transpose I_F 6j           -> -0.1581, +0.0000, +0.0000
    transpose I_Th 6j          -> +0.0000, -0.0000, +0.0000
    (5.173) phase on BRA F1    -> -0.4000, +0.1714, +0.5714

  dJ = +-1 axial hyperfine at I_Th = 0 vs B&C (9.51), units A_par^F:
     J=2<-1, F=1.5: recoupled -0.433013 | B&C(9.51) -0.433013
     J=3<-2, F=2.5: recoupled -0.471405 | B&C(9.51) -0.471405
     J=4<-3, F=3.5: recoupled -0.484123 | B&C(9.51) -0.484123
     max |dev| over J=2..4, Om=+-1, all F : 5.551e-17

  dJ = +-1 at I_Th = 5/2 (J=2<-1), units A_par^F, F1 diagonal:
     F1=1.5 F=1.0: +0.295804
     F1=1.5 F=2.0: -0.177482
     F1=2.5 F=2.0: +0.252982
     F1=2.5 F=3.0: -0.180702
     F1=3.5 F=3.0: +0.185577
     F1=3.5 F=4.0: -0.144338

  c_I^F at I_Th = 0 vs B&C (8.20) closed form [F(F+1)-I(I+1)-J(J+1)]/2:
     max |dev| = 2.220e-16
```

**Table 1 agrees with [TH] §3.4 to all five digits** (+0.75000, +0.41667, +0.29167, +0.22500), and independently with the §2.4 closed form and with the projection theorem.

**Table 2 agrees with [TH] §4.5 to all four digits** (−0.4000, +0.1714, +0.5714). The **sign inversion at F₁ = 3/2** relative to the I_Th = 0 value of +0.7500 survives three independent routes — the full 6j recoupling, the projection-theorem closed form above, and [TH]'s earlier independent implementation. Its mechanism is visible in the closed form: the bracket `J(J+1) + F₁(F₁+1) − I_Th(I_Th+1)` is `2 + 3.75 − 8.75 = −3` at (J, F₁) = (1, 3/2), i.e. **J is anti-aligned with F₁ there**, so the ¹⁹F doublet ordering flips inside that F₁ manifold. At `A∥^F = −20.1 MHz` the three splittings are +8.0, −3.4 and −11.5 MHz. This is the sharpest falsifiable ²²⁹ThF⁺ prediction in this document.

**The falsification block shows the tables can fail.** Corrupting either 6j's column order sends Table 2 to `−0.1581, 0, 0` or to `0, 0, 0`. In particular the I_Th 6j corruption — the one the §9.1 `I_Th = 0` check is blind to — is caught here, loudly.

**ΔJ = ±1 reproduces B&C (9.51) exactly at I_Th = 0**, max deviation 5.551 × 10⁻¹⁷ over J = 2–4, Ω = ±1, all F. So the axial `T¹(n̂)·T¹(I_F)` form is the operator (9.51) is the matrix element of, and `heff`'s existing `hyperfine_A_par_dJ1` sits inside the master gate.

**The (5.173) phase, settled by an independent route.** Because Hermiticity is blind to it (both readings give `max |M − Mᵀ| = 0`), the element was rebuilt from scratch in the fully **decoupled** product basis `|J,m_J⟩|I_Th,m₁⟩|I_F,m₂⟩`, with the coupled states assembled from Clebsch–Gordan coefficients and `T¹(J)·T¹(I_F) = J_z I_z + ½(J₊I₋ + J₋I₊)` written out as an explicit matrix — a route that uses **neither (5.173) nor (5.174)**. Printed verbatim:

```
S9.3 (independent) decoupled-basis rebuild of the (5.173)+(5.174) chain
  J = 1, A = T1(J), I_Th = 5/2:  <F1'|O|F1> / [J(J+1)]
     F1'   F1    F  | decoupled  | (5.173) ket-phase | bra-phase
      1.5  1.5  1.0 | +0.250000  |    +0.250000      | +0.250000
      1.5  1.5  2.0 | -0.150000  |    -0.150000      | -0.150000
      1.5  2.5  2.0 | -0.374166  |    -0.374166      | +0.374166
      2.5  1.5  2.0 | -0.374166  |    -0.374166      | +0.374166
      2.5  2.5  2.0 | -0.100000  |    -0.100000      | -0.100000
      2.5  2.5  3.0 | +0.071429  |    +0.071429      | +0.071429
      2.5  3.5  3.0 | -0.319438  |    -0.319438      | +0.319438
      3.5  2.5  3.0 | -0.319438  |    -0.319438      | +0.319438
      3.5  3.5  3.0 | -0.321429  |    -0.321429      | -0.321429
      3.5  3.5  4.0 | +0.250000  |    +0.250000      | +0.250000
```

**The ket-F₁ phase — B&C's literal reading — is right; the bra-F₁ variant is wrong on every `ΔF₁ = ±1` element.** Code the phase as `(−1)^{F₁(ket) + F + I_F}`.

**Sizes.** The `ΔF₁ = ±1` ¹⁹F elements at I_Th = 5/2, J = 1 are 0.32–0.37 × `A∥^F` ≈ 6.4–7.5 MHz against F₁ spacings of 190–2670 MHz ([TH] §4.1), so F₁ remains a good quantum number and this sign moves levels only at the ≲ 25 kHz level. It is nonetheless a sign in a matrix element that check V1 cannot see, which is why it is pinned here rather than left to the implementer.

---

### 9.4 The quadrupole, and the eQq₂ normalisation bridge

#### 9.4.1 The B&C matrix elements, verbatim

**B&C Eq. (9.52), PDF p. 636 / book p. 604** (B&C's primes are on the ket):

```
⟨η, Λ; S, Σ; J, Ω, I, F, M_F| H_Q |η′, Λ′; S, Σ; J′, Ω′, I, F, M_F⟩
  = −(1/2) eQ Σ_q (−1)^{J′+I+F+J−Ω} {(2J+1)(2J′+1)}^{1/2}
    × { J′  I  F ;  I  J  2 } ( J  2  J′ ; −Ω  q  Ω′ ) ( I  2  I ; −I  0  I )^{−1}
    × ⟨η, Λ| T²_q(∇E) |η′, Λ′⟩
```

**B&C Eq. (9.53), PDF p. 637 / book p. 605**, the q = 0 specialisation:

```
  = (e q₀ Q / 4) (−1)^{J′+I+F+J−Ω} {(2J+1)(2J′+1)}^{1/2}
    × { J′  I  F ;  I  J  2 } ( J  2  J′ ; −Ω  0  Ω ) ( I  2  I ; −I  0  I )^{−1} ,

where q₀ is the negative of the electric field gradient, and eq₀Q is the quadrupole
coupling constant. The 3-j symbol indicates that matrix elements with J′ = J, J ± 1
and J ± 2 are non-zero, but the diagonal elements are, of course, the most significant.
```

With `I → I_Th`, `F → F₁` (justified by §9.2), and with the primes moved onto the bra so that the molecule-frame component reads `q = Ω′ − Ω` as everywhere else in §9, these are the two Th quadrupole terms `heff` needs.

**Two implications read straight off the 3j `(J 2 J′; −Ω q Ω′)`:**

- `ΔJ = 0, ±1, ±2` (B&C say so in the sentence quoted above).
- `q = Ω_bra − Ω_ket`, so inside the Ω = ±1 block of ³Δ₁ there are exactly two terms: **q = 0, ΔΩ = 0** (constant `eq₀Q`) and **q = ∓2, ΔΩ = ±2** (constant `eq₂Q`). The second is parity-even, diagonal in J, F₁, F, m_F, and sits in the same matrix position as the Ω-doubling operator ([TH] §3.2, §4.4).
- **Both vanish identically for I ≤ ½**, because `(I 2 I; −I 0 I)` in the denominator requires the triangle `(I, 2, I)`, i.e. `I ≥ 1`. So ²²⁷ThF⁺ (I_Th = ½, §9.6) has no Th quadrupole at all, exactly as ¹⁹F has none (§2.11). `[derived]`

**Comparing (9.52) at q = 0 with (9.53) fixes B&C's own definition of the constant** `[derived]`:

```
−(1/2) eQ ⟨η,Λ| T²₀(∇E) |η,Λ⟩ = e q₀ Q / 4      ⟹      e q₀ Q = −2 eQ ⟨η,Λ| T²₀(∇E) |η,Λ⟩
```

which is the algebraic content of B&C's sentence "q₀ is the negative of the electric field gradient".

#### 9.4.2 The Casimir cross-check, and the factor −1 that is easy to lose

[TH] §3.4 reports that (9.53) at Ω = 0, J′ = J reproduces the textbook Casimir function with a **uniform ratio of exactly −1**. Re-run here independently `[derived]`; printed verbatim:

```
S9.4  quadrupole: B&C (9.53) vs the Casimir function
     (J, I) |  F  |  (9.53)/eq0Q  |  Casimir  |  ratio
     (1,1) | 0.0 |  -0.5000000   | +0.5000000 | -1.000000
     (1,1) | 1.0 |  +0.2500000   | -0.2500000 | -1.000000
     (1,1) | 2.0 |  -0.0500000   | +0.0500000 | -1.000000
     (2,1) | 1.0 |  -0.2500000   | +0.2500000 | -1.000000
     (2,1) | 2.0 |  +0.2500000   | -0.2500000 | -1.000000
     (2,1) | 3.0 |  -0.0714286   | +0.0714286 | -1.000000
     (2,1.5) | 0.5 |  -0.2500000   | +0.2500000 | -1.000000
     (2,1.5) | 1.5 |  +0.0000000   | +0.0000000 | +nan
     (2,1.5) | 2.5 |  +0.1785714   | -0.1785714 | -1.000000
     (2,1.5) | 3.5 |  -0.0714286   | +0.0714286 | -1.000000
     (3,2.5) | 0.5 |  -0.2000000   | +0.2000000 | -1.000000
     (3,2.5) | 1.5 |  -0.1100000   | +0.1100000 | -1.000000
     (3,2.5) | 2.5 |  +0.0066667   | -0.0066667 | -1.000000
     (3,2.5) | 3.5 |  +0.1000000   | -0.1000000 | -1.000000
     (3,2.5) | 4.5 |  +0.1000000   | -0.1000000 | -1.000000
     (3,2.5) | 5.5 |  -0.0833333   | +0.0833333 | -1.000000
```

Casimir function used: `[¾C(C+1) − I(I+1)J(J+1)] / [2I(2I−1)(2J−1)(2J+3)]`, `C = F(F+1) − I(I+1) − J(J+1)`. The single `nan` row is `0/0` (both sides vanish at (J, I, F) = (2, 3/2, 3/2)), not a failure. **The uniform −1 *is* B&C's q₀ convention; code it explicitly or the sign of every quadrupole splitting flips.**

**Third-source corroboration of the angular skeleton** (not of the normalisation): Skripnikov, Petrov, Titov & Flambaum, arXiv:1408.5368, Eq. (5), give the ThO ³Δ₁ MQM shift as `δ(J,F) = (−1)^{Ω+I+F+1} C(J,F) W_M M` with `C(J,F) = [(2J+1)/2] (J 2 J; −Ω 0 Ω)/(I 2 I; −I 0 I) {J I F; I J 2}` — the identical rank-2 case-(c) structure, differing only by an overall constant and by an overall sign consistent with B&C's q₀ convention ([TH] §3.4). Three sources, one angular structure.

#### 9.4.3 OPEN-16 — the sign of A∥(Th) is settled, and it is not a convention

Recorded here because Tasks 4 and 5 need the ruling and must not re-open it. `docs/lit/lookup-apar-th-sign-convention.md` (committed `0b3e5fa`) audited the Skripnikov & Titov 2015 (A∥ = −4163 μ_Th/μ_N MHz) versus Denis et al. 2015 (+1833 MHz) disagreement. Its findings: both groups do state their molecular-axis convention elsewhere in the same paper and the axes are opposite (Skripnikov ζ from Th to F, Denis F→Th); **A∥ as both define it is invariant under a consistent axis reversal** `[derived there]`; and the two groups **agree** on the sign of the analogous HfF⁺ constant. **The disagreement is therefore real physics, not a convention fork**, and the audit recommends `A∥(²²⁹Th, ThF⁺ X ³Δ₁) = −1.51(11) GHz`. **Arian ruled 2026-09-05 that the package default is `a_par_th_sign = 'negative'`.** The flag records which calculation you trust, not which convention you work in — unlike `n_hat`, which is a genuine convention switch — and the docstring must say so. Do not re-derive this; see the audit for the four-step chain and the one email that would close it.

#### 9.4.4 OPEN-17 — the eQq₂ normalisation bridge: NOT resolved, escalated

**What Petrov 2018 prints.** From `docs/lit/petrov2018-CP-violation-HfFplus.pdf`, page 2 (Eq. 19, rendered and read as an image because the text layer mangles it) and page 3 (Eqs. 22–23):

```
Ĥ_hfs = … + −e² Σ_q (−1)^q Q̂²_q(I¹) Σ_i √(2π/5) Y_{2q}(θ_{1i}, φ_{1i}) / r_{1i}³        (19)

eQq₀ = 2 eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₀(θ_{1i}, φ_{1i}) / r_{1i}³ |³Δ₁⟩                        (22)

eQq₂ = 2 √6 eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₂(θ_{1i}, φ_{1i}) / r_{1i}³ |³Δ₋₁⟩                    (23)

with  Q = 2 ⟨U^Hf_{I¹I¹}| Q̂²₀(I¹) |U^Hf_{I¹I¹}⟩
```

**What B&C print** (PDF p. 164 / book p. 132, Ch. 4, verbatim):

```
H_Q = −e T²(∇E) · T²(Q),                                                              (4.30)

T²(∇E) = −(1/4πε₀) Σ_i (e_i / R_i³) C²(θ_i, φ_i),                                     (4.31)

e T²(Q) = e Σ_p R_p² C²(θ_p, φ_p),                                                    (4.32)
```

**The algebra, as far as the printed equations take it** `[derived]`. Racah's tensor is `C²_q = √(4π/5) Y_{2q}`, so Petrov's electronic operator is `𝒱_q ≡ Σ_i √(2π/5) Y_{2q}/r_i³ = (1/√2) Σ_i C²_q/r_i³`. With electron charges `e_i = −e` in (4.31), and in atomic units,

```
T²_q(∇E) = + Σ_i C²_q(i) / r_i³ = √2 𝒱_q .
```

Both Q definitions agree (`Q = 2⟨T²₀(Q)⟩` in each), so, using `eq₀Q = −2eQ⟨T²₀(∇E)⟩` from §9.4.1 and the same definition extended to `q = ±2`,

```
eq₀Q(B&C)  = −√2 · eQq₀(Petrov 2018)
eq₂Q(B&C)  = −√2 · eQq₂(Petrov 2018) / √6  =  −eQq₂(Petrov 2018) / √3
```

**Why this is a candidate, not the answer.** Two things are not pinned by the printed text:

1. **The scalar-product pairing in Eq. (19).** As printed, *both* tensor indices are `q` — `Σ_q (−1)^q Q̂²_q 𝒱_q` — which is not a scalar. A scalar product requires `Σ_q (−1)^q A_q B_{−q}`. The derivation above assumes the standard pairing (the electronic factor carrying `−q`). Under that reading Eq. (19) is **1/√2 of B&C's (4.30)**, which is exactly what a systematic `√(2π/5)` for `√(4π/5)` would produce, and Eqs. (19), (22), (23) are then mutually consistent but sit √2 below the conventional normalisation.
2. **Whether `√(2π/5)` is intended.** If it is a typo for `√(4π/5)` (i.e. Racah's `C²_q`), then `eq₀Q(B&C) = −eQq₀(Petrov)` — which is the implicit assumption behind [TH] §4.3's use of Petrov's −2100 MHz — and `eq₂Q(B&C) = −eQq₂(Petrov)/√6`. **The two readings differ by exactly √2 and the printed equations do not discriminate.**

**Attempt to discriminate from Petrov's own internal consistency** `[derived]`. Petrov Eqs. (24)–(25) give `eQq₂ = 483 w Q ⟨1/r³⟩_{5d}` MHz with `w ≈ G∥ + 0.002319`, and at `G∥ = 0.011768` this reproduces his quoted 110 MHz (`483 × 0.014087 × 3.365 × 4.86 = 111.3` MHz ✓). Reconstructing the prefactor 483 from Eq. (23) with the model Petrov states (spin–orbit admixture, weight `w`, of a Π state with leading configuration |5s5dπ|, so the active matrix element is `⟨5d, λ=+1| 𝒱₂ |5d, λ=−1⟩` times `⟨1/r³⟩_{5d}`):

```
  Petrov 2018 Eq.(24) prefactor probe  [eQq2 = 483 w Q <1/r^3>_5d MHz]
    sqrt(2 pi/5)  as printed : <d+1|Y22|d-1> ang = -0.220728, prefactor = -284.8 MHz
    sqrt(4 pi/5)  = C^2_q    : <d+1|Y22|d-1> ang = -0.220728, prefactor = -402.8 MHz
    (Petrov's own numbers: w = 0.014087, Q = 3.365 b, <1/r^3> = 4.86 -> 111.3 MHz vs the 110 MHz quoted)
```

(Angular factor `∫Y*_{21}Y_{22}Y_{2,−1}dΩ`; conversion 234.97 MHz per unit `Q[barn]·⟨1/r³⟩[a.u.]`.) **Neither reading reproduces 483**, so this probe is **inconclusive** — the single-orbital, `w`-as-weight model is too crude to separate a √2 from the two-electron and spin–orbit detail Petrov does not print. It is reported because it was run and because it rules out the lazy conclusion that the factor is obvious.

**What *is* pinned, and is usable** `[derived]`. Petrov's Eqs. (22) and (23) carry the **same** electronic-operator normalisation `√(2π/5) Y_{2q}`; only the leading factor differs (2 versus 2√6). So whatever the absolute bridge turns out to be, the **relative** one is fixed:

```
eq₂Q(B&C) / eq₀Q(B&C)  =  [ eQq₂(Petrov) / √6 ] / eQq₀(Petrov)
```

i.e. **if a published Petrov-style `eQq₀` is entered into (9.53) as-is (which is what [TH] §4.3 and §4.4 do), then consistency requires entering `eQq₂ / √6` into (9.52) at `q = ±2`.** That relation follows from the printed equations alone and does not depend on either open question above.

**Ruling for the code, per the brief's stop-work condition.**

- `eqq2_norm = 'bc_9p52_q2'` is the **only implemented normalisation**: the package computes in B&C's normalisation, where the parameter is `eq₂Q` defined by `−2eQ⟨η,Λ|T²_{±2}(∇E)|η,Λ′⟩`.
- The `'petrov2018_eq23'` converter must **raise `NotImplementedError`**, naming OPEN-17 and quoting the two candidate factors (`−1/√3` and `−1/√6`) and the reason they cannot be separated. **A guessed factor is a stop-work condition; none is guessed here.**
- Every `eQq₂`-derived number in [TH] §4.4 (the ~35–80 MHz Ω = +1 ↔ Ω = −1 element at J = 1) inherits a √2 and a sign caveat and is order-of-magnitude only. So does [TH] §4.3's eQq₀ estimate, which took Petrov's −2100 MHz directly into (9.53).

**OPEN-17 (escalated, unchanged in status, sharpened).** Two questions for the authors, either of which settles it in one sentence: (i) in PRA 98, 042502 (2018) Eq. (19), is the electronic factor `Y_{2,−q}` (standard scalar product)? (ii) is the `√(2π/5)` in Eqs. (19), (22), (23) intended, or is it `√(4π/5) = C²_q`? Absent an answer, the package refuses the conversion rather than picking a branch.

---

### 9.5 The rank-K two-photon operator

**The starting operator** ([2γ] §3.1; Cossel PhD thesis (Colorado, 2014) Eqs. (6.29), (6.34), pp. 213, 219). Adiabatic elimination of a far-detuned intermediate manifold `{|i⟩}` gives, between X-state levels,

```
T_eff = Σ_i ( d·ε₂* ) |i⟩⟨i| ( d·ε₁ ) / Δ_i ,        Δ_i = E_i − E_g − ħω₁
```

with amplitude `A_{f←g} = ⟨f|T_eff|g⟩` and line strength `|A|²` — the sum over `p₁, p₂` taken **before** squaring (Cossel Eq. 6.34; this is the two-photon analogue of `heff`'s gate B8, and Cossel's Fig. 6.18 shows a real measured cancellation from it).

#### 9.5.1 The rank decomposition and the closure form

**B&C Eq. (5.141), PDF p. 198 / book p. 166, verbatim:**

```
T^K_p(A₁, B₁) = (−1)^{k₁−k₂+p} (2K+1)^{1/2} Σ_{p₁ p₂} ( k₁ k₂ K ; p₁ p₂ −p )
                 × T^{k₁}_{p₁}(A₁) T^{k₂}_{p₂}(B₁)
```

**B&C Eq. (5.142), same page, verbatim:**

```
⟨η, j‖T^K(A₁, B₁)‖η′, j′⟩ = (2K+1)^{1/2} (−1)^{K+j+j′} Σ_{η″ j″}
      { k₁  k₂  K ;  j′  j  j″ } ⟨η, j‖T^{k₁}(A₁)‖η″, j″⟩ ⟨η″, j″‖T^{k₂}(B₁)‖η′, j′⟩
```

With `k₁ = k₂ = 1` these are *literally* the two-photon operator: the `Σ_{η″ j″}` is the sum over intermediate states, and the 6j `{1 1 K; j′ j j″}` is what a resolved-intermediate calculation would otherwise carry numerically. Two rank-1 operators couple to `K = 0, 1, 2` and nothing else.

**(1) The closure form** `[derived]`. When `Δ_i → Δ` is common it comes out of the sum, and (5.142) can be read backwards: the intermediate sum **is** the reduced element of a single rank-K operator,

```
Σ_{η″ j″} { 1  1  K ;  j′  j  j″ } ⟨η,j‖T¹(d)‖η″,j″⟩ ⟨η″,j″‖T¹(d)‖η′,j′⟩
      = (−1)^{K+j+j′} (2K+1)^{−1/2} ⟨η,j‖T^K(d, d)‖η′,j′⟩
```

so that

```
⟨η,j‖α^K‖η′,j′⟩ ≡ (1/Δ) ⟨η,j‖T^K(d, d)‖η′,j′⟩
```

is **one scalar per (K, ΔΩ) channel times parameter-free geometry** ([2γ] §3.4). This is the Placzek-type polarisability picture; `heff` never sums over intermediates in this form.

**(2) Spectator reduction over both nuclear spins.** `α^K` acts on the electronic–rotational part only, so §9.1's chain applies **unchanged with `k → K`**:

```
⟨J′,Ω′,F₁′,F′,m′_F| α^K_P |J,Ω,F₁,F,m_F⟩
  = (−1)^{F′−m′_F} ( F′  K  F ; −m′_F  P  m_F )                                   ← (5.172)
  × (−1)^{F + F₁′ + K + I_F} √((2F′+1)(2F+1)) { F₁  F  I_F ;  F′  F₁′  K }        ← (5.174), I_F spectator
  × (−1)^{F₁ + J′ + K + I_Th} √((2F₁′+1)(2F₁+1)) { J  F₁  I_Th ;  F₁′  J′  K }    ← (5.174), I_Th spectator
  × (−1)^{J′−Ω′} √((2J′+1)(2J+1)) ( J′  K  J ; −Ω′  q  Ω ) ⟨η′‖α^K‖η⟩             ← (5.186)
```

`q = Ω′ − Ω`, `Δm_F = P`. **This is the same four lines as §9.1 with `k` replaced by `K` and `p` by `P`** — the one-photon `{J F I; F′ J′ 1}` becomes `{J F I; F′ J′ K}`, exactly as [2γ] §3.2 says. Validity condition: the spectator reduction is exact **only when the intermediate hyperfine structure is unresolved**; if `Δ_i` depends on the intermediate `F′`, the `F′` sum cannot be factored out and the reduction fails.

**(3) The lab contraction and the polarisation dyad** `[derived from B&C (5.141)]`. B&C (5.141) with `k₁ = k₂ = 1` is exactly the Clebsch–Gordan coupling, since `⟨k₁p₁ k₂p₂|KP⟩ = (−1)^{k₁−k₂+P}(2K+1)^{1/2}(k₁ k₂ K; p₁ p₂ −P)`. So the polarisation dyad is

```
(ε₁ ⊗ ε₂)^K_P = (−1)^P (2K+1)^{1/2} Σ_{p₁ p₂} ( 1  1  K ; p₁  p₂  −P ) ε₁^{p₁} ε₂^{p₂}
              = Σ_{p₁ p₂} ⟨1 p₁ 1 p₂ | K P⟩ ε₁^{p₁} ε₂^{p₂} ,          P = p₁ + p₂
```

with the lab spherical components of a Jones vector `(ε_x, ε_y, ε_z)` in Condon–Shortley phase (as fixed in [HAM] §2):

```
ε_{+1} = −(ε_x + i ε_y)/√2 ,      ε_0 = ε_z ,      ε_{−1} = +(ε_x − i ε_y)/√2
```

and the full contraction

```
T_eff = Σ_{K=0,1,2} Σ_P (−1)^P (ε₁ ⊗ ε₂)^K_{−P} α^K_P .
```

`K = 0` is the scalar `ε₁·ε₂`; `K = 1` is the antisymmetric part, `∝ ε₁ × ε₂`, non-zero only for non-parallel or elliptical polarisations; `K = 2` is the symmetric traceless part.

#### 9.5.2 Which (K, ΔΩ) channels exist — the algebra, not the assertion

Within X ³Δ₁ both `|Ω|` are 1, so `q = Ω′ − Ω ∈ {0, ±2}`. The molecule-frame factor is the 3j `(J′ K J; −Ω′ q Ω)`, and **a 3j vanishes identically unless each projection satisfies `|m| ≤ j`** — here `|q| ≤ K`. Since `k₁ = k₂ = 1` bounds `K ≤ 2`:

```
ΔΩ = 0   : K = 0, 1, 2      (|q| = 0 ≤ K always)
ΔΩ = ±2  : K = 2 only       (|q| = 2 requires K ≥ 2, and K ≤ 2)
```

Verified numerically over `J, J′ = 1, 2, 3` at `Ω′ = −1, Ω = +1` (so `q = −2`); printed verbatim:

```
  molecule-frame 3j (J' K J; -Om' q Om) at Om'=-1, Om=+1 (q = -2):
     K = 0: max |3j| = 0.000000
     K = 1: max |3j| = 0.000000
     K = 2: max |3j| = 0.447214
```

This is the algebraic content of [2γ] §3.3's table (Ω = +1 → −1 requires an Ω = 0 intermediate, because each E1 leg has `|q_i| ≤ 1` and `ΔΩ = q₁ + q₂`).

#### 9.5.3 OPEN-21 — does K = 1 survive? Answer: no, not in exact closure

**The algebra** `[derived]`. `K = 1` is precisely the antisymmetric part of the dyad: the `K = 1` Clebsch–Gordan coefficients `⟨1 p₁ 1 p₂|1 P⟩` are antisymmetric under `p₁ ↔ p₂`, so `α^1 ∝ ½(d_a 𝒫 d_b − d_b 𝒫 d_a)` with `𝒫 = Σ_i |i⟩⟨i|` the intermediate projector. Now take exact closure: for `|g⟩` in X, `d|g⟩` lies **entirely** in the opposite-parity space, so if `𝒫` projects onto that whole space it acts as the identity on `d|g⟩`, and

```
α^1 ∝ P_X ½ ( d_a d_b − d_b d_a ) P_X = P_X ½ [d_a, d_b] P_X = 0
```

because `d = −e Σ_i r_i` is a vector operator whose Cartesian components commute. **`K = 1` is identically zero in exact closure**, confirming [SPEC-v2] §3.2's inference.

**The order at which it reappears** `[derived]`. Writing `1/Δ_i = (1/Δ)(1 − δ_i/Δ + …)` with `δ_i = E_i − Ē`, the zeroth-order term is the closure term above and vanishes; the first-order term is `−(1/Δ²) P_X d_a (H − Ē) d_b P_X`, whose antisymmetric part `−(1/2Δ²) P_X (d_a H d_b − d_b H d_a) P_X` does **not** vanish. So `K = 1` is suppressed by one power of **(intermediate splitting)/(detuning)** — equivalently, it is the well-known antisymmetric Raman tensor, which exists only when the two time-orderings' denominators differ (resolved intermediates, or `ω₁ ≠ ω₂`).

**Numerical verification, with both outcomes reachable.** Model: the sphere-harmonic space `{|l,m⟩ : l ≤ 7}`, with `d_p = C¹_p` — a genuine vector operator whose Cartesian components commute (multiplication by a function on the sphere) and which is parity-odd, connecting even `l` to odd `l` only. X = even `l ≤ 6`; intermediates = odd `l ≤ 7`, which is **complete** for `d` acting on X. The dyad is decomposed into ranks with B&C (5.141). Printed verbatim:

```
S9.5  does K = 1 survive exact closure?
  Delta =    20.0 : ||K=0|| 2.6243e-02  ||K=1|| 4.0177e-03  ||K=2|| 1.6598e-02   K1/K0 = 1.531e-01
  Delta =   200.0 : ||K=0|| 2.8582e-03  ||K=1|| 1.4411e-04  ||K=2|| 2.0587e-03   K1/K0 = 5.042e-02
  Delta =  2000.0 : ||K=0|| 2.8839e-04  ||K=1|| 2.0330e-06  ||K=2|| 2.2563e-04   K1/K0 = 7.050e-03
  Delta = 20000.0 : ||K=0|| 2.8865e-05  ||K=1|| 2.1122e-08  ||K=2|| 2.2867e-05   K1/K0 = 7.318e-04
  max |P_X (d_a P_odd d_b - d_b P_odd d_a) P_X| = 3.331e-16   (unrestricted [d_a,d_b] at the l<=L edge: 4.667e-01)
  EXACT CLOSURE (common denominator): ||K=0|| 5.7735e-01  ||K=1|| 2.2204e-16  ||K=2|| 4.5803e-01
  ENERGY-WEIGHTED (Delta = 20)     : ||K=0|| 2.6243e-02  ||K=1|| 4.0177e-03  ||K=2|| 1.6598e-02
```

Three things this shows. (i) With a **common** denominator, `||K=1|| = 2.2 × 10⁻¹⁶` — zero to machine precision — while `K = 0` and `K = 2` are of order 0.5, so the test is not vacuous: the same code returns a non-zero `K = 1` in the weighted case. (ii) `K1/K0` falls as `1/Δ` asymptotically (7.05 × 10⁻³ → 7.32 × 10⁻⁴ for a factor-10 increase in Δ), confirming the "first order in `δ/Δ`" statement. (iii) The commutator is zero **on the subspace that matters** (3.3 × 10⁻¹⁶); the 0.467 is the unrestricted commutator at the `l ≤ L` truncation edge, which never enters because the projectors exclude it.

**Ruling for Task 7: register `K ∈ {0, 2}` only.** There is no `alpha_K1_*` parameter in the closure operator, and none should be invented. A `K = 1` term becomes meaningful only if the *resolved* form ([2γ] §3.5, §4.2) is ever implemented, where it enters at `O(δ/Δ)`; at that point it would need its own `placeholder` α with the order recorded here.

#### 9.5.4 Selection rules, as data

From the lab Wigner–Eckart 3j `(F′ K F; −m′_F P m_F)` and `K ≤ 2` `[derived]`:

| rule | value | source |
|---|---|---|
| `\|ΔF\| ≤ K` | so `ΔF ∈ {0, ±1, ±2}` | B&C (5.172) triangle |
| `Δm_F = P`, `\|P\| ≤ K` | so `Δm_F ∈ {0, ±1, ±2}`, **never ±3** | B&C (5.172) projection; [2γ] §3.2 |
| `ΔΩ ∈ {0, ±2}` | `±2` at `K = 2` only | §9.5.2 |
| `ΔJ`: `\|ΔJ\| ≤ K` | `ΔJ ∈ {0, ±1, ±2}` | 3j `(J′ K J; …)` |
| parity | **even** — at zero field it does not connect e to f | [2γ] §3.3, `[derived]` there: `P d P† = −d`, twice |
| under σ± only (JILA) | `p₁, p₂ ∈ {±1}` ⟹ `Δm_F ∈ {0, ±2}` | Ng thesis p. 102 fn. 4; [2γ] §3.3 probe |

The narrowing to `{0, ±2}` is the operational statement: Ng's target `|J=1, F=3/2, m_F=+3/2⟩ → |m_F=+1/2⟩` is `Δm_F = −1` and genuinely out of reach with σ± only, exactly as he says; `m_F = +3/2 → −1/2` is `Δm_F = −2` and **is** reachable with a σ⁻σ⁻ pair ([2γ] §3.3, `[derived]` there, not claimed to be JILA's intent).

#### 9.5.5 The validity condition, in one sentence the notebook can quote

> The closure form of the two-photon operator requires a detuning large compared with the intermediate rotational structure, `Δ ≫ 2B ≈ 7 GHz` for ThF⁺, whereas the JILA experiments run at 0.16–1.5 GHz — so this is the right *operator shape* and the wrong *limit for the current experiment*, and the α's could later be generated by a resolved sum once the intermediate ladder and its 0⁺/0⁻ labels are settled.

([2γ] §3.4: Cossel thesis p. 218 for the 160 MHz → ≈1.5 GHz detunings; Gresh 2016 Table 2 for `B ≈ 0.23 cm⁻¹ ≈ 6.9 GHz`; [2γ] §4.2 for "A can generate B's parameters, B cannot generate A's spectra". The intermediate *hyperfine* structure is a different story and being hyperfine-unresolved is plausible; being rotationally unresolved is not, at JILA detunings.)

---

### 9.6 ²²⁷ThF⁺: the A∥(Th) placeholder from the Schmidt moment

**Why a placeholder is needed at all.** ²²⁷Th has **no measured or estimated magnetic dipole moment anywhere** — `docs/lit/lookup-227th-nuclear-moment.md` queried Stone INDC(NDS)-0794 (2019) p. 42, the earlier PSI-hosted Stone compilation p. 149, and the IAEA NDS live nuclear-moments database, and all three jump from ²²⁷Ac straight to ²²⁹Th with no A = 227 row at Z = 90; four targeted searches returned no numeric value. The ENSDF ground-state assignment is **(1/2⁺)** — parenthesised, i.e. tentative — with 9.3 keV (5/2⁺) and 24.38 keV (3/2⁺) above it, the pattern of a K = 1/2 rotational band ([TH] §1.4).

**Arian's ruling (2026-09-05): use the Schmidt single-particle moment, not `μ(²²⁷) = μ(²²⁹)`.**

**Which orbital a 1/2⁺ odd neutron implies** `[derived]`. Parity `(−1)^ℓ = +1` forces even ℓ; `j = 1/2` with even ℓ forces `ℓ = 0` (an ℓ = 2 neutron gives j = 3/2 or 5/2). So the odd neutron (Z = 90, N = 137) sits in an **s₁/₂ orbital** — in the N > 126 shell, `4s₁/₂` — and this is the **`j = ℓ + ½`** Schmidt case.

**The Schmidt formula** (Schmidt, Z. Phys. **106**, 358 (1937); standard modern statement in the shell-model review arXiv:2201.13247 §1, fetched 2026-09-05: *"the single particle magnetic moments, commonly called the Schmidt moments, for an odd neutron are μ = μ_n for j = ℓ + 1/2, and μ = −j/(j+1) μ_n for j = ℓ − 1/2"*, with `μ_n = −1.913` in nuclear magnetons):

```
j = ℓ + ½ :   μ = ( j − ½ ) g_ℓ + ½ g_s          (in μ_N)
j = ℓ − ½ :   μ = [ j/(j+1) ] [ ( j + 3/2 ) g_ℓ − ½ g_s ]
```

with, for a **neutron**, `g_ℓ = 0` and `g_s = −3.826`.

**The number** `[derived]`. With `ℓ = 0`, `j = ½ = ℓ + ½`:

```
μ_Schmidt(²²⁷Th) = ½ g_s = ½ (−3.826) = −1.913 μ_N        (= the free-neutron moment, as it must be
                                                            for an s₁/₂ neutron)
g_N = μ / I = −1.913 / (1/2) = −3.826

A∥(²²⁷Th, ThF⁺ X ³Δ₁) = G_el × g_N = (−10 408 MHz) × (−3.826) = +39 821 MHz ≈ +39.8 GHz
```

with `G_el = A∥/g_N = −10 408 MHz`, the isotope-free electronic factor from Skripnikov & Titov 2015 as extracted in [TH] §2.2 (and the value that reproduces `A∥(²²⁹Th) = −10 408 × 0.1464 = −1524 MHz`).

**Status: `placeholder`, and it must be labelled as one everywhere it appears.** Two independent reasons, both real:

1. **Schmidt values for deformed actinides are typically wrong by a factor ~2.** Observed odd-A moments lie between the Schmidt lines, not on them; the standard remedy is a quenched effective spin g-factor (`g_s^eff ≈ 0.6 g_s^free` is the usual order), which alone would move `A∥` from +39.8 GHz to ~+24 GHz.
2. **A K = 1/2 band decouples.** ²²⁷Th's level pattern is a K = 1/2 rotational band, and the moment of a K = 1/2 band carries a **decoupling-parameter** term (Bohr & Mottelson) that is absent from the spherical single-particle estimate and can change the magnitude and, in principle, the sign. On top of that the `(1/2⁺)` spin assignment is itself tentative, and the 5/2⁺ level sits only 9.3 keV away.

**The alternative, recorded for comparison** `[derived]`: taking `μ(²²⁷) = μ(²²⁹) = 0.366 μ_N` with `I = 1/2` gives `g_N = 0.732` and `A∥ = −10 408 × 0.732 = −7619 MHz ≈ −7.62 GHz` — same order as the Schmidt value in magnitude but **opposite in sign** and 5× smaller. The two candidates disagree in sign, which is the honest measure of how unconstrained this is.

**Structural consequence that does not depend on which number is used** `[derived]`. At `I_Th = 1/2` the Th hyperfine has the same algebraic shape as the ¹⁹F structure already coded (§9.2), with `F₁ = J ± 1/2` and splitting `A∥^Th (2J+1)/[2J(J+1)]`, i.e. `0.75 A∥^Th` at J = 1. At the Schmidt value that is **+29.9 GHz at J = 1, against `4B = 29.1 GHz`** — the Th hyperfine splitting is as large as the whole J = 1 → 2 rotational interval. The ΔJ = ±1 element of B&C (9.51) at (J = 2 ← 1, F₁ = 3/2, Ω = 1) is `−0.433013 × A∥^Th = −17.2 GHz`, i.e. **0.59 × the rotational spacing**: J is not even approximately a good quantum number for ²²⁷ThF⁺ at this A∥, second-order perturbation theory is meaningless, and the J_max truncation must be re-tested for this isotopologue specifically (OPEN-22). Even at the −7.62 GHz alternative the ΔJ = ±1 element is 3.30 GHz (0.113 × the rotational spacing), i.e. ~340× the 8.7–9.7 MHz ²³²ThF⁺ elements of §2.5. **²²⁷ThF⁺ is not a small perturbation of the v1 model, whichever placeholder is used.**

**Escalated as OPEN-20** (the tentative spin and the absent moment). What would close it: a measurement of `μ(²²⁷Th)`, or a deformed-shell-model / DFT calculation with the K = 1/2 decoupling parameter included. The one plausibly relevant unread source is Kälber et al., Z. Phys. A **334**, 103 (1989) (collinear laser spectroscopy of stored Th⁺ ions across ²²⁷–²³²Th), paywalled in the session that produced `lookup-227th-nuclear-moment.md`; its abstract synopsis reports hyperfine analysis for ²²⁹Th only, but that is not a verified absence.
