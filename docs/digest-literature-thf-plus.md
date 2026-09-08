# Digest: effective Hamiltonian of ThF⁺ X ³Δ₁ (and its HfF⁺ ³Δ₁ lineage)

Assembled 2026-09-04 from the sources and page ranges listed in §2.
Conventions used in this file: **[conv.]** = unit conversion I performed (arithmetic shown); **[inferred]** = my inference, not stated by the source; everything else is a direct report with a cite. HfF⁺ values are always labelled as such and are never carried into the ThF⁺ table.

---

## 1. Bottom line and confidence

**Code this Hamiltonian** (Ng thesis, Appendix C, Eq. C.1, p. 318 of the PDF / p. 298 of the thesis):

    H = H_hyperfine + H_Ω-doub + H_rot + H_Stark + H_Zeeman + H_offset(v,J) + H_eEDM

in the **product basis |J, Ω = ±1, F, m_F⟩** (electronic Ω × rotational Wigner-D × ¹⁹F nuclear spin), truncated to J = 1 (12 states) or J = 1,2 (32 states). Every term is given in closed form with explicit 3j/6j matrix elements. Ng's Appendix C is the single best reference to code from: it is ThF⁺-specific, it consolidates the scattered JILA formalism, it states the frames explicitly, and its numbers are the measured ones. Second reference, for the *derivation* and for everything involving the rotating trap field (Berry phase, avoided crossings Δ^u/Δ^ℓ, δg_F(E)): Leanhardt et al. 2011 J. Mol. Spectrosc. 270, 1, Sec. IV (Eqs. 10–43, 64–67) — but that paper is written for HfF⁺ parameters. Third, for an independent numerical implementation with an explicit operator-level Hamiltonian: Petrov & Skripnikov arXiv:2503.02840 (ThF⁺) plus Petrov, Skripnikov & Titov PRA 96, 022508 / arXiv:1704.06631 (HfF⁺, Eqs. 5–8), which give Ĥ_rot, Ĥ_hfs, Ĥ_ext as operators rather than as reduced matrix elements.

**Parameter completeness for ²³²Th¹⁹F⁺ X ³Δ₁ (v = 0): good.** Measured: A∥, ω_ef, d_mf, |g_F=3/2|, δg_F/g_F, B₀, D₀, ω_e, ω_eχ_e, α_e, T₀ of the neighbouring states, v = 1 lifetime, T₁. Computed only: E_eff, W_T,P, G∥, and the E-field dependence of g^u, g^ℓ. **Not found anywhere: the sign of g_F is not measured** (theory forces g_F < 0, G∥ > 0 — §5.7); the ¹⁹F quadrupole/dipolar hyperfine terms beyond A∥ are never quoted (they are dropped, §3.1); and there is no published D-value or centrifugal correction to ω_ef.

Confidence: high on the term list, the basis, and the measured constants (multiple independent cross-checks close numerically to <1%, §7). Medium on two sign/factor conventions that differ between sources (H_eEDM factor of 2 and the direction of n̂) — both are documented in §5 and §7 and must be pinned in code, not guessed.

---

## 2. Search scope and coverage

**Library search.** Title and attachment-name queries run on 2026-09-04:
  `*ThF*`, `*HfF*`, `*Gresh*`, `*Cairncross*`, `*Leanhardt*`, `*Cornell*`, `*Skripnikov*`, `*Petrov*`, `*Roussy*`, `*Loh*`, `*Zhou*`, `*Denis*`, `*Fleig*`, `*eEDM*`, `*electron electric dipole*`, `*thorium*`, `*Ng -*`, `*velocity modulation*`, `*trapped molecular*`, `*Omega-doubl*`, `*Grau*`, `*Meyer*`, `*Vutha*`, `*g factor*`, `*g-factor*`, `*spin-rotational*`, `*Titov*`, `*Mosyagin*`, `*thesis*`.
The broader `*thesis*` query found author-title attachments missed by author-name
queries. Exact-title searches covered ThF spectroscopy, velocity-modulation
spectroscopy, coherence, electric-field-dependent g factors, and the HfF⁺
spin-rotational Hamiltonian.

**Online search.** arXiv queries covered the exact ThF⁺ spectroscopy title,
ThF⁺ hyperfine and Ω-doubling, ThF g factors, Skripnikov/Petrov ThF⁺ work,
and Petrov's Zeeman papers. Hit counts and exclusions below preserve the search
boundary.

**Sources obtained as full text:**

| Source | Route | Note |
|---|---|---|
| Ng et al., PRA 105, 022823 (2022) | arXiv:2202.01346 (the JILA-hosted `PhysRevA.105.022823.pdf` URL in the brief returns an HTML error page, 34 kB, not a PDF) | read in full |
| Ng PhD thesis (JILA), 325 pp. | Zotero `QX2C5ZA4/KiaBoonNgThesis.pdf` | App. B (p. 316–317), App. C (p. 318–325), Ch. 2.4 (p. 35–40), Ch. 4.1 (p. 76–96) read in full |
| Gresh et al., J. Mol. Spectrosc. 319, 1 (2016) | arXiv:1509.03682 | read in full incl. Tables 1–3 (rendered) |
| Gresh PhD thesis (JILA), 145 pp. | Zotero `YM4P3CV6/dgresh_thesis.pdf` | TOC only |
| Leanhardt et al., J. Mol. Spectrosc. 270, 1 (2011) | arXiv:1008.2997 | Sec. III–IV read in full (pp. 12–25 of the arXiv PDF) |
| Petrov & Skripnikov, arXiv:2503.02840 (2025) | arXiv | read in full, pp. 3–4 rendered |
| Petrov, Skripnikov & Titov, arXiv:2302.02856 (PRA 107, 062814) | arXiv | read in full |
| Petrov, Skripnikov & Titov, PRA 96, 022508 (2017) | arXiv:1704.06631 | Theory section read in full |
| Skripnikov & Titov, PRA 91, 042504 (2015) | arXiv:1503.01001 | Tables I–II + PT-odd definitions read |
| Denis et al., New J. Phys. 17, 043005 (2015) | Zotero `KC5JC4TW` | abstract, Tables 5/10/11, conclusions |
| Cairncross et al., PRL 119, 153001 (2017) + Supplement | Zotero `E7RPWZE2` (18 pp., supplement included) | scanned; supplement is data collection/systematics, **not** the Hamiltonian |
| Cairncross PhD thesis (JILA), 260 pp. | Zotero `JKS9NAEE` | Ch. 2 pp. 36–50 extracted |
| Roussy et al., Science 381, 46 (2023) | arXiv:2212.11841 + Zotero `9VBXBDEZ` | skimmed only |
| Zhou et al., PRL 124, 053201 (2020) | Zotero `R553MJBH` | skimmed only |
| Loh, Grau, Stutz theses (JILA) | Zotero | copied, not read |

**Not obtained:** Loh et al., Science 342, 1220 (2013) and its supplement (paywalled; no arXiv version located — searched `au:Petrov AND ti:"Zeeman interaction"` and the ThF⁺ title queries above; the Loh thesis in Zotero is the substitute and was not read). Caldwell et al., PRA 108, 012804 (2023) (systematics; cited by Petrov 2025 as [51]). Petrov, PRA 108, 062804 (2023) (cited as [23] for the S₂/S₃ formulas). Skripnikov JCP 147, 021101 (2017) — this is HfF⁺ E_eff, not ThF⁺, so it was deliberately not pursued. Fleig & Nayak: the Zotero copy (`TCMGSTQ9`) is **ThO**, not ThF⁺; no Fleig ThF⁺ E_eff paper was located separately from Denis et al. 2015 (on which Fleig is a co-author). The Tier-2 item "Rotational splittings in diatomic molecules of interest to searches for new physics (2025)" was not pursued (§9).

---

## 3. Hamiltonian form by source

### 3.1 Ng PhD thesis, Appendix C, PDF pp. 318–325 — **the ThF⁺ statement**

Frames (App. C.1, p. 318): two frames. (i) The quantization-axis frame **rotates with E_rot**; F, m_F, J are defined there. (ii) The molecule frame, quantization axis along the internuclear axis **"pointing towards thorium"**; Ω is defined there. (Note: Ch. 2.4, p. 35, defines Ω = J_a·n̂ with n̂ "pointing towards the cation" — that is the *neutral ThF* chapter, so the two statements are consistent only for ThF⁺ where Th is the cation end.)

Basis (App. C.3.1, p. 323): **{|J, Ω, F, m_F⟩}**; 12 states for J = 1 alone, 32 states for J = 1 & 2. Time-dependent propagation done with QuTiP.

**H = H_hyperfine + H_Ω-doub + H_rot + H_Stark + H_Zeeman + H_offset(v,J) + H_eEDM**  (Eq. C.1, p. 318)

1. **Hyperfine** (Eq. C.2, p. 319): `H_hf = A∥ (F² − I² − J²) / (2 J²)`, from the projection theorem applied to E ∼ (I·n)(J·n). Only ¹⁹F contributes (²³²Th has I = 0). Eigenvalue: `A∥ [F(F+1) − I(I+1) − J(J+1)] / (2 J(J+1))` [inferred, elementary]. Consequences stated in the source: ΔE_hf(J=1)/ΔE_hf(J=2) = 9/5 (p. 319); J = 1 splitting = (3/4)|A∥|, J = 2 splitting = (5/12)|A∥| (Ng 2022 Fig. 2, p. 3 — verified by rendering that page). **No dipolar, contact, or quadrupole hyperfine term appears** — the whole ¹⁹F hyperfine structure is one constant A∥.
2. **Ω-doubling** (Eq. C.3, p. 319): `H_Ω-doub = ⊕_{J≥Ω} [(−1)^J / 2] ℏ ω_ef Ω_x^{(J)}` with
   `Ω_x^{(J)} = [J(J+1)/2] ( |J,Ω=+1⟩⟨J,Ω=−1| + |J,Ω=−1⟩⟨J,Ω=+1| )`.
   Parity states: `|Ω = ±1⟩ = (|+⟩ ± |−⟩)/√2`. The `(−1)^J` "accounts for the alternating parity with J". Splitting at level J = `ω_ef · J(J+1)/2` [inferred from the operator; verified against Ng 2022 Fig. 2, which draws ω_ef at J = 1 and 3ω_ef at J = 2, and against Gresh's k, §7.1]. ω_ef is *defined* as the J = 1 splitting.
3. **Rotating-frame ("rotation") term** (Eq. C.4, p. 320): `H_rot = ℏ ω_rot F_x`, with **z along E_rot and x anti-parallel to the (counter-clockwise) rotation vector**. Couples Δm_F = ±1 within the same F, J, Ω. "This is the part of the Hamiltonian that is responsible for Berry's phase."
4. **Stark** (Eq. C.5, p. 320): `H_Stark = −d_mf n̂·E`, written out as
   `−d_mf Σ_{p∈{0,±1}} E_p (−1)^{F+J′+1+I} {J F I; F′ J′ 1} × (−1)^{F′−m′_F} √((2F+1)(2F′+1)) (F′ 1 F; −m′_F p m_F) × (−1)^{J′−Ω′} √((2J+1)(2J′+1)) (J′ 1 J; −Ω′ q=0 Ω)`
   (round brackets = 3j, curly = 6j). q = 0 because the dipole lies along the internuclear axis. Mixes different J and F at fixed Ω and m_F. Derivation pointer: Brown & Carrington p. 167.
5. **Zeeman** (Eq. C.6, p. 321): `H_Zeeman = −G∥ μ_B (J·n̂)(n̂·B) − g_N μ_N I·B`, with the same 3j/6j expansion for (J·n̂)(n̂·B) as for the Stark term, plus a standard I·B expansion. The relation to the hyperfine g-factor is given explicitly (p. 321):
   `g_F = −G∥ [F(F+1) + J(J+1) − 3/4] / (2F(F+1)J(J+1)) + g_N (μ_N/μ_B) [F(F+1) − J(J+1) + 3/4] / (2F(F+1))`
   which "evaluates to (1/3)(−G∥ + g_N μ_N/μ_B) for F = 3/2, J = 1". Derivation pointer: Brown & Carrington p. 606.
6. **Offset** (Eqs. C.7a–d, p. 322): `H_offset = G(v) + E(v,J)`, `G(v) = ℏω_e(v+½) − ℏω_eχ_e(v+½)²`, `E(v,J) = B_e J(J+1) − D[J(J+1)]² − α_e(v+½)J(J+1) = B_v J(J+1) − D[J(J+1)]²`. Note the rotational energy is `B J(J+1)` with **no −Ω² term** — the Hund's-case-(c) convention of Ch. 2.4 Eq. 2.2a (p. 37), `E(J) = B J(J+1)`, with the constant absorbed into the band origin. This is why the measured J = 1 → 2 interval is quoted as exactly `4B` (Ch. 4.1, p. 78: "we are able to extract 4B = 29.09733(4) GHz (see Equation 2.2a)").
7. **eEDM** (Eq. C.8, p. 322): `H_eEDM = −d_e E_eff Ω/|Ω|`, E_eff ≈ 35 GV/cm. **Footnote 3 on that page**: "There is an extra factor of 1/2 in the cited paper [Leanhardt 2011], which uses a different convention for how E_eff is defined in papers that calculate the value. The factor of 1/2 is removed for consistency with the cited value of E_eff." Also stated: Stark shifts go as m_F Ω while the eEDM term goes as Ω, so the states degenerate under Stark (|m_F=+3/2, Ω=+1⟩ and |m_F=−3/2, Ω=−1⟩) have opposite eEDM shifts.
8. **Two-level effective Hamiltonian** (App. C.4, pp. 324–325): defines `g = (g^u + g^ℓ)/2`, **`δg = (g^u − g^ℓ)/2`**, non-reversing field B_nr, Berry frequency `3 f_rot(S̃⟨α⟩ + α_nr)`, avoided crossings `Δ = (Δ^u + Δ^ℓ)/2`, `Δ_D = (Δ^u − Δ^ℓ)/2`, doublet switch D̃ with `Δ^{u,ℓ} = Δ + D̃ Δ_D`.

### 3.2 Leanhardt et al. 2011, J. Mol. Spectrosc. 270, 1 (arXiv:1008.2997) — the derivation, for ³Δ₁ ions generally

Basis/approach (Sec. IV A, p. 13): effective Hamiltonian within the ³Δ₁ manifold, "as elaborated by Brown and Carrington", citing Brown and Nelis for the ³Δ effective Hamiltonian.

`H_struct = H_elec + H_vib + H_SO + H_tum + H_SS + H_SR + H_HFS + H_LD`  (Eq. 10, p. 13), with (Eqs. 11–16, p. 13)

- `H_SO = A Λ Σ`
- `H_tum = B_e (J − S)² − D (J − S)⁴`
- `H_SS = (2/3) λ (3Σ² − S²)`
- `H_SR = γ_SR (J − S)·S`
- `H_HFS = a I_z L_z + b_F I·S + (c/3)(3 I_z S_z − I·S) + ½ e_Δ (J₊I₊S₊² + J₋I₋S₋²)`
- `H_LD = ½ (o_Δ + 3p_Δ + 6q_Δ)(S₊²J₊² + S₋²J₋²)`

Stated simplifications (p. 13–14): the first four terms "primarily describe an overall shift of the ³Δ₁ J-level, and can be ignored in evaluating energy differences in the states we care about". The e_Δ hyperfine term "is expected to be even smaller than the already small Λ-doublet splitting itself, however, and will be ignored". For Λ-doubling, only the Ω = +1 ↔ Ω = −1 term is kept; `|o_Δ + 3p_Δ + 6q_Δ| = õ_Δ ≈ Σ C_{Π,Σ,Π′} A²B_e² / [(E_Δ−E_Π)(E_Δ−E_Σ)(E_Δ−E_Π′)]` (Eq. 17, p. 13), and **"we shall express the energy difference in parity levels for the J = 1 as ω_ef = 4õ_Δ"** (p. 14). Off-diagonal-in-Ω couplings are suppressed by ~10⁻⁶ (hyperfine/spin-orbit) and dropped (p. 14).

Net statement (p. 14): "the basic molecular structure of interest to the ³Δ₁, J = 1 state is governed by two constants: the hyperfine splitting E_hf (given by 3A∥/4 for J = 1, I = 1/2) and the Λ-doublet splitting ω_ef."

External fields (Sec. IV B, pp. 14–15):
- `H_Stark = −d⃗_mf·E⃗` (Eq. 18); in the fully polarized limit d_mf E ≫ ω_ef, `E_Stark = −m_F Ω γ_F d_mf E` (Eq. 20) with **`γ_F = [J(J+1) + F(F+1) − I(I+1)] / (2F(F+1)J(J+1))`** (Eq. 21), γ_{F=3/2} = 1/3, γ_{F=1/2} = 2/3. Sign statement (p. 14): "The electric field therefore raises the energy of the states with m_F Ω < 0 (denoted 'upper' states, superscript u), and lowers the energy of states with m_F Ω > 0 ('lower')."
- `H_Zeeman = −μ⃗·B⃗` (Eq. 19); leading Ω-preserving terms:
  `H_Zeeman = ( γ_F [ ((g_L+g_r)Λ + (g_S+g_r)Σ) Ω − g_r J(J+1) ] − g_I κ_F ) m_F μ_B B`  (Eq. 22, p. 14)
  with `κ_F = [F(F+1) + I(I+1) − J(J+1)] / (2F(F+1))`. Simplified: `g_{F=3/2} = γ_{F=3/2}(g_L Λ + g_S Σ)Ω` (Eq. 23) and `g_{F=1/2} = 2 g_{F=3/2}` (Eq. 24).
- `H_EDM = −d⃗_e·E⃗_eff = d_e E_eff σ⃗₁·n̂` (Eq. 25) with **n̂ "pointing from the more negative atom to the more positive one; in our case from the fluorine or hydrogen to thorium, platinum, or hafnium"** and "E_eff as positive if it is anti-parallel to n̂"; `E_EDM = −d_e E_eff Ω/(2|Ω|)` (Eq. 26, p. 15). **The 1/2 that Ng removes.**
- Linearized level energies (Eq. 28, p. 15):
  `E^{u/ℓ}_nr(F, m_F, Ω; E, B) ≈ (1/3)(F(F+1) − 11/4) E_hf − m_F Ω γ_F d_mf E + m_F g^{u/ℓ}_F μ_B B − (d_e E_eff/2|Ω|)Ω`
  valid for ω_ef ≪ d_mf E ≪ E_hf, d_mf E ≫ g_F μ_B B∥, B⊥ = 0.
- State labels (Eq. 29, p. 15): |a⟩ = |F=3/2, m=+3/2, Ω=−1⟩, |b⟩ = |3/2, −3/2, +1⟩, |c⟩ = |3/2, +3/2, +1⟩, |d⟩ = |3/2, −3/2, −1⟩; `W^u = E_a − E_b = 3g^u_F μ_B B + d_e E_eff`, `W^ℓ = E_c − E_d = 3g^ℓ_F μ_B B − d_e E_eff` (Eq. 30).
- Strong-field F-mixing (Eqs. 31–32, pp. 15–16): closed forms for E_nr(F̃∼3/2 or 1/2, m_FΩ = ±1/2) when d_mf E is comparable to E_hf, including a ∓3ω_ef²/(4 d_mf E) term.

Rotating field, large-angle/dressed-state limit (Sec. IV D, pp. 17–18):
- `H_dressed = H₀ − d⃗_mf·E⃗_rot + H_rot`, `H_rot = −ω_rot ( cos(θ) F_z − sin(θ) F_x )` (Eqs. 36–37, p. 17). At θ = π/2 this reduces to Ng's `+ℏω_rot F_x` [inferred].
- Geometric phase: `2π(1 − cosθ) m_F`; with α = π/2 − θ, `E_geo = −m_F ω_rot sin α ≈ −m_F ω_rot α` (Eqs. 38–39, p. 18).
- Third-order avoided crossing between |a⟩ and |b⟩: `Δ ∼ ω_ef (ω_rot/d_mf E_rot)^{2m_F}` (Eq. 41) with the numerically evaluated prefactor `Δ^{u/ℓ} ≈ 170 ω_ef (ω_rot / d_mf E_rot)³` for m_F = ±3/2 (Eqs. 42–43, p. 18).

δg_F contributions (Sec. IV G, pp. 24–25):
- Zero-field: the Ω → −Ω Zeeman terms omitted from Eq. 22, written by Brown et al. and Nelis et al. as
  `H_ZeemanDist = −½ g_rS μ_B (B₊J₋S₊S₋ + B₋J₊S₋S₊)`, `H_ZeemanDoub = +½ g′_rS μ_B (B₊J₊S₊² + B₋J₋S₋²)` (Eq. 64, p. 24); H_ZeemanDoub is parity-dependent, `|g′_rS| ≈ ω_ef/(2B_e)` (Eq. 65). ~10⁻⁶ for HfF⁺.
- E-field induced (the dominant one in the operating regime):
  `δg_F(E_rot) = Σ_{J′,F′} [d_mf E_rot / (B_e (J+1))] (g_F/γ_F) Ω [F,F′,J,J′]² (F 1 F′; −m_F 0 m_F)² (J 1 J′; −Ω 0 Ω)² {F′ J′ I; J F 1}²`  (Eq. 66, p. 24)
  evaluating for the J = 1, F = 3/2 case to **`δg_{F=3/2}/g_{F=3/2} = 9 d_mf E_rot / (40 B_e)`** (Eq. 67, p. 24). Same scaling applies to γ_F: δγ_F/γ_F ≈ 10⁻⁴.

### 3.3 Petrov & Skripnikov, arXiv:2503.02840 (2025) — ThF⁺, numerical-diagonalization form

Basis (Eq. 1, p. 2): `Ψ_Ω · θ^J_{M,Ω}(α,β) · U^F_{M_I}`, `θ^J_{M,Ω} = √((2J+1)/4π) D^J_{M,Ω}(α,β,γ=0)`; M_I = ±1/2; **M_F = M_I + M** (projection on lab ẑ) which, they warn in the companion HfF⁺ paper (arXiv:2302.02856, p. 3), "is not equal to m_F", the projection on the *rotating* field direction.

`Ĥ_mol = Ĥ_el + Ĥ_rot + Ĥ_hfs + Ĥ_ext` (Eq. 2, p. 2). The explicit operators are in Petrov, Skripnikov & Titov PRA 96, 022508 / arXiv:1704.06631, Eqs. 5–8 (p. 2 of that PDF):
- `Ĥ_rot = B_rot J² − 2 B_rot (J·J⃗ᵉ)`  (Eq. 6)
- `Ĥ_hfs = g_F I · Σ_i (α_i × r_i)/r_i³`  (Eq. 7), α = Dirac matrices, g_F = 5.25773 the ¹⁹F nuclear g-factor
- `Ĥ_ext = μ_B (L⃗ᵉ − g_S S⃗ᵉ)·B − g_F (μ_N/μ_B) I·B − D⃗·E`  (Eq. 8), g_S = −2.0023
- Zeeman convention: `E_Zeeman = −g μ_B B M_F` (Eq. 1); the adiabatic zero-hyperfine g-factor is
  `g = −G∥ [F(F+1)+J(J+1)−3/4]/(2F(F+1)J(J+1)) + g_F (μ_N/μ_B) [F(F+1)−J(J+1)+3/4]/(2F(F+1))` (Eqs. 2–3) — **identical to Ng thesis Eq. C.6's companion formula.**
- `G∥ = ⟨³Δ₁| L̂ᵉ_n̂ − g_S Ŝᵉ_n̂ |³Δ₁⟩` (arXiv:2302.02856 Eq. 10, p. 3) — the definition of G∥.

ThF⁺ electronic basis actually used (arXiv:2503.02840 Eq. 3, p. 2): ³Δ₁ (T₀ = 0), ¹Σ⁺ (314 cm⁻¹), 1³Δ₂ (1052 cm⁻¹), 2 Ω=2 (T_v = 5806), 3 Ω=2 (T_v = 6402), ³Π₀₋ (3044), ³Π₀₊ (3395). Off-diagonal electronic matrix elements Eqs. 4–15 (p. 2), e.g. `G^{(1)}_⊥ = ⟨³Δ₁|L̂ᵉ₊ − g_S Ŝᵉ₊|¹Σ⁺⟩ = 0.48` and `Δ^{(1)} = 2B_rot⟨³Δ₁|Ĵᵉ₊|¹Σ⁺⟩ = 0.250 cm⁻¹`. **The Ω-doubling is not a fitted constant in this approach**: matrix element (4) "provides the main contribution to the Ω-doubling effect. It was adjusted to reproduce the experimental value of 5.29 MHz … we increased the calculated value by 8.8%" (p. 3). Δg parametrization near a working field (Eq. 18, p. 4): `Δg = Δg₀ + Δg₁ E`, with **Δg ≡ g^u − g^ℓ**.

Rotating fields (arXiv:2302.02856 Eqs. 7–8, p. 3): `E_rot(t) = E_rot(x̂ cos ω_rot t + R̃ ŷ sin ω_rot t)`, same for B_rot(t); R̃ = ±1 sets the rotation sense, ω⃗_rot = R̃ ω_rot ẑ; ω_rot and E_rot always positive. Handled by transformation to the rotating frame.

### 3.4 Skripnikov & Titov, PRA 91, 042504 (2015) (arXiv:1503.01001) — the PT-odd operators for ThF⁺

- `W_d = (1/Ω) ⟨Ψ| Σ_i H_d(i)/d_e |Ψ⟩` (Eq. 1, p. 2), with `H_d = 2d_e [[0,0],[0,σ·E]]` (Eq. 2), E the inner molecular field, **Ω = ⟨Ψ|J·n|Ψ⟩ with n the unit vector along the molecular axis directed from Th to F** (Ω = +1 for ³Δ₁). **`E_eff = W_d |Ω|`** (p. 2).
- Scalar–pseudoscalar electron–nucleon: `H_{T,P} = i (G_F/√2) Z k_{T,P} γ⁰γ⁵ ρ_N(r)` (Eq. 3, p. 2); the tabulated coefficient is `W_{T,P}`, in kHz, so that the energy shift is `W_{T,P} k_{T,P}` [inferred from the table units].

### 3.5 Cairncross et al., PRL 119, 153001 (2017) + Supplemental Material

Ng 2022 (p. 1) points the reader to "the Supplementary Material of Ref. 14" for the effective Hamiltonian. **That is not what the supplement contains.** The 12-page supplement (`cairncross2017-PRL-HfFplus-eEDM.pdf`, pp. 7–18) covers switch-state timing, data processing, the eight-channel parity decomposition (Eq. S1), blinding, cuts, frequency-channel modelling and systematics. No structural Hamiltonian is written there. The Hamiltonian statement to use is Leanhardt 2011 or the Ng thesis appendix.

---

## 4. Parameter table for ThF⁺ X ³Δ₁

Isotopologue throughout: **²³²Th¹⁹F⁺** (²³²Th: I = 0; ¹⁹F: I = 1/2), v = 0 unless noted. Exp. = measured; Th. = computed.

| Symbol | Value ± uncertainty (as printed) | Units as printed | In MHz / cm⁻¹ | Source (paper, table/eq., page) | Exp/Th | Notes |
|---|---|---|---|---|---|---|
| A∥ (¹⁹F magnetic hyperfine) | −20.1(1) | MHz (quoted as A∥/2π) | −20.1(1) MHz | Ng 2022 Table I, p. 5; Ng thesis Table 4.1, p. 87 | Exp | Negative ⇒ F = 1/2 lies **above** F = 3/2 in J = 1 |
| A∥ (theory) | −21.5 | MHz | — | Ng 2022 Table I, p. 5 (X2CAMF-CCSD(T), L. Cheng) | Th | no uncertainty quoted |
| E_hf(J=1) = 3A∥/4 | −15.1 | MHz (2π×) | −15.075 MHz [conv.: ¾ × 20.1] | Ng thesis Table B.2, p. 317 | Exp (derived) | "hyperfine splitting" |
| ω_ef (Ω-doubling, J = 1) | 5.29(5) | MHz (quoted as ω_ef/2π) | 5.29(5) MHz | Ng 2022 Table I, p. 5 | Exp | defined as the **J = 1** parity splitting |
| ω_ef (earlier) | 5.21(4) | MHz | 5.21(4) MHz | Ng 2022 Table I, p. 5, quoting Gresh 2016 | Exp | = 2k, see next row |
| k (Λ/Ω-doubling, ∝ J(J+1)) | 0.869(7) × 10⁻⁴ | cm⁻¹ | 2.605(21) MHz [conv.: ×29979.2458 MHz/cm⁻¹] | Gresh 2016 §3.1 p. 5 and Table 1 p. 10 (Ω=0⁺[10.47]←³Δ₁ (0,0), column k″) | Exp | ω_ef = 2k → 5.210(42) MHz [conv.], matching the 5.21(4) quoted above |
| k, other bands | −0.887(6), −0.90(1), −0.876(3), −0.882(6), −0.882(5) ×10⁻⁴ (Ω=0⁻ bands); +0.874(13), +0.892(12), +0.884(14), +0.891(13) ×10⁻⁴ (Ω=0⁺ bands) | cm⁻¹ | — | Gresh 2016 Table 1, p. 10 | Exp | sign is relative to the s = ±1 (e/f) labelling of the *upper* state, not physical |
| 4B (J = 1 → 2 interval) | 29.09733(4) | GHz | — | Ng 2022 §II B p. 2; Ng thesis §4.1.1 p. 78 | Exp | quoted explicitly as "4B" |
| B₀ | — | — | 7.2743325(10) GHz = 0.242647 cm⁻¹ [conv.: 29.09733/4] | derived from the row above | Exp (derived) | matches the next two rows |
| B₀ | 0.24264(3) | cm⁻¹ | 7274.1(9) MHz [conv.] | Gresh 2016 §3.1 p. 5 and Table 1 p. 10 | Exp | 95% CI |
| B_rot (as used in the ThF⁺ ab initio model) | 0.243 | cm⁻¹ | — | Petrov & Skripnikov arXiv:2503.02840, p. 3 | (input) | cites Ng 2022 |
| B_e | 0.24311(7) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | |
| α_e | 1.00(4) × 10⁻³ | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | B_v = B_e − α_e(v+½) |
| D₀ (centrifugal) | 1.30(4) × 10⁻⁷ | cm⁻¹ | 3.90(12) kHz [conv.] | Gresh 2016 Table 1, p. 10, column D″ for Ω=0⁺[10.47]←³Δ₁ (0,0) | Exp | other bands give 1.27–1.40 ×10⁻⁷ |
| ω_e | 656.96(1) | cm⁻¹ | 19.6934 THz [conv.] | Gresh 2016 Table 2, p. 11 (and §3.1 p. 5) | Exp | computed without a Morse assumption |
| ω_eχ_e | 1.920(3) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | |
| v = 0 → 1 spacing | 653 | cm⁻¹ | — | Ng 2022 Fig. 6, p. 5 | Exp | = ω_e − 2ω_eχ_e = 653.1 [conv., consistency check] |
| R_e | 3.75 (CCSD(T)); 3.74(4) (exp., Barker) | a₀ | 1.984 Å [conv.] | Skripnikov & Titov 2015 Table I, p. 7 | Th/Exp | Denis 2015 works at R = 3.779 a₀ |
| d_mf (molecule-frame dipole) | 3.37(9) | D | 1.70 MHz/(V/cm) [Ng thesis Table B.2 p. 317; my check: 3.37 × 0.50341 MHz/(V/cm)/D = 1.696 ✓] | Ng 2022 Table I, p. 5 | Exp | referenced to the **centre of mass**, see §5.4 |
| d_mf (theory, this-work CCSD(T)) | 3.46 | D | — | Ng 2022 Table I, p. 5 | Th | |
| d_mf (theory, Skripnikov & Titov) | 2.74 (w.r.t. **Th nucleus**) | D | 3.46 D w.r.t. centre of mass [conv., mine: +e·r(Th→c.m.) = 4.80320 D/Å × 0.1499 Å = 0.720 D; 2.74+0.72 = 3.46] | Skripnikov & Titov 2015 Table II, p. 8, and its footnote "*The dipole moment is calculated with respect to Th nucleus" | Th | reproduces the "3.46 [19]" that Ng 2022 quotes — the conversion is confirmed by that agreement |
| d_mf (theory, Denis) | 4.03 (4.029 for ³Δ₁) | D | — | Denis 2015 abstract p. 1 and Table 10; origin stated as the centre of mass (p. 1455 of the extract, table caption) | Th | ~16% above experiment |
| d(d_mf)/dr | 7(2) exp; 6.66 theory | D/Å | — | Ng 2022 §III p. 6 and §II D p. 5 | Exp/Th | sets the v = 1 lifetime |
| |g_{F=3/2}| | 0.0149(3) | — | — | Ng 2022 Table I, p. 5 | Exp | **sign not measured**; Ng thesis Table B.2 p. 317 adopts −0.0149 |
| g^u, g^ℓ (F = 3/2, |M_F| = 3/2) | ≈ −1.485 × 10⁻² and −1.494 × 10⁻² at 50–100 V/cm | — | — | Petrov & Skripnikov arXiv:2503.02840 Fig. 2, p. 4 (read from the rendered figure) | Th | **negative**, in the convention E_Zeeman = −gμ_B B M_F |
| G∥ | 0.048(2) if g_F < 0; −0.042(2) if g_F > 0 | — | — | Ng 2022 §II E, p. 5; Ng thesis p. 87 | Exp (derived from |g_F|) | rotational contribution neglected |
| G∥ (theory) | 0.034 (Skripnikov & Titov 2015); 0.035 (Ng 2022/L. Cheng); 0.047 (Petrov 2025, chosen to reproduce g = 0.0149 with non-adiabatic mixing) | — | — | Skripnikov & Titov 2015 Table II p. 8; Ng 2022 §II E p. 5; arXiv:2503.02840 p. 3 | Th | all positive ⇒ g_F < 0 |
| δg_{F=3/2} (paper convention) | |δg| = 0.0003(3) | — | — | Ng 2022 Table I, p. 5 | Exp | "difference in magnetic g-factors between the upper and lower doublets" |
| δg_F/g_F at E_rot = 60 V/cm | −0.00255(6) exp; −0.00223 from the 32-level model | — | — | Ng thesis §4.1.3 p. 85 and Table 4.1 p. 87 | Exp/Th | thesis convention δg = (g^u−g^ℓ)/2 (§3.1 item 8) |
| Δg₀, Δg₁ (Δg = g^u − g^ℓ = Δg₀ + Δg₁E) | at E = 60 V/cm: 10⁷Δg₀ = 233.8, 10⁷Δg₁ = 8.7 cm/V (full table 40–150 V/cm) | — | — | Petrov & Skripnikov arXiv:2503.02840 Table I, p. 5 | Th | ⇒ Δg(60 V/cm) = 7.56 × 10⁻⁵ [conv.] |
| Δg at E = 0 | 2.3 × 10⁻⁴ | — | — | arXiv:2503.02840 p. 4 | Th | claimed accuracy "not worse than 10%" (p. 3) |
| E_eff | 37.3 | GV/cm | — | Skripnikov & Titov 2015 Table II, p. 8 (FINAL) | Th | ±7% stated (p. 7) |
| E_eff | 35.2 | GV/cm | — | Denis 2015 abstract, p. 1 | Th | |
| E_eff (value adopted by JILA) | ≈ 35 | GV/cm | — | Ng thesis Eq. C.8 discussion, p. 322 | Th | |
| W_{T,P} (scalar–pseudoscalar, k_{T,P} coefficient) | 50 | kHz | — | Skripnikov & Titov 2015 Table II, p. 8 | Th | ±7% |
| W_{T,P} | 48.4 | kHz | — | Denis 2015 abstract, p. 1 | Th | |
| W_M (nuclear MQM) | 0.88 × 10³³ | Hz/(e·cm²) | — | Skripnikov & Titov 2015 Table II, p. 8 | Th | **irrelevant for ²³²Th (I = 0)**; applies to ²²⁹ThF⁺ |
| A∥(²²⁹Th) | −4163 | (μ_Th/μ_N)·MHz | — | Skripnikov & Titov 2015 Table II, p. 8 | Th | **Th, not F.** Not to be confused with the ¹⁹F A∥ |
| A∥(²²⁹Th, I = 5/2) | 1833 | MHz | — | Denis 2015 abstract p. 1, conclusions p. 1653 | Th | same caveat; "awaits confirmation" |
| d_e E_eff | 2π × 0.851 μHz for d_e = 10⁻³¹ e·cm | — | — | Ng thesis Table B.2, p. 317 | derived | my check with E_eff = 35 GV/cm gives 0.846 μHz [conv.] |
| Radiative lifetime of X ³Δ₁ | ∞ (it is the electronic **ground** state) | — | — | Ng 2022 §III p. 6, quoting Gresh 2016 | Exp | Gresh 2016 §3.1 p. 6 establishes the ordering from the Ω=0⁺[10.47] band positions |
| τ(X ³Δ₁, v = 1 → 0) | 0.16(11) | s | — | Ng 2022 §III p. 6 | Exp | v = 2 → 1: 0.08(6) s |
| τ(a ¹Σ⁺, v = 0 → X) | ≈ 6 | s | — | Ng 2022 §III p. 7 | Exp | |
| T₁ (blackbody-limited) | 3.2 s (300 K), 13 (200 K), 21 (180 K), 51 (150 K), 160 (120 K), 2100 (77 K) | s | — | Ng 2022 Table II, p. 7 | Exp+model | combined v = 1 and a¹Σ⁺ channels |
| T₀(a ¹Σ⁺) | 314.282(7) | cm⁻¹ | 9.4224 THz [conv.] | Gresh 2016 Table 2, p. 11 | Exp | T_e = 314.0(2); Ng 2022 Fig. 6 quotes 314 |
| T₀(1³Δ₂) | 1052.5(1.0) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 (from Barker et al.) | Exp | B_e = 0.24342(2), ω_e = 657.38(22) |
| T₀(³Δ₃) | 3150(30) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | |
| T₀(³Π₀₋), T₀(³Π₀₊) | 3044, 3395 | cm⁻¹ | — | arXiv:2503.02840 Eq. 3, p. 2 | Th | vertical, R = 3.75 a₀ |
| T₀(Ω = 0⁺[10.47]) | 10471.889(7) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | couples to **both** X ³Δ₁ and a ¹Σ⁺ (Gresh Fig. 4, p. 6) |
| T₀(Ω = 0⁻) | 14589.09(2) | cm⁻¹ | — | Gresh 2016 Table 2, p. 11 | Exp | the optical-pumping state; Ng 2022 says "approximately 14600 cm⁻¹" (p. 2) |
| ¹⁹F nuclear g-factor | 5.25773 | — | — | Petrov et al. arXiv:1704.06631, Eq. 3 discussion, p. 1 | (constant) | enters g_F via g_N μ_N/μ_B = 2.8634 × 10⁻³ [conv.] |
| B_v with v-dependence beyond α_e | **not found** (searched: Gresh 2016 Tables 1–3, Ng thesis Ch. 2 and App. C, Ng 2022) | | | | | only B_e, α_e, D₀ are reported |
| Sign of g_F | **not found** (searched: Ng 2022 §II E p. 5, Ng thesis §4.1.4 p. 87, arXiv:2503.02840 p. 4 — all state the experiment is insensitive to it) | | | | | theory implies g_F < 0 [inferred, §5.7] |
| ¹⁹F dipolar / contact hyperfine constants (a, b_F, c, e_Δ separately) | **not found** for ThF⁺ (searched: Ng 2022, Ng thesis App. C, Gresh 2016, arXiv:2503.02840) | | | | | the model uses only A∥; Leanhardt 2011 p. 13 states e_Δ is dropped |
| Second-order Zeeman / g_r (rotational g-factor) for ThF⁺ | **not found** (Ng 2022 §II D p. 5 says "It might be of interest to compute the rotational g-factor of ThF⁺" and notes it is ~6% of the total in ThO, citing Petrov et al. PRA 89, 062505) | | | | | |

**HfF⁺ values encountered — do not transfer.** ³Δ₁ T_e = 976.930 cm⁻¹; A∥(¹⁹F) = −62.0 MHz; g(J=1) = 0.00306; G∥ = 0.011768; D∥ = −1.53(2) a.u. w.r.t. centre of mass; B_rot = 0.2989 cm⁻¹; E_eff = 22.5(0.9)–24 GV/cm; Stark splitting at E_rot = 58 V/cm ≈ 114 MHz. (arXiv:2302.02856 Eqs. 9–13 and Table I, pp. 3–4; arXiv:1704.06631 Eq. 9, p. 2.)

---

## 5. Conventions to preserve

1. **Two frames, and m_F is defined in the rotating one.** Quantization axis = instantaneous direction of E_rot; F, m_F, J live there. Ω lives in the molecule frame. (Ng thesis App. C.1, p. 318.) Petrov's papers use M_F = M_I + M with M the projection on the *lab* ẑ, and explicitly warn "M_F is not equal to m_F" (arXiv:2302.02856, p. 3). If a code mixes the two it will silently mislabel every Zeeman sublevel.
2. **Direction of n̂ is not agreed between the two lineages.** JILA: n̂ points **from F to Th** ("from the more negative atom to the more positive one; in our case from the fluorine … to thorium", Leanhardt 2011 p. 15; "pointing towards thorium", Ng thesis p. 318). Petrov/Skripnikov: n̂ is "the molecular axis directed **from Th to F**" (arXiv:2503.02840 p. 3; Skripnikov & Titov 2015 p. 2). Ω, the signed d, and the sign convention for E_eff all flip with this choice. Pick one and assert it in the code.
3. **Signed d.** In the Petrov convention (n̂: Th→F) the body-fixed dipole is negative: D = −1.33 a.u. (arXiv:2503.02840 p. 3, but see §7.2 on the printed digit) ≡ −3.38 D [conv.]. In the JILA convention d_mf = +3.37 D and H_Stark = −d_mf n̂·E with d_mf > 0 (Ng thesis Eq. C.5, p. 320).
4. **Origin of the dipole moment matters — ThF⁺ is an ion.** The quoted d_mf is referenced to the **centre of (nuclear) mass** (arXiv:2503.02840 p. 3, explicitly; arXiv:2302.02856 Table I caption for HfF⁺, p. 4). Skripnikov & Titov 2015 tabulate 2.74 D w.r.t. the **Th nucleus** (Table II footnote, p. 8) — a 0.72 D offset [conv., mine]. Getting this wrong is a 20% error in every Stark shift.
5. **Ω-doubling constant ω_ef is the J = 1 splitting**, and the level-J splitting is ω_ef·J(J+1)/2 (Ng thesis Eq. C.3, p. 319; ω_ef at J = 1, 3ω_ef at J = 2 in Ng 2022 Fig. 2, p. 3). Equivalently ω_ef = 2k with Gresh's spectroscopic k (Gresh 2016 Table 1 caption, p. 10: "k and k_D, proportional to J(J+1) and J²(J+1)²"). In the microscopic language, ω_ef = 4õ_Δ where õ_Δ = |o_Δ + 3p_Δ + 6q_Δ| (Leanhardt 2011 p. 14). Three different constants, three different factors — do not interchange them.
6. **Parity of the Ω-doublet alternates with J.** The `(−1)^J` prefactor in Ng thesis Eq. C.3 (p. 319) encodes this; |Ω = ±1⟩ = (|+⟩ ± |−⟩)/√2. In Ng 2022 Fig. 2 (verified from the rendered page): black = positive parity, grey = negative; at J = 1 the negative-parity component is the **upper** one of each doublet, at J = 2 the positive-parity component is upper.
7. **Sign of g_F.** The experiment measures |g_{F=3/2}| = 0.0149(3) only (Ng 2022 §II E, p. 5). Ng gives the two branches: G∥ = −0.042(2) if g_F > 0, +0.048(2) if g_F < 0, and notes the latter is close to theory (0.034, 0.035). Petrov & Skripnikov's calculation shows g^{u,ℓ} ≈ −1.49 × 10⁻² directly (Fig. 2, p. 4, verified from the render). **So g_{F=3/2} < 0** [inferred, but supported by both theory routes]; Ng's own thesis Table B.2 (p. 317) writes −0.0149 with the note "Sign not measured". My arithmetic check: (1/3)(−0.048 + 5.25773/1836.15) = −0.01505 [conv.], matching the measured 0.0149(3) to 1%.
8. **Zeeman shift convention:** `E_Zeeman = −g μ_B B M_F` (arXiv:2503.02840 Eq. 16, p. 4; arXiv:1704.06631 Eq. 1, p. 1; "This definition matches the ones in [Ng 2022] and [Petrov et al. 2017]"). Ng thesis Eq. C.6 (p. 321) is the same, written as an operator. Note this is the *opposite* overall sign from Leanhardt's Eq. 22, which writes `E = (…) m_F μ_B B` with the g-factor absorbed differently — Leanhardt's g_{F=3/2} = γ_F(g_LΛ + g_SΣ)Ω has no leading minus.
9. **Upper vs lower doublet.** Two equivalent statements: (a) "the electric field raises the energy of the states with m_F Ω < 0 (denoted 'upper')" (Leanhardt 2011 p. 14); (b) "For upper (lower) Stark doublet n is parallel (antiparallel) to the external rotating electric field", with n̂ from metal to F (arXiv:2302.02856 Fig. 1 caption, p. 2). Ng thesis (p. 323) just says "we call the pair that is more (less) energetic the 'upper doublet' ('lower doublet')". The experimentalists' doublet switch is D̃ = +1 for **lower** (arXiv:2302.02856, p. 2).
10. **δg vs Δg differ by a factor 2 between the thesis and the theory papers.** Ng thesis App. C.4 (p. 324): `δg = (g^u − g^ℓ)/2`. Petrov & Skripnikov (Eq. 18, p. 4): `Δg = g^u − g^ℓ`. See §7.3 for the numerical check that confirms this is a real factor-2 convention split, not an error.
11. **The eEDM term differs by a factor 2 between Leanhardt and everyone since.** Leanhardt Eq. 26 (p. 15): `E_EDM = −d_e E_eff Ω/(2|Ω|)`. Ng thesis Eq. C.8 (p. 322): `H_eEDM = −d_e E_eff Ω/|Ω|`, with an explicit footnote saying the 1/2 was removed for consistency with the convention under which E_eff is *calculated*. The observable relation everyone now uses is `f^{BD} = 2 d_e E_eff` (arXiv:2302.02856 Eq. 3, p. 2; arXiv:2503.02840 p. 4).
12. **Rotational energy is B J(J+1), with no −Ω² term** (Ng thesis Eq. 2.2a p. 37 and Eq. C.7c p. 322; Hund's case (c), J ≥ Ω). This is why "4B = 29.09733(4) GHz" is a literal statement about the J = 1 → 2 interval (Ng thesis p. 78) and why B₀ from Gresh matches it to 5 × 10⁻⁵ relative [conv., §7.1].
13. **Hyperfine sign convention.** A∥ < 0 with `H_hf = A∥(F² − I² − J²)/(2J²)` puts F = 1/2 **above** F = 3/2 in J = 1 (Ng thesis Eq. C.2 p. 319; drawn that way in Ng 2022 Fig. 2, p. 3).
14. **No explicit statement of divergence from Brown & Carrington or PGopher** appears in any source read. Ng thesis cites B&C p. 167 (Stark) and p. 606 (Zeeman) as the derivation source for its own expressions; Leanhardt 2011 (p. 13) says its H_struct follows "Brown and Carrington" and Brown/Nelis's ³Δ effective Hamiltonian. Gresh 2016 uses the standard e/f parity convention `s = +1 for e-symmetry, s = −1 for f-symmetry` in Eq. 1 (p. 3).

---

## 6. Modeling lineage and the best reference to code from

**Leanhardt, Bohn, Loh, Maletinsky, Meyer, Sinclair, Stutz & Cornell, J. Mol. Spectrosc. 270, 1 (2011)** is the origin. It takes Brown and Nelis's ³Δ effective Hamiltonian, restricts it to the |Ω| = 1, J = 1 manifold, argues away everything except the hyperfine splitting E_hf = 3A∥/4 and the Ω-doublet ω_ef (Sec. IV A), adds Stark and Zeeman in the F basis with the Landé-like factors γ_F and κ_F (Eqs. 20–22), adds the eEDM shift (Eqs. 25–26), and then — the part nobody else had — works out the physics of the *rotating* bias field: the dressed-state Hamiltonian (Eq. 36), the geometric/Berry energy (Eqs. 38–39), the third-order avoided crossing Δ ≈ 170 ω_ef(ω_rot/d_mfE_rot)³ (Eq. 42), and the E-field-induced g-factor difference δg_F/g_F = 9 d_mf E_rot/(40 B_e) (Eq. 67). Everything the JILA experiment does downstream is a refinement of that document. Its parameters, however, are HfF⁺.

**Petrov, Skripnikov & Titov** built the parallel, numerical lineage: PRA 96, 022508 (2017) / arXiv:1704.06631 states the operator-level Ĥ_mol = Ĥ_el + Ĥ_rot + Ĥ_hfs + Ĥ_ext (Eqs. 5–8) in the Ψ_Ω·D^J_{M,Ω}·U^F_{M_I} basis and computes g^u(E), g^ℓ(E) for HfF⁺ by diagonalization over a small set of low-lying electronic states, with the Ω-doubling not fitted but generated by ⟨³Δ₁|L̂₊ − g_SŜ₊|Ω=0⟩ matrix elements tuned to the measured ω_ef. arXiv:2302.02856 (2023) validated the whole chain against the Roussy 2023 HfF⁺ data (all four f-channels), and arXiv:2503.02840 (2025) ported it to ThF⁺ with the ThF⁺ electronic manifold (³Δ₁, ¹Σ⁺, 1–3 Ω=2, ³Π₀±) and produced the E-dependent Δg₀, Δg₁ table.

**Extensions after 2011, in order:** (i) hyperfine mixing between J levels folded into δg and into the avoided crossing (Leanhardt Sec. IV G already; then numerically); (ii) the E-field-dependent g-factor, first analytic (Eq. 67) then ab initio (Petrov 2017 for HfF⁺, 2025 for ThF⁺); (iii) non-adiabatic mixing with other electronic states, which changes the G∥ needed to reproduce a given g (arXiv:2503.02840 p. 3: G∥ = 0.047 vs Ng's adiabatic 0.048); (iv) systematics parametrizations S₁, S₂, S₃ tying f^B/f^{BD} to g^u, g^ℓ, Δg₀, Δg₁ (arXiv:2503.02840 Eqs. 17–20, p. 4); (v) the ThF⁺-specific consolidation in the Ng thesis.

**Single best reference to code from: Ng thesis, Appendix C, pp. 318–325.** It is the only place where the complete term list, the basis, the frames, the 3j/6j matrix elements, the 12/32-state truncation, the parameter values, and the two-level reduction are stated together for ThF⁺. Pair it with Appendix B Table B.2 (p. 317) for the numerical parameter set, Ng 2022 Table I (p. 5) for the measured constants with uncertainties, and Leanhardt Sec. IV D–G for anything involving the rotating field that the appendix compresses into one line.

---

## 7. Substantive disagreements between sources

1. **Internal cross-checks that PASS** (worth stating because they validate the conventions above):
   - ω_ef = 2k: Gresh's k = 0.869(7) × 10⁻⁴ cm⁻¹ → 5.210(42) MHz [conv.], exactly the "5.21(4)" that Ng 2022 attributes to Gresh.
   - 4B: 29.09733(4)/4 = 7.274333 GHz = 0.2426466 cm⁻¹ [conv.] vs Gresh's B₀ = 0.24264(3) cm⁻¹ — agreement to 3 × 10⁻⁵ relative.
   - Hyperfine ratios: the operator A∥(F²−I²−J²)/(2J²) gives 3A∥/4 at J = 1 and 5A∥/12 at J = 2, ratio 9/5, matching both Ng 2022 Fig. 2 and the thesis text (p. 319).
   - G∥ → g_F: (1/3)(−0.048 + g_Nμ_N/μ_B) = −0.01505 [conv.] vs measured |g| = 0.0149(3).
   - d_mf origin shift: 2.74 D (Th nucleus) + 0.720 D [conv.] = 3.46 D, exactly the value Ng 2022 attributes to Skripnikov & Titov.
   - d_eE_eff: 10⁻³¹ e·cm × 35 GV/cm = 0.846 μHz [conv.] vs the thesis's 0.851 μHz (they used 35.2).
2. **Petrov & Skripnikov's printed body-fixed dipole for ThF⁺ appears to be a typo.** The paper prints "We used D = −0.133 a.u. (with respect to center of nuclear mass; the molecular axis directed from Th to F) and A∥ = −20.1 MHz from Ref. [11]" (arXiv:2503.02840 p. 3 — I rendered the page and read it as an image to be sure). −0.133 a.u. = −0.338 D [conv.], eight times smaller than the d_mf = 3.37(9) D they claim to be taking from the same reference. **−1.33 a.u. = −3.381 D** [conv.] matches Ng exactly, and matches the companion HfF⁺ paper's format (D∥ = −1.53(2) a.u., arXiv:2302.02856 Table I, p. 4). Treat the intended value as −1.33 a.u.; the calculation almost certainly used that. **[inferred]**
3. **δg vs Δg: a genuine factor-2 convention split, and the two are consistent once it is applied.** Ng thesis App. C.4 (p. 324) defines δg = (g^u − g^ℓ)/2 and reports δg/g = −0.00255(6) at E_rot = 60 V/cm (Table 4.1, p. 87). Petrov & Skripnikov Table I (p. 5) at E = 60 V/cm gives Δg = Δg₀ + Δg₁E = 2.338 × 10⁻⁵ + 8.7 × 10⁻⁷ × 60 = 7.56 × 10⁻⁵ [conv.]. Halving: 3.78 × 10⁻⁵. Ng: |δg| = 0.00255 × 0.0149 = 3.80 × 10⁻⁵ [conv.]. **Agreement to 0.5%.** Signs also agree: g < 0 and δg/g < 0 ⇒ δg > 0 ⇒ g^u > g^ℓ, which is what Petrov's Fig. 2 shows. The Ng 2022 *paper*, by contrast, defines δg_{F=3/2} in words as "the difference in magnetic g-factors between the upper and lower doublets" (Table I caption, p. 5) and reports 0.0003(3); Petrov quotes this as "|Δg| = 3(3) × 10⁻⁴" and compares it to their zero-field Δg = 2.3 × 10⁻⁴. So the paper and the thesis use the symbol δg with different factors. **Pin this in code with a comment.**
4. **eEDM Hamiltonian factor of 2** — Leanhardt Eq. 26 vs Ng Eq. C.8, described in §5.11. Ng's footnote is the resolution; the modern convention (E_eff as defined by the calculators, H = −d_e E_eff Ω/|Ω|, f^{BD} = 2d_eE_eff) is the one to code.
5. **E_eff: 35.2 GV/cm (Denis 2015) vs 37.3 GV/cm (Skripnikov & Titov 2015), a 6% spread**, both with ~7% claimed uncertainty (Skripnikov p. 7). Ng thesis adopts ≈35, Petrov 2025 adopts 37.3. Both supersede the older Meyer & Bohn estimate of 90 GV/cm, which Skripnikov calls "more than twice overestimated" (p. 7) and Denis calls "more than 60% smaller than" (conclusions).
6. **d_mf theory: 4.03 D (Denis) vs 3.46 D (Skripnikov & Titov, converted to the c.m. origin) vs 3.37(9) D measured.** Denis is 16% above experiment; Skripnikov is 3% above. Ng 2022 (p. 5) says its own 3.46 D "is in good agreement with calculations from previous work [Denis, Skripnikov]" without commenting on the Denis discrepancy.
7. **G∥: 0.034 (Skripnikov & Titov) vs 0.035 (Cheng, in Ng 2022) vs 0.047 (Petrov 2025) vs 0.048(2) (derived from the measurement).** Petrov attributes their higher value to including non-adiabatic interactions between ³Δ₁ and other electronic states (p. 3); Ng notes "we do not have a systematic estimate for the error in the theoretical value of G∥" (p. 5). The 30% gap between the two ab initio values is unresolved in the literature read.
8. **δg/g model discrepancy at 60 V/cm: measured −0.00255(6), 32-level model −0.00223** — a 15% gap that Ng attributes to "coupling between the X ³Δ₁ and the ³Δ₂ states that is not accounted for in the 32-level Hamiltonian; we saw a similar correction in HfF⁺" (thesis p. 85). Petrov's ThF⁺ model *does* include 1³Δ₂ and lands at 7.56 × 10⁻⁵ ⇒ δg/g = −0.00254 [conv.], i.e. essentially on the measurement. That is direct evidence that **³Δ₂ must be in the electronic basis** if the code is to reproduce the differential g-factor.
9. **Ground-state ordering.** Gresh 2016 (§3.1, p. 6) establishes experimentally that ³Δ₁ is the ground state and a ¹Σ⁺ lies 314.282(7) cm⁻¹ above it, overturning Barker et al. 2012 (who put ¹Σ⁺ lowest, with Ω = 1 higher by 65 cm⁻¹ — Denis 2015 p. 208 of the extract). Denis 2015 independently found Ω = 1 more than 300 cm⁻¹ below Ω = 0⁺ (abstract), and Skripnikov & Titov agree. Gresh's Table 3 shows Barker's T_e assignments shifted by exactly this reversal.

---

## 8. Evidence gaps and open primary-source questions

1. **Sign of g_{F=3/2}: never measured.** Everything in §5.7 is inference from the ab initio G∥ > 0. If the code exposes a signed g_F, make the sign a switch.
2. **No published ¹⁹F hyperfine structure beyond A∥.** No dipolar (c), contact (b_F), or Ω-doubling-coupled (e_Δ) hyperfine constants for ThF⁺ exist in any source read. Leanhardt 2011 (p. 13) justifies dropping e_Δ; nobody justifies dropping the anisotropic parts explicitly for ThF⁺ — the projection-theorem form of Eq. C.2 simply assumes the axial term dominates. Searched: Ng 2022 (full), Ng thesis App. C and Ch. 4, Gresh 2016 (full), arXiv:2503.02840 (full).
3. **No rotational g-factor (g_r) for ThF⁺.** Ng 2022 (p. 5) flags it as ~6% of the total in a comparable molecule (ThO, Petrov et al. PRA 89, 062505) and says it "might be of interest to compute". The measured |g_F| therefore contains an unseparated g_r contribution, and the G∥ extracted "neglecting the rotational contribution" (p. 5) carries that systematic.
4. **No centrifugal-distortion correction to ω_ef** (no k_D for the X ³Δ₁ state: Gresh Table 1 lists k_D = 0 for every ³Δ₁ lower state). Fine for J = 1,2; a problem if the code sweeps high J.
5. **No published D-value/higher-order rotational terms beyond D₀ for X ³Δ₁**, and no v-dependence of D.
6. **The value of Δ^u and Δ^ℓ (avoided crossings) for ThF⁺ is model-output, not measurement.** Ng thesis Table B.2 (p. 317) lists Δ^u = 2π × 0.511 Hz and Δ^ℓ = 2π × 1.43 Hz for E_rot = 60 V/cm, ω_rot = 2π × 147.5 kHz — obtained by diagonalizing the 12-level Hamiltonian (Fig. 4.8 caption, p. 85: "The energy difference at the avoided crossing is fixed by our ab initio calculations"). The Leanhardt closed form (Eq. 42) uses a numerical prefactor of 170 whose "form within perturbation theory is rather complicated" (p. 18). A kernel test that reproduces both would be a strong Hamiltonian-agnostic check.
7. **The Ng 2022 reference to "the Supplementary Material of Ref. 14" for the effective Hamiltonian is wrong or at least misleading** (§3.5) — worth telling Arian so he doesn't hunt for it.
8. **Loh et al., Science 342, 1220 (2013) supplement was not obtained.** It is cited by Ng 2022 as the source of the Ramsey/g-factor measurement procedure ("Following the procedure reported in Ref. 34", p. 4) and may contain a further Hamiltonian statement for HfF⁺. The Loh thesis (`notes/lit/loh-thesis-JILA.pdf`) is on disk as a substitute and was not read.
9. **Questions only Arian (or the JILA authors) can settle:** which n̂ convention he wants the code to adopt (§5.2); whether δg should follow the thesis (halved) or the theory-paper (unhalved) definition (§5.10); whether the eEDM term should carry the Leanhardt 1/2 (§5.11); and whether ³Δ₂ belongs in the electronic basis from the start (§7.8 says yes if δg matters).

---

## 9. Screening table and exclusions

| Source | Why included | System / regime | Method | Relevant result | Important limitation |
|---|---|---|---|---|---|
| Ng et al., PRA 105, 022823 (2022) (arXiv:2202.01346) | Tier 1; the only measurement of A∥, d_mf, g_F in ThF⁺ | ²³²Th¹⁹F⁺ X ³Δ₁, v = 0, J = 1,2; E ≤ ~60 V/cm | µwave + Ramsey spectroscopy on trapped ions; X2CAMF-CCSD(T) for comparison | Table I (p. 5): A∥ = −20.1(1), ω_ef = 5.29(5) MHz, d_mf = 3.37(9) D, |g_{3/2}| = 0.0149(3), |δg| = 0.0003(3); Table II blackbody T₁ | Does not write the Hamiltonian; the pointer it gives (Cairncross supplement) does not contain it either. g_F sign unmeasured |
| Ng PhD thesis (JILA), 325 pp. | The fullest ThF⁺ Hamiltonian statement | same | consolidation + QuTiP diagonalization | App. C (pp. 318–325): full H, basis, 3j/6j matrix elements, two-level reduction; App. B Table B.2 (p. 317): the numerical parameter set | δg convention differs from the paper (§7.3); H_eEDM footnote silently reverses Leanhardt |
| Leanhardt et al., JMS 270, 1 (2011) (arXiv:1008.2997) | Tier 1; the origin of the ³Δ₁ rotating-field formalism | HfF⁺ (parameters), ³Δ₁ generally (formalism) | effective-Hamiltonian derivation from Brown/Nelis + perturbation theory | Eqs. 10–43, 64–67: term list, γ_F, κ_F, Stark/Zeeman, Berry phase, Δ ≈ 170ω_ef(ω_rot/d_mfE)³, δg_F/g_F = 9d_mfE/(40B_e) | HfF⁺ numbers; eEDM term carries a factor 1/2 later dropped |
| Gresh et al., JMS 319, 1 (2016) (arXiv:1509.03682) | Tier 1; the rotational/vibrational constants and the state ordering | ThF⁺, all low-lying electronic states, J up to ~100 | frequency-comb + cw velocity-modulation spectroscopy | Tables 1–3: B₀, B_e, D₀, ω_e, ω_eχ_e, α_e, k, T₀ of ¹Σ⁺/³Δ₂/³Δ₃/Ω=0± ; ³Δ₁ established as the ground state | No hyperfine, no Ω-doublet resolved in the ion (k is a fit parameter, not a µwave measurement); strong high-J perturbations |
| Petrov & Skripnikov, arXiv:2503.02840 (2025) | Tier 1; the E-dependent g factors | ²³²ThF⁺, J = 1, F = 3/2, |M_F| = 3/2, E = 0–200 V/cm | FS-CCSD electronic matrix elements + numerical diagonalization | Fig. 2 (g^u, g^ℓ ≈ −1.49 × 10⁻²), Table I (Δg₀, Δg₁, S₁–S₃) | Printed D = −0.133 a.u. appears to be a typo (§7.2); no uncertainties on individual g values |
| Petrov, Skripnikov & Titov, arXiv:2302.02856 (2023) | Tier 1 (brief's "modern form of the PT-odd operators") | ¹⁸⁰Hf¹⁹F⁺ | same numerical framework, validated against Roussy 2023 | Eqs. 5–8 basis and Ĥ_mol; Eqs. 7–8 rotating fields; Fig. 1 the upper/lower ↔ n̂ convention; Table I D∥ = −1.53(2) a.u. (c.m.) | HfF⁺ only; the "revisited PT-odd" content is mostly channel-level validation, not new operators |
| Petrov, Skripnikov & Titov, PRA 96, 022508 / arXiv:1704.06631 (2017) | The reference that arXiv:2503.02840 cites for Ĥ_ext | HfF⁺ ³Δ₁ | ab initio + diagonalization | Eqs. 1–8: the operator-level Ĥ_rot, Ĥ_hfs, Ĥ_ext, the g-factor formula, g_N(¹⁹F) = 5.25773 | HfF⁺ parameters |
| Skripnikov & Titov, PRA 91, 042504 (2015) (arXiv:1503.01001) | Tier 1 theory; E_eff, W_{T,P}, G∥, d for ThF⁺ | ThF⁺ ³Δ₁ | 2c-CCSD(T), 38 correlated electrons | Table II: E_eff = 37.3 GV/cm, W_{T,P} = 50 kHz, W_M = 0.88 × 10³³, A∥(Th) = −4163 µ_Th/µ_N MHz, G∥ = 0.034, d = 2.74 D (Th nucleus); Table I: R_e, ω_e | Dipole origin is the Th nucleus, not c.m.; A∥ is the **Th** constant, not ¹⁹F; 7% stated uncertainty |
| Denis et al., NJP 17, 043005 (2015) | Tier 1 theory; independent E_eff | ThF⁺, low-lying states | MRCC / GASCI, DIRAC | E_eff = 35.2 GV/cm, W_{P,T} = 48.4 kHz, A∥(²²⁹Th, I = 5/2) = 1833 MHz, d = 4.03 D (c.m.); ³Δ₁ below ¹Σ⁺ by >300 cm⁻¹ | d is 16% above the later measurement; no ¹⁹F hyperfine |
| Cairncross et al., PRL 119, 153001 (2017) + Suppl. | Cited by Ng 2022 as the Hamiltonian source | HfF⁺ | measurement | Eq. S1 eight-channel parity decomposition; systematics treatment | **Does not contain the effective Hamiltonian** (§3.5) |
| Cairncross PhD thesis (JILA) | Possible fuller HfF⁺ Hamiltonian | HfF⁺ | — | Ch. 2 "Theoretical calculations in HfF⁺," pp. 36–50; App. "Molecular Data," p. 258 | Not read in detail — HfF⁺, superseded for this scope by the Ng thesis |
| Roussy et al., Science 381, 46 (2023) (arXiv:2212.11841) | Tier 1 list; the current eEDM bound | HfF⁺ | measurement | |d_e| < 4.1 × 10⁻³⁰ e·cm | Skimmed; HfF⁺; supplement not separately obtained |
| Zhou et al., PRL 124, 053201 (2020) | Tier 1 list; coherence | HfF⁺ | measurement | second-scale coherence, QPN limit | Skimmed; HfF⁺ |

**Screened and excluded, one line each:**
- Barker et al., JCP 136, 104305 (2012) (Zotero `KMIGLPGB` is the 2011 HfF⁺ Communication, not the 2012 ThF paper): superseded on the state ordering by Gresh 2016 and Denis 2015; its B_e/ω_e values are reproduced in Gresh Table 3.
- Meyer & Bohn, PRA 78, 010502 (2008) (Zotero `CDEXSY3K`): its E_eff = 90 GV/cm for ThF⁺ is repudiated by both modern calculations (§7.5); of historical interest only.
- Petrov et al. 2018, "Evaluation of CP violation in HfF⁺" (Zotero `2CWMNAD7`): HfF⁺ MQM/CP analysis, no ThF⁺ parameter.
- Baturo et al. 2021, "Electric-field-dependent g factor…" (Zotero `28HBEWB9`): predecessor of arXiv:2503.02840 for a different species; not opened after the 2025 ThF⁺ paper was obtained.
- Fleig & Nayak 2014 (Zotero `TCMGSTQ9`): **ThO**, not ThF⁺, despite the brief's expectation.
- Kozlov et al. 1987 / 1992 (Zotero `6VF9I64N`, `VNCBEH4A`): PbF/HgF P,T-odd spin-rotational Hamiltonians; foundational but superseded by Skripnikov's operator definitions for our species.
- Zhou et al. 2019 JMS (neutral ThF): the neutral constants are in Ng thesis Table 2.1 (ω_e = 601.00(2) cm⁻¹, B_e = 0.2339(2) cm⁻¹, r_e = 2.026(3) Å for X ²Δ_{3/2}) — not needed for the ion, since every ion constant was measured directly.
- "Rotational splittings in diatomic molecules of interest to searches for new physics" (2025, ResearchGate): not pursued; the Ω-doubling scaling question it would answer is settled directly by Gresh's k (∝ J(J+1)) and Ng's operator (§5.5).
- Grau, Loh, and Stutz theses: not read; all are HfF⁺-era sources superseded by the Ng thesis for this ThF⁺ scope.
