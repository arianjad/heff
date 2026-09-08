# Digest: Jadbabaie 2025 Springer thesis — the effective-Hamiltonian program

Source: Arian Jadbabaie, *Measuring Fundamental Symmetry Violation in Polyatomic
Molecules*, Springer Theses (2025), 344 PDF pages. This dated source review uses
private text extracts; the published repository includes citations and analysis.

**Page convention.** `p.NNN` = PDF page index, 1-based, as `fitz` counts. **Printed page = PDF page − 19** (verified on PDF pp. 20 / 54 / 61 / 288 / 301 → printed 1 / 35 / 42 / 271 / 284). Where useful I give `p.NNN (pr. NNN)`. The PDF's own bookmark/TOC page numbers are PDF indices, not printed numbers.

Everything below is read directly from the PDF unless tagged **[inferred]** or **[from notes file]**. Quotes are ≤15 words. Where the thesis is silent I say so explicitly rather than filling from Brown & Carrington.

**One extraction artifact you must know about**: the text layer interleaves the rendered glyphs with the LaTeX source (alt-text), so a line often reads `H = B ⃗J 2` followed by `H = B \vecarrow{J}^2`. This is *helpful* — the LaTeX is the ground truth for symbols — but it doubles the token cost of every page. Also, `fi` is a ligature (U+FB01) in this PDF, so a byte-oriented `grep` for `hyperfine` returns **zero** while `hyperﬁne` returns 176. Any absence claim about a word containing `fi`/`ﬁ` must be re-run with the ligature.

---

## 1. Section map

| Topic | Chapter / section | PDF pages | Printed pages |
|---|---|---|---|
| **PT-odd moments, form factors, Schiff's theorem** | §1.2.3, §1.3.1 | 26–33 | 7–14 |
| **PT-odd interaction Hamiltonian form** `H_PTV = ξ W_elec (J_odd·n̂)`; eEDM / NSM / NMQM operators | §1.3.3 (Eq. 1.10, 1.11) | 34–36 | 15–17 |
| Nuclear (octupole-deformation) enhancements | §1.3.4 | 36–37 | 17–18 |
| SMEFT / high-energy mapping, `d_e`–`C_S` plane | §1.3.5 | 37–41 | 18–22 |
| **Parity doublets, molecular orientation control, `⟨n̂⟩ ∝ ⟨MK⟩`, the 4-state ±M/±MK picture** | §1.4.2 (Eq. 1.15, Figs. 1.3–1.4) | 43–45 | 24–26 |
| Quantum projection noise scaling | §1.4.3 (Eq. 1.16–1.18) | 46–47 | 27–28 |
| **Angular momentum, coupling, 3j** | §2.1.1 (Eq. 2.1) | 55–57 | 36–38 |
| **Rotations, Euler angles, D-matrices, spherical tensors, Wigner-Eckart** | §2.1.2 (Eq. 2.2–2.6) | 57–59 | 38–40 |
| Atomic states, Russell-Saunders vs jj | §2.1.3 | 59–61 | 40–42 |
| Ligand-field origin of Λ, ±Λ degeneracy | §2.1.4 | 61–62 | 42–43 |
| BO separation, energy scales | §2.1.5 (Eq. 2.8–2.11) | 62–64 | 43–45 |
| **Molecule frame, anomalous commutation, lab↔molecule tensor transformation, symmetric-top wavefunction** | §2.1.6 (Eq. 2.12–2.19) | 64–67 | 45–48 |
| **Vibration: normal modes, ℓ, `\|v₂,ℓ⟩`, q± ladder operators, δ_ℓ phase choice, Coriolis ζ** | §2.1.7 (Eq. 2.20–2.27) | 67–71 | 48–52 |
| **Hund's cases (a)/(b βJ)/(b βS)/(c); Table 2.1 angular momenta; Table 2.2 basis states; vibronic term symbol `^{2S+1}K_P`; case (a)↔(b) transform Eq. 2.30; (b βJ)↔(b βS) Eq. 2.31** | §2.1.8 | 71–76 | 52–57 |
| **Effective-Hamiltonian principle (EFT analogy, fitted parameters)** | §2.2.1 | 76–77 | 57–58 |
| **Effective-Hamiltonian details: H_Rot, H_SO in case (a), the three operator classes, tracing out L⊥ and G⊥, "no effective Stark terms"** | §2.2.2 (Eq. 2.32–2.33, Fig. 2.1) | 77–82 | 58–63 |
| **Parity doubling: symmetrized states Eq. 2.34, phase `p`, mechanism, Ω=3/2 suppression, parity-dependent Zeeman, Λ vs ℓ phase Eq. 2.35–2.36** | §2.2.3 | 82–84 | 63–65 |
| **Renner-Teller: V₁₁/V₂₂, effective H_RT Eq. 2.39, ε from force constants Eq. 2.40, ε₁/ε₂/g_K Eq. 2.41–2.44** | §2.2.4 | 84–87 | 65–68 |
| Cryogenic buffer-gas source, targets, chemistry | Ch. 3 | 90–140 | 71–121 |
| Decay rates / branching ratios / cross sections (spectra machinery) | §3.2.4.1 (Eq. 3.1 ff.) | 107–112 | 88–93 |
| **X̃(010) effective Hamiltonian (Eq. 4.4 operator, Eq. 4.5 spherical tensor); case (b) parity phase; l-doubling sign convention; how the spectrum is predicted** | §4.3.2 | 148–150 | 129–131 |
| Field-free fit, Table 4.1 parameters, q_G ↔ Coriolis Eq. 4.6 | §4.3.3.1 | 150–153 | 131–134 |
| **Stark and Zeeman: effective Zeeman Hamiltonians Eq. 4.7a/b, `H_E = −D⃗_mol·E⃗`, fitted g_S and D_mol, `n̂·Ẑ ∝ M_N ℓ`** | §4.3.3.2 | 154–156 | 135–137 |
| **TDMs, intensity borrowing, case-(a) TDM matrix element Eq. 4.15, basis-change Eq. 4.12–4.14, interference** | §4.3.3.3 | 156–162 | 137–143 |
| **Renner-Teller / K-resonance modeling: parity-symmetrized case (a) basis Eq. 4.16–4.18, H_RT Eq. 4.19, Table 4.2 (4×4 K-resonance matrix), H_Λ Eq. 4.20, phases Eq. 4.21–4.22, parity-splitting estimate Eq. 4.23** | §4.4.2 | 164–171 | 145–152 |
| **Fitting procedure, parameter fixing, correlations, F-test confidence intervals, Table 4.4; deperturbation; 2×2 reduced form Eq. 4.24** | §4.4.5 | 175–182 | 156–163 |
| Bend-to-bend transitions, saturation | §4.4.6 | 182–185 | 163–166 |
| Transition dipole moments (measured) | §4.4.7 | 185–189 | 166–170 |
| 4f states, unassigned bands, reassignment | §4.5 | 189–197 | 170–178 |
| Rabi/Ramsey theory | §5.1 | 204–210 | 185–191 |
| **ACME / JILA schemes; ³Δ₁ in ThO and HfF⁺ (the only HfF⁺ mention)** | §5.2 | 210–213 | 191–194 |
| CPT, dark states, unresolved hyperfine | §5.3–5.4 | 213–218 | 194–199 |
| **Lindblad master-equation simulation in QuTiP; how H and TDMs are assembled** | §5.6.3.2 | 240–244 | 221–225 |
| Two-photon detuned Raman, X̃(010) hyperfine fit | §5.6.5 | 247–251 | 228–232 |
| **Zero g-factor states: full effective Hamiltonian Eq. 5.29a–h (rot + SR + ℓ-doubling + hyperfine + Zeeman + Stark + ODT), 2D (E_Z,B_Z) grid diagonalization, M_F-block diagonalization, adiabatic state ordering** | §5.7.4 | 263–271 | 244–252 |
| Transverse-field Schrieffer-Wolff treatment Eq. 5.30; imperfect field reversal | §5.7.4.1–5.7.4.2 | 268–271 | 249–252 |
| Ramsey with zero-g states; `H = g_S μ_B B_Z M_S − d_e E_eff Σ` Eq. 5.31 | §5.7.5 | 271–274 | 252–255 |
| **App A.1 Basis states (Eq. A.1–A.8): case (a), case (b), (b βS) coupling** | App A.1 | 288–290 | 269–271 |
| **App A.2 Parity / time-reversal / phase conventions (Eq. A.10–A.20) — authoritative** | App A.2 | 290–294 | 271–275 |
| **App A.2.3 Electronic (Λ) parity doubling, derivation of the −1 phase (Eq. A.21–A.26)** | App A.2.3 | 294–296 | 275–277 |
| **App A.2.4 Vibrational (ℓ) parity doubling, derivation of the +1 phase (Eq. A.27–A.30)** | App A.2.4 | 296–297 | 277–278 |
| **App A.3.1 N² vs R² formalisms (origin offset −K²)** | App A.3.1 | 297–298 | 278–279 |
| App A.3.2 Pure precession | App A.3.2 | 298 | 279 |
| **App A.3.3 Hamiltonian transformations (Van Vleck / contact / Schrieffer-Wolff, Eq. A.31–A.33)** | App A.3.3 | 298–299 | 279–280 |
| **App A.4 Recipe for evaluating matrix elements (Eq. A.34–A.47): 6 numbered steps, 9j/6j/spectator/dot-product/same-system-composite/reduced-ME/D-matrix** | App A.4 | 299–305 | 280–286 |
| **App A.5 Sample matrix elements: A.5.1 no hyperfine (γ_G, p_G, q_G, dipole, spin, Λ-doubling, TDM Eq. A.48–A.54); A.5.2 with hyperfine (b_F, c, Zeeman, Stark Eq. A.55–A.57)** | App A.5 | 305–310 | 286–291 |
| Hanle magnetic-field calibration | App B | 311–313 | 292–294 |
| Bending angle from parity doubling / from hyperfine | App C | 314–315 | 295–296 |
| Line list (science-state lines) | App D | 316–317 | 297–298 |
| Polarization moments | App E | 318–319 | 299–300 |
| MQM science chamber designs; **appendix reference list** | App F | 320–336 | 301–317 |

---

## 2. The effective-Hamiltonian program in the thesis's words

**What an effective Hamiltonian is.** The goal is "to reduce the effects of electronic, vibrational, and rotational interactions to a single, finite Hamiltonian, expressed within a basis of good quantum numbers" (p.76–77, pr. 57–58). He draws the explicit analogy to EFT: the fitted parameters "are like the Wilson coefficients in EFTs" and "encode all of the complicated physics that we have 'integrated out'" (p.77, pr. 58). The price is that the model "is only as accurate as the residuals of the fit" (p.77).

**The formal machinery.** Degenerate perturbation theory on `H = H⁽⁰⁾ + V`, equivalent "up to third order with the contact transformation approach, also referred to as Van Vleck transformations" or Schrieffer-Wolff (p.298, pr. 279). Second-order matrix element (Eq. A.33, p.298–299):
`⟨ψ⁽⁰⁾,i| H̃₂ |ψ⁽⁰⁾,j⟩ = Σ_{ψ⁽ᵅ⁾≠ψ⁽⁰⁾} Σ_k ⟨ψ⁽⁰⁾,i|V|ψ⁽ᵅ⁾,k⟩⟨ψ⁽ᵅ⁾,k|V|ψ⁽⁰⁾,j⟩ / (E⁽⁰⁾ − E⁽ᵅ⁾)`.
Crucially for a term-list code: "for `V = Σ_m V_m`, at second order we must consider `⟨V_m⟩⟨V_n⟩` for all m and n" (p.299) — i.e. **cross terms between distinct perturbations**, not just squares.

**What gets traced out.** The interactions he wants to eliminate are those containing `L⊥ = L_{x,y}` — they appear in `H_Rot = B(J⃗−L⃗−S⃗)²` and `H_SO = A(L⃗·S⃗)` (Eq. 2.32–2.33, p.79, pr. 60). Operators are sorted into three classes (p.79–80): (i) those acting within `ψ⁽⁰⁾` (`J⃗²`, `J₊S₋`) — computed with angular-momentum algebra; (ii) those identical for all states in `ψ⁽⁰⁾` (`S⃗²`) — absorbed into the **origin**; (iii) those containing `L±` that leave the manifold — reappear at second order as effective terms. He notes the origin/diagonal split "is somewhat arbitrary" — `L_zS_z` can go either way depending on whether the basis spans multiple Ω (p.79). The same procedure applies to `G⊥` for polyatomic bending (p.81, pr. 62).

**A structural constraint the code should enforce**: any effective term must respect rotation and parity, so in free field it is "written as scalar products of possibly many angular momenta", and Wigner-Eckart requires an operator involving `S` to satisfy `2S ≥ k` (p.78, pr. 59). Hence no spin-spin term for a single unpaired electron (p.78), but spin-spin *is* expected "with bending modes of polyatomics, or spin-spin interactions in triplet systems" (p.82, pr. 63).

**A deep structural claim worth encoding as a test.** "Curiously, there are no effective Stark terms in the effective Hamiltonian" (p.81, pr. 62) — because such a term would couple the P-odd/T-even `n̂` to the P-even/T-odd `L⃗`, which is exactly the PT-violating signal being searched for. This is a symmetry invariant a generic engine can assert.

**Which basis.** "Since we are most interested in interactions internal to the molecule frame, it is easier to use case (a) representations" (p.78, pr. 59). In practice he builds each electronic/vibronic state in its natural case — X̃(010) in case (b), Ã²Π_{1/2} in case (a) — then converts to case (a) with Eq. 2.30 whenever he needs a TDM (p.75, p.159).

**How he diagonalizes.** Separately per state, then combines. See §5 below.

---

## 3. Term catalogue

All equations verified against the PDF text layer; entries marked ✅rendered were additionally checked against a PNG render of the page.

### 3a. Free-field / internal structure

| Term | Symbol(s), units | Eq. no. | PDF p. (pr.) | Basis | Notes |
|---|---|---|---|---|---|
| Rotation, R² formalism, case (a) | `B` (also "bare" `A`) | 2.32 | 79 (60) | case (a) | `H_Rot = B(J⃗−L⃗−S⃗)²`, expanded to `B(J²+L²+S²−2J_z(L_z+S_z)−2L_zS_z −J₊L₋−J₋L₊−J₊S₋−J₋S₊−L₊S₋−L₋S₊)`. Triatomic generalization `R⃗ = J⃗−L⃗−G⃗_ℓ−S⃗` (p.79) ✅rendered |
| Spin-orbit | `A` | 2.33 | 79 (60) | case (a) | `H_SO = A(L_zS_z + ½(L₊S₋+L₋S₊))` ✅rendered |
| Rigid rotor / symmetric top | `B`; `A`,`B`,`C` | 2.15–2.17 | 66 (47) | — | Prolate `H = B(J⃗²−J_z²) + A J_z²`; convention `I_r` (a,b,c ↔ z,x,y); only prolate tops considered |
| Rotation (bending mode) | `B` | 4.4 | 149 (130) | case (b) | `B(N⃗² − ℓ²)` ✅rendered |
| Spin-rotation (bending-corrected) | `γ` (MHz) | 4.4, 4.5 | 149 (130) | case (b) | `γ(N⃗·S⃗ − N_zS_z)`. **The `−γ N_zS_z` subtraction is deliberate** and "crucial for accurate description of low-N spectra" (p.149); for a linear molecule `N_z=0` so `N⃗·S⃗` implicitly only carries `N_xS_x + N_yS_y` ✅rendered |
| Axial spin-rotation | `γ_G` (a.k.a. `γ′`) | 4.4 | 149 (130) | case (b) | `+γ_G N_zS_z`. 1st order from magnetic dipole (negligible for Yb); 3rd order via vibronic × spin-orbit mixing with Π states (p.149, p.152) |
| Parity-dependent spin-rotation (ℓ-type, p-type) | `p_G` | 4.4 | 149 (130) | case (b) | `+ (p_G/2)(N₊S₊e^{−i2φ} + N₋S₋e^{+i2φ})` ✅rendered |
| Rotational ℓ-type doubling | `q_G` | 4.4 | 149 (130) | case (b) | `− (q_G/2)(N₊²e^{−i2φ} + N₋²e^{+i2φ})`. **Note the sign asymmetry: `+p_G/2`, `−q_G/2`** ✅rendered |
| Same, spherical-tensor form | — | 4.5 | 149 (130) | case (b) | `H_X̃ = T₀ + B(N²−ℓ²) + γ(N·S − T¹_{q=0}(N)T¹_{q=0}(S)) + γ_G T¹_{q=0}(N)T¹_{q=0}(S) + Σ_{q=±1} e^{−2iqφ}(p_G T²_{2q}(N,S) − q_G T²_{2q}(N,N))` ✅rendered |
| Electronic Λ-doubling | `p_e`, `q_e` (cm⁻¹) | 4.20 | 169 (150) | case (a) | `H_Λ = ½(p_e+2q_e)(J₊S₊e^{−2iθ} + J₋S₋e^{2iθ}) − (q_e/2)(J₊²e^{−2iθ} + J₋²e^{2iθ})` = `(p_e+2q_e)Σ_{q=±1}e^{−2iqθ}T²_{2q}(J,S) − q_e Σ_{q=±1}e^{−2iqθ}T²_{2q}(J,J)` ✅rendered. Selection rules `Δℓ=0, ΔΛ=±2, ΔΣ=0,∓1, ΔP=±2,±1` |
| Vibrational ℓ-doubling (case a) | `p_G`, `q_G` | (from 4.20) | 169 (150) | case (a) | Obtained from Eq. 4.20 by `p_e→p_G, q_e→q_G, θ→φ`. **The thesis literally prints "θ → θ" here — a typo; the bending azimuthal angle is φ.** Selection rules `Δℓ=±2, ΔΛ=0, ΔΣ=0,∓1, ΔP=±2,±1` ✅rendered |
| Renner-Teller (effective) | `εω₂`, `g_K`, `εω_{2,D}` (cm⁻¹; ε unitless) | 2.39 = 4.19 | 85 (66), 166 (147) | case (a) vibronic | `H_RT = ½εω₂(q₊²e^{−2iθ} + q₋²e^{2iθ}) + g_K(G_z+L_z)L_z + ½εω_{2,D}(q₊²e^{−2iθ}+q₋²e^{2iθ})N⃗²`, with `G_z=ℓ`, `L_z=Λ`. Selection rules `ΔΛ = −Δℓ = ±2, Δv₂ = 0,±2, ΔK=0, ΔΣ=0, ΔP=0` (p.167). The centrifugal `εω_{2,D}` term is "included for completeness" and not used |
| RT microscopic origin | `V₁₁` (dipolar), `V₂₂` (quadrupolar) | 2.37–2.38, 4.8 | 84 (65), 158 (140) | — | `H_RT = V₁₁ q₂cos(θ−φ) + V₂₂ q₂²cos²(θ−φ)+…`; exponential form shows `L±q∓` (dipolar) and `L±²q∓²` (quadrupolar). Eq. 4.8 gives the Herzberg-Teller form `H_RT = (V₁₁/2)(L₊q₋e^{i(θ−φ)} + L₋q₊e^{−i(θ−φ)})` |
| RT parameter from force constants | `ε` | 2.40 | 86 (67) | — | `ε = (k′−k″)/(k′+k″)`; `ε<0` in the M-OH molecules ⇒ the reflection-symmetric electronic state is lower |
| RT effective-parameter formulas | `ε₁`, `ε₂`, `g_K` | 2.41–2.44 | 86 (67) | — | `εω₂ = (ε₁+ε₂)ω₂`; `ε₂ω₂ = ⟨η\|V₂₂\|η⟩`; `ε₁ω₂ = −Σ_{η⊂Σ}(−1)^s \|⟨η\|V₁₁\|η′⟩\|²/(2ΔE) (1 + (ω₂/ΔE)²)`; `g_K = (ω₂/4)Σ_{η⊂Σ,Δ}(−1)^p \|⟨η\|V₁₁\|η′⟩\|²/(ΔE)²`. `s=0` for Σ⁺, `s=1` for Σ⁻; `p=0` for Σ, `p=1` for Δ. **Gauyacq & Jungen (Ref. [28]) disagree**, writing `ε = (ε₁+ε₂)(1+ε₁)⁻¹` (p.86) |
| K-resonance matrix (the deliverable table) | `T₀, ω, A, B, γ, g_K, ε, p_e, q_e, p_G, q_G` | **Table 4.2** | 168 (149) | case (a), parity-symmetrized | 4×4 over `{²Δ_{5/2}, ²Δ_{3/2}, κ²Σ_{1/2}, μ²Σ_{−1/2}}`; **2×2 for J=1/2, 3×3 for J=3/2, 4×4 for J≥5/2**; symmetric; block diagonal in J without hyperfine. Shorthands `z = (J+½)²−1 = J(J+1)−3/4`, `B* = B − γ/2`, `ω = ω₂`, `ε_corr = ε(1 + ((εω)²/4)·(8ω²−6A²)/(4ω²−A²)²)`. Upper/lower signs = overall parity; parity phase `p = J−S−ℓ` ✅rendered |
| Reduced 2×2 K-resonance form | — | 4.24 | 181 (162) | rotated basis `\|²Δ_{3/2}⟩ ± (−1)^p\|μ²Σ_{1/2}⟩` | `H_eff = E₀ + BJ(J+1) + [[∓(p_e+2q_e)/2 · z^{1/2}, g_K−B−(εω)²/(2(2ω+A))],[…, ±(p_e+2q_e)/2 · z^{1/2}]]`, with `E₀ = T₀ + ω − A/2 + g_K + B/4 − γ − (εω)²/(2ω+A)` |
| Vibrational energy (diatomic) | `ω_e, ω_ex_e, ω_ey_e` | 2.20 | 68 (49) | — | standard Dunham-style expansion |
| `B_v` from vibration | `B_e, α_i` | 2.21, 2.27 | 68 (49), 71 (52) | — | `B_v = B_e − Σ_i α_i(v_i + d_i/2) + …`, `d_s=1`, `d_t=2` |
| Polyatomic vibrational energy incl. `g_{kk′}ℓ_kℓ_{k′}` | `ω_i, g_{kk′}` | 2.26 | 71 (52) | — | `g_{kk′}ℓ_kℓ_{k′}` generates the anharmonic ℓ-splitting |
| `q_G` ↔ Coriolis constants | `ζ_{2n}`, `ω_n` | 4.6 | 152 (133) | — | `q_G = −(v₂+1)(B²/ω₂)(1 + Σ_{n=1,3} ζ²_{2n}·4ω₂²/(ω_n²−ω₂²))`; with `Σ_i ζ²_{ij}=1` for linear molecules (p.71) |
| Curl relation | — | text | 169 (150) | — | `q ~ p B/A` in the unique-perturber limit; used to justify dropping `q_e`, `q_G` |

### 3b. Hyperfine

| Term | Symbol | Eq. no. | PDF p. (pr.) | Basis | Notes |
|---|---|---|---|---|---|
| Fermi contact | `b_F` (MHz) | 5.29e; A.55 | 266 (247); 308 (289) | (b βJ) | `b_F I⃗·S⃗`. Full 6j-factorized matrix element given in A.55 |
| Dipolar (axial) | `c` (MHz) | 5.29e; A.56 | 266 (247); 309 (290) | (b βJ) | `(c/3)(3I_zS_z − I⃗·S⃗)` in Eq. 5.29e; the matrix element A.56 is written for `T²_{q=0}(I,S)` with a `−√(5/3)` prefactor and a 9j |
| (b βJ) ↔ (b βS) transform | — | 2.31 | 76 (57) | — | 6j-based; used for odd-Yb isotopologues where `G_Yb = S⃗+I⃗_Yb` couples first |
| Multi-spin coupling order | — | text | 76 (57) | — | ¹⁷¹/¹⁷³YbOH: `G_Yb = S+I_Yb`, then `F₁ = N+G_Yb`, then `F = F₁+I_H` |

**Absences to record.** The thesis gives **no nuclear electric quadrupole (`eQq`) operator or matrix element** and **no nuclear spin-rotation term**. Confirmed by `grep -ciE "eQq" _full.txt` → 5, all of which are false positives from the LaTeX macro `\coloneqq` (contains "eqq"); and `grep -ciE "nuclear spin.rotation"` → 0 (no `fi` ligature in that phrase, so the count is trustworthy). The nuclear electric quadrupole appears only *physically*, as the reason ¹⁷³YbOH couples in case (b βS): "the strong electric quadrupole interaction between the Yb-centered electron … and the non-spherical Yb nucleus" (p.132, pr. 113). It also gives no `a`, `d` (orbital / dipolar-perpendicular) hyperfine parameters, and no spin-spin term (explicitly deferred to triplets, p.82).

### 3c. Field terms

| Term | Symbol | Eq. no. | PDF p. (pr.) | Notes |
|---|---|---|---|---|
| Stark, ground and excited | `D_mol` (D) | text; 5.29g | 155 (136); 266 (247) | `H_E = −D⃗_mol·E⃗`; in Eq. 5.29g `H_Stark = −D_Z E_Z`. Measured `D_mol = 2.16(1) D = 1.09 h MHz/(V/cm)` for YbOH X̃(010) |
| **Dipole sign convention** | — | footnotes | 34 (15) fn.12; 156 (137) fn.4 | "We take `n̂` to point along `D⃗_mol`, that is pointing from `− → +`. In M-OH molecules, `D⃗_mol` points from O to M in the physics convention." And: "The molecule `ẑ` axis and dipole moment `D_mol` both point from O to Yb." |
| Zeeman, ²Σ ground state | `g_S` | 4.7a | 154 (135) | `H^Zee_X = g_S μ_B S_Z B_Z` — **isotropic electron spin only**. Fitted `g_S = 2.07(2)` for X̃(010) |
| Zeeman, ²Π_{1/2} excited state | `g_S′, g_L, g_l′` | 4.7b | 154 (135) | `H^Zee_A = g_S′ μ_B S_Z B_Z + g_L L_Z B_Z + g_l′ μ_B (e^{−2iθ}S₊B₊ + e^{2iθ}S₋B₋)`. Fixed to `g_S′=1.860`, `g_L=1.0`, `g_l′=−0.724`. **The third term is the anisotropic / parity-dependent piece** |
| Zeeman, zero-g modeling | `g_S` | 5.29f | 266 (247) | `g_S μ_B B_Z S_Z`; **nuclear and rotational Zeeman explicitly ignored**, estimated to contribute at `≤10⁻³ μ_B` (p.267, pr. 248) |
| Transverse Zeeman | — | text; 5.30 | 268–270 | Adds `B_X S_X + B_Y S_Y` (`ΔM_F = ±1`); second-order effect via Schrieffer-Wolff, Eq. 5.30 |
| AC Stark / optical dipole trap | `α(ω)` | 5.29h | 266 (247) | `H_ODT = −D⃗·E⃗_ODT`, `E⃗_ODT = E₀/2(ε̂_ODT e^{−iωt} + c.c.)`; 1064 nm dynamic polarizabilities taken from the literature. Scalar/vector/tensor decomposition referenced but not re-derived |
| Anisotropic-g estimate | — | text | 154 (135) | Curl-type relations estimate anisotropic spin interactions at `6×10⁻³ μ_B`, below the resolution of that experiment |

### 3d. PT-odd

| Term | Symbol | Eq. no. | PDF p. (pr.) | Notes |
|---|---|---|---|---|
| Generic PT-odd interaction | `ξ_PTV`, `W_elec` | **1.10** | 35 (16) | `H_PTV = ξ_PTV W_elec (J⃗_odd · n̂)`. `W_elec` "must be calculated from electronic structure theory" and its accuracy gauged against hyperfine parameters |
| eEDM | `d_e`, `E_eff` | 1.10; 5.31 | 35 (16); 271 (252) | `J⃗_odd = S⃗`. In the Ramsey Hamiltonian: `H = g_S μ_B B_Z S⃗·Ẑ − d_e E_eff S⃗·n̂ = g_S μ_B B_Z M_S − d_e E_eff Σ` |
| NSM | — | text | 35 (16) | `J⃗_odd = I⃗` |
| NMQM | `W_M`, `M` | **1.11** | 36 (17) | `H_MQM = W_M M J⃗_odd·n̂ = −(M/(2I(2I−1))) S⃗·T⃡·n̂ = (M/(2I(2I−1)))√(20/3) T¹(S,T²(I,I))·n̂`, with `T_ij = {I_i,I_j} − (2/3)δ_ij I(I+1) = 2T²(I,I)` |
| eEDM sensitivity observable | `Σ = S⃗·n̂` | 5.31 | 271 (252) | Computed as an expectation value of eigenvectors, not as a Hamiltonian term |
| Effective g-factor observable | `g_eff` | Fig. 5.28 caption; 271 (252) | 269, 272 | `g_eff = g_S μ_B (⟨M_S⟩_{M=+1} − ⟨M_S⟩_{M=−1})`, i.e. a **differential** between the ±M pair, computed numerically |
| Numeric `W_M` | — | text | 286 (267) | `W_M ≈ −1.07×10³³ Hz/(e cm²)/c` for ¹⁷³YbOH; `ΔE = h W_M M × 0.2` |

**Absence.** The symbols `W_d` and `W_s` do not appear; enumerating every `W_<sub>` token in the whole PDF (`grep -oE "W_[A-Za-z{]…"`) returns only `W_elec`, `W_M`, `W_{eg}`, `W_{ge}`, `W_{ji}` (the last three are optical transition frequencies, not PT-odd constants). There is no `W_P` / NSD-PV κ-parameter treatment: `grep -ci "NSD"` → 0. Parity violation is mentioned as motivation (8 hits, e.g. ¹⁷¹YbOH, p.162) but no P-odd/T-even Hamiltonian term is written down.

---

## 4. Conventions the code must preserve

1. **Rotation convention (active, Brown & Carrington).** `J_X, J_Y, J_Z` generate *active* rotations; `R_n̂(θ) = e^{−iθJ_n̂}`. Euler angles `ω = (φ,θ,χ)` with `R(ω) = R_Z(φ)R_Y(θ)R_Z(χ)`, applied right to left; `R⁻¹(φ,θ,χ) = R(−χ,−θ,−φ)`. **Eq. 2.2, p.57–58 (pr. 38–39).** The thesis flags that this is "opposite of the convention followed by Ref. [19]" — Budker/Kimball/DeMille (footnote 2, p.57).
2. **Wigner-Eckart.** `⟨η,J,M|T^k_p(A)|η′,J′,M′⟩ = (−1)^{J−M} (J k J′; −M p M′) ⟨η,J||T^k(A)||η′,J′⟩`. **Eq. 2.6, p.59 (pr. 40)** — boxed in the thesis "as we will use this Eq. 2.6 over and over again."
3. **Condon-Shortley phase everywhere.** "`J_X` and `J_x` are real and positive, known as the Condon and Shortley phase, which takes `δ_P = 0`. We use this phase choice everywhere in this thesis." **p.292 (pr. 275)** ✅rendered.
4. **Symmetric-top wavefunction.** `Ψ_JKM = ⟨ω|JKM⟩ = √((2J+1)/8π²) D^{(J)}_{M,K}(ω)*`, normalized so `∫dω Ψ*Ψ = 1`. **Eq. 2.18, p.67 (pr. 48).**
5. **Anomalous commutation.** `[J_a,J_b] = −iJ_c ε_abc` in the molecule frame. **Applies to `J` and `N` and anything containing them; `L`, `S`, `I`, `G_ℓ` are NOT anomalous** (p.65, pr. 46; Table 2.1 p.73; and restated at p.289 and p.299). Remedy: always rotate anomalous operators to the lab frame before evaluating, using Eq. 2.13–2.14 / A.34.
6. **Lab↔molecule tensor transformation.** `T^k_p(A) = Σ_q D^{(k)}_{p,q}(ω)* T^k_q(A)` and `T^k_q(A) = Σ_p (−1)^{p−q} D^{(k)}_{−p,−q}(ω)* T^k_p(A) = (−1)^q D^{(k)}_{.,−q}(ω)* · T^k(A)`. **Eq. 2.13–2.14 (p.65), restated Eq. A.34 (p.300, pr. 281).** The `D^{(k)}_{.,−q}` notation means "reduced in the lab frame", the dot product summed over lab components `p`.
7. **Parity operator.** `E* = σ_xz R_y(π)` — the Brown & Carrington / Hirota choice, "This is the convention we use in this thesis". The alternative `E* = σ_yz R_x(π)` (Bunker & Jensen, Zare) is given "for completeness" and **not used**. **p.291 (pr. 274)** ✅rendered. Euler-angle action: `(φ,θ,ξ) → (π+φ, π−θ, π−ξ)`; molecule-frame function action: `f(x,y,z) → f(x,−y,z)`.
8. **Component parity phases.** `R_y(π)|J,P,M⟩ = (−1)^{J−P}e^{−2iPδ_P}|J,−P,M⟩` (Eq. A.11); `σ_xz|v,ℓ⟩ = e^{−2iℓδ_ℓ}|v,−ℓ⟩` (Eq. A.12); `σ_xz|S,Σ⟩ = (−1)^{S−Σ}|S,−Σ⟩` (Eq. A.13); `σ_xz|Λ⟩ = (−1)^{s+Λ}|−Λ⟩` with `s=1` for Σ⁻ and `s=0` otherwise (Eq. A.14). **pp. 291–292 (pr. 274–275)** ✅rendered. *(Note: the printed Eq. A.14 second line reads `σ_yz|S,Σ⟩ = (−1)^s|S,−Λ⟩` — a mismatched ket, a typesetting error in the published thesis.)*
9. **Composite parity phase — the one the code must hard-code.** `E* Ψ_{Λ,ℓ,Σ,P} = (−1)^{J−P}(−1)^{S−Σ}(−1)^{Λ+s} e^{−2iℓδ_ℓ} Ψ_{−Λ,−ℓ,−Σ,−P} = (−1)^{J−S−ℓ+s} Ψ_{−Λ,−ℓ,−Σ,−P}`. **Eq. A.15, p.292 (pr. 275)** ✅rendered. With `δ_ℓ = 0`, `P = Σ+ℓ+Λ`, and `S = |Σ|`. So **case (a): `p = J − S − ℓ` (+`s`); case (b): `p = N − ℓ`** (p.82, pr. 63; p.148, pr. 129; p.165, pr. 146). Consequence he flags explicitly: "this choice means the definition of parity changes upon exciting odd number of ℓ quanta" (p.293, pr. 276) — `Δℓ=±1` interactions mix symmetric and anti-symmetric parity superpositions. Ref. [24] (in App A) proposes `δ_ℓ = π/2` to remove this; **he does not use it**, because it would modify the `q±` matrix elements.
10. **Vibrational phase `δ_ℓ = 0`** (following Brown and Hirota), so `Ψ_{v₂,|ℓ|} = Ψ_{v₂,−|ℓ|}` and `E*|v,ℓ⟩ = |v,−ℓ⟩`. `δ_v = 0` follows from requiring real positive `q±` matrix elements. **p.292–294 (pr. 275–277).**
11. **The two azimuthal-operator phases — the single most consequential sign in the code.**
    `⟨Λ = ±1| e^{±2iφ_e} |Λ′ = ∓1⟩ = −1` (Eq. A.21, p.294 / Eq. 2.35, p.83 / Eq. 4.21, p.170)
    `⟨ℓ = ±1| e^{±2iφ_n} |ℓ′ = ∓1⟩ = +1` (Eq. A.30, p.297 / Eq. 2.36, p.83 / Eq. 4.22, p.170)
    **They differ in sign.** Derived from scratch in App A.2.3 and A.2.4 (the `Θ_{L,−|Λ|} = (−1)^Λ Θ_{L,|Λ|}` property of associated Legendre functions vs. the `ℓ`-independence of the radial `Ψ_{v,ℓ}(q)`). The Λ convention is Mulliken & Christy, reiterated by Brown and by Brown & Carrington (p.294, pr. 275).
12. **Λ-doubling sign consequence.** In this convention "a positive `q_e` electronic Λ-doubling parameter in a ¹Π state corresponds to the `(−1)^J` parity level lying above the `(−1)^{J+1}` parity level" (p.294, pr. 275). When written in case (a), **the `J±S±` terms carry a positive prefactor and the `J±J±` terms a negative prefactor** — matching Eq. 4.20. In the YbOH Ã state `p_e + 2q_e` is negative, so the `−` parity state lies below `+`.
13. **Time reversal.** `T|J,M⟩ = (−1)^M|J,−M⟩`, `T²|J,M⟩ = (−1)^{2M}|J,M⟩` (Eq. A.16–A.17); `T|J,P,M⟩ = (−1)^{M−P}e^{2iJη_J}|J,−P,−M⟩` with **`η_J = 0` always**, following B&C (Eq. A.18); `T|v,ℓ⟩ = e^{2ivδ_v}|v,−ℓ⟩` with `δ_v = 0` (Eq. A.19); `T T^k_p(A) T = (−1)^p T^k_{−p}(A)` (Eq. A.20). The half-integer phase choice is `η_T = i^{2J}` (so the integer form is recovered); the alternative `η_T = 1` gives `(−1)^{J−M}` and differs only by `(−1)^J`. **pp. 293–294 (pr. 276–277).**
14. **Coupling order matters and must be tracked.** "As long as we are always consistent with which angular momenta is `J₁` and which is `J₂`, the order of coupling does not matter for the end result" — but the *spectator-theorem phase factors differ* between the `A₁`-on-`J₁` and `A₂`-on-`J₂` forms (Eq. A.39 vs A.40), and "This distinction is important when combining multiple terms that act on different components", e.g. Zeeman in hyperfine-coupled states. **p.302–303 (pr. 283–284).** Also stated at p.308 (pr. 289): "we must be consistent with the order of coupling `N` and `S` to form `J`."
15. **Case (a) ↔ case (b) transformation.** `|N,K;(N,S)J,M⟩ = Σ_{Σ,P} (−1)^{N−S+P} √(2N+1) (J S N; P −Σ −K) |S,Σ⟩|J,P,M⟩`, `P = Λ+Σ+ℓ`, `K = Λ+ℓ`. **Eq. 2.30 (p.75, pr. 56), restated Eq. A.7 (p.289) and Eq. 4.11 (p.159).** The thesis says explicitly: **"We note the formula in Ref. [13] has a typo"** — i.e. the Brown & Carrington version is wrong; his is from Brown's original paper and matches Hirota. He uses it "extensively to calculate transition dipole moments … by expressing all states in Hund's case (a)."
16. **R² formalism, not N².** He uses R² (p.298, pr. 279). The two differ in electronic origin by `−(L_z + G_{ℓz})² = −K²`, and consequently in their centrifugal-distortion parameters. `R⃗ = N_x x̂ + N_y ŷ + (N_z − K)ẑ` connects them. Reason given: the YbOH origin-band paper used R², "with matrix elements taken from the Appendix of Ref. [7]" — which, per the appendix reference list (p.335), is **Brown, Kopp, Malmberg & Rydh, Phys. Scripta 17, 55 (1978), on AlF hyperfine** — and he flags **a typo in that paper's `H_sr` matrix element: "Σ² − S(S+1)" should be "ΩΣ − S(S+1)" for a diatomic** (footnote 3, p.298).
17. **Reduced matrix element of the D-matrix.** `⟨J,P||D^{(k)}_{.,q}||J′,P′⟩ = (−1)^{J−P}√((2J+1)(2J′+1)) (J k J′; −P q P′)`, "reduced in the lab frame" (B&C eq. 5.186). Case (b) form by `J→N, P→K`. **Eq. A.47, p.305 (pr. 286).**
18. **Reduced matrix elements of angular momenta.** `⟨J||T¹(J)||J′⟩ = δ_{J,J′}√(J(J+1)(2J+1))` (Eq. A.45); `⟨J||T²(J)||J′⟩ = δ_{J,J′}(2J−1)J/(√6 (J 2 J; J 0 J))` (Eq. A.46); `⟨j||1||j′⟩ = δ_{j,j′}√(2j+1)` (footnote 4, p.302). **pp. 304–305.**
19. **Composite-tensor coupling identity used for MQM and hyperfine.** `T²(C)·T²(A,B) = −√(5/3) T¹(A)·T¹(C²,B)` (B&C eq. 8.459). **Eq. A.43, p.304 (pr. 285).**
20. **Three further "the textbook is wrong" flags the code should heed.** (a) Hirota tab. 2.4 item 2 "has a typo in the 3j symbol's lower row" (p.300, pr. 281). (b) "We caution the reader from using Brown and Carrington, eq. 5.177, which seems to be missing the extra factors from the spectator theorem" (p.303, pr. 284). (c) The Cu Cl₂ / Brown references for the K-resonance matrix both have "minor typos in the expression for `ε_corr`", corrected in Table 4.2 (p.169, pr. 150). Also: in Ref. [16] `±` means `e/f` parity, while in Ref. [20] `±` means overall parity (p.167, pr. 148) — a live trap when importing published matrices.
21. **`e/f` parity definition.** `e` levels have `P = (−1)^{J−S−l}`, `f` levels have `P = −(−1)^{J−S−l}` (Fig. 4.14 caption, p.181, pr. 162).

**PGopher: not mentioned anywhere.** `grep -ciE "pgopher" _full.txt` → **0** (no `fi` ligature in the word, so the zero is trustworthy; control query `Wigner-Eckart` on the same file returns 25). So the thesis makes **no statement** about where his conventions differ from PGopher. Differences from Brown & Carrington *are* stated, and are listed in items 15, 16, 20 above.

---

## 5. Computational procedure

**Software the thesis names.** Only three, none of them an effective-Hamiltonian package:
- **QuTiP** for the Lindblad master equation (p.240, pr. 221; footnote 13 gives `https://qutip.org/`).
- **AtomicDensityMatrix** package in **Mathematica**, for angular-momentum probability surfaces (Fig. 5.3 caption, p.211, pr. 192).
- **LightTools** for fluorescence-collection ray tracing (p.323, pr. 304) — not structure.

**The thesis never names the code used to build and diagonalize the effective Hamiltonians.** There is no repository link, no "our code", no software section. Searched `our code|custom code|software|open.source|repositor|jupyter|notebook` over the whole text: the only hits are the Springer copyright boilerplate, the LightTools mention, and an Arduino mention in the CV.

**How levels are computed (field-free, spectra).** §4.3.2, p.150 (pr. 131):
1. Build and **diagonalize the ground and excited effective Hamiltonians *separately***.
2. Truncate: `N″ = 6` for X̃(010), `J′ = 15/2` for Ã; **include the `P = 3/2` manifold when diagonalizing Ã** (i.e. do not truncate to the Ω=1/2 block).
3. Convert all eigenvectors to **Hund's case (a)** (Eq. 2.30 / 4.11).
4. Compute TDM matrix elements (Eq. 4.15 / A.54).
5. For transitions with non-zero TDM, line position = difference of excited and ground eigenvalues.

**How field maps are computed (Stark/Zeeman spectroscopy).** §4.3.3.2, p.154 (pr. 135): "we fix the field-free parameters and diagonalize the **combined Stark, Zeeman, and field-free Hamiltonian**" — i.e. one full diagonalization of the summed Hamiltonian, not perturbation theory. Free field parameters (`g_S`, `D_mol`) then obtained by least-squares on observed-minus-predicted line positions.

**How 2D field maps and state tracking are done (the closest thing to a scan engine).** §5.7.4, p.267 (pr. 248):
- "By diagonalizing `H_eff` over a **grid of `(E_Z, B_Z)` values**, we can obtain 2D plots of g-factors and eEDM sensitivities."
- "Using the Z-symmetry of the Hamiltonian, we **separately diagonalize each `M_F` block to avoid degeneracies at `B_Z = 0`**." — a concrete, reusable trick.
- "**Continuous 2D surfaces for eigenvalues and eigenvectors are obtained by ordering eigenstates at each value of `(E,B)` according to their adiabatically correlated free field state.**" — this is his state-tracking/labelling rule, stated once and only here. It is *adiabatic correlation to the zero-field state*, not eigenvector-overlap-to-previous-step or Hungarian matching. He also labels the finite-field states "in terms of their adiabatically correlated zero-field quantum numbers `|N,J′,F,M⟩`" (p.267).
- Observables extracted from eigenvectors: `⟨M_S⟩`, `⟨Σ⟩`, `⟨M_Nℓ⟩`, and the **differential** `g_eff = g_S μ_B(⟨M_S⟩_{M=+1} − ⟨M_S⟩_{M=−1})` and `⟨Σ⟩_{M=+1} − ⟨Σ⟩_{M=−1}`.

**The full effective Hamiltonian assembled for that calculation (Eq. 5.29a–h, p.266, pr. 247)** — this is the cleanest single "term list" in the thesis and the best template for a generic engine:
```
H_eff  = H_Rot + H_SR + H_ℓ + H_Hyp + H_Zeeman + H_Stark + H_ODT
H_Rot  = B(N⃗² − ℓ²)
H_SR   = γ(N⃗·S⃗ − N_z S_z)
H_ℓ    = −q_ℓ (N₊² e^{−i2φ} + N₋² e^{+i2φ})
H_Hyp  = b_F I⃗·S⃗ + (c/3)(3 I_z S_z − I⃗·S⃗)
H_Zeeman = g_S μ_B B_Z S_Z
H_Stark  = −D_Z E_Z
H_ODT    = −D⃗·E⃗_ODT
```
with `z`,`±` molecule-frame and `Z` lab-frame, φ the nuclear bending coordinate. Note `p_G` is *absent* from this list even though it was measured in Ch. 4 — an inconsistency worth asking about (Open Question 5).

**How a fit is run.** §4.4.5, p.177 (pr. 158):
1. Reduce the parameter count using physical arguments (fix `p_e` to the (000) value; set `q_e = q_G = 0` justified by the Curl relation `q_e ≈ p_e B/A ≈ 2×10⁻⁴ p_e`; fix `A = 1350 cm⁻¹`, `T₀`, and `γ = 0` because "the effects of γ are largely indistinguishable from B and origin offsets").
2. Seed initial values from a homologous molecule (CaOH Ã(010) for `ε`, `g_K`) or a sibling state (YbOH Ã(000) for `B`, `γ`).
3. **Nelder-Mead** to get a first parameter set, then use those as initial values for **Levenberg-Marquardt** non-linear least squares.
4. Report residual standard deviation (0.0012 cm⁻¹ ≈ 36 MHz) against experimental error.

**How spectra/TDMs are computed.** Case (a) E1 TDM, Eq. 4.15 = Eq. A.54 (p.161 / p.308): `δ_{Σ,Σ′}δ_{ℓ,ℓ′} × (−1)^{J−M}(J 1 J′; −M p M′) × √((2J+1)(2J′+1))(−1)^{J−M} × Σ_q (J 1 J′; −P q P′) δ_{Λ,Λ′+q} × ⟨Λ||T¹_q(d)||Λ′⟩`. Selection-rule mapping: `q = ±1` drives `ΔΛ = ±1` (perpendicular bands), `q = 0` drives `ΔΛ = 0` (parallel bands) (p.157, pr. 138). **Amplitudes from different intensity-borrowing admixtures are summed and then squared** — "the TDM is squared after the sum, allowing TDMs from different states to interfere" (p.161, pr. 142). Intensity model uses a fixed rotational temperature (2 K).

**Dynamics.** Lindblad/OBE in QuTiP: `dρ/dt = −(i/ħ)[H,ρ] + Σ_ij γ_ij D[c_ij]ρ`, `D[c]ρ = cρc† − ½{c†c,ρ}`, `c_ij = |i⟩⟨j|`, `γ_ij` the partial width from Eq. 3.1 (p.240–241, pr. 221–222). **The Hamiltonian for the dynamics is assembled by putting "the eigenvalues of direct diagonalization on the diagonals, and the computed TDMs for the off diagonals"** (p.241) — i.e. the structure code feeds the dynamics code. Rotating-frame transformation `R = e^{iξt}` with diagonal `ξ` (Eq. 5.22 ff., p.243).

---

## 6. Parameters and uncertainty handling

- **Provenance.** Parameters come from fits to his own spectra, from published sibling states, or from homologous molecules. Explicit examples: X̃(010) `B, γ, γ_G, q_G, p_G` fitted here while Ã-state parameters were "held fixed … to the values given in Ref. [94]" (p.151, pr. 132). Ã(010) fit fixes `A`, `T₀`, `p_e`, `q_e`, `q_G`, `γ` and floats `ω₂, ε, B, g_K, p_G` (p.177, pr. 158). Zero-g modeling takes field-free parameters from the literature but re-fits `b_F, c` to observed lines (p.266, pr. 247).
- **Uncertainty treatment — the part that matters for the new code's test design.** For the Ã(010) fit he reports "strong correlations (>90%) between `ω₂`, `ε`, and `g_K`, moderate correlation (~60%) of those parameters with `p_G`, and little correlation (~15%) with `B`" (p.177, pr. 158). Then, explicitly: **"While the strong correlations could be removed by fixing `ω₂`, we decided not to, as this results in parameter uncertainties that were unreasonably small."** Instead he plots **2-D confidence intervals obtained from F-tests** (Fig. 4.12, p.178) and validates the fit by checking that "all of the parameter intervals correspond to physically reasonable values."
- **Propagation.** `εω₂` errors "obtained by multiplication and propagation of errors"; parentheses in Table 4.4 are **2σ errors (68% confidence interval)** as the table footnote states (p.179, pr. 160) — note the internal inconsistency between "2σ" and "68%".
- **Running the same model with different parameter sets.** He does exactly this, deliberately: **"Finally, we can set `H_K = 0` in the model, essentially turning off the K resonance"** to produce a deperturbation diagram (p.179, pr. 160; Fig. 4.14). He also does the `Λ → ℓ` substitution to approximate X̃(010)/Ã(010) as a ²Π_{1/2} state for the master-equation model, flagging that it "is not strictly accurate" (p.241, pr. 222). Both are term-toggling / parameter-substitution studies, not snapshot comparisons — good support for a design where terms are switchable and parameter sets are swappable.
- **Cross-checks against independent estimates rather than pinned numbers.** The emergent spin-rotation `γ ~ 0.3 cm⁻¹` from the K-resonance is compared to `γ ~ 0.38 cm⁻¹` estimated from B&C eq. 7.122 with the pure-precession value `L⊥² ~ λ(λ+1) = 2` (p.179–180, pr. 160–161). `q_G` is inverted to a Coriolis constant `ζ₂₁ ≈ 0.137` and compared to CaOH (0.1969) and SrOH (0.179) (p.152, pr. 133). Parity splitting is estimated analytically as `ΔE± ≈ 4Bεω₂/A ~ 700 MHz` and compared to the observed ~500 MHz (Eq. 4.23, p.171; p.175). **This is his verification style: order-of-magnitude closed forms and cross-species trends, not byte-identical reproduction.** It matches the constraint in the shared brief.

---

## 7. Polyatomic extensions

- **Basis.** `|ψ_elec⟩ ⊗ |ψ_vib⟩ ⊗ |ψ_rot⟩` (Eq. A.1). Vibrational: `|v₁, v₂^ℓ, v₃⟩ = |v₁⟩⊗|v₂,ℓ⟩⊗|v₃⟩` (Eq. A.3, p.288). Bending state `|v₂,ℓ⟩ = (1/√(2π)) e^{iℓφ} Ψ_{v₂,ℓ}(q)` with `Ψ_{v,ℓ}(q) = (−1)^{(v+|ℓ|)/2} N_{v,ℓ} q^{|ℓ|} e^{−q²/2} L^{|ℓ|}_{(v+|ℓ|)/2}(q²)` (Eq. 2.22–2.23, p.70; Eq. A.27–A.28, p.296).
- **Projection algebra.** `K = Λ + ℓ = N⃗·n̂`; `P = Λ + Σ + ℓ = J⃗·n̂`. Vibronic term symbol `^{2S+1}K_P` with `K = Σ, Π, Δ, Φ, …`, the polyatomic analogue of `^{2S+1}Λ_Ω`; the `±` superscript applies only to `K=0` (Σ) vibronic states (Eq. 2.29, p.74, pr. 55).
- **ℓ values.** `v₂=1 → ℓ=1`; `v₂=2 → ℓ=0,2`; `v₂=3 → ℓ=1,3`; `v₂=4 → ℓ=0,2,4` (p.69, pr. 50).
- **Ladder operators.** `q±|v₂,ℓ⟩ = √((v₂+2±ℓ)/2)|v₂+1,ℓ±1⟩ + √((v₂∓ℓ)/2)|v₂−1,ℓ±1⟩` (Eq. 2.24, p.70) — note `q₊` **must** raise `ℓ` but can raise *or* lower `v₂`. True `Δv₂=±1, Δℓ=±1` ladders require momentum: `F^{±(±)} = q_{(±)} ∓ i p_{(±)}`; Ref. [48] (Yamada) defines `R^{±(±)}` with an extra `∓i`, `F = ∓iR`. Matrix elements from Di Lauro & Mills, all in the `δ_ℓ=0` convention (p.70, pr. 51).
- **Vibrational angular momentum from Coriolis constants.** `G⃗_ℓ = Σ_i Σ_j ζ_{ij} Q_i P_j`, with `ζ_{ij} = −ζ_{ji}` and, for linear molecules, `Σ_i ζ²_{ij} = 1` (Eq. 2.25, p.71, pr. 52).
- **`K`-labelling rule for higher `v₂`** — a real implementation trap he calls out (p.166, pr. 147). For `v₂=1` all `|K|>0` states are "unique": `|K| = |ℓ| + |Λ|`. For `v₂=3` you can have `|K| = |ℓ₂| − |Λ|` with `|ℓ₂|=3` *and* `|K| = |ℓ₁| + |Λ|` with `|ℓ₁|=1`, both giving `|K|=2`. His disambiguation rule: **write `ℓ₁` states with `Λ=+1` as the first ket of the parity superposition (the one not multiplied by the parity phase), and `ℓ₂` states with `Λ=−1` as the first ket.**
- **8-level bookkeeping in `v₂=1`.** Three projections (Λ, ℓ, Σ) × two signs = 8 levels; his convention is that the first ket in the superposition always has `Λ=+1`, and the rest follows from `K` (Eq. 4.16–4.18, p.165–166). He notes explicitly that writing μ as `²Σ_{−1/2}` "makes it clear the `(−1)^{p_a}` phase factor is on the `P>0` ket, unlike the case in `²Δ`" — sign bookkeeping that a code must reproduce.
- **K-resonance.** Both Λ-type and ℓ-type doubling mix `²Σ` and `²Δ`; "The full K-resonance Hamiltonian, `H_K = H_Λ + H_ℓ`, has many off-diagonal couplings, **requiring full diagonalization** of the effective Hamiltonian" (p.170, pr. 151). Regime: `A ≫ |εω₂|` groups vibronic states by `|Ω|`; when `A > ω₂` (BrCN⁺, ICN⁺) spin-orbit must enter at zeroth order. Severe-mixing limit gives `⟨|K|⟩ ≈ 1` and a level pattern that "look[s] like … a case (b) ²Π state" — the thesis's headline claim: "**To our knowledge, this is the first observation of a K-resonance severe enough that it changes the energy level pattern from case (a) to case (b)**" (p.182, pr. 163).
- **Parity for polyatomic bending states.** Covered by conventions items 9–11 above; `p = J−S−ℓ` (case a) / `p = N−ℓ` (case b). The `(−1)^ℓ` factor means "the action of the parity operator on a singly excited bending mode is similar to that of a `Σ⁻` electronic state" (p.148, pr. 129). Basis-change phase factors `(−1)^{P−1/2}` are inserted **by hand** in Eq. 4.12–4.14 "explicitly included to preserve parity" and "only valid for `P=1/2`" (footnote 5, p.160).
- **Orientation observable for a bending state.** `n̂·Ẑ = (N⃗·Z⃗)(N⃗·n̂)/(N(N+1)) ∝ M_N ℓ`; free-field `⟨M_Nℓ⟩ = 0`; fully mixed doublets give `2N+1` orientations along `M_Nℓ/(N(N+1))` (p.156, pr. 137). Asymptotic eEDM sensitivity in case (b): `|Σ| → S/(N(N+1)) = 1/4` for `N=1` (Fig. 5.27 caption, p.268).

**Spin statistics: the thesis does not cover it.** `grep -ciE` over the full text for `spin statistic` → 0, `nuclear statistic` → 0, `identical nuclei` → 0, `symmetrization` → 0 (none of these contain a `fi` ligature, so the zeros are trustworthy; control `Wigner-Eckart` → 25 on the same file). "homonuclear" appears twice, both in the §1.3.4 analogy between quadrupole-deformed nuclei and homonuclear molecules (p.36). The molecules studied (YbOH, CaOH) have no equivalent nuclei, so the topic never arises. A generic engine that needs spin statistics will get **no guidance from this thesis**.

---

## 8. Ions / target-species search result

**Queries run** (on `notes/thesis-text/_full.txt`, the complete 344-page text; and on a page-tagged copy `_num.txt` for context):

```
grep -ciE "ThF"              -> 0
grep -ciE "HfF"              -> 1
grep -ciE "3\s*Delta"        -> 1   (a bibliography entry, "H3Delta1")
grep -c   "varDelta_1"       -> 5   (the real in-text ³Δ₁ mentions; 6 lines total with the above)
grep -niE "varDelta_1|ThO"   -> 8 lines
grep -ciE "Omega.?doubl"     -> 1
grep -ciE "molecular ion"    -> 5
grep -ciE "JILA"             -> 5
grep -ciE "Cornell"          -> 1   (a Cornell University Press citation, not Eric Cornell)
grep -ciE "co-?magnetometer" -> 0
grep -ciE "NaF"              -> 0
grep -ciE "cation"           -> 97  (all false positives: "cation" inside "applications", "identification", …; the real uses are H₂⁺ on p.61 and the generic "cation" for H₂⁺)
```
Control: the same query shape on a known-populated target (`Wigner-Eckart` → 25, `hyperﬁne` → 176) returns hits, so the zeros are real absences, not a broken query.

**Result: ThF⁺ is never mentioned in the thesis. Not once.** There is also no NaF, no Ω-doublet-as-internal-co-magnetometer discussion by that name, and no molecular-ion effective Hamiltonian anywhere.

**Every hit, read in context:**

1. **p.210 (pr. 191), §5.2 — the only substantive ³Δ₁ / HfF⁺ passage in the thesis.** Quoting the key sentence fragments (≤15 words each): "The `³Δ₁` states in ThO and HfF⁺, for example, both provide two important features". Feature 1: "the states have g-factor cancellation arising from the anti-alignment of Λ and Σ, resulting in `g ~ μ_B/50`". Feature 2: "The states have small parity-doublets from the **Ω-doubling**, with the parity splitting measured to be ≈400 kHz in the lowest `J=1` state of ThO." Consequence given: "The parity doublets allow[] for full polarization in low fields, as well as control over the sign of the EDM interaction (Σ) without changing external fields." He then describes the ACME sequence: molecules "begin in the `X ¹Σ` ground state and are first optically pumped into the `H ³Δ₁` state, similar to how we have to populate the science state in metal hydroxides."
   **This is the entire ³Δ₁ treatment. There is no ³Δ₁ effective Hamiltonian, no Ω-doubling matrix element, no `g_eff` formula for a ³Δ₁, no S=1/Λ=2 matrix element anywhere in the thesis.**
2. **p.211 (pr. 192)** — "the suppression of the g-factor in the `H ³Δ₁` state aids in this aspect", in the context of the `|X⟩/|Y⟩` aligned superposition being maximally B-sensitive.
3. **p.212 (pr. 193)** — "we soon ran into difficulties … that are absent in ThO and other similar diatomics", the unresolved-hyperfine motivation for pivoting away from CPT. Also: "The JILA experiment performs state preparation in an entirely different way, essentially using a rotating electric field to implement `π/2`-pulses". He explicitly declines to elaborate: "We do not discuss their scheme further, as it leverages aspects unique to their experiment."
4. **p.270 (pr. 251), §5.7.4.2** — "The electric field dependence of `g_eff` can mimic an eEDM signal when combined with other systematic effects, very much like in `³Δ₁` molecules." This is the only place ³Δ₁ physics informs his own modeling, and it is a systematics analogy, not a structure calculation.
5. **p.42 (pr. 23)** — "Already, molecular ion experiments obtain seconds long coherence times …, though the number of ions in their traps are limited by Coulomb repulsion." The only statement in the thesis about molecular ions as a platform.
6. **p.40 (pr. 21)** — combining the JILA and ACME constraints to bound the `(d_e, C_S)` plane: `|d_e| < 2.1×10⁻²⁹ e cm` and `|C_S| < 1.9×10⁻⁹`.
7. **p.203 (pr. 184), p.234 (pr. 215)** — abstract-level mentions of JILA/ACME as the two schemes reviewed.
8. Remaining hits are bibliography entries (Grau's HfF⁺ thesis, Vutha's ThO magnetic/electric dipole moments paper, Edvinsson & Lagerqvist ThO rotational analysis, PRL 119 153001).

**Bottom line for the coordinator.** For the ThF⁺ X ³Δ₁ target, this thesis supplies the *framework* (case (a) basis, parity symmetrization, contact transformation, matrix-element recipe, Λ-doubling operator forms, field terms, diagonalization/labelling procedure) but **zero species-specific content and zero triplet/Δ-state matrix elements**. The `S=1` case is explicitly out of scope in one place: he notes he need not consider spin-spin "when working with a state with a single unpaired valence electron spin" but that this changes for "spin-spin interactions in triplet systems" (p.78, p.82). All his worked matrix elements assume `S=1/2`. For `Λ=2` he does point at the right literature: Brown, Cheung & Merer, "Λ-Type doubling parameters for molecules in Δ electronic states", J. Mol. Spectrosc. 124, 464 (1987) — Ref. [18] of Ch. 2 (p.87, pr. 68) — and Brown & Merer 1979 for Λ-doubling in ³Π and higher multiplicity (Ref. [16], p.87). He also states the general suppression argument that applies directly to `Ω=1` in a `³Δ`: connecting `P = 3/2 → −3/2` needs "two effective Hamiltonian terms, effectively a four operator interaction", giving "a much smaller parity splitting" (p.83, pr. 64).

---

## 9. Discrepancies with the notes file

I checked ~30 load-bearing claims in `~/.claude/notes/about-arian-physics.md` against the PDF. **The note is accurate on essentially everything I checked.** Verified correct: the X(010) parameter table (T₀, B, γ, γ_G, q_G, p_G all match Table 4.1, p.152); the Ã(010) parameter table (all ten entries match Table 4.4, p.179); 39 transitions / 6.1 MHz residual (p.151); `D_mol = 2.16(1) D` and `g_S = 2.07(2)` (p.155); `⟨Σ⟩ = 0.40` at `E = 101 V/cm` (p.162); mixing coefficients `(0.28, −0.49, 0.83)` and the 69/24/7% split (p.161); `b_F = 4.07(18)`, `c = 3.49(38)` for X(010) vs `4.80(18)` for X(000) (p.250); zero-g crossings at 59.6 and ≈64 V/cm (p.272); YbF `γ = −13.4 MHz`, crossing at 866 V/cm (p.265); 30 ms coherence time (p.273); the parity phases `p = J−S−ℓ` / `p = N−ℓ`; `E* = σ_xz R_y(π)`; `δ_ℓ = 0`, `δ_v = 0`; the anomalous-commutation rule and which momenta are exempt; the Λ vs ℓ azimuthal phase signs; the R² choice and the `−K²` origin offset; the H_RT form; the Curl relation `q ~ pB/A`; the Λ/ℓ/RT selection rules; the case-(a)↔(b) transform; `\|d_e\| < 4.1×10⁻³⁰ e·cm`.

Six items to flag, all minor:

1. **[correction] "R² … matches Steimle's YbOH paper" — the matrix-element source is a different paper.** The thesis says the YbOH origin-band paper used the R² form "with matrix elements taken from the Appendix of Ref. [7]" (p.298), and the appendix reference list (p.335) makes Ref. [7] = **Brown, Kopp, Malmberg & Rydh, Phys. Scripta 17, 55 (1978), on AlF**. So: the *choice* of R² follows the YbOH origin-band paper; the *matrix elements* come from the AlF paper. The note also omits that he flags a typo in that appendix (`Σ² − S(S+1)` should read `ΩΣ − S(S+1)`).
2. **[correction] "Sears-resonance interaction `q±S∓`" attribution.** The note reads as if the Sears resonance is a Renner-Teller variant. The thesis (p.164, pr. 145) says it is a *vibronic spin-orbit* term coupling `ℓ` and `Σ` (mixing `Ω=1/2` and `Ω=3/2`), derived in Refs. [72, 84], "most prominent … when spin-orbit and vibrational splittings are approximately equal", and named a Sears resonance in the effective-Hamiltonian picture. The GeCH/BrCN⁺/ICN⁺ examples are correct, but BrCN⁺/ICN⁺ are cited at p.167 for the `A > ω₂` regime, not as Sears-resonance examples.
3. **[precision] "he uses the R² formalism" is stated flatly; the thesis also gives the N² alternative and the conversion.** `R⃗ = N_x x̂ + N_y ŷ + (N_z − K)ẑ`; conversion formulae for centrifugal/distortion parameters are deferred to B&C §7.5.3 rather than reproduced (p.298).
4. **[precision] "Modeled X(010) as ²Π_{1/2} with substituted Λ→ℓ".** Correct, but the note omits that the thesis calls the approximation out as wrong in principle: the E1 transition "should couple via `T¹_{q=±1}(d)` components", and the substitution makes it parallel (`q=0`) instead (p.241, pr. 222).
5. **[precision] Table 4.4 error bars.** The note reproduces the values without the footnote; the thesis footnote says "Parentheses represent 2σ errors (68% confidence interval)" (p.179) — internally inconsistent (2σ ≠ 68%). Worth not copying either interpretation forward without asking.
6. **[gap, not error] The note's "canonical references" section lists Brown & Carrington eq. numbers "5.169, 5.174, 5.94, 7.122, eq. 8.459"** — all of these do appear in the thesis (App A.4 items 1–5 and p.180), so the note is right, but it omits the three places where the thesis says a textbook equation is **wrong**: B&C's case-(a)↔(b) transform (p.75), B&C eq. 5.177 (p.303), and Hirota tab. 2.4 item 2 (p.300). Those are the highest-value items in the note for code purposes and should be promoted.

Not a discrepancy but worth recording: the note's file path for the PDF (`~/Downloads/Arian_Jadbabaie_Springer_Thesis.pdf`) is not the Zotero path used here; content matches.

---

## 10. Open questions for Arian

1. **What code produced the Ch. 4 and Ch. 5 diagonalizations?** The thesis names QuTiP, Mathematica/AtomicDensityMatrix and LightTools but never the structure code. Is it one of the two repos in scope, a Caltech-era private script, or Pilgram's/Augenbraun's code?
2. **State tracking.** The only tracking rule in the thesis is "order eigenstates … according to their adiabatically correlated free field state" (p.267). Was that implemented as eigenvector overlap against the zero-field basis, step-to-step overlap along the sweep, or by hand for a small manifold? Does it survive a genuine avoided crossing (the K-resonance in Ã(010) is exactly the case that would break naive energy ordering)?
3. **`s` in the parity phase.** Eq. A.15 carries `(−1)^{J−S−ℓ+s}` with `s=1` for `Σ⁻`, but §2.2.3 and §4.4.2 both write `p = J−S−ℓ` with no `s`. Is `s` intended to be folded in always, or only when a `Σ⁻` component is present in the basis? This matters for a code that must handle `Σ⁻` and bending states simultaneously.
4. **The `θ → θ` in the `H_ℓ` replacement rule (p.169).** I read this as a typo for `θ → φ`. Confirm, since it determines which azimuthal phase (`−1` for Λ, `+1` for ℓ) attaches to the `p_G`/`q_G` operators.
5. **Why does Eq. 5.29 omit `p_G`** (and `γ_G`) when Ch. 4 measured both to be non-zero in YbOH X̃(010)? Was that a deliberate simplification for CaOH (where they may be negligible), or an inherited Hamiltonian? A generic engine needs to know whether `p_G` belongs in the default bending-mode term list.
6. **`Ω`/`Σ` sign conventions for `S ≥ 1`.** All worked matrix elements assume `S = 1/2`. For a `³Δ₁` (`Λ=2, Σ=−1`), does the parity phase `(−1)^{J−S−ℓ+s}` and the "first ket has `Λ=+1`" bookkeeping rule carry over unchanged, or does the `|Σ| = S` simplification used in deriving Eq. A.15 (`S − 2Σ = −S`) break when `Σ ≠ ±S`?
7. **Which Λ-doubling operator form for `Λ = 2`?** The thesis points at Brown, Cheung & Merer 1987 for Δ states but never writes the operator. Does he want the code to implement the `ΔΛ = ±4` (or two-step `ΔΛ=±2` composed) doubling for `³Δ`, or to treat Ω-doubling in `³Δ₁` phenomenologically as a single fitted splitting?
8. **Units convention.** The thesis mixes cm⁻¹ (Ch. 4 vibronic fits) and MHz (Ch. 4 X-state fits, Ch. 5) freely, with `1.09 h MHz/(V/cm)` for the dipole. What should the engine's canonical internal unit be, and should parameters carry units as metadata?
9. **Basis truncation policy.** `N″=6`, `J′=15/2`, "include the `P=3/2` manifold". Were those convergence-tested, or chosen by eye? A generic engine should have a convergence check; is there a precedent he trusts?
10. **Hyperfine scope.** No `eQq`, no nuclear spin-rotation, no `a`/`d` parameters anywhere in the thesis. For ThF⁺ (¹⁹F, `I=1/2`) `eQq` is irrelevant, but for the "other diatomic and polyatomic" targets it will not be. Is extending the hyperfine term list part of scope-one, or deferred?

---

## 11. Method log

**Extraction.** `PyMuPDF (fitz)` in conda env `claude-code`, via a script at the scratchpad path (not written into any repo). Worked first try; `pypdf`/`pdfplumber` never needed. Two scripts: `extract.py` (per-section text + `_full.txt`) and `render.py` (8 page PNGs).

**Text extracts written** to `notes/thesis-text/`:
`_full.txt` (all 344 pages, `=== PDFPAGE N ===` markers), `_num.txt` (same, with `lineno|pNNN|` prefixes for grep-with-page-number), `front.txt` (1–19), `ch1-intro.txt` (20–53), `ch2-molecules.txt` (54–89), `ch3-production.txt` (90–140), `ch4-yboh.txt` (141–202), `ch5-spm.txt` (203–281), `ch6-concl.txt` (282–287), `appA-matrix-elements.txt` (288–310), `appB-F.txt` (311–336), plus focused slices `sub-ch1-ptodd.txt` (31–48), `sub-ch4-modeling-432.txt` (148–153), `sub-ch4-starkzeeman.txt` (154–163), `sub-ch4-RT-442.txt` (164–172), `sub-ch4-fit-445.txt` (175–184), `sub-ch4-tdm-447.txt` (185–190), `sub-ch5-edm-zerog.txt` (209–215), `sub-ch5-zerog-574.txt` (263–275), `sub-ch5-masterEq.txt` (239–245).

**Read in full** (every line): Chapter 2 in its entirety, PDF pp. 54–89 (pr. 35–70) including the reference list; Appendix A in its entirety, pp. 288–310 (pr. 269–291); §4.3.2–4.3.4 pp. 148–163; §4.4.1–4.4.6 pp. 164–184; §1.3.1–1.3.4 pp. 31–37; §1.4.2–1.4.3 pp. 43–47; §5.7.4–5.7.5 pp. 265–273; §5.6.3.2 pp. 240–243; §5.2–5.3 pp. 210–213.

**Skimmed / grep-targeted only**: Chapter 3 (production, pp. 90–140) — read only p.132 (isotopologue/`G_Yb` coupling); §4.4.7 and §4.5 (pp. 185–197); Chapter 5 §5.1, §5.4–5.6.6 outside the passages above; Chapter 6; Appendices B–F except the App F reference list on p.335.

**Rendered to PNG and viewed with the Read tool (4 of a 25-page budget used)**: p.292 (Eq. A.12–A.15, the parity phases), p.149 (Eq. 4.4–4.5, the X̃(010) Hamiltonian), p.168 (Table 4.2, the K-resonance matrix), p.169 (Eq. 4.20, `H_Λ`). Four more were rendered but not needed after the LaTeX alt-text proved reliable: pp. 75, 79, 266, 305. PNGs live in the scratchpad `pages/` directory.

**Searches run** (all on `_full.txt` / `_num.txt`, whole document): the item-7 species battery (`ThF|HfF|3\s*Delta|varDelta_1|Δ1|Omega.?doubl|Ω.?doubl|molecular ion|cation|JILA|Cornell|co-?magnetometer|NaF|science state`); software (`pgopher|qutip|mathematica|python|scipy|numpy|matlab|julia|github|sympy|our code|custom code|software|repositor|jupyter|notebook`); numerics (`Nelder|Levenberg|least.squares|F-test|diagonaliz`); typo/disagreement flags (`typo|disagree|differs from|error in Ref|incorrect`); hyperfine and PT-odd coverage (`eQq|quadrupole|nuclear spin.rotation|spin.spin|W_[A-Za-z{]|Welec|anapole|NSD|C_S|Eeff`); spin statistics (`spin statistic|nuclear statistic|identical nuclei|symmetrization|homonuclear|Pauli`); parameter values (`b_F|4\.07|3\.49|coherence time|30 ms`).

**One methodological correction I made mid-run, worth passing on**: my first `grep -ci "hyperfine"` returned **0**, which would have been a badly wrong absence claim. The PDF uses the `ﬁ` ligature (U+FB01), and a byte-oriented `.` does not match a 3-byte UTF-8 character, so `hyper.ne` also returns 0. `grep -c "hyperﬁne"` returns 176. I re-validated every absence claim in this digest against that failure mode: none of the words in my zero-result queries (`ThF`, `PGopher`, `eQq`, `spin statistic`, `identical nuclei`, `symmetrization`, `nuclear spin-rotation`, `W_d`, `W_s`, `g_r`, `g_N`, `NSD`, `NaF`) contains `fi` or `fl`, and I ran a positive control (`Wigner-Eckart` → 25) on the same file with the same query shape.

**Could not access / did not do**: nothing was blocked. The `pdf-mcp` and `zotero` MCP servers both failed to connect this session (HTTP 401), as the brief anticipated, so all PDF work went through PyMuPDF. The PDF and everything under `~/Zotero` were read-only throughout; no files written outside `notes/thesis-text/`, `notes/digest-thesis-effective-hamiltonian.md`, and the session scratchpad.

**Approximate effort**: ~35 tool calls, ~330k tokens of PDF content passed through context (of which ~150k read in full), roughly 45 minutes wall clock.
