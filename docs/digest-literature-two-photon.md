# Effective two-photon (2 × E1) transitions within X ³Δ₁ — literature digest

Assembled 2026-09-05 for `heff`. Purpose: give the tutorial notebook a physically
defensible **effective two-photon (Raman-type) spectrum** within X ³Δ₁, alongside
the existing one-photon E1 spectra.

**Reading rules for this document.** Every factual claim carries an inline cite to
a page/equation I read at source in this session. Anything I derived is tagged
`[derived]`. Anything I could not verify is tagged **UNVERIFIED** and is *not*
used to support a downstream conclusion. Where a matrix element or a number is
simply not in the literature I read, the row says **GAP** and section 5 records
the exact query I ran. I have not guessed a single matrix element.

**Sources read at source in this session** (local PDFs under `docs/lit/`, text
extracted with PyMuPDF, tables rendered as images where the PDF text layer
scrambles them):
**Cossel PhD thesis** ("Techniques in molecular spectroscopy: from broad
bandwidth to high resolution", Colorado 2014, downloaded from the public JILA
Cornell-group thesis list, PDF pp. 179, 184, 218–226, 234–242) — **the single
most important source for this task**;
Cairncross PhD thesis (JILA/Colorado, "Searching for time-reversal symmetry
violation with molecular ions", PDF pp. 34, 50–54, 70–73, 94–96, 124–128,
149–153); Ng PhD thesis (JILA/Colorado, PDF pp. 55, 60, 70–71, 118–126);
Ng et al. PRA **105**, 022823 (2022) §II and Fig. 6; Gresh et al. JMS **319**, 1
(2016) Table 2 (rendered); Denis et al. NJP **17**, 043005 (2015) Tables 3, 9, 10
(Table 9 rendered as an image — its PDF text layer is rotated and scrambled);
Petrov & Skripnikov arXiv:2503.02840; Leanhardt et al. JMS **270**, 1 (2011)
= arXiv:1008.2997; Brown & Carrington §5.5.4–5.5.6 and Appendix 5.1 (local
extract `docs/lit/bc-pages/bc_p195-210_ch5-matrixelements.txt`, Eqs. 5.141,
5.142, 5.165–5.179); Bray & Hochstrasser, Mol. Phys. **31**, 1199 (1976),
**abstract only** — body paywalled.

---

## 0. The headline, before the detail

1. **No published experiment drives a two-photon transition with *both* endpoints
   inside X ³Δ₁.** Every Raman/two-photon step in the HfF⁺ and ThF⁺ eEDM
   programmes is either (a) an *inter-electronic* transfer ¹Σ⁺ → ³Π₀₊ → ³Δ₁
   (HfF⁺ Gen 1), or (b) a two-photon *dissociation* (REMPD) out of ³Δ₁ for
   readout. Transfer within ³Δ₁ — Ω-doublet, hyperfine, ΔJ, Δm_F — is done with
   **microwaves, incoherent optical pumping, and E_rot ramps**, not optical
   Raman.
2. **Ng's thesis explicitly proposes the within-X two-photon Raman that Arian
   wants, and says the spectroscopy to enable it has not been done** (Ng thesis
   p. 102, "We would need to do some spectroscopy to make this happen"). That is
   the state of the art: a proposal, no measured rate, no chosen intermediate.
3. **The formalism, however, is fully written down and experimentally
   validated** — in **Cossel's PhD thesis §6.3.5, Eqs. (6.34)–(6.38)**, whose
   final states *are* the X ³Δ₁ J = 1 Stark/hyperfine manifold. It is the exact
   `m1·m2·m3` decomposition `heff.elements_c.dipole_geometry` already uses, with
   the molecule-frame `q` left free instead of frozen at 0, and with the
   amplitudes summed over intermediate sublevels and polarisations **before**
   squaring — the two-photon form of `heff`'s gate-B8 convention. Cossel also
   documents a *measured* destructive interference between two σ pathways
   (Fig. 6.18, p. 224). §3 builds on this.
4. **Transition dipoles to Ω = 0 intermediates exist.** ThF⁺, *ab initio* only
   (Denis 2015 Table 9): Ω = 0 states at 6344 / 6528 / 6747 cm⁻¹ with
   |⟨i|D|³Δ₁⟩| = 0.455 / 0.571 / 0.391 D. HfF⁺, calculated (Cossel Table 5.6):
   ³Π₀₋ and ³Π₀₊ both 0.27 e·a₀, plus a higher pair at 0.36 / 0.39 e·a₀ — and
   ³Π₀|d|³Δ₁ ≈ **0.27(3) e·a₀ measured** (Cossel p. 218). These are exactly the
   ΔΩ = ±2 pathway Arian names.
5. **The biggest gap is the ThF⁺ intermediate ladder itself**, not the
   formalism: Denis (Ω = 0 at 6344–6747 cm⁻¹) and Petrov (³Π₀∓ at 3044/3395
   cm⁻¹) disagree irreconcilably, nothing is measured between 3150 and
   10 472 cm⁻¹, no ThF⁺ excited-state lifetime or hyperfine constant is
   published, and Denis does not label 0⁺ vs 0⁻ — which is the sign that decides
   whether two ΔΩ = ±2 amplitudes add or cancel. §5.
6. **A structural result worth knowing before plotting** `[derived, numerically
   checked this session]`: the 2 × E1 operator is **parity-even**, so at zero
   field it connects e→e and f→f and never e↔f; and with only σ± available in
   the E_rot plane the reach is **Δm_F ∈ {0, ±2}**, never ±1 or ±3. §3.3.

---

## 1. Experimental schemes: what the JILA experiments actually do

Two-photon / Raman steps are marked **2γ**; everything else is listed because it
is what the experiments use *instead* of a within-X Raman, and the notebook
should not misrepresent that.

| # | Molecule / generation | Transition | Type | Intermediate state | Wavelengths / detuning | Polarisation geometry | Source |
|---|---|---|---|---|---|---|---|
| 0 | HfF⁺ (original **proposal**) | ¹Σ₀(J=0) → ³Δ₁(J=1), the two \|m_F\|=3/2 sublevels | **2γ** off-resonant stimulated Raman, as *both* Ramsey π/2 pulses | **¹,³Π₁** — "off-resonant from the intermediate ¹,³Π₁ states" | not specified | not specified | Leanhardt et al. JMS 270, 1 (2011) = arXiv:1008.2997 p. 8 |
| 1 | HfF⁺ Gen 1 | ¹Σ⁺(v=0, J=0) → ³Δ₁(v=0, J=1, F=3/2) | **2γ** stimulated Raman, driven adiabatically by the ions' own oscillating Doppler shift | **³Π₀₊** (v = 1, J = 1) | pump ≈ 899.6665 nm ("Stella", S), Stokes ≈ 986.4175 nm ("Toptica", T); one-photon detuning **Δ = 160 MHz initially, raised to ≈ 1.5 GHz** to suppress spontaneous emission; second laser detuning set by an AOM to select the **upper or lower Stark doublet**. Rabi rates **Ω₁₂ = 2π×6 MHz, Ω₂₃ = 2π×19 MHz, two-photon Ω₁₃ = 2π×210 kHz** | two cw lasers **co-propagating along Ẑ** (Λ-type, so the two-photon transition is Doppler-free), **both linearly polarised at a relative angle θ**; transfer ∝ **sin²θ** — zero at θ = 0, maximum at θ = π/2, confirmed experimentally | Cairncross thesis pp. 61, 63, 115; **Cossel thesis pp. 205–224, esp. §6.3.3–6.3.5** |
| 1b | HfF⁺ Gen 1 | as #1, **used as a spectroscopic probe** | **2γ** | ³Π₀₊ | second laser scanned through two-photon resonance with each hyperfine/Stark level of ³Δ₁ J = 1 | as #1 | Cossel thesis §6.3.5, Fig. 6.17 — this is the "Raman spectrum of ³Δ₁(v=0, J=1)" Cairncross p. 42 points at: **a two-photon spectrum whose final states are the X ³Δ₁ J = 1 Stark manifold**, fitted to extract d_mf = 1.401(5) e·a₀ and A∥ = 62.2(1) MHz |
| 2 | HfF⁺ Gen 1 | same as #1 | **2γ** attempted **STIRAP** and **Raman ARP** | ³Π₀₊(v=1, J=1) — chosen because it has "moderately strong coupling to both ¹Σ⁺(v=0,J=0) and ³Δ₁(v=0,J=1)" | Doppler width ≈ 40 MHz for the one-photon NIR steps but only **≈ 5 MHz for the ¹Σ⁺→³Π₀₊→³Δ₁ two-photon (Raman) transition**; ARP swept the Stokes laser across the Doppler width, sweep durations ~0.1–1 ms | as #1 | Cairncross thesis pp. 114–118 |
| 3 | HfF⁺ Gen 1 | m_F-selective **depletion** out of ³Δ₁(J=1, F=3/2) | **1γ**, not Raman | ³Σ⁻₀₊ (Q(1) line) | resonant | **circularly polarised**, strobed synchronously with E_rot so k is ∥ or ∦ E_rot ⇒ drives σ± to F′ = 3/2 | Cairncross thesis p. 63 |
| 4 | HfF⁺ Gen 1 | π/2 pulse: \|m_F = +3/2⟩ → (\|+3/2⟩+\|−3/2⟩)/√2, i.e. **Δm_F = 3 within J = 1** | **not photonic at all** — E_rot amplitude ramp | none | E_rot reduced briefly; rotation-induced coupling Δ_u/l between m_F = ±3/2 does the transfer in ≈ 1 ms | — | Cairncross thesis p. 63 |
| 5 | HfF⁺ Gen 1 | readout: ³Δ₁ → dissociation | **2γ** REMPD (bound → bound → continuum) | \|Ω_n\| = 0, 1 or 2 depending on line; best doublet contrast with \|Ω_n\| = 2, R(1) | pulsed UV | ε₁, ε₂ chosen per Table 6.1: (σ±, π), (σ+, σ+), (σ+, σ−), (π, π) all tabulated with closed-form I(θ,φ) | Cairncross thesis pp. 140–142, Table 6.1 |
| 6 | ThF⁺ Gen 3 | ³Δ₁(v=0) → X ³Δ₁(v=0, J=1, F=3/2, m_F=−3/2, Ω=1) | **1γ optical pumping + microwaves**, *no* Raman | **Ω = 0⁻ at ≈ 14 600 cm⁻¹** | cw 685 nm pumping, 717 nm repump; **microwaves at 29 GHz (J=1↔2) and 43 GHz (J=2↔3)** | not stated in the paper text I read | Ng 2022 §II and Fig. 1; Ω=0⁻ energy from Gresh 2016 (14 589.09 cm⁻¹) |
| 7 | ThF⁺ (proposed, **not done**) | \|J=1, F=3/2, m_F=+3/2⟩ → \|J=1, F=3/2, m_F=+1/2⟩ — **within X ³Δ₁** | **2γ** Raman / STIRAP | "One would need to **find** a good intermediate state" | — | Constraint: JILA optical access puts all lasers and microwaves **in the plane of the rotating E_rot**, so only **σ± or x̂ polarisations are available, never π**. Hence a naive σ⁻ + π two-photon path is impossible; Ng's alternative is "a two-photon Raman process, **using photons of opposite helicities**" | Ng thesis p. 102 and its footnote 4 |
| 8 | ThF⁺ (proposed, **not done**) | \|J=1, m_F=1/2⟩ ↔ \|J=2, m_F=1/2⟩ (the ≈ 29 GHz CP-FIT rotational transition) | **2γ** Raman as an alternative to π-polarised microwaves | not chosen | — | opposite-helicity photon pair (σ⁺σ⁻ or two π, the latter unavailable) | Ng thesis p. 102 |

**The load-bearing sentence for this whole exercise**, Ng thesis p. 102, verbatim:
"An alternative to the π-polarized microwaves is to use a two-photon Raman
process, using photons of opposite helicities. **We would need to do some
spectroscopy to make this happen.**" So: the scheme Arian wants to plot is a
documented JILA *intention* with no measured intermediate state, no measured
detuning, and no published rate.

**Δm_F reach, and why helicity matters** `[derived from the geometry constraint in Ng p. 102]`:
with only σ± available in the E_rot plane, a two-photon pair spans
Δm_F ∈ {−2, 0, +2}; adding π (not available at JILA, available in a general
calculation) opens Δm_F = ±1. The Δm_F = ±3 transfer that the eEDM π/2 pulse
performs is therefore **not** reachable by 2 × E1 at all — it is done by the
E_rot ramp (row 4). Worth saying in the notebook so the two-photon panel is not
mistaken for the π/2 pulse.

**Which Δ-Ω a given intermediate gives** `[derived]` — each E1 step carries a
molecule-frame component q_i = Ω_i − Ω (the 3j `(J' 1 J; −Ω' q Ω)` is non-zero
only there, `heff.elements_c.dipole_geometry` docstring; Cairncross thesis
p. 141 states the same rule for REMPD, "we can simply set q = Ω′ − Ω"), and
|q_i| ≤ 1 because each photon is rank 1:

| Ω_i of intermediate | q₁ (from Ω=+1) | q₂ (to Ω′) | ΔΩ = q₁+q₂ reachable |
|---|---|---|---|
| 0 | −1 | ±1 | **0 and −2** ⇒ this is the **only** route to Ω = +1 → −1 |
| 1 | 0 | 0 | 0 only |
| 2 | +1 | −1 | 0 only |

So the Ω-doublet-flipping (ΔΩ = ±2) two-photon operator **requires an Ω = 0
intermediate**, confirming Arian's statement. Ω = 1 and Ω = 2 intermediates
contribute only to the ΔΩ = 0 (orientation-preserving) part.

---

## 2. The ThF⁺ excited-state ladder

### 2.1 Measured (Gresh et al. JMS 319, 1 (2016), Table 2, rendered and read)

T₀ relative to X ³Δ₁ v = 0; all in cm⁻¹; quoted errors 95 %.

| State | T₀ (cm⁻¹) | B_e (cm⁻¹) | ω_e (cm⁻¹) | Ω | ΔΩ pathway from Ω=±1 |
|---|---|---|---|---|---|
| X ³Δ₁ | 0 | 0.24311(7) | 656.96(1) | 1 | — |
| a ¹Σ⁺ | 314.282(7) | 0.24601(8) | 657.90(32) | 0 | 0, ±2 (but E1 to ³Δ₁ is spin- and Λ-forbidden, see below) |
| ³Δ₂ | 1052.5(1.0) | 0.24342(2) | 657.38(22) | 2 | 0 |
| ³Δ₃ | 3150(30) | 0.24388(7) | 659.31(29) | 3 | none (needs q = 2) |
| Ω = 0⁺[10.47] | 10 471.889(7) | 0.23583(5) | 626.67(56) | 0 | 0, ±2 |
| Π₁[12.99] | 12 994.3(1.0) | 0.23049(2) | 580.33(10) | 1 | 0 |
| Π₁[14.22] | 14 223.6(1.0) | 0.23101(4) | — | 1 | 0 |
| Ω = 0⁺[14.29] | 14 295.967(8) | 0.23187(5) | — | 0 | 0, ±2 |
| Ω = 0⁻ | 14 589.09(2) | 0.22947(3) | 593.467(8) | 0 | 0, ±2 |
| ³Π₂ | 16 682(30) | 0.22858(7) | 582.04(75) | 2 | 0 |

The Ω = 0⁻ at 14 589.09 cm⁻¹ is the state Ng 2022 uses for optical pumping
("approximately 14600 cm⁻¹", Ng 2022 §II, 685 nm ⇒ 14 599 cm⁻¹ `[derived]`).

**GAP: no lifetime is quoted for any of these states in the sources I read.**
Gresh 2016 Table 2 lists T₀, T_e, B_e, ω_e, ω_eχ_e, α_e only — no Γ, no τ. See §5.

### 2.2 Ab initio (Denis et al. NJP 17, 043005 (2015), Table 9)

TZ′ basis, CI model *IV*^CI, R = 3.779 a₀, origin at the centre of mass. Diagonal
entries are static dipole moments ⟨D_z⟩ (negative, ≈ −4 D for the low states);
off-diagonals are |⟨i|**D**|j⟩| in Debye. Column reproduced here: the transition
dipole **to X ³Δ₁**.

| State (Denis label) | T_v (cm⁻¹) | Ω | \|⟨i\|D\|³Δ₁⟩\| (D) | ⟨D_z⟩ (D) | ΔΩ it enables |
|---|---|---|---|---|---|
| ¹Σ₀⁺ | 274 | 0 | **0.012** | −4.004 | 0, ±2 (but see note) |
| ³Δ₁ | 0 | 1 | (−4.075 diag) | −4.075 | — |
| ³Δ₂ | 724 | 2 | 0.070 | −4.022 | 0 |
| ³Δ₃ | 2198 | 3 | 0.000 | −4.075 | none |
| ¹Σ₀(³Π₀) | 6344 | **0** | **0.455** | −3.752 | **0, ±2** |
| ³Π₀ | 6528 | **0** | **0.571** | −2.116 | **0, ±2** |
| ¹,³Π₁(³Σ₁) | 6639 | 1 | 0.142 | −2.375 | 0 |
| ³Π₀(¹Σ₀) | 6747 | **0** | **0.391** | −2.717 | **0, ±2** |
| ¹,³Δ₂(³Π₂) | 7008 | 2 | 0.473 | −2.734 | 0 |
| ³Σ₁ | 7490 | 1 | 0.069 | −4.463 | 0 |
| ¹,³Π₁ | 7918 | 1 | 0.052 | −2.708 | 0 |
| ³Φ₂(³Π₂) | 8245 | 2 | **1.338** | −2.271 | 0 |

Denis's own commentary: the largest element,
|⟨³Φ₂|**D**|³Δ₁⟩| = 1.34 D, "is spin-allowed (ΔS = 0) and also orbital angular
momentum allowed (here ΔΛ = ±1). In addition, ΔΩ = ±1 is also satisfied", and
"very small transition moments are typically found for spin-forbidden
transitions" (Denis 2015 §3.2.3). The 0.012 D for ¹Σ₀⁺ ↔ ³Δ₁ is the quantitative
version of Ng 2022's statement that "the transition between X ³Δ₁ and a ¹Σ⁺ is
forbidden by selection rules" (Ng 2022 Fig. 6 caption).

**Three caveats on this table, all load-bearing.**
- Denis's T_v are *vertical* energies from a CI calculation and disagree with the
  measured T₀: ¹Σ⁺ at 274 vs 314.282(7) measured, ³Δ₂ at 724 vs 1052.5(1.0)
  measured, ³Δ₃ at 2198 vs 3150(30) measured (Gresh 2016 Table 2). Use measured
  T₀ for energies and Denis only for the *dipoles*.
- Denis labels the Ω = 0 states as ¹Σ₀(³Π₀), ³Π₀, ³Π₀(¹Σ₀) **without stating the
  reflection parity 0⁺ vs 0⁻**. The 0⁺/0⁻ distinction is exactly what decides the
  sign of each state's contribution to the Ω-doubling (see §2.3 of
  `docs/thf-plus-x3delta1-effective-hamiltonian.md`, σ_i = +1 for 0⁺, −1 for 0⁻)
  and it may well matter for the sign of the interfering two-photon amplitudes.
  **GAP.**
- The Denis Ω = 0 states sit at 6344–6747 cm⁻¹. **Petrov & Skripnikov
  arXiv:2503.02840 place their Ω = 0 perturbers far lower**: a ¹Σ⁺ at T₀ = 314,
  ³Π₀₋ at T₀ = 3044, ³Π₀₊ at T₀ = 3395 cm⁻¹ (arXiv:2503.02840 lines 179–184 of
  the local text extract). These two ab initio sets are **not consistent** on
  where the low-lying Ω = 0 manifold is, and neither ³Π₀± at 3044/3395 appears in
  Gresh's measured list. Any notebook number that depends on the intermediate
  energy must say which set it used.

### 2.3 HfF⁺ transition dipoles, for contrast

**Cossel thesis Table 5.6** (calculated molecule-frame dipole and transition
dipole matrix elements in HfF⁺, in e·a₀; "–" means < 0.01, "?" means non-zero but
not calculated), the columns that matter here:

| State | ↔ ¹Σ⁺ | ↔ ³Δ₁ | ↔ ³Δ₂ |
|---|---|---|---|
| ¹Σ⁺ (diagonal) | −1.2 | | |
| ³Δ₁ (diagonal) | 0.02 | −1.4 | |
| ³Π₀₋ | 0 | **0.27** | – |
| ³Π₀₊ | 0.15 | **0.27** | – |
| ³Π₁ | 0.25 | 0.03 | 0.22 |
| ³Π₂ | – | 0.03 | 0.01 |
| ¹Π₁ | 0.55 | 0.03 | 0.18 |
| ³Σ₀₊ | 0.23 | 0.08 | – |
| ³Φ₂ | – | **0.48** | 0.10 |
| ³Π₀₋ (second block) | – | **0.36** | – |
| ³Π₀₊ (second block) | 0.12 | **0.39** | – |
| ³Σ₀₋ | – | 0.17 | – |

The ³Σ⁻₀₊ row (0.23 / 0.08) is the pair Cairncross quotes at p. 118. **Unlike
ThF⁺, one of these is measured**: from the transfer data Cossel extracts
³Π₀|d|³Δ₁ ≈ **0.27(3) e·a₀** — "matches the theory" — while ³Π₀|d|¹Σ⁺ comes out
"low by about a factor of 2" against Table 5.6 (Cossel thesis p. 218). Note that
this table **does** carry the 0⁺/0⁻ labels that Denis's ThF⁺ table lacks, and it
gives ³Π₀₋ and ³Π₀₊ *equal* couplings to ³Δ₁ (0.27 each) in the lower block —
i.e. near-cancelling contributions to any quantity that weights them with
opposite sign.

---

## 3. The effective two-photon operator

### 3.0 What is actually written down in the literature, and by whom

There is **one** source I found that writes the two-photon transfer strength for
this exact family of states: **Cossel PhD thesis (Colorado, 2014), §6.3.5,
Eqs. (6.34)–(6.38)**, for HfF⁺ ¹Σ⁺ → ³Π₀₊ → ³Δ₁(J = 1). It is worth quoting the
shape verbatim because it is the template for `heff`:

> "Because the single-photon detuning is large compared to the single-photon
> transition width, we assume that the total rate of transfer W_f into a
> particular Stark level of ³Δ₁ is proportional to the modulus squared of the
> two-photon transition dipole moment,
> W_f ∝ |Σ_{i,m} Σ_{p_T,p_S} ⟨f|T⁽¹⁾_{p_T}|m⟩⟨m|S⁽¹⁾_{p_S}|i⟩|²." (Eq. 6.34)

and each one-photon factor is expanded (Eqs. 6.36–6.37) as the **same
three-factor product `heff.elements_c.dipole_geometry` already computes**:

```
⟨η′; I; J′,Ω′; F′,m′_F| T⁽¹⁾_p |η; I; J,Ω; F,m_F⟩ = m1 · m2 · m3
  m1 = (−1)^{F′−m′_F} ( F′ 1 F ; −m′_F  p  m_F )
  m2 = (−1)^{F+J′+I+1} √((2F′+1)(2F+1)) { J  F  I ; F′ J′ 1 }
  m3 = Σ_q (−1)^{J′−Ω′} √((2J′+1)(2J+1)) ( J′ 1 J ; −Ω′  q  Ω ) ⟨η′‖T⁽¹⁾‖η⟩
```
(Cossel thesis Eq. 6.37.) `heff`'s `dipole_geometry` is exactly `m1·m2·m3` with
**q frozen to 0** and the reduced element factored out; Cossel's `m3` keeps the
sum over q. **The generalisation `heff` needs is one line: let q = Ω′ − Ω.**

Cossel also gives the lab↔rotating-frame link and the polarisation dependence:
`S⁽¹⁾_p = Σ_{p′} D⁽¹⁾_{pp′}(ω_rot t)* S⁽¹⁾_{p′}`, `T⁽¹⁾_p` the same with
`ω_rot t + θ` (Eq. 6.35), from which he finds **W_f ∝ sin²θ**, θ the angle
between two linear polarisations — matching the measurement (no transfer at
θ = 0, maximum at θ = π/2; Cossel thesis p. 219).

**Interference is real and measured here.** Cossel's Fig. 6.18 (p. 224) draws the
pathways with amplitudes 1/√3, 1/√3, 1/√6, −1/√6 and states: "the two-σ pathways
(all blue lines) **cancel because of the signs of the Wigner 3j coefficients**,
so the transition to m_J = 0 is forbidden" (p. 223). That is the two-photon
analogue of `heff`'s gate B8 (sum the amplitudes, *then* square) and it is not a
formal nicety — it changes which lines exist.

Two more references, **abstract-level only** (I read the abstracts, not the
papers):
- **Bray & Hochstrasser, Mol. Phys. 31, 1199 (1976)**, "Two-photon absorption by
  rotating diatomic molecules": explicit two-photon cross sections in linear and
  circular polarisation for diatomics, covering "all branches of **ΔΩ = 0, ±1,
  ±2** transitions"; σ_λλ/σ_cc = 2/3 for every branch except the Q branch of a
  ΔΩ = 0 transition, where it becomes J-dependent and "may range from infinity to
  4/1 (for high J) for different relative contributions to the absorption cross
  section of intermediate states having either the same or different electronic
  symmetry" (abstract, tandfonline.com/doi/abs/10.1080/00268977600100931). That
  last clause is the statement that the **interference between intermediate
  symmetries is experimentally accessible in the Q branch** — directly relevant
  to a ΔΩ = 0 two-photon spectrum within X ³Δ₁.
- **Bonin & McIlrath, JOSA B 1, 52 (1984)**, "Two-photon electric-dipole
  selection rules": derives an expression for two-photon E1 selection rules for
  two photons of unequal frequency. Optica's full text is paywalled and the
  abstract page returned a login wall in this session, so I have this at
  search-summary level only. **UNVERIFIED as to its detailed content.**

I did **not** find a published Raman/two-photon Rabi rate for transitions
*within* ³Δ₁ in any source. **GAP.**

### 3.1 The operator

Adiabatic elimination of a far-detuned intermediate manifold {|i⟩} gives the
standard second-order effective operator between X-state levels |g⟩ and |f⟩:

```
T_eff = Σ_i  ( d·ε₂* ) |i⟩⟨i| ( d·ε₁ ) / Δ_i ,        Δ_i = E_i − E_g − ħω₁
```

with the amplitude A_{f←g} = ⟨f|T_eff|g⟩ and the line strength |A|². The
two-photon Rabi rate in the same normalisation is **Ω₁₃ = Ω₁₂Ω₂₃/(2Δ₀)**
(Cossel thesis Eq. 6.29 and surrounding text, p. 213). Cossel's Eq. (6.34) is
this expression with a **common** Δ pulled outside the sum, which is legitimate
exactly when the intermediate structure is unresolved (§3.4).

Writing each scalar product in spherical components,
`d·ε = Σ_p (−1)^p ε_{−p} T¹_p(d)`, the amplitude becomes a double sum over the
two lab-frame photon components p₁, p₂ weighted by the polarisation components
of ε₁ and ε₂ — i.e. **the same `weights`-dict machinery `heff.spectra` already
has, applied to a pair of indices rather than one**.

### 3.2 Decomposition into ranks K = 0, 1, 2

Two rank-1 operators acting on the *same* part of a coupled scheme (here the
rotational–electronic part, with I the spectator) couple into irreducible
tensors of rank K = 0, 1, 2 by **B&C Eq. (5.141)** (local extract
`bc_p195-210_ch5-matrixelements.txt`), verbatim:

```
T^K_p(A₁,B₁) = (−1)^{k₁−k₂+p} (2K+1)^{1/2} Σ_{p₁p₂} ( k₁ k₂ K ; p₁ p₂ −p )
                 × T^{k₁}_{p₁}(A₁) T^{k₂}_{p₂}(B₁)
```

and its reduced matrix element is **B&C Eq. (5.142)**, verbatim:

```
⟨η,j‖T^K(A₁,B₁)‖η′,j′⟩ = (2K+1)^{1/2} (−1)^{K+j+j′} Σ_{η″ j″}
      { k₁ k₂ K ; j′ j j″ } ⟨η,j‖T^{k₁}(A₁)‖η″,j″⟩ ⟨η″,j″‖T^{k₂}(B₁)‖η′,j′⟩
```

With k₁ = k₂ = 1 this is *literally* the two-photon operator: the sum over
`η″, j″` is the sum over intermediate states, and the 6j
`{1 1 K; j′ j j″}` is what a resolved-intermediate calculation would otherwise
carry numerically. Correspondingly, the **polarisation dyad** ε₁ ⊗ ε₂ decomposes
into the same three ranks: K = 0 (the scalar ε₁·ε₂), K = 1 (the antisymmetric
part, ∝ ε₁ × ε₂ — non-zero only for non-parallel or elliptical polarisations),
K = 2 (the symmetric traceless part). The full contraction is
`T_eff = Σ_K Σ_P (−1)^P (ε₁⊗ε₂)^K_{−P} α^K_P`, with α^K the effective rank-K
molecular two-photon operator. `[derived from B&C 5.141 + the standard
tensor-contraction identity B&C 5.167; the explicit contraction line is
algebra I did not check numerically — treat the *ranks* as cited and the
*assembled contraction* as derived]`

**Lab-frame Wigner–Eckart structure.** α^K_P is a rank-K lab-frame tensor, so
B&C Eq. (5.172) gives `⟨F′,m′_F|α^K_P|F,m_F⟩ = (−1)^{F′−m′_F} (F′ K F; −m′_F P m_F)
⟨F′‖α^K‖F⟩`, hence the selection rules **|ΔF| ≤ K** and **Δm_F = P, |P| ≤ K**.
With K ≤ 2 the two-photon operator reaches ΔF ∈ {0,±1,±2} and Δm_F ∈ {0,±1,±2},
and never Δm_F = ±3.

**Spectator theorem over I.** Both dipole factors act on the electronic–rotational
part only; I is a spectator in each. B&C Eq. (5.174) (operator on the first part
of the coupled scheme, the form `heff` already uses for the one-photon dipole)
therefore applies to the *composite* rank-K operator too, giving one 6j
`{J F I; F′ J′ K}` in place of the one-photon `{J F I; F′ J′ 1}`. This is exact
**only when the intermediate hyperfine structure is unresolved** — if Δ_i depends
on F′, the F′ sum cannot be factored out and the spectator reduction fails.
`[derived; each ingredient equation is cited, the composite is not numerically
checked]`

### 3.3 Molecule-frame structure and the ΔΩ selection rule

Each E1 step carries a molecule-frame component `q_i` fixed by the Ω's through
the 3j `(J′ 1 J; −Ω′ q Ω)`; since each photon is rank 1, `|q_i| ≤ 1`
(Cossel Eq. 6.37; Cairncross thesis p. 141: "we can simply set q = Ω′ − Ω").
Hence **ΔΩ = q₁ + q₂ ∈ {0, ±1, ±2}**, and, within X ³Δ₁ where the initial and
final |Ω| are both 1, the reachable ΔΩ is 0 or ±2 with

| intermediate Ω_i | ΔΩ reachable from Ω = +1 |
|---|---|
| **0** | **0 and −2** (i.e. Ω = +1 → −1 requires an Ω = 0 intermediate) |
| ±1 | 0 only |
| ±2 | 0 only |

**Verified numerically this session** with a scratch probe built on `heff.wigner`
and `heff.conventions` (script kept out of the repo, at
`…/scratchpad/twophoton_probe.py`), summing coherently over the intermediate
manifold for X J = 1,2, Ω = ±1, I = ½:

```
leg 1 q values (X -> Omega=0): [-1.0, 1.0]
leg 2 q values (Omega=0 -> X): [-1.0, 1.0]
Omega=0 intermediate: DeltaOmega present = [-2.0, 0.0, +2.0]   (every polarisation pair)
Omega=1 intermediate: DeltaOmega present = [0.0]               (both pairs tried)
```

The test is non-vacuous: the Ω = 1 case returns the *other* answer.

**Parity: the two-photon operator is parity-EVEN, so at zero field it does not
connect e to f.** Proof: P d P† = −d, and at zero field P commutes with H so the
parity partner |i′⟩ = P|i⟩ has the same Δ_i; therefore
`P T_eff P† = (−1)² Σ_{i′} d|i′⟩⟨i′|d / Δ_{i′} = T_eff`. `[derived]` Numerically,
the same probe gives `max|[T,P]| ≤ 5.6 × 10⁻¹⁷` for every polarisation pair
tried, with `max|T| ≈ 0.6`, using `heff.conventions.parity_operator`. So the
ΔΩ = ±2 channel is **not** an "Ω-doublet flip" in the parity sense: in the
signed-Ω basis it maps |+1⟩ ↔ |−1⟩, but on the parity eigenstates it acts as
`|e⟩ → +|e⟩`, `|f⟩ → −|f⟩` — the same structure as the Ω-doubling term itself
(§2.3 of `docs/thf-plus-x3delta1-effective-hamiltonian.md`), which is why `heff`'s
V8 gate finds parity commuting with it. **At non-zero E_rot parity is no longer
good** (V8: parity does not commute with the Stark term), and the restriction
lifts — which is the regime the JILA experiment actually runs in. The notebook
should say this rather than advertise a two-photon "parity flip".

**Δm_F reach under the JILA polarisation constraint.** Same probe:

```
(p1,p2)=(+1,-1)  Delta m_F = [0]
(p1,p2)=(+1,+1)  Delta m_F = [+2]
(p1,p2)=(-1,-1)  Delta m_F = [-2]
```

So with only σ± available in the E_rot plane (Ng thesis p. 102 footnote 4), the
reachable set is **Δm_F ∈ {0, ±2}**. Two consequences:
- Ng's stated target |J=1, F=3/2, m_F=+3/2⟩ → |…, m_F=+1/2⟩ is **Δm_F = −1** and
  therefore genuinely out of reach with σ± only, exactly as he says.
- But |m_F = +3/2⟩ → |m_F = −1/2⟩ is Δm_F = −2 and **is** reachable with a
  σ⁻σ⁻ pair. Since the CP-FIT proposal Ng analyses uses the m_F = ±1/2 pair
  symmetrically (Ng thesis pp. 100–104), a σ⁻σ⁻ two-photon step landing on
  m_F = −1/2 may serve the same purpose. `[derived this session; Ng does not
  discuss this and I am not claiming it is what JILA intends]`

### 3.4 Unresolved form (Δ ≫ intermediate B and hyperfine)

When the detuning is large compared with the intermediate rotational and
hyperfine spacings, Δ_i → Δ is common and comes out of the sum:

```
A_{f←g} = (1/Δ) Σ_i ⟨f|d·ε₂*|i⟩⟨i|d·ε₁|g⟩
```

and the closure over the intermediate manifold can be done analytically: with a
common denominator, the `η″, j″` sum in **B&C Eq. (5.142)** collapses to the
reduced element of a single rank-K operator, so

```
T_eff = (1/Δ) Σ_{K=0,1,2} Σ_P (−1)^P (ε₁⊗ε₂)^K_{−P} α^K_P
```

with α^K a **parameter-free geometry × one scalar per (K, ΔΩ)**. In the
molecule frame the surviving components are ΔΩ = 0 (two independent scalars,
"parallel" and "perpendicular" polarisability, from Ω_i = ±1, ±2 and the ΔΩ = 0
part of Ω_i = 0) and ΔΩ = ±2 (one scalar, from Ω_i = 0 only). This is the
Placzek-type polarisability picture; **the standard textbook treatment is
Long, "The Raman Effect" (2002), which I did not have access to — UNVERIFIED as a
citation, named only so the reader knows where the conventional derivation lives.**
Cossel's Eq. (6.34) is the operational version of this limit (his justification,
p. 219: "the single-photon detuning is large compared to the single-photon
transition width"), and Bray & Hochstrasser's closed-form diatomic cross sections
(abstract) are the ΔΩ = 0, ±1, ±2 rotational line strengths in the same limit.

**Caveat that decides whether this limit is usable for ThF⁺.** Cossel's HfF⁺
transfer ran at a one-photon detuning of 160 MHz, later increased to ≈ 1.5 GHz to
suppress spontaneous emission (Cossel thesis p. 218). The ThF⁺ intermediate
rotational constant is B ≈ 0.23 cm⁻¹ ≈ **6.9 GHz** (Gresh 2016 Table 2, e.g.
0.23187(5) cm⁻¹ for Ω = 0⁺[14.29]) `[derived, 1 cm⁻¹ = 29 979.2458 MHz]`. So at
a JILA-style detuning of ~1.5 GHz **the intermediate rotational structure is
resolved, not unresolved** — Δ ≪ B. The unresolved form is the wrong limit at
those detunings; it becomes right only for Δ ≫ ~7 GHz. The intermediate
*hyperfine* structure is a different story (X ³Δ₁ has A∥ = −20.1 MHz; no excited-
state A∥ is measured for ThF⁺ — **GAP**), so hyperfine-unresolved is plausible
while rotationally-unresolved is not.

### 3.5 Resolved form (sum over intermediate J′, F′, Ω_i)

Keeping every intermediate level explicit:

```
A_{f←g} = Σ_{Ω_i, J′, F′, m′_F}
     [ m1·m2·m3 ]_{f←i}(p₂, q₂) ε₂(p₂) · [ m1·m2·m3 ]_{i←g}(p₁, q₁) ε₁(p₁)
     / ( E_i(J′,F′) − E_g − ħω₁ )
```

summed over p₁, p₂ **before** squaring (Cossel Eq. 6.34), with the m1/m2/m3
factors of Cossel Eq. (6.37) and q_i = Ω_i − Ω. Two structural notes:
- **Ω_i must be summed with sign.** Leanhardt's original proposal routes the
  Raman through the ¹,³Π₁ states (Leanhardt 2011 = arXiv:1008.2997 p. 8:
  "a two-photon, stimulated Raman pulse, off-resonant from the intermediate
  ¹,³Π₁ states, which will coherently transfer population from the ¹Σ₀, J = 0
  ground state to the two |m_F| = 3/2 magnetic sublevels of the ³Δ₁, J = 1
  level"), and reaching both signs of final Ω requires both signs of Ω_i.
  `[derived]`
- **E_i must include the intermediate's own Ω-doubling if the two 0⁺/0⁻
  components are within Δ of each other**, because their contributions enter
  with opposite relative sign — the same σ_i = ±1 structure that produces ω_ef in
  the second-order picture (Petrov & Skripnikov arXiv:2503.02840; summarised in
  §2.3 of `docs/thf-plus-x3delta1-effective-hamiltonian.md`). Denis 2015 does not
  label 0⁺ vs 0⁻ for its Ω = 0 states (§2.2), so this sign is currently
  **unavailable for ThF⁺. GAP.**

---

## 4. What the package needs

### 4.1 Parameters

| parameter | meaning | best available value | status |
|---|---|---|---|
| `d_X_i` | reduced E1 dipole ⟨η_i‖T¹(d)‖X ³Δ₁⟩ per intermediate | ThF⁺: Denis 2015 Table 9 column "³Δ₁" — 0.455 / 0.571 / 0.391 D for the three Ω = 0 states. HfF⁺: Cossel thesis Table 5.6 — ³Π₀₋ 0.27, ³Π₀₊ 0.27, ³Σ₀₊ 0.08, ³Φ₂ 0.48 e·a₀, plus a second block at 0.36 / 0.39 / 0.17 e·a₀ | ThF⁺ ab initio only; HfF⁺ ³Π₀|d|³Δ₁ ≈ **0.27(3) e·a₀ measured** (Cossel p. 218) |
| `T_i` | intermediate term energy | ThF⁺ measured: Gresh 2016 Table 2 (§2.1). ThF⁺ ab initio: Denis Table 9, Petrov 2025 (mutually inconsistent, §2.2) | measured for the ≥10 000 cm⁻¹ states; contested below that |
| `Omega_i` | signed Ω of the intermediate, **and its 0⁺/0⁻ reflection parity** | Ω from Gresh/Denis labels | 0⁺/0⁻ missing for Denis's low Ω = 0 states — **GAP** |
| `B_i` | intermediate rotational constant | Gresh 2016 Table 2, B_e ≈ 0.229–0.236 cm⁻¹ for the measured excited states | measured |
| `A_par_i` | intermediate hyperfine | — | **GAP, never measured for a ThF⁺ excited state** |
| `Delta` | one-photon detuning | HfF⁺ practice: 160 MHz → 1.5 GHz (Cossel p. 218) | a user knob, not a molecular constant |
| `eps1`, `eps2` | polarisation Jones vectors | Cairncross thesis p. 141 represents laser fields exactly this way ("E = (X̂ + iŶ)/√2 represents a circularly polarised laser propagating along Ẑ") | user input |
| `theta` | relative linear-polarisation angle | Cossel Eq. 6.35, W_f ∝ sin²θ | user input |
| `Gamma_i` | intermediate linewidth, for the scattering-rate caveat | — | **GAP: no lifetime published for any ThF⁺ excited state in the sources read (§2.1)** |

The E1 legs' **overall scale is a single multiplicative constant** — as with
`heff.spectra`, the geometry can be returned dimensionless and the caller
supplies `d_X_i` products and `1/Δ`. That keeps the "no un-cited constants in
the kernel" property `heff` already has.

### 4.2 Two code shapes

**Option A — generalise `dipole_geometry` to arbitrary q, then take the
second-order product numerically.**

- `dipole_geometry(bra, ket, I, p)` currently returns 0 when `bra["Om"] != ket["Om"]`
  and hard-codes `q = 0` in `w3j(J2, 1, J1, -Om, 0, Om)`. The generalisation is
  `q = Om_bra - Om_ket`, guarded by `|q| ≤ 1`, with the existing three-factor
  body otherwise untouched — Cossel Eq. (6.37) is the citation, and it is the
  same equation `heff` already cites for q = 0.
- Then a `heff.spectra`-shaped function builds `D₂ @ diag(1/Δ_i) @ D₁`, sums over
  (p₁, p₂) with polarisation weights, and squares — the exact two-photon analogue
  of `_strengths_from_matrices`, and literally Cossel Eq. (6.34).
- *For*: no new tensor algebra; reuses the kernel, the blocking, and the
  sum-then-square convention (gate B8) unchanged; **exact in the resolved case**,
  which §3.4 shows is the physically relevant one for ThF⁺ at realistic
  detunings; interference between pathways (Cossel's cancelling σ pathways) comes
  out for free; one new public function plus a one-line kernel change. I built
  and ran precisely this shape as a scratch probe in <60 lines (§3.3), so the
  cost estimate is measured, not guessed.
- *Against*: needs an intermediate manifold — a second `StateSpec` with its own
  Ω_i and J range, its own energies, and (if you want it right) its own
  Hamiltonian. Introduces Δ as a knob that is not a molecular constant. Cost
  scales as dim_g × dim_i × dim_f. Does **not** fit the `H = Σ c_k M_k`
  term-matrix architecture, because the amplitude depends on Δ non-linearly.

**Option B — a direct rank-K effective operator within X only.**

- Register three parameter-free geometry matrices per (K, ΔΩ) channel — K = 0, 1, 2
  crossed with ΔΩ ∈ {0, ±2} — built from B&C Eqs. (5.141)/(5.142) with the
  spectator 6j `{J F I; F′ J′ K}`, and let the user supply the effective
  two-photon polarisabilities α^K_{ΔΩ} as parameters.
- *For*: fits `heff`'s architecture exactly — a rank-K operator's matrix is
  parameter-free geometry times a coefficient, so it is `M_k` in
  `H = Σ c_k M_k`, gets the `tensordot` fast path, gets Hellmann–Feynman exact
  derivatives, and gets the `@term` selection-rule gate for free. No intermediate
  basis, no Δ. The polarisation enters only as the K-dependent weights
  `(ε₁⊗ε₂)^K_P`, which is a clean and small user-facing API.
- *Against*: valid **only in the unresolved limit**, which §3.4 argues is *not*
  the ThF⁺ regime at 0.1–2 GHz detunings; the α^K_{ΔΩ} are then free parameters
  with no measured value for ThF⁺ (they would have to be generated by Option A,
  or fitted); and it hides exactly the intermediate-symmetry interference that
  Bray & Hochstrasser identify as the experimentally interesting content of the
  ΔΩ = 0 Q branch.

No decision taken here. The one observation worth carrying into that decision:
**A can generate B's parameters, B cannot generate A's spectra.**

---

## 5. Gaps, with the query I actually ran

Every row names the exact search so the negative is auditable rather than an
assertion.

| # | Gap | Query run | Result |
|---|---|---|---|
| 1 | Published Raman/two-photon **rate or line strength for a transition within X ³Δ₁** (either molecule) | PyMuPDF regex `raman\|two-photon\|two photon` over all of `ng-thesis-JILA.pdf`, `cairncross-thesis-JILA.pdf`, `gresh-thesis-JILA.pdf`, `cossel.pdf`, plus `grep -rniE "raman\|two-photon\|two photon\|2-photon\|stimulated"` over every `docs/lit/*.txt` | Nothing. Every hit is inter-electronic transfer, REMPD/REMPI, or fibre-optic Raman. Ng thesis p. 102 proposes it and says the spectroscopy is undone. **Positive control**: the same query shape found 12 real Raman hits in `cairncross-thesis-JILA.pdf` and the whole of Cossel §6.3, so the query is not silently mis-shaped |
| 2 | **0⁺ vs 0⁻ reflection parity** of Denis's low-lying Ω = 0 ThF⁺ states | Read Denis Table 9 as a rendered image (its PDF text layer is rotated and scrambled) and read §3.2.3 in full | Denis labels them ¹Σ₀(³Π₀), ³Π₀, ³Π₀(¹Σ₀) with no ± superscript. Needed for the relative sign of interfering ΔΩ = ±2 amplitudes (§3.5) |
| 3 | **Which ab initio ladder is right** below 10 000 cm⁻¹ | Compared Denis 2015 Table 9 (Ω = 0 at 6344/6528/6747 cm⁻¹) against Petrov & Skripnikov arXiv:2503.02840 (³Π₀₋ at 3044, ³Π₀₊ at 3395) against Gresh 2016 Table 2 (nothing measured between 3150 and 10 472) | Irreconcilable from the sources read. Any two-photon detuning quoted for ThF⁺ inherits this |
| 4 | **Excited-state lifetimes / linewidths** for ThF⁺ | Gresh 2016 Table 2 read in full (rendered); `grep -niE "lifetime"` over `ng2022*.txt` | Gresh Table 2 carries T₀, T_e, B_e, ω_e, ω_eχ_e, α_e only. Ng 2022's lifetimes are all for X ³Δ₁ v = 0/1 and a ¹Σ⁺, not for the optical intermediates |
| 5 | **Excited-state hyperfine A∥** for ThF⁺ | same sweep as #4 | Nothing. Blocks any resolved-hyperfine intermediate sum |
| 6 | **Bonin & McIlrath 1984** full text | `ctx_fetch_and_index` on `opg.optica.org/josab/abstract.cfm?uri=josab-1-1-52` and `ui.adsabs.harvard.edu/abs/1984JOSAB...1...52B/abstract` | Optica returned a login wall; ADS returned HTTP 405. Held at search-summary level and marked UNVERIFIED in §3.0 |
| 7 | **Bray & Hochstrasser 1976** full text | `ctx_fetch_and_index` on `tandfonline.com/doi/abs/10.1080/00268977600100931`, and on the 1999 follow-up `…/00268979909482813` | Abstract retrieved and quoted in §3.0; body paywalled. The 1999 "Interference in two-photon rotational line strengths of diatomic molecules" (Mol. Phys. 97, 1) returned empty content — **it is the single most on-target unread reference for the interference question** and is worth an institutional-access retrieval |
| 8 | **Zare, "Angular Momentum" ch. 5** | not attempted — no local copy and no open source | The B&C extracts (5.141, 5.142, 5.165–5.178) cover the same algebra and are local, so this is a redundancy gap, not a blocking one |
| 9 | A **measured** X ³Δ₁ ↔ Ω=0 transition dipole for **ThF⁺** | Denis 2015 §3.2.3 and Tables 9/10 read in full; Gresh 2016 Table 2 | Only ab initio. HfF⁺ has the one measured number (0.27(3) e·a₀, Cossel p. 218); ThF⁺ has none |

**Tool note for the record.** `pdf-mcp` and the `zotero` MCP server both failed to
connect this session (ConnectionRefused and HTTP 401 respectively), so every PDF
above was read with PyMuPDF (1.27.2.2) out of the `claude-code` conda env, and
`C:/Users/Arian/Zotero/storage` was not searched — it exists (5543 entries) and
may hold the two paywalled Molecular Physics papers.

**One file was downloaded**, Cossel's thesis, from
`https://www.colorado.edu/jila/media/696` (13.3 MB, linked from the public JILA
Cornell-group thesis list). It lives in the session scratchpad, not in the repo.

