# The sign of A∥(²²⁹Th) in ThF⁺ X ³Δ₁ — convention audit

Resolves gap **G4** of `docs/digest-literature-th-hyperfine.md` §2.2 (Skripnikov &
Titov 2015 print A∥ = −4163 (μ_Th/μ_N) MHz; Denis et al. 2015 print A∥ = +1833 MHz;
magnitudes agree to 2.2 % after rescaling to a common μ, signs disagree).

Assembled 2026-09-05. Tags: `[verbatim]` quoted from a source read this session,
`[derived]` algebra run here, `[inferred]` an argument, not a citation.

**Headline.** Both groups **do** state their molecular-axis convention — not in the
passage that defines A∥, but elsewhere in the same paper. The axes are **opposite**
(Skripnikov: ζ from Th to F; Denis: z with F at negative z, i.e. F→Th). That
opposition does **not** explain the sign difference: A∥ as both groups define it is
invariant under a consistent axis reversal `[derived]`, and the two groups are shown
below to **agree** on the sign of the analogous constant in HfF⁺. The disagreement is
therefore **real**, not conventional, and the evidence favours **A∥(²²⁹Th) < 0**.

---

## 1. Each group's stated convention

### 1.1 Skripnikov & Titov, PRA 91, 042504 (2015) = arXiv:1503.01001

Local text: `docs/lit/skripnikov2015-arXiv1503.01001-ThFplus-theory.txt`.

The definition of A∥ is their Eq. (10), p. 2, and it carries no convention statement:

> `[verbatim]` "To obtain A∥ on ²²⁹Th in the ²²⁹ThF⁺ theoretically, one can evaluate
> the following matrix element:
> A∥ = (μ_Th / IΩ) ⟨Ψ| Σ_i (α_i × r_i / r_i³)_ζ |Ψ⟩,  (10)
> where μ_Th is magnetic moment of an isotope of ²²⁹Th nucleus having spin I. In the
> present paper we do not consider fluorine nuclear spin."

The convention for ζ and Ω is fixed **one page earlier**, in the sentence that
introduces Eq. (1)–(2) for W_d, and it governs every ζ-projected quantity in the
paper (Eqs. 1, 4, 9, 10, 11):

> `[verbatim]` (p. 2, immediately after Eq. 1) "…and Ω = ⟨Ψ|J · n|Ψ⟩, J is the total
> electronic momentum, n is the unit vector **along the molecular axis ζ directed
> from Th to F (Ω = 1 for the considered ³Δ₁ state of ThF⁺)**"

**So: n̂ = ζ = Th→F, Ω = +1.** This is the *opposite* axis direction to `heff`'s
`n_hat = F_to_Th` (README "Conventions").

Their Table II (p. 8) reports, for the ³Δ₁ state of ThF⁺, with signs carried
explicitly: d = 2.74 D (origin at the Th nucleus), E_eff = 37.3 GV/cm, W_T,P = 50 kHz,
W_M = 0.88 ×10³³ Hz/(e·cm²), **A∥ = −4163** (units μ_Th/μ_N · MHz), G∥ = 0.034; and for
ThO, **A∥ = −2949**. The minus sign is printed in the table, not an artefact.

### 1.2 Denis, Nørby, Jensen, Gomes, Nayak, Knecht & Fleig, NJP 17, 043005 (2015)

Local text: `docs/lit/denis2015-NJP-ThFplus-theory.txt`.

Their Eq. (2), §2.1, is the same operator with subscript `z` instead of `ζ`
(text extraction is mangled; the reconstructed form is
`A∥ = (μ_Th /(I Ω)) ⟨ψ| Σ_i (α_i × r_i / r_i³)_z |ψ⟩`), introduced by:

> `[verbatim]` "The parallel magnetic hyperfine interaction constant A∥ is defined as
> the z projection of the expectation value of the corresponding perturbative
> Hamiltonian in Dirac theory"

No convention statement accompanies it. The axis **is** pinned, in the caption of
their Table 10 (molecular dipole moments, p. 13):

> `[verbatim]` "…using the TZ basis set and the CI model. The origin is at the center
> of mass, and the internuclear distance is R = 3.779 a₀ **(F nucleus at z⃗ = z e_z
> with z < 0)**."

**So: z = F→Th, Ω = +1** — the same axis direction as `heff`'s `n_hat = F_to_Th`, and
opposite to Skripnikov's.

Signs they carry: their E_eff column header is literally **"−E_eff (GV cm⁻¹)"** with
entries 35.2 etc., i.e. E_eff = **−35.2** GV/cm (Tables 5–7, 11 and the abstract:
`[verbatim]` "we obtain an effective electric field of E_eff = −35.2 GV cm⁻¹"). Their
A∥ column header is plain "A∥ (MHz)" with entries **+1833**, and the abstract reads
`[verbatim]` "a magnetic hyperfine interaction constant of A∥ = 1833 MHz for ²²⁹Th
(I = 5/2)". W_P,T = +48.4 kHz.

### 1.3 The same group's later paper, where the axis is stated up front

Fleig (sole author; senior author of Denis 2015), *P,T-Odd and Magnetic Hyperfine
Interaction Constants and Excited-State Lifetime for HfF⁺*, arXiv:1706.02893,
submitted 9 Jun 2017. Retrieved 2026-09-05 via ar5iv. (The arXiv record carries no
`Journal ref:` field, so no journal citation is asserted here; cite the preprint.)

Definition, his Eq. (3), §II.2 — identical in structure to Denis Eq. (2), with the
1/(2 c I m_p Ω) prefactor written out:

> `[verbatim]` "A_{||}(K) = μ_K[μ_N]/(2 c I m_p Ω) ⟨Ψ_Ω| Σ_{i=1}^{n} (α⃗_i × r⃗_{iK} /
> r_{iK}³)_z |Ψ_Ω⟩"

Axis, from his Table 1 caption:

> `[verbatim]` "All values are given for the **Hf nucleus at the origin of the
> reference frame and the F nucleus at −3.4384 a.u. on the z axis.** The molecule-frame
> dipole moment corresponds to an origin located at the center of mass."

Same convention as Denis 2015 (heavy atom at origin, F at negative z ⇒ z = F→Hf), and
the same as `heff`'s `n_hat = F_to_Th`.

### 1.4 Why the axis opposition cannot by itself explain the sign

`[derived]`. Let ê = Th→F and let |a⟩ be the physical state with ⟨a|J·ê|a⟩ = +1.

- Skripnikov's convention: ζ = ê, Ω = +1 ⇒ the state is |a⟩, and
  A_S = (μ/I)(1/+1) ⟨a| V_ê |a⟩ with V ≡ Σ_i (α_i × r_i)/r_i³.
- Denis/Fleig's convention: z = −ê. The state with ⟨J·z⟩ = +1 is the mirror partner
  |b⟩ = σ_v|a⟩, which has ⟨b|J·ê|b⟩ = −1. So
  A_D = (μ/I)(1/+1) ⟨b| V_z |b⟩ = −(μ/I) ⟨b| V_ê |b⟩.
- V is an axial vector, so its axial component is odd under σ_v:
  ⟨b|V_ê|b⟩ = −⟨a|V_ê|a⟩.
- Hence A_D = +(μ/I)⟨a|V_ê|a⟩ = **A_S**.

The two flips cancel. The same conclusion follows from the effective-Hamiltonian side:
in Brown & Carrington Eq. (9.50) the constant is A∥ = {aΛ + (b_F + (2/3)c)Σ}Ω, which is
*quadratic* in the axis reversal (Λ, Σ and Ω all flip together) and is by construction
the coefficient of the observable pattern [F(F+1) − J(J+1) − I(I+1)]/[2J(J+1)].
**A∥ is convention-independent.** So is E_eff up to the separate and well-known
d_e-sign fork, which is why the E_eff signs (+37.3 vs −35.2) are *not* usable as a
cross-check and are set aside below.

---

## 2. Calibration cases

The test that discriminates "convention offset between the two groups" from "real
disagreement": find a system where both groups compute the same A∥, and/or where a
computed A∥ can be compared with a measured one. HfF⁺ X ³Δ₁ supplies both. It is the
right anchor — same ³Δ₁ term, same (n)d¹(n+1)s¹ heavy-metal configuration, same ¹⁹F
partner.

| # | system / nucleus | who | computed A∥ | measured A∥ | verdict |
|---|---|---|---|---|---|
| C1 | ¹⁸⁰Hf¹⁹F⁺, **¹⁹F** | Fleig 2017 (Denis-group convention, z = F→Hf) | **−43.0 MHz** | **−62.0(2) MHz**, Cairncross et al. 2017 | **signs agree** (magnitude 31 % low, and he says so) |
| C2 | ¹⁷⁷Hf¹⁹F⁺, **¹⁷⁷Hf** (g_N = +0.2267) | Fleig 2017 | **−1328 MHz** | none | — |
| C3 | ¹⁷⁷Hf¹⁹F⁺, **¹⁷⁷Hf** | Skripnikov group (Ref. 15 of Petrov 2018) | **−1429 MHz** | none | agrees with C2 in sign, 7 % in magnitude |
| C4 | ¹⁷⁹Hf¹⁹F⁺, **¹⁷⁹Hf** (g_N = −0.1424) | Skripnikov group | **+898 MHz** | none | sign tracks g_N, as it must |
| C5 | ¹⁸⁰Hf¹⁹F⁺, **¹⁹F** | Petrov, Skripnikov & Titov 2023 (arXiv:2302.02856) adopt the measured value into their spin-rotational Hamiltonian | — | **−62.0 MHz** | their H_hfs takes A∥ negative |
| C6 | ²³²Th¹⁹F⁺, **¹⁹F** | — | — | **−20.1(1) MHz**, Ng et al. 2022 Table I | the value `heff` already carries |
| C7 | ²²⁹ThO H ³Δ₁, **²²⁹Th** | Skripnikov & Titov 2015 Table II | **−2949** (μ_Th/μ_N) MHz | none | same sign as their ThF⁺ |

Sources for the calibration rows:

- **C1, C2** — Fleig 2017 Table 1, row `vTZ/MR12-CISD(20)`: D = 4.19 D, E_eff = −22.7
  GV/cm, W_S = 20.0 kHz, G∥ = 0.0127, τ = 2.7 s, **A∥(¹⁷⁷Hf) = −1328 MHz**,
  **A∥(¹⁹F) = −43.0 MHz**. And §III.2, verbatim:
  > `[verbatim]` "The magnetic hyperfine interaction constant for the ¹⁹F nucleus
  > agrees qualitatively with the value of −62.0(2) [MHz] measured in reference
  > [Cairncross:2017fip]. However, the present electronic-structure model has not been
  > designed with a focus on properties depending on spin density in the vicinity of
  > the fluorine nucleus. On the other hand, the calculated hyperfine interaction
  > constant for the ¹⁷⁷Hf nucleus allows for an assessment of the accuracy of the
  > molecular wavefunction … At present, the author is not aware of a measurement of
  > A∥(¹⁷⁷Hf) in the Ω = 1 state of HfF⁺ to compare with."

  This is the decisive row. Fleig computes A∥ from the *same* Eq. (2)/(3) definition
  used in Denis 2015, in the *same* axis convention as Denis 2015, gets a **negative**
  number for ¹⁹F, and compares it directly with the **negative** measured value as an
  agreement. His pipeline therefore reproduces measured A∥ signs with no relabelling.

- **C3, C4** — Petrov, Skripnikov, Titov & Flambaum, PRA 98, 042502 (2018), §III,
  local `docs/lit/petrov2018-CP-violation-HfFplus.txt`:
  > `[verbatim]` "For the ¹⁷⁷Hf and ¹⁷⁹Hf isotopes I₁ = 7/2, g_Hf = 0.2267, Q = 3.365 b
  > and I₁ = 9/2, g_Hf = −0.1424, Q = 3.793 b, respectively. The magnetic dipole
  > hyperfine structure constant A^F∥ = −62.0 MHz was measured in Ref. [1]. The
  > magnetic dipole hyperfine structure constants A^Hf∥ = −1429 MHz and A^Hf∥ = 898 MHz
  > for ¹⁷⁷Hf¹⁹F⁺ and ¹⁷⁹Hf¹⁹F⁺, respectively, were calculated in Ref. [15]."

  Their Fig. 1 plots MQM shifts against A^Hf∥ with a vertical line at −1429 MHz, so the
  sign is load-bearing in their own numerics, not decorative.

- **C5** — `docs/lit/2023-arXiv2302.02856-revisited-PT-odd-HfFplus.txt`:
  > `[verbatim]` "…except for the hyperfine structure constant A∥ = −62.0 MHz measured
  > in Ref. [9]"

- **C6** — Ng et al., PRA 105, 022823 (2022) Table I; already the `A_par` entry in
  `heff/params.py` (−20.1 MHz, `status="measured"`).

**Result of the calibration.** On the one system where both groups compute the same
constant (C2 vs C3), they agree: **−1328 vs −1429 MHz**, same sign, 7 % apart. There
is no systematic sign offset between the Skripnikov/Petrov convention and the
Fleig/Denis convention. And the Fleig/Denis convention is anchored to experiment by
C1. The premise of gap G4 — that the ThF⁺ discrepancy might be a convention artefact —
**fails its own test**.

---

## 3. Which sign is physically right

Four independent lines, all pointing negative.

**(a) The ab initio decomposition, stated by Skripnikov himself.** For the ThF⁺ ³Δ₁
state he reports the W-reduced density-matrix differences (his Eq. 33, p. 9):

> `[verbatim]` "Δ(e7s½,½ , e7s½,½) − Δ(e7s½,−½ , e7s½,−½) = −0.99
>  Δ(f7p½,½ , f7p½,½) − Δ(f7p½,−½ , f7p½,−½) = −0.47
>  Δ(f6d3/2,3/2 , f6d3/2,3/2) − Δ(f6d3/2,−3/2 , f6d3/2,−3/2) = 0.88"

and

> `[verbatim]` "The leading matrix element of the HFS operator in the basis of
> reference functions is between the e7s functions."

Read that off: the 6d₃/₂ electron sits at m = +3/2 (positive difference) and the 7s₁/₂
electron at m = −1/2 (negative difference), summing to Ω = +1. The dominant hyperfine
contributor is the 7s electron, and it is *anti*-aligned with Ω. A positive nuclear
moment then gives a negative A∥. `[inferred, from the source's own numbers]`

**(b) The same statement in effective-Hamiltonian language.** For ³Δ₁, Λ = +2, S = 1,
Ω = +1 ⇒ **Σ = −1**. In B&C Eq. (9.50), A∥ = {aΛ + (b_F + (2/3)c)Σ}Ω = 2a −
(b_F + (2/3)c). For a σ(7s) electron with μ_Th > 0 the Fermi-contact constant b_F is
large and positive, and it enters with the minus sign because the spin is anti-aligned
in the Ω = 1 component of a ³Δ. It dominates the δ-orbital term 2a. Hence A∥ < 0.
`[inferred]` This is the argument §2.2 of the digest sketched and declined to code; it
is here corroborated by (a), which is a source statement rather than an argument.

**(c) Isovalent consistency.** ThF⁺ (Th²⁺ 7s¹6d¹), ThO (same effective Th state) and
HfF⁺ (Hf²⁺ 6s¹5d¹) are the same ³Δ₁ σδ problem. Every published heavy-nucleus A∥ for
this family with a positive g_N is negative: ThF⁺ −4163, ThO −2949 (Skripnikov &
Titov 2015 Table II), ¹⁷⁷HfF⁺ −1429 (Skripnikov group) and −1328 (Fleig 2017); and the
one case with a negative g_N, ¹⁷⁹HfF⁺, flips to +898, exactly as A∥ ∝ g_N requires.
Denis 2015's +1833 for ²²⁹ThF⁺ (g_N > 0) is the sole outlier in the family, including
against a paper from its own senior author.

**(d) The ¹⁹F constants.** Measured A∥(¹⁹F) is negative in both ions (−62.0 MHz in
HfF⁺, −20.1 MHz in ThF⁺) with μ(¹⁹F) = +2.628 μ_N, the same Σ = −1 mechanism. Not a
proof for the heavy nucleus, but it removes any suspicion that "negative A∥ with
positive μ" is anomalous in these states.

**On the ²²⁹Th moment.** μ(²²⁹Th) = +0.366(6) μ_N (Porsev, Safronova & Kozlov 2021),
positive — digest §1.2. So a negative A∥ means a negative *electronic* factor, which is
what (a)–(c) establish. Nothing in the Bohr–Weisskopf correction changes a sign; it is
a few-percent reduction of the magnitude.

---

## 4. Conclusion and the recommended sign for `heff`

**The disagreement is real, not conventional.**

1. Both conventions are now known: Skripnikov ζ = Th→F, Ω = +1; Denis/Fleig z = F→Th,
   Ω = +1. They are opposite.
2. That opposition is *not* a sign difference in A∥: the quantity is invariant under a
   consistent axis reversal (§1.4, two independent derivations).
3. The calibration confirms this empirically — the two groups agree, sign and 7 % in
   magnitude, on A∥(¹⁷⁷Hf) in HfF⁺ (§2, C2 vs C3) — and the Denis/Fleig convention is
   independently anchored to a measured sign (C1).
4. Therefore Denis et al. 2015's **+1833 MHz is an outlier**: either a sign slip in
   that paper, or a magnitude quoted without its sign. I cannot distinguish those two
   from the text; their abstract and every table entry read as a genuine `+`, and they
   *did* carry the minus sign for E_eff (via a "−E_eff" column header), which argues
   against a blanket magnitude convention.
5. The physics (§3) independently requires A∥ < 0 for a ³Δ₁ σδ state with μ > 0.

### Recommended value

**A∥(²²⁹Th, ThF⁺ X ³Δ₁) = −1.51(11) GHz**, i.e. **negative**. The magnitude is the
digest §2.2 rescaling to μ = 0.366(6) μ_N, mean of the two calculations
(−1524 / −1491 MHz), with the spread and the authors' ~7 % ab initio uncertainty in
quadrature. Status: `[ab initio]`, sign now `[resolved]` rather than `UNVERIFIED`.

### It does not depend on `heff`'s `n_hat`

`heff/elements_c.py::hyperfine_A_par` implements B&C (9.50) as the bare Casimir
kernel — verbatim from the source:

```python
J, F, I = float(ket["J"]), float(ket["F"]), float(ctx.I)
return (F * (F + 1.0) - I * (I + 1.0) - J * (J + 1.0)) / (2.0 * J * (J + 1.0))
```

Neither Ω nor n̂ appears. The `A_par` parameter multiplying it is therefore the
convention-independent observable coefficient — the same object as the measured
−20.1 MHz already in `heff/params.py`. So:

- **Enter A_par(²²⁹Th) = −1510 MHz** (or −1524 if you prefer Skripnikov's value alone,
  rescaled), `status="ab initio"`, source Skripnikov & Titov 2015 Table II rescaled to
  μ = 0.366 μ_N, with a note that Denis et al. 2015 print the opposite sign and that
  this document resolves it.
- **Switching `n_hat` between `F_to_Th` and `Th_to_F` must not change it.** That is
  worth one test: `A_par` is invariant under the convention flip, unlike `d_mf` and the
  Zeeman/PT-odd terms, whose signs are convention-carrying. `hyperfine_A_par_dJ1` does
  divide by the signed Ω — correct, since B&C's brace is A∥/Ω — so the Ω = +1 and
  Ω = −1 blocks get opposite-signed ΔJ = ±1 elements; that is physics (they are parity
  partners), not a convention leak, and it is unaffected by the choice of n̂.

### Residual uncertainty and the single check that would close it

The chain is: (i) A∥ is convention-invariant `[derived]`; (ii) the groups agree on
HfF⁺ `[cited]`; (iii) Fleig's convention reproduces a measured sign `[cited]`;
(iv) the σ(7s) anti-alignment forces A∥ < 0 `[inferred from Skripnikov's own density
matrices]`. Step (iv) is the only inference, and it is corroborated, not load-bearing
alone.

What is *not* settled: **why** Denis et al. print `+`. I did not find any later paper
by that group restating the ThF⁺ A∥, nor an erratum. Searches run this session:
arXiv API `au:"Fleig_T" AND abs:"hyperfine"` → 0 results (the 1706.02893 paper is
indexed under the `abs:"HfF"` query, not that author-string form);
`abs:"HfF" AND abs:"hyperfine"` → 12 results, the relevant ones being 1706.02893,
1704.06631, 2308.12832; `abs:"ThF" AND abs:"hyperfine"` → 3 results
(2503.02840, 2202.01346, 1503.01001) — **no post-2015 Fleig-group ThF⁺ hyperfine
paper exists**. Not consulted: the published NJP version's supplementary material;
pdf-mcp and the Zotero MCP server both failed to connect this session, so no
PDF-native text search was possible.

**The single check that would close it definitively:** email M. Denis / T. Fleig (or
S. Knecht) asking whether the 1833 MHz in NJP 17, 043005 Table 6 is signed, and if so
how it is reconciled with A∥(¹⁷⁷Hf) = −1328 MHz in Fleig, arXiv:1706.02893 Table 1,
computed from the same Eq. (2)/(3) in the same axis frame. One sentence of
reply settles it. Absent that reply, the recommendation above stands on the
calibration and the physics, and I would publish ²²⁹ThF⁺ level orderings with it while
flagging the one-sentence provenance.
