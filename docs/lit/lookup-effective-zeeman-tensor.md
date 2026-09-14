# The effective Zeeman G tensor for ²³²ThF⁺ X ³Δ₁: sources, derivation, and the `heff` operator

The physics behind `heff`'s three-parameter body-frame Zeeman tensor (`zeeman_Gzz`,
`zeeman_Gxx`, `zeeman_Gyy`; parameters `G_zz`, `G_xx`, `G_yy`), with every constant
traced to its source and the ΔΩ = ±2 operator derived and checked numerically. Source for
the resolution of `OPEN-7` and `OPEN-10` of `docs/open-questions.md`.

Tags: `[verbatim]` quoted from a source read at the time of writing, `[derived]` algebra
run here, `[inferred]` an argument, not a citation.

Local source texts (extracted from ar5iv; untracked per `.gitignore`, like every other
`docs/lit` source text; regenerate from the arXiv IDs):
- `docs/lit/petrov2025-arXiv2503.02840-ThFplus-gfactors.txt`
- `docs/lit/petrov2014-arXiv1404.4024-ThO-zeeman.txt`
- `docs/lit/petrov2017-arXiv1704.06631-HfFplus-zeeman.txt`

A typeset companion with a level diagram of the perturbers, the Van Vleck reduction, the
exact mapping to Petrov's g(J) form, and the ThF⁺/ThO/HfF⁺ checks:
[2026-09-14-thf-zeeman-second-order-summary.pdf](2026-09-14-thf-zeeman-second-order-summary.pdf). Its G∥ = 0.04680 is Petrov's
G₀ − C̄₀ + C₂/2; the package carries G_zz = 0.046532 from the measured sum (§5.5).

**Summary.** Through second order in the Coriolis × Zeeman cross terms with the nearby
Ω = 0±, 2 electronic states, the ³Δ₁ g-factor spectrum is exactly that of one
J-independent body-frame tensor G = diag(G_xx, G_yy, G_zz). Petrov & Skripnikov's published
ThF⁺ matrix elements fix G⊥ = (G_xx+G_yy)/2 and G_Δ = (G_xx−G_yy)/2; the single measured
number, |g(J=1, F=3/2)| = 0.0149(3), fixes only G_zz + G⊥. The e level (parity +(−1)^J)
carries the larger |g| and lies below f.

---

## 1. The second-order structure is Petrov 2014, verbatim

Petrov 2014 (ThO), Eqs. (12)–(13), p. 4:

> `[verbatim]` g_e(J) = −G∥/[J(J+1)] + G⊥⁽²⁾Δ⁽²⁾/[T_e(H³Δ₁) − T_e(A³Π₀₊)]
> + G⊥⁽¹⁾Δ⁽¹⁾/[T_e(H³Δ₁) − T_e(Q³Δ₂)] · (J+2)(J−1)/[2J(J+1)]   (12)
>
> g_f(J) = −G∥/[J(J+1)] + G⊥⁽²⁾Δ⁽²⁾/[T_e(H³Δ₁) − T_e(³Π₀₋)]
> + G⊥⁽¹⁾Δ⁽¹⁾/[T_e(H³Δ₁) − T_e(Q³Δ₂)] · (J+2)(J−1)/[2J(J+1)]   (13)

with G∥ = (1/Ω)⟨H³Δ₁|L̂ᵉ_n − g_S Ŝᵉ_n|H³Δ₁⟩ (Eq. 4) and g_S = −2.0023, so the magnetic
operator is L + |g_S|S. Write G₀ ≡ G∥, C₀± ≡ Σ_{Ω=0±} G⊥ₙΔₙ/(E_X − Eₙ),
C₂ ≡ Σ_{Ω=2} G⊥ₙΔₙ/(E_X − Eₙ). The generalisation from ThO (one state per Ω class) to
ThF⁺ (two at Ω=0⁺, one at Ω=0⁻, three at Ω=2) is a sum, nothing more `[inferred]`.

**e/f assignment.** Σ⁺ and 0⁺ states have only e levels, Σ⁻ and 0⁻ only f (Brown et al.
1975, Table I note, p. 502; B&C PDF p. 283), and the e/f label depends only on J and
parity, not on Λ, S, Ω or coupling case (Brown 1975, p. 501). So the ³Δ₁ e level
(parity +(−1)^J) is the one that mixes with the Ω = 0⁺ states, g_e carries C₀₊, and
G_Δ = (C₀₊ − C₀₋)/2 adds to the e-level g-factor. This is convention-free.

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

Through second order the ThO/ThF⁺ g-factor spectrum **is** the spectrum of a
J-independent rank-2 body-frame tensor, and a separate leading rotational g_R is
redundant (it is absorbed into G⊥). §5.2 shows the same holds for the full operator,
ΔJ = ±1 and ΔΩ = ±2 elements included, not only for the diagonal.

### 1.2 Which level is upper

The second-order algebra fixes which parity component gets which correction, not the
energy ordering. The ordering is not measured: Ng's thesis, Fig. 1.4 caption (PDF p. 30)
`[verbatim]`: "We assume that the positive parity state has the lower eigenenergy in
this schematic. We did not determine the energy ordering of the parity states in our
experiments." Ng 2022 Fig. 2 and thesis Fig. 4.2 are that schematic.

Gresh 2016 fits E = ν₀ + (B′ − s′k′/2)J′(J′+1) − … with s = +1 for e (Eq. 1, p. 3) and
finds k″ = −0.876…−0.90 × 10⁻⁴ cm⁻¹ on the Ω=0⁻ ← ³Δ₁ bands but +0.869…+0.892 × 10⁻⁴ on the
Ω=0⁺ ← ³Δ₁ bands (Table 1, p. 10); "the sign is changed, which is expected if the upper
state is of a different overall parity" (p. 5). A property of ³Δ₁ alone cannot flip with
the upper state, so the sign is label bookkeeping relative to it. The [10.47] upper state
is 0⁺ by its Q-branch-free combination band with ¹Σ⁺ (p. 5); on that band k″ > 0 puts e
lower if s = +1 marks the P/R-connected component `[inferred]`.

Petrov's model decides it: the ¹Σ⁺ state 314 cm⁻¹ above ³Δ₁ mixes with the e levels and
pushes them down, and the same mixing gives them the larger |g|. In the §5.2 model the
ratio (Zeeman Ω-flip element)/(Ω-doubling Ω-flip element) is −1.92 per cm⁻¹ for both
reflection signs, so upper − lower = g_f − g_e > 0, the sign Petrov prints. `heff` adopts
**e below f at every J** (per Arian, 2026-09-14): `omega_doubling` carries
+ω_ef J(J+1)/4 off-diagonal, and Δg = g^u − g^l = −2G_Δ · (2/3) = +2.74 × 10⁻⁴ at
J = 1, F = 3/2.

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
to both reproduces the sums below; the raw values give C₀₊ = −5.473×10⁻⁴ instead.

| contribution | value |
|---|---|
| ¹Σ⁺ (Ω=0⁺, ×1.088²) | −4.5239 × 10⁻⁴ |
| ³Π₀₊ (Ω=0⁺) | −1.6511 × 10⁻⁴ |
| **C₀₊** | **−6.17497 × 10⁻⁴** |
| **C₀₋** = ³Π₀₋ | **−2.06337 × 10⁻⁴** |
| 1³Δ₂ | −1.16835 × 10⁻³ |
| 2Ω=2 | −1.1678 × 10⁻⁶ |
| 3Ω=2 | −6.2029 × 10⁻⁵ |
| **C₂** | **−1.231543 × 10⁻³** |

### 2.1 G∥ = 0.047 is a fitted parameter, not an ab initio one

> `[verbatim]` (p. 3) "We also used G∥ = 0.047 to reproduce the g = 0.0149 value for the
> g-factor of the F = 3/2 state. Our value is slightly different from G∥ = 0.048 estimated
> in Ref. [Ng et al. 2022] since we account for nonadiabatic interactions between ³Δ₁ and
> other states in the basis set."

The independent theory content is C₀± and C₂, the nonadiabatic corrections. G₀ absorbs
whatever remains to land on the measured 0.0149, and it is printed to two figures:
±0.0005 in G∥ moves g(J=1, F=3/2) by ±1.1 %. Comparing a reconstructed J = 1 value to
0.0149 therefore tests nothing `[inferred]`; the one measured number constrains only
G∥ + G⊥ (§4.1).

---

## 3. Petrov's Δg as a function of field

> `[verbatim]` (§IV Results, p. 4) "The calculated difference Δg = gᵘ − gˡ = 2.3 × 10⁻⁴
> between the g-factors of the upper (gᵘ) and lower (gˡ) levels of the Ω doublets **for a
> zero electric field** is in agreement with the experimental value |Δg| = 3(3) × 10⁻⁴
> [Ng et al. 2022] (the sign of Δg was not determined in Ref. [Ng et al. 2022])."

Table I is indexed by electric field. Each row is a local linearisation, Eq. (18)
`[verbatim]` "Within a small area of the value of static electric field, the g factor
difference can be presented as Δg = Δg₀ + Δg₁E", with Δg₁ in cm/V: at E = 40 V/cm,
10⁷Δg₀ = 349.4 and 10⁷Δg₁ = 6.3; at E = 60 V/cm, 233.8 and 8.7. So Δg(60 V/cm) =
7.56 × 10⁻⁵ and Δg₀ is a per-field coefficient, not an E = 0 intercept. Δg falls steeply
off zero field to a minimum of about 7 × 10⁻⁵ near 25 V/cm (Fig. 2) then rises slowly,
the Stark-suppression shape Petrov 2017 reports for HfF⁺ (minimum Δg = 3 × 10⁻⁶ at
7 V/cm). Fig. 2's own E → 0 endpoints read g^u ≈ −1.478 × 10⁻², g^l ≈ −1.507 × 10⁻²
`[derived, digitized]`, a split of about 2.8 × 10⁻⁴ against the 2.3 × 10⁻⁴ in the text.

---

## 4. The second-order tensor and how far it misses

Substituting §2's sums (G₀ = 0.047, I(¹⁹F) = 1/2, g_N μ_N/μ_B = 2.86345 × 10⁻³):

```
G∥ = 0.0467961      G⊥ = 1.02769e-3      G_Δ = -2.05580e-4
G ≈ diag(8.221e-4, 1.2333e-3, 4.6796e-2)      G⁽²⁾₀ = 3.737e-2
```

| quantity | 2nd-order tensor | Petrov, full 7-state | measured (Ng 2022) |
|---|---|---|---|
| g(J=1, F=3/2) | −0.0149868 | −0.0149 (fit target) | 0.0149(3), sign undetermined |
| Δg(E=0), F=3/2 | 2.741 × 10⁻⁴ | 2.3 × 10⁻⁴ | \|Δg\| = 3(3) × 10⁻⁴ |

- **J = 1 g factor: 0.58 % from Petrov's fit target**, inside the two-figure rounding of
  his G∥ (§2.1). Not a measure of truncation.
- **Δg: 19 % high.** This is the number that probes the second-order truncation, and the
  measurement cannot adjudicate it (100 % error bar).
- **Independent magnitude check** `[derived]`: Leanhardt 2011 Eq. (65) gives
  |g′_rS| ≈ ω_ef/(2B_e) for the parity-dependent term. For ThF⁺,
  5.29 MHz / (2 × 0.243 cm⁻¹ × 29979.2458) = 3.63 × 10⁻⁴, against 2|G_Δ| = 4.11 × 10⁻⁴
  (rovibronic, before F projection), 13 % apart, same ¹Σ⁺ mixing mechanism.

### 4.1 ThF⁺ data cannot separate G∥ from G⊥

At J = 1 the tensor collapses: `g_J(1) = −(G∥ + G⊥)/2`. Ng's single measurement therefore
constrains only the sum, to ±1.9 %. Splitting the two needs J ≥ 2:

```
g_J(1) − g_J(2) = −(G∥ − G⊥)/3 = −0.015256      (64 % of g_J(1))
```

and no J ≥ 2 g-factor of ThF⁺ has been measured (Ng 2022 §II D: "it might be of interest
to compute" the rotational g-factor). ThO is where multi-J g factors exist, and where
Petrov 2014 validated the mechanism (his Table 2: pure case-(a) 1/[J(J+1)] scaling is
"badly violated", and the Ω = 0,2 mixing restores agreement). Fitting to measured shifts
gives one equation for three unknowns, so G⊥ and G_Δ come from §2 and the measurement
sets G∥ + G⊥.

---

## 5. The operator in `heff`

### 5.0 Ranks: lab versus molecule frame

A rotated tensor component carries one rank for both indices. B&C (5.143), PDF p. 199:
T^k_p(A) = Σ_q D^{(k)}_{pq}(ω)* T^k_q(A), and (5.146):

```
⟨J,Ω,M| D^{(k)}_{pq}(ω)* |J′,Ω′,M′⟩ = (−1)^{M−Ω} [(2J+1)(2J′+1)]^{1/2} (J k J′; −Ω q Ω′)(J k J′; −M p M′)
```

The same k sits in the body 3j and the lab 3j; of the reduced form (5.148) B&C say it
"is reduced as far as orientation in the space-fixed coordinate system is concerned, but
not for the molecule-fixed axes". `axial_geometry(k, …)` in `elements_c2.py` is exactly
this element (its first line is the lab 3j with k), so its k = 2 channel is a lab rank-2
operator: quadrupolar in m_F (the J=2, Ω-flip element at M = 1 and M = 2 has ratio −2)
with ΔJ, ΔF up to ±2. A Zeeman term is linear in B, lab rank 1, and is never that
element. Hirota (2.3.34), p. 37, is the identical transform.

Body components of J are not a spherical tensor: B&C (5.152)–(5.153), p. 200, [J_x, J_y]
= −iJ_z, and the would-be T_{q=±1}(J) "do not" satisfy the tensor definition (Hirota
(2.1.17–19), p. 17, the anomalous sign; Jadbabaie thesis p. 299, "For operators involving
anomalous commutation, such as J or N, we evaluate the operator in the lab frame").
Brown and Howard's method, B&C (5.155)–(5.162), pp. 201–202: keep T¹_p(J) in the lab
frame and insert closure over J″; the result (5.162) is
(−1)^{J−Ω}(J 1 J; −Ω q Ω′)[J(J+1)(2J+1)]^{1/2}δ_JJ′. The parity-dependent Zeeman term in
B&C's complete Zeeman Hamiltonian is (9.70) term (vii), p. 652,

```
−g_{e′r} μ_B B_Z Σ_{q=±1} Σ_p exp(−2iqφ) (−1)^p D^{(1)}_{−p,−q}(ω)* T¹_p(J−S) D^{(1)}_{0,−q}(ω)*
```

two rank-1 D matrices supplying the body index −q twice (ΔΩ = ±2), one lab index from
B_Z, the lab-frame T¹(J) contracted over p; its matrix element (9.71), p. 653, is a
closure sum over Ω″ with (J 1 J; −Ω −q Ω″)[J(J+1)(2J+1)]^{1/2}(J 1 J′; −Ω″ −q Ω′),
averaged over the two orderings. The S in J − S is absorbed into the effective constant
in a case-(c) Ω = ±1 manifold, as Leanhardt 2011 Eq. (64)'s `+½ g′_rS μ_B (B₊J₊S₊² +
B₋J₋S₋²)` absorbs it into g′_rS.

### 5.1 Decomposition and the three terms

`B·G·J` with G body-diagonal decomposes into three covariant pieces `[derived]`:

```
B·G·J = (G∥ − G⊥)(B·n̂)(J·n̂) + G⊥ (B·J) + G_Δ [(B·êx)(J·êx) − (B·êy)(J·êy)]
```

whose diagonal elements are Ω²m/[J(J+1)], m, and ±(the doublet term), i.e. exactly
`g(J) = −G⊥ − (G∥−G⊥)/[J(J+1)] ± G_Δ`. `heff` carries the same operator by Cartesian body
components: `zeeman_Gzz` = G_zz (B·n̂)(J·n̂), `zeeman_Gxx` = G_xx B_xJ_x = ½G_xx(P + C),
`zeeman_Gyy` = G_yy B_yJ_y = ½G_yy(P − C), with P = B·J − (B·n̂)(J·n̂) and
C = B_xJ_x − B_yJ_y; G∥ = G_zz, G⊥ = (G_xx+G_yy)/2, G_Δ = (G_xx−G_yy)/2.

| piece | operator | in `heff` |
|---|---|---|
| A | (J·n̂)(n̂·B), B&C (9.57) | `zeeman_Gzz`, coefficient `G_zz`; J-space factor = `dipole_geometry`'s default `inner` |
| P | B·J − (B·n̂)(J·n̂): the T¹₀(J) piece of B&C (9.60) minus A | the `inner_J − Ω·inner_axial` half of `zeeman_Gxx`/`zeeman_Gyy`; lab rank 1, ΔJ=0,±1, ΔΩ=0, ΔF=0,±1 |
| C | B&C (9.70) term (vii) / (9.71), §5.0 | the `inner_flip` half, entering `zeeman_Gxx` with + and `zeeman_Gyy` with −; lab rank 1, ΔJ=0,±1, ΔΩ=±2, ΔF=0,±1 |

Every piece is a lab rank-1 operator acting on J, so the nuclear-spin recoupling is the
same spectator chain (B&C 5.172, 5.174; two spectators in `axial_geometry`) for all
three, and only the J-space factor differs.

### 5.2 Term C, the ΔΩ = ±2 piece `[derived]`

In body spherical components:

```
G_Δ (B_x J_x − B_y J_y)  =  G_Δ (B^b_{+1} J^b_{+1} + B^b_{−1} J^b_{−1})
B^b_{±1} = D¹_{0,±1} B = −D¹*_{0,∓1} B        (body from lab, B&C 5.144; D¹_{0,±1} = −D¹*_{0,∓1})
J^b_{±1} = ∓ J^b_± /√2,   J^b_± |J,Ω⟩ = √[J(J+1) − Ω(Ω∓1)] |J, Ω∓1⟩   (anomalous)
⇒ C = (1/√2) [ D¹*_{0,−1} J^b_+  −  D¹*_{0,+1} J^b_− ]
```

The relative minus is forced by Hermiticity ((D¹*₀,₋₁)† = −D¹*₀,₊₁); with a relative plus
the operator is anti-Hermitian. The two factors in each product commute
([J^b_+, D¹*₀,₋₁] ∝ D¹*₀,₋₂ = 0), so ladder-on-ket and ladder-on-bra orderings agree.

Matrix elements in the signed-Ω primitive basis, no nuclear spin, B ∥ ẑ:

```
⟨J′, −1, M| C |J, +1, M⟩ = +√2 · √[J(J+1)] · ⟨J′, −1, M| D¹*_{0,−1} |J, 0, M⟩
⟨J′, +1, M| C |J, −1, M⟩ = −√2 · √[J(J+1)] · ⟨J′, +1, M| D¹*_{0,+1} |J, 0, M⟩
⟨J′Ω′M| D¹*_{0q} |JΩM⟩ = (−1)^{J′−M} (J′ 1 J; −M 0 M) (−1)^{J′−Ω′} √[(2J+1)(2J′+1)] (J′ 1 J; −Ω′ q Ω)
```

i.e. `dipole_geometry`'s J-space line evaluated with the intermediate Ω = 0 in the ket
slot and q = ∓1, times the ladder factor √[J(J+1)], times √2, with the channel sign. The
overall √2 makes ⟨J, ±| C |J, ±⟩ = ±M exactly and J-independently, so that
`g(J) = … ± G_Δ` with G_Δ as defined in §1.1. Under `conventions.parity_operator` (phase
(−1)^{J−1} for ³Δ) the e level gets ⟨e|C|e⟩ = −M, hence g_e = ḡ + G_Δ.

Hyperfine recoupling: C is a product of a lab rank-1 tensor (D¹*₀q) and a lab scalar
(J^b_±), so it is lab rank 1 and the k = 1 spectator chain applies to its reduced element
unchanged. Rules: ΔJ = 0, ±1; ΔΩ = ±2; ΔF = 0, ±1 (ΔF1 = 0, ±1 in the two-spin backend);
Δm_F = 0. F = 3/2 projection factor at J = 1 is 2/3.

### 5.3 Numerical check of the operator form

`docs/lit/2026-09-14-termC-second-order-model.py` (run with
`conda run -n structure python …`; documentation, not a test). Model: electronic states
X(Ω = ±1) and one intermediate n (Ω = 0 or ±2) at T_n above X; V = H_Cor + H_Z with
H_Cor = −B_rot(J^b₊J^e₋ + J^b₋J^e₊), H_Z = μ_B B Σ_q D¹*₀q G^b_q, Petrov's G⊥, Δ as the
electronic elements, a reflection sign s = ±1 relating the Ω < 0 side. Second order with
electronic denominators (Petrov's approximation): Z_eff = P H_Z P − Σ_n [P H_Cor Q_n H_Z P
+ h.c.]/T_n on the X manifold, J = 1…5, M = 1. Fit (G∥−G⊥, G⊥, G_Δ) by least squares over
all elements against {A = Ω D¹*₀₀, B = M·1, C as in §5.2}; report the max residual.

| intermediate | s | G∥−G⊥ | G⊥ | G_Δ | residual {A,B,C} | residual {A,B,k=2} | residual {A,B} |
|---|---|---|---|---|---|---|---|
| Ω=0, ¹Σ⁺ inputs | +1 | 0.047000 | 2.2619e-4 | +2.2619e-4 | **0.0** | 4.9e-4 | 5.5e-4 |
| Ω=0 | −1 | 0.047000 | 2.2619e-4 | −2.2619e-4 | **0.0** | 4.9e-4 | 5.5e-4 |
| Ω=2, 1³Δ₂ inputs | ±1 | 0.045832 | 5.8417e-4 | 0 | **0.0** | 0.0 | 0.0 |

Against §1.1: Ω=0 case, C₀ = G⊥Δ/(−T) = −4.5239e-4 → G⊥ = −C₀/2 = 2.2619e-4 ✓,
G∥ = G₀ − C₀/2 = 0.047226 ✓, G_Δ = C₀/2 ✓. Ω=2 case, C₂ = −1.16835e-3 → G⊥ = −C₂/2 =
5.8417e-4 ✓, G∥ − G⊥ = G₀ + C₂ = 0.045832 ✓. Diagonal g(J) in parity eigenstates equals
Eqs. (12)–(13) to 1e-9 for J = 1–3. The constant tensor with this C reproduces the
second-order operator exactly, off-diagonals included; the anti-Hermitian C gives residual
5.5e-4, identical to having no C, and the k = 2 element gives 4.9e-4.

Ordering (Ω=0, both s): the field-free second-order doubling from the same state gives a
J = 1 split of 7.06 MHz (¹Σ⁺ alone, before the 0± cancellation that lands on 5.29 MHz) and
g_upper − g_lower = +4.52e-4 in both cases: the level mixed with the 0⁺ state is lower and
carries the larger |g| (§1.2).

### 5.4 Parameters

Cartesian body-diagonal triple, in `heff/params.py` and all three isotopologue blocks of
`heff/models/thf_plus.toml`:

```
G_zz = 0.04756 − G⊥ = 0.046532        measured-sum anchor (|g(J=1,F=3/2)| = 0.0149 exact), status derived
G_xx = G⊥ + G_Δ = 8.2211e-4           estimate, from the §2 sums
G_yy = G⊥ − G_Δ = 1.23327e-3          estimate, from the §2 sums
```

G⊥ and G_Δ are second-order estimates (§4: about 20 % on G_Δ). Convention:
G_Δ = (G_xx − G_yy)/2 adds to the e-level g-factor; x is the transverse axis for which
that holds. The earlier axial-only fit G∥ = 0.04756 and Ng's 0.048 are the same J = 1
datum as Petrov's 0.047.

### 5.5 Closed form checked by the test suite

At E = 0, with κ_F = [F(F+1) − J(J+1) + I(I+1)]/[2F(F+1)]:

```
g_J(e) = −G⊥ − (G_zz − G⊥)/[J(J+1)] + G_Δ
g_J(f) = −G⊥ − (G_zz − G⊥)/[J(J+1)] − G_Δ
g_F    = g_J · [F(F+1) + J(J+1) − I(I+1)]/[2F(F+1)] + g_N (μ_N/μ_B) κ_F
```

### 5.6 Open

- **Third order.** Needs excited-to-excited Coriolis and magnetic matrix elements; Petrov
  publishes only ³Δ₁-to-X elements, so this is blocked on new ab initio input.
- **Seven-state benchmark.** Exact diagonalization of Petrov's published Hamiltonian for
  J = 1…5 is buildable from §2's table alone (add Ω=0±, 2 kets plus the Coriolis and
  electronic-Zeeman couplings) and is the only way to measure the truncation without a
  J ≥ 2 measurement. It would supersede the three-constant model rather than test it.
- **`D = −0.133 a.u.`** is what arXiv:2503.02840 p. 3 prints `[verbatim]`; it is 10× too
  small for ThF⁺ (−1.33 a.u. ≡ −3.38 D is the physical value), a typo in the paper.
