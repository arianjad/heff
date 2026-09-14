# Effective Zeeman G tensor for ThF⁺ X ³Δ₁ — implementation brief (v2)

Status: implemented 2026-09-14, all four milestones landed, suite green. Supersedes the v1 brief
of the same date; the audit that corrected it is
`docs/superpowers/reports/2026-09-14-zeeman-tensor-handoff-audit.md` (read it first, it
has the derivation, the citations and the numerical check).
Closes `OPEN-10` of `docs/open-questions.md`.

## Decisions (per Arian, 2026-09-14)

1. **Doublet ordering follows Petrov's model: e (parity +(−1)^J) lies BELOW f at every J.**
   Ng never measured the ordering (Ng thesis Fig. 1.4 caption, PDF p.30: "We did not
   determine the energy ordering of the parity states"). Gresh 2016 Table 1 gives
   k″ > 0 on the Ω=0⁺ ← ³Δ₁ bands, which with Eq. (1) (s = +1 for e) puts e lower, if
   s labels the P/R-connected component. The current `omega_doubling` sign (e upper)
   must be flipped, with V4 and the docs updated to say why.
2. **Parameters are the Cartesian body-frame diagonal tensor `G_xx`, `G_yy`, `G_zz`**
   (x, y ⊥ n̂; z ∥ n̂). `G_par` is renamed `G_zz`. Comment the connection everywhere the
   values are set: `G_par ≡ G_zz`, `G_perp ≡ (G_xx + G_yy)/2`, `G_Delta ≡ (G_xx − G_yy)/2`.
3. `G_zz` keeps the one measured constraint: `G_zz = 0.04756 − G_perp`. `G_perp` and
   `G_Delta` come from Petrov's second-order sums and are ESTIMATES (say so in the notes).
4. No new gates beyond the existing suite plus one closed-form check. The second-order
   model that verified term C is documented in the audit report, not landed as a test.

## The operator

```
H_Z = μ_B B·G·J,   G = diag(G_xx, G_yy, G_zz)  (molecule frame)
B·G·J = G_zz (B·n̂)(J·n̂) + G_xx B_x J_x + G_yy B_y J_y
B_x J_x = ½ [P + C],   B_y J_y = ½ [P − C]
P ≡ B·J − (B·n̂)(J·n̂)           lab T¹₀(J) minus the axial piece       B&C (9.60), (9.57) PDF p.638
C ≡ B_x J_x − B_y J_y           ΔΩ = ±2, parity-dependent                B&C (9.70) term (vii), (9.71) PDF pp.652–653
```

Body components of J have anomalous commutation (B&C 5.152–5.153, PDF p.200), so C is
evaluated Brown–Howard style (B&C 5.155–5.162, pp.201–202): lab-frame T¹(J), closure over
the intermediate Ω. In body spherical components, with B ∥ ẑ_lab:

```
C = (1/√2) [ D¹*₀,₋₁ J^b₊  −  D¹*₀,₊₁ J^b₋ ] × (normalisation, below)
J^b_± |J,Ω⟩ = √[J(J+1) − Ω(Ω∓1)] |J,Ω∓1⟩   (anomalous: raising lowers Ω)
```

The relative minus is Hermiticity ((D¹*₀,₋₁)† = −D¹*₀,₊₁); the two factors in each product
commute, so ordering is immaterial.

### Matrix elements, one spin, basis |J, Ω, I, F, m_F⟩, B ∥ z (p = 0)

Every term is a lab rank-1 operator acting on J, so the spectator recoupling is the SAME
`a · b` as `elements_c.dipole_geometry` (B&C 5.174, 5.172), and only the J-space factor
(`dipole_geometry`'s `c` line) changes. Factor it out as `inner(J', Ω', J, Ω)`:

```
inner_axial(J',Ω',J,Ω) = δ_ΩΩ' (−1)^{J'−Ω} √[(2J+1)(2J'+1)] (J' 1 J; −Ω 0 Ω)         [today's c line]
inner_J(J',Ω',J,Ω)     = δ_JJ' δ_ΩΩ' √[J(J+1)(2J+1)]                                  [B&C 9.60]
inner_flip(J',Ω',J,Ω)  = −sign(q) √2 √[J(J+1) − Ω(Ω+s)] (−1)^{J'−Ω'} √[(2J+1)(2J'+1)] (J' 1 J; −Ω' s Ω+s)
                          with q = Ω'−Ω ∈ {+2, −2}, s = q/2; zero otherwise
```

Term matrices (all real, Hermitian, Δm_F = 0):

```
zeeman_Gzz : μ_B · zeeman_sign · Ω · [a b inner_axial]          rules dJ 0,±1; dOm 0; dF 0,±1   (= today's zeeman_Gpar)
zeeman_Gxx : μ_B · ½ · [ (a b inner_J − Ω a b inner_axial) + a b inner_flip ]   rules dJ 0,±1; dOm 0,±2; dF 0,±1
zeeman_Gyy : μ_B · ½ · [ (a b inner_J − Ω a b inner_axial) − a b inner_flip ]   same rules
```

Two-spin backend: `elements_c2.axial_geometry` lines 1–3 stay (k = 1); line 4 is the
`inner`. `zeeman_Gxx/Gyy` there declare dF1 0,±1 as `stark_z` does. At I_Th = 0 both
backends must agree to round-off (the existing V16 master reduction covers this once the
new names are in `V2_TO_V1`).

### Normalisation and sign (what `inner_flip`'s prefactor guarantees)

With `inner_flip` as written, ⟨J,∓1,M| C |J,±1,M⟩ = +M for every J (J-independent), and
under `conventions.parity_operator` (phase (−1)^{J−1} for ³Δ) the e level
(parity +(−1)^J) has ⟨e|C|e⟩ = −M, the f level +M. Hence with E = −g μ_B B M:

```
g_J(e) = −G_perp − (G_zz − G_perp)/[J(J+1)] + G_Delta
g_J(f) = −G_perp − (G_zz − G_perp)/[J(J+1)] − G_Delta,   G_Delta = (G_xx − G_yy)/2
g_F     = g_J · [F(F+1) + J(J+1) − I(I+1)] / [2F(F+1)]  + g_N (μ_N/μ_B) κ_F        (κ_F as in V5 today)
```

The e level is the one that mixes with the Ω = 0⁺ states (¹Σ⁺, ³Π₀₊), so it carries
Petrov's C₀₊ and G_Delta = (C₀₊ − C₀₋)/2 < 0: |g_e| > |g_f|. This is the convention to
comment next to the parameters: **G_Delta adds to the e-level g-factor.** With decision 1
(e lower) the prediction is Δg = g^u − g^l = g_f − g_e = −2G_Delta·(2/3) = +2.74e-4 at
J = 1, F = 3/2, same sign as Petrov's +2.3e-4.

## Parameters (all three isotopologue blocks of `heff/models/thf_plus.toml`, and `params.py`)

| symbol | value | status | note |
|---|---|---|---|
| `G_zz` | 0.046532 | derived | = 0.04756 − G_perp; keeps \|g_{F=3/2}\| = 0.0149 (Ng 2022) exact; ≡ G_par; uncertainty 0.002 as before |
| `G_xx` | 8.2211e-4 | estimate | = G_perp + G_Delta; G_perp = 1.02769e-3, G_Delta = −2.0558e-4 from Petrov & Skripnikov arXiv:2503.02840 Eqs. (3)–(15) second-order sums, `docs/lit/lookup-effective-zeeman-tensor.md` §2, §4 — an ESTIMATE (second-order PT, ~20 % on G_Delta) |
| `G_yy` | 1.23327e-3 | estimate | = G_perp − G_Delta; same source |

229/227 blocks: transferred copies with the existing transfer note pattern (add the three
names to the transfer list in `thf_v2`).

## Ordering flip

`elements_c.omega_doubling` returns `+J(J+1)/4` (was minus); the c2 delegate follows.
Rewrite its cite: e (parity +(−1)^J) below f at every J, per Petrov 2025's model (¹Σ⁺ at
314 cm⁻¹ pushes e down), consistent with Gresh 2016 Table 1 k″ > 0 on the 0⁺ bands; Ng's
schematics assume an ordering they did not measure. Splitting law unchanged. Update
`tests/test_elements_c_fieldfree.py::test_V4_...` (upper parity is −(−1)^J), the
`conventions.py` docstrings, `docs/open-questions.md` OPEN-2 and OQ-A, and
`docs/thf-plus-x3delta1-effective-hamiltonian.md` §2.3.

## Checks

- Existing suite, with `G_par` → `G_zz` renamed and V5's closed form replaced by the
  e/f closed form above (levels are now split by parity, so group by (J, F, parity)).
- One closed-form test at E = 0: `observe.g_factors` reproduces g_F(e/f) above for
  J = 1…4 at the shipped parameters, to 1e-10. Nothing else.

## Files

`heff/elements_c.py`, `heff/elements_c2.py`, `heff/params.py`, `heff/models/thf_plus.toml`,
`heff/backends/case_c.py` (unit table), `heff/conventions.py` (docstrings), `heff/terms.py`
(docstring), `heff/spec.py` (docstring), `tests/*` that name `G_par`/`zeeman_Gpar`,
`tests/test_elements_c2_reduction.py::V2_TO_V1`, `notebooks/build_tutorial.py:419`,
`docs/open-questions.md` (OPEN-2, OPEN-7, OPEN-10), `docs/thf-plus-x3delta1-effective-hamiltonian.md`
§2.3 and §2.8. Leave `results/*.json` (frozen outputs) and the `.ipynb` outputs alone
unless `build_tutorial.py` regenerates them.

## Working rules

Commit at each green milestone: (i) rename + ordering flip, suite green; (ii) `inner`
refactor in both geometries, suite green; (iii) `zeeman_Gxx/Gyy` + params + closed-form
check green; (iv) docs. Docs describe the code as it is now, no changelog prose. Do not
add gates. Stop and report rather than guess on any sign the closed form cannot pin.

## Progress

All four milestones landed; 372 passed, 2 skipped, plus the two opt-in
`HEFF_RUN_LITERATURE=1` comparisons.

| milestone | commit |
|---|---|
| (i) rename + ordering flip | `8b76964` |
| (ii) `inner` refactor | `55a9db1` |
| (iii) `zeeman_Gxx/Gyy`, parameters, closed-form check | `ad7289d` |
| (iv) docs | this commit |

**The sign check passed as predicted.** `⟨J,∓1,M|C|J,±1,M⟩ = +M` exactly and
J-independently for J = 1–5, Hermitian to zero. On the J = 1, F = 3/2 stretched block the
parity eigenvector of parity `+(−1)^J` — the e level — has `⟨C⟩ = −M`, so
`g_e = g_f + 2 G_Delta`, `G_Delta < 0`, `|g_e| > |g_f|`, and
`Δg = g^u − g^ℓ = +2.74 × 10⁻⁴` at J = 1, F = 3/2.

**Three things the brief did not anticipate.**

1. `tests/test_track.py`'s global-gauge check reconstructed the ground-reference component
   with plain `argmax`, not the tie-robust rule `apply_gauge` uses. At E = 0 the ground
   state is a parity eigenstate, so its two Ω components are equal in magnitude to
   rounding and the two rules pick different ones. The old ordering left both positive, so
   the mismatch was invisible; the flip exposed it. The test now uses the same rule.
2. `test_delta_g_collapses_when_the_stark_J_mixing_is_removed` needs a second knob. With
   ΔJ = ±1 Stark killed, δg/g no longer collapses to zero — `G_Delta` is a second source
   of a differential g-factor, leaving 4.8 × 10⁻⁴ at 60 V/cm. Setting `G_xx = G_yy` on top
   of that does collapse it, so the fabricated failure is still reachable.
3. The opt-in tier-D V7 comparison split in two. Ng's 32-level model is axial-only, so its
   −0.00223 is now compared against this package with `G_Delta` switched off (−0.00213).
   The full model gives −0.00261, which is 2.4 % from the measured −0.00255(6) — closer
   than the old axial model's 13 %.

**Left open.** `notebooks/ThF_plus_Isotopologues.ipynb` still carries the stored output
line "upper component is e at every J (OQ-A, closed)". Its builder
`notebooks/build_isotopologues.py` is corrected, but regenerating the notebook would
discard the stored outputs, so that was not done here. `OQ-A` has no entry in
`docs/open-questions.md`; the citations that pointed at it now point at `OPEN-2`.
