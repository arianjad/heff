# The ThF⁺ X ³Δ₁ effective Hamiltonian in `heff`: form and values

One page. The derivations, citations and conventions are in
[the Hamiltonian reference](thf-plus-x3delta1-effective-hamiltonian.md); the Zeeman tensor's
sources are in [docs/lit/lookup-effective-zeeman-tensor.md](lit/lookup-effective-zeeman-tensor.md);
the limits of every input are in [open-questions.md](open-questions.md). Every value below
is the shipped default in `heff/models/thf_plus.toml`.

## Basis and conventions

- ²³²Th¹⁹F⁺, X ³Δ₁ (Λ = 2, S = 1, Ω = ±1), v = 0, Hund's case (c) kets |J, Ω, I_F = ½; F, m_F⟩,
  J = 1…4 by default (96 states), grouped into signed-m_F blocks for collinear fields.
- Parity: E*|J, Ω⟩ = (−1)^{J−1}|J, −Ω⟩. The e level has parity +(−1)^J (Brown 1975) and
  lies **below** f at every J.
- n̂ points from F to Th (JILA). Zeeman energies are E = −g μ_B B m_F. Δg ≡ g^u − g^l.
- All matrices are in MHz; E in V/cm, B in G; μ_B = 1.3996245 MHz/G, μ_N = μ_B/1836.15267;
  1 D = 0.5034118 MHz/(V/cm).

## The Hamiltonian

```
H = B₀ J(J+1) − D₀ [J(J+1)]²                                      rotation
  + H_Ωd                                                          Ω-doubling
  + A∥ [F(F+1) − I(I+1) − J(J+1)] / [2J(J+1)]  + (ΔJ = ±1 partner)  ¹⁹F axial hyperfine
  + c_I [F(F+1) − I(I+1) − J(J+1)] / 2                            ¹⁹F spin–rotation
  − d_mf n̂·E                                                      Stark
  + μ_B B·G·J − g_N μ_N I·B,   G = diag(G_xx, G_yy, G_zz) (molecule frame)   Zeeman
  − (d_e E_eff + W_TP k_TP) Ω/|Ω|                                 PT-odd (opt-in)
```

Term by term:

| term | operator and matrix element | rules |
|---|---|---|
| rotation | B₀ J(J+1); no −Ω² term (case (c), absorbed in the origin) | diagonal |
| centrifugal | −D₀ [J(J+1)]² | diagonal |
| Ω-doubling | ⟨J, −Ω\|H\|J, +Ω⟩ = +ω_ef J(J+1)/4, so the e/f splitting is ω_ef J(J+1)/2 with e lower | ΔΩ = ±2, ΔJ = 0 |
| ¹⁹F hyperfine | A∥ [F(F+1) − I(I+1) − J(J+1)]/[2J(J+1)] (B&C 9.50), plus the ΔJ = ±1 element of B&C 9.51 with the same A∥ | ΔJ = 0, ±1 |
| ¹⁹F spin–rotation | c_I [F(F+1) − I(I+1) − J(J+1)]/2 (B&C 8.7, 8.20) | diagonal |
| Stark | −d_mf n̂·E; polarised limit −m_F Ω γ_F d_mf E with γ_F = [J(J+1)+F(F+1)−I(I+1)]/[2F(F+1)J(J+1)] | ΔJ = 0, ±1; ΔΩ = 0; ΔF = 0, ±1; parity-odd |
| Zeeman, axial | G_zz μ_B (B·n̂)(J·n̂): Ω × the Stark geometry | ΔJ = 0, ±1; ΔΩ = 0; ΔF = 0, ±1 |
| Zeeman, perpendicular | G_xx μ_B B_xJ_x + G_yy μ_B B_yJ_y = G⊥ μ_B [B·J − (B·n̂)(J·n̂)] + G_Δ μ_B (B_xJ_x − B_yJ_y), with G⊥ = (G_xx+G_yy)/2 and G_Δ = (G_xx−G_yy)/2; the G_Δ piece is the ΔΩ = ±2 term of B&C (9.70)(vii) | ΔJ = 0, ±1; ΔΩ = 0, ±2; ΔF = 0, ±1 |
| nuclear Zeeman | −g_N μ_N I·B (B&C 9.59) | ΔJ = 0; ΔF = 0, ±1 |
| eEDM / scalar–pseudoscalar | −(d_e E_eff + W_TP k_TP) Ω/\|Ω\|; the observable is f^BD = 2 d_e E_eff | diagonal, parity-odd |

g-factors at zero field (the closed form the test suite checks):

```
g_J(e) = −G⊥ − (G_zz − G⊥)/[J(J+1)] + G_Δ
g_J(f) = −G⊥ − (G_zz − G⊥)/[J(J+1)] − G_Δ
g_F    = g_J [F(F+1) + J(J+1) − I(I+1)]/[2F(F+1)] + g_N (μ_N/μ_B) [F(F+1) − J(J+1) + I(I+1)]/[2F(F+1)]
```

giving g_F(e) = −0.015036, g_F(f) = −0.014762 at J = 1, F = 3/2: mean 0.0149 (Ng's
measurement, by construction) and Δg = g^u − g^l = +2.74 × 10⁻⁴ (Petrov 2025 prints +2.3 × 10⁻⁴).

## Values, ²³²Th¹⁹F⁺

| symbol | value | status | source |
|---|---|---|---|
| B₀ | 7274.3325(1) MHz | derived | Ng 2022, J = 1→2 interval 29.09733(4) GHz ÷ 4 |
| D₀ | 3.897(12) kHz | measured | Gresh 2016 Table 1 |
| ω_ef | 5.29(5) MHz | measured | Ng 2022 Table I |
| A∥ | −20.1(1) MHz | measured | Ng 2022 Table I (may absorb 40–120 kHz of c_I) |
| c_I | 20 kHz | estimate | CsF analogue scaled by B (factor ~3) |
| d_mf | 3.37(9) D, centre-of-mass origin | measured | Ng 2022 Table I |
| G_zz | 0.046532(2000) | derived | 0.04756 − G⊥; the measured sum G_zz + G⊥ fixes \|g_{F=3/2}\| = 0.0149 |
| G_xx | 8.2211 × 10⁻⁴ | estimate | G⊥ + G_Δ, Petrov & Skripnikov 2025 second-order sums |
| G_yy | 1.23327 × 10⁻³ | estimate | G⊥ − G_Δ, same; G⊥ = 1.02769 × 10⁻³, G_Δ = −2.0558 × 10⁻⁴ (~20 %) |
| g_N(¹⁹F) | 5.25773 | held fixed | nuclear g-factor |
| E_eff | 35.0(2.5) GV/cm | ab initio | JILA-adopted; Skripnikov & Titov 37.3, Denis 35.2 |
| W_TP | 50.0(3.5) kHz | ab initio | Skripnikov & Titov 2015 |
| d_e, k_TP | 0 | held fixed | set nonzero to switch the PT-odd block on |

## Odd-thorium isotopologues

For ²²⁹Th (I = 5/2) and ²²⁷Th (I = 1/2) the coupling is F₁ = J + I_Th, F = F₁ + I_F, the
¹⁹F terms above are recoupled as outer-spin operators, the Zeeman tensor and all ²³²ThF⁺
constants are transferred unchanged (status `estimate`), and the Th nucleus adds:

| symbol | ²²⁹Th | ²²⁷Th | status |
|---|---|---|---|
| A∥(Th) | −1510 MHz | +39821 MHz | ²²⁹: rescaled ab initio, sign chosen negative (OPEN-16); ²²⁷: Schmidt-model placeholder (OPEN-20) |
| c_I(Th) | 0 | 0 | no value anywhere (OPEN-19) |
| eQq₀(Th) | −2600 MHz | none (I = ½) | placeholder from an HfF⁺ anchor (OPEN-18) |
| eQq₂(Th) | +300 MHz | none | placeholder; normalisation bridge unresolved (OPEN-17) |
| g_N(Th) | 0.1464(24) | −3.826 | ²²⁹: μ = 0.366(6) μ_N; ²²⁷: Schmidt placeholder |

## What the model does not contain

Third-order Zeeman terms, hyperfine-dependent Ω-doubling (e_Δ), the anisotropic ¹⁹F
hyperfine constants, and the ³Δ₂ admixture that shifts δg/g at 60 V/cm by about 15 %.
See [open-questions.md](open-questions.md).
