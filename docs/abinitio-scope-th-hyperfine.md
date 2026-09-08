# Ab initio scope: what calculation would give us the missing Th constants in ThF⁺

Companion to `docs/digest-literature-th-hyperfine.md` ([DIGEST]) and
`docs/thf-plus-x3delta1-effective-hamiltonian.md` ([HAM]). [DIGEST] §5 lists five
gaps with no published value anywhere: `eQq₀(Th)`, `eQq₂(Th)`, `c_I(Th)`,
`c_I(¹⁹F)`, and the *sign* of `A∥(Th)`. This 2026-09-05 source review scopes
calculations to replace the exploratory parameter estimates; consult
[current limitations](open-questions.md) for the implemented model's status.

This document answers: what each constant is as a property operator, which code
computes it and how, what it would cost, who already has the pipeline, and what
to actually do.

**Tags.** `[measured]`, `[ab initio]`, `[derived]` (arithmetic/algebra run here
from cited inputs), `[inference]` (my reasoning, not in any source),
`UNVERIFIED` (I could not check it this session).

**Not consulted this session** — declare so absence is not read as evidence:
pdf-mcp and the zotero MCP server both failed to connect. **`diracprogram.org`
was down** (returned a French "Site en maintenance" page, retrieved 2026-09-05),
so every DIRAC manual quotation below comes from the manual *source* in the
DIRAC git repository, `https://gitlab.com/dirac/dirac/-/raw/master/doc/manual/`
(`properties.rst`, `hamiltonian.rst`, `wave_function.rst`,
`one_electron_operators.rst`), retrieved 2026-09-05 — that is the master branch,
not a tagged release, so keyword availability in a specific installed DIRAC
version must be re-checked. The CFOUR manual pages for `PROP` were 404 /
unreachable; the CFOUR EFG statement below is correspondingly weaker.

---

## 0. Executive summary

| constant | operator | best route | status of that route |
|---|---|---|---|
| `A∥(Th)`, `A∥(F)` | electron–nucleus magnetic dipole (relativistic, one operator) | X2CAMF-CCSD(T) in CFOUR, or 4c-DC-CCSD(T) in DIRAC, or Skripnikov's 2-step GRECP | **routine** — done for ThF⁺ twice already |
| `eQq₀(Th)` | EFG at Th, `q₀ = φ⁽²⁾_zz/e` | 4c/X2C-CCSD(T) EFG expectation value, DIRAC `.EFG`/`.NQCC` or CFOUR | **routine** — exactly what Petrov 2018 did for HfF⁺ |
| `eQq₂(Th)` | *off-diagonal* rank-2 element between Ω = +1 and Ω = −1 | not a keyword anywhere; needs KRCI/MRCI with an off-diagonal property, or finite-field symmetry breaking | **research project** (§5) |
| `c_I(Th)`, `c_I(F)` | nuclear spin-rotation tensor, `⟨⟨ĥ^hfs; ĥ^BO⟩⟩` linear response | DIRAC `**PROPERTIES .SPIN-ROTATION` | keyword exists; **open-shell Ω = 1 applicability UNVERIFIED** |
| sign of `A∥(Th)` | — | compute `A∥(Th)` and `A∥(¹⁹F)` in the *same* run and anchor on the measured F value | **free rider** on the `A∥` calculation (§4d) |

**One calculation, if only one:** a single X2CAMF-CCSD(T) (CFOUR) or
4c-DC-CCSD(T) (DIRAC) run at R(Th–F) = 3.75 a₀ that reports, together, in one
stated axis convention: `A∥(²²⁹Th)`, `A∥(¹⁹F)`, the EFG `q₀` at Th in atomic
units, and `G∥`. Three of the four are checkable against measurement or against
prior theory, and it kills gaps G2 (partly) and G4 outright.

---

## 1. The constants as property operators

### 1.1 Magnetic hyperfine `A∥` — the electron–nuclear spin–orbit operator

Both ThF⁺ papers use the *same* defining expression, the axial projection of the
relativistic hyperfine operator divided by IΩ:

- **Skripnikov & Titov**, PRA 91, 042504 (2015) = arXiv:1503.01001, their
  Eq. (10) (verbatim from `docs/lit/skripnikov2015-…txt`):
  `A∥ = (μ_Th/(IΩ)) ⟨Ψ| Σ_i (α_i × r_i / r_i³)_ζ |Ψ⟩`
- **Denis, Nørby, Jensen, Gomes, Nayak, Knecht, Fleig**, NJP 17, 043005 (2015),
  their Eq. (2): the identical expression, described as "the z projection of the
  expectation value of the corresponding perturbative Hamiltonian in Dirac
  theory".
- **Petrov et al.**, PRA 98, 042502 (2018) = arXiv:1809.06701, Eqs. (20)–(21),
  the same operator written with the nuclear g-factor pulled out:
  `A^Hf∥ = g_Hf ⟨³Δ₁| Σ_i (α_i × r_1i/r_1i³)_ζ |³Δ₁⟩`.

In a four- or two-component framework this is **one** operator — there is no
separate Fermi-contact / dipolar / orbital decomposition to assemble; those are
the nonrelativistic limits of `α × r/r³`. DIRAC writes the same thing as
`ĥ^hfs_K = −Σ_i m_K · B̂^el_K(i)`, with
`B̂^el_K(i) = −(1/4πε₀c²)(r_iK × ecα)/r_iK³`
(DIRAC manual `properties.rst`, `.NMR` and `.SPIN-ROTATION` entries, retrieved
2026-09-05). The Fermi-contact piece is separately available as integrals with
label `'FC NAMab'` for finite-field work (`one_electron_operators.rst`).

**Two nuclear-side corrections that matter at the few-percent level**, both of
which the digest's ±7 % error bar silently contains:

1. **Bohr–Weisskopf (finite magnetization distribution).** Skripnikov, JCP 153,
   114114 (2020) = arXiv:2008.01520, abstract: "accurate prediction of the
   hyperfine structure of the RaF molecule requires to take into account the
   finite magnetization distribution inside the radium nucleus", and its
   magnitude "depends on the model of the nuclear magnetization distribution
   which is usually not well known". Th (Z = 90) is heavier than Ra (Z = 88), so
   this is not negligible for ²²⁹ThF⁺. Porsev, Safronova & Kozlov
   (arXiv:2107.14723) already flagged that the 2013 μ(²²⁹Th) error bar "did not
   include uncertainty due to the magnetization distribution" ([DIGEST] §1.2).
2. **Which μ.** `A∥ ∝ g_N = μ/I` exactly ([DIGEST] §2.2, and stated by Petrov
   2018 §III: "The ratios … correspond to the ratios for the nuclear g-factors
   and the quadrupole moments"). So the calculation should report the
   **isotope-independent electronic factor**, not `A∥` for one assumed μ.

**Accuracy achieved on analogous systems, against experiment:**

| system | method | computed | measured | source |
|---|---|---|---|---|
| ThF⁺ X³Δ₁, `A∥(¹⁹F)` | X2CAMF-CCSD(T)/ANO-RCC-unc, CFOUR | **−21.5 MHz** | **−20.1(1) MHz** (7 %) | Ng et al., PRA 105, 022823 (2022) = arXiv:2202.01346, §II D and Table I |
| ThF⁺ X³Δ₁, `d_mf` | same run | 3.46 D | 3.37(9) D (2.7 %) | ibid., Table I |
| HfF⁺ ³Δ₁, `A∥(¹⁹F)` | 4c MR12-CISD(20), TZ | "agrees qualitatively" | −62.0(2) MHz | Fleig, PRA 96, 040502 (2017) = arXiv:1706.02893, Results §; measurement Cairncross et al. PRL 2017 |
| BaF X and excited, `A∥(¹⁹F)` | 4c FSCC, DIRAC, full uncertainty budget | agrees within theoretical uncertainties | recent expt | NL-eEDM Collaboration (Denis, Haase, Mooij, …, Borschevsky et al.), PRA 105, 052811 (2022), abstract |
| ²²⁹Th³⁺ (atomic ion) | DCB + QED relativistic CC | used to extract μ = 0.359(9), Q = 2.95(7) e·b | Campbell et al. PRL 106, 223001 (2011) | Li, Qiao, Tang & Shi, PRA 104, 062808 (2021), abstract |

The Ng 2022 row is the one that matters: **it is ThF⁺, in the X³Δ₁ state, with
the ¹⁹F hyperfine constant, and it is right to 7 % with the right sign.** That
is the empirical accuracy floor for any A∥ prediction in this molecule.

Fleig's own caveat is worth carrying (arXiv:1706.02893, verbatim): "the present
electronic-structure model has not been designed with a focus on properties
depending on spin density in the vicinity of the fluorine nucleus" — a CI model
optimised for E_eff on the heavy centre is not automatically good at the light
nucleus, and vice versa.

### 1.2 Electric quadrupole `eQq₀` — the axial EFG at Th

Petrov 2018 Eq. (22):
`eQq₀ = 2eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₀(θ_1i,φ_1i)/r_1i³ |³Δ₁⟩`.

This is a plain diagonal expectation value of a one-electron rank-2 operator on
the Th centre, i.e. exactly the electric field gradient. DIRAC's `.EFG` computes
the electronic contribution

`φ^[2]el_αβ(R_K) = (−e/4πε₀) ⟨(3 r_Kα r_Kβ − δ_αβ r_K²)/r_K⁵⟩`

plus the point-charge nuclear term, and `.NQCC` converts it
(`properties.rst`, verbatim):

> NQCC [in MHz] = 234.9647 × Q [in b] × q [in atomic units E_h/e a₀²]

with `q = φ^[2]_zz/e` in the principal axis system. For a linear molecule the
principal axis is the internuclear axis and the asymmetry η = 0 by symmetry, so
`.NQCC` at the Th centre *is* `eQq₀` up to the sign convention (§4e).

**Sign-convention warning, and it is a live one.** B&C's convention is that
"q₀ is the negative of the electric field gradient" ([DIGEST] §3.2, §3.4 — a
verified factor of −1 against the Casimir function). Separately, Derevianko,
Perera, Krośnicki, Nalikowski, Morgan & Veryazov, arXiv:2601.07098 (abstract),
explicitly "address persistent differences in EFG sign conventions and tensor
definitions employed in the literature and in widely used quantum chemistry
codes". So there are at least three conventions in play (B&C, Petrov's Eq. 22,
whatever the code prints) and they must be reconciled explicitly, not assumed.

**Level of theory used for the analogous system.** Petrov 2018 §IV, verbatim:

> To compute eQq₀ in the ³Δ₁ state of HfF⁺ we have performed relativistic
> coupled cluster calculations within the Dirac-Coulomb Hamiltonian using the
> DIRAC12 code. In all the calculations the Hf-F internuclear distance in the ³Δ₁
> state was set to 3.41 a₀. In the main calculation all 80 electrons of HfF⁺ were
> included in the correlation treatment within the coupled cluster calculations
> with single, double, and perturbative triple amplitudes, CCSD(T), using the
> uncontracted Dyall CVTZ basis set for Hf and the ccpVTZ basis set for F. We
> have also applied the correction on the basis set expansion up to the
> uncontracted Dyall AEQZ basis set for Hf and the aug-ccpVQZ for F. In the
> calculation 1s–3d core electrons of Hf were excluded from the correlation
> treatment within the CCSD(T) method. Accounting for the perturbative triple
> cluster amplitudes contributes ≈88 MHz in eQq₀(¹⁷⁷HfF⁺).

They quote no uncertainty on the resulting −2100 MHz. The triples contribution
of 88 MHz is 4 % of the total, which is a weak *lower* bound on the correlation
sensitivity `[inference]`.

**Modern accuracy benchmark for molecular EFG.** Fabbro, Brandejs & Saue,
arXiv:2602.00855 (abstract), determine Q(⁸⁷Sr) from the measured NQCCs of SrO
and SrS plus computed EFGs: X2C molecular-mean-field Hamiltonian, CCSD(T)/CCSD-T
/CCSD-T̃, Gaunt included, even-tempered optimised quadruple-ζ basis, vibrational
corrections, arriving at Q = 0.33666 ± 0.00258 b — a **0.8 % total uncertainty**,
about 10 % away from the accepted standard value. That is the current
state-of-the-art for this exact class of property; it also shows how much
machinery (iterative triples with a distributed tensor library) goes into 1 %.

For an *estimate* rather than a benchmark, `[inference]` a plain
CCSD(T)/uncontracted-CVTZ EFG in the Petrov 2018 style should be good to
10–20 %, which is already a factor 5–10 better than [DIGEST] §4.3's
"−2.6(1.0) GHz" HfF⁺-anchored guess.

### 1.3 Nuclear spin-rotation `c_I`

`H_nsr = c_I T¹(J)·T¹(I)` (B&C Eq. 8.7, [DIGEST] §3.2d). The property-theory
object is the nuclear spin-rotation tensor `M_K`, a **second-order (linear
response)** quantity, not an expectation value: it is the cross-derivative of the
energy with respect to the nuclear magnetic moment and the molecular rotation.
DIRAC's `**PROPERTIES .SPIN-ROTATION` entry (verbatim from `properties.rst`):

> Evaluate nuclear spin-rotation constants: linear response, expectation value
> and nuclear contributions, see `Aucar_JCP2012` and `AucarChap2019`.
> `M_{K;μν} = M^nuc_{K;μν} + M^elec_{K;μν}` with
> `M^elec_{K;μν} = −(ℏ²/h) ∂²/∂I_{K;μ}∂L_ν ⟨⟨ĥ^hfs_K; ĥ^BO⟩⟩₀`

where `ĥ^BO = −ω·Ĵ_e`, `ω = L·I⁻¹` is the first-order correction to the
Born–Oppenheimer approximation. So the operator pair is *hyperfine × rotational*
— the same `ĥ^hfs` as §1.1, contracted against the rotational perturbation. That
is the relativistic generalisation of the paramagnetic + diamagnetic Ramsey
decomposition.

**On the Flygare relation.** The nonrelativistic Ramsey–Flygare mapping between
the nuclear shielding tensor σ and the spin-rotation tensor M (W. H. Flygare,
*Chem. Rev.* **74**, 653 (1974), "Magnetic interactions in molecules and an
analysis of molecular electronic charge distribution from magnetic parameters" —
title and year from a web search, 2026-09-05; **I did not open the paper**, so
treat the exact page as UNVERIFIED) is the standard shortcut for estimating
`c_I` from a shielding calculation. **Do not use it at Th.** The relativistic
treatment is what DIRAC implements, and the modern literature is explicit that
the nonrelativistic mapping fails: see Aucar, Colombo Jofré & Aucar,
arXiv:2208.04788, "Relativistic relationship between nuclear-spin-dependent
parity-violating NMR shielding and nuclear spin-rotation tensors" (2022), and
the review chapter DIRAC cites as `AucarChap2019` ("Relativistic Theory of
Nuclear Spin-Rotation Tensor", Springer). At Z = 90 the Flygare relation is a
back-of-envelope, not a prediction `[inference]`.

**Nobody has computed `c_I` for either nucleus of ThF⁺** ([DIGEST] §5, G3). I
found no molecular NSR calculation for any actinide diatomic in this session's
searches (arXiv `abs:"nuclear spin-rotation"` → 13 results, none actinide;
`abs:"spin-rotation" AND abs:"relativistic" AND abs:"shielding"` → 2 results,
both Aucar-group methodology papers).

**Caveat on applicability.** DIRAC's `.SPIN-ROTATION` is a linear-response
property built on a reference wavefunction. Whether it is usable for an
open-shell Ω = 1 state (as opposed to closed-shell diamagnetic molecules, which
is the NSR literature's home turf) is **UNVERIFIED**. This is the single
weakest-supported item in this document. For a ³Δ₁ state the "nuclear
spin-rotation" constant that appears in B&C Eq. (8.7) also picks up an effective
contribution from the electronic-orbital motion, and I have not checked that the
two definitions coincide.

### 1.4 `eQq₂` — the ΔΩ = 2 quadrupole constant

Petrov 2018 Eq. (23):
`eQq₂ = 2√6 eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₂(θ_1i,φ_1i)/r_1i³ |³Δ₋₁⟩`.

Note the ket: **Ω = −1**. This is not an expectation value; it is an off-diagonal
matrix element between the two degenerate components of the ³Δ₁ manifold. That
single fact is why no code has a keyword for it. Full treatment in §5.

Petrov 2018 §IV, verbatim on why it is nonzero at all:

> As follows from Eq. (23), eQq₂ has no nonzero matrix elements within a main
> nonrelativistic term ³Δ₁. The main contribution to eQq₂ is due to the
> spin-orbit admixture of the [Π] state with the leading configuration |5s5dπ|
> composed of 5s and 5d atomic orbitals of Hf. Then one can obtain
> eQq₂ = 483 w Q ⟨1/r³⟩_5d MHz, (24)
> where w is the weight of the admixture … and ⟨1/r³⟩_5d = 4.86 a.u. as obtained
> from the Hartree-Fock-Dirac calculations for Hf⁺. The value of w can be
> estimated from the sensitivity to it of the body fixed g-factor
> G∥ = (1/Ω)⟨³Δ₁| L̂_e n̂ − g_S Ŝ_e n̂ |³Δ₁⟩ ≈ 2 − 2.002319 + w (25)

`[derived]` consistency check on the prefactor: 483 / 234.9647 = **2.056**, where
234.9647 is DIRAC's universal Q[b]·q[a.u.] → MHz conversion (§1.2). The residual
≈ 2 is plausibly the angular factor in Eq. (23), which would make **483 a pure
unit-conversion constant, transferable to Th unchanged**. That supports
[DIGEST] §4.4's scaling. Marked UNVERIFIED — Petrov does not say.

---

## 2. Which code does what

### 2.1 DIRAC (four-component and X2C)

Retrieved from the manual source in the DIRAC git repo, 2026-09-05
(`diracprogram.org` itself was down).

The complete `**PROPERTIES` keyword list on master is:

```
PRINT, ABUNDANCIES, RKBIMP, NOPCTR, RDCCDM, DIPOLE, QUADRUPOLE, EFG, NQCC,
POLARIZABILITY, FIRST ORDER HYPERPOLARIZABILITY, TWO-PHOTON, NMR, SHIELDING,
MAGNET, ROTG, SPIN-SPIN COUPLING, DSO, NSTDIAMAGNETIC, ESR, OPTROT, VERDET,
MOLGRD, PVC, PVCSHI, PVCSR, PVCEFG, RHONUC, EFFDEN, SPIN-ROTATION
```

and `**WAVE FUNCTION` offers:

```
SCF, RESOLVE, COSCI, MP2, MVO, MP2 NO, RELCCSD, RELADC, POLPRP, DIRRCI, LUCITA,
EXACC, CASPT2, REORDER MO, ORBROT, POST SCF REORDER MO, PHCOEF, KRCI, KRMCSCF,
LAPLCE
```

Mapping to our constants:

| constant | DIRAC route | notes |
|---|---|---|
| `eQq₀` | `**PROPERTIES .EFG` or `.NQCC` | expectation value; **correlated** via `.RDCCDM`, which "activates the reading of the file CCDENS obtained from a previous CC calculation with either `.RELCCSD` or `.EXACC`" (verbatim). So the CC density, not just HF, feeds the EFG. |
| `A∥` | no predefined keyword; user-defined `**HAMILTONIAN .OPERATOR` with the `α × r/r³` components, or finite-field with `'FC NAMab'`-style integrals | The manual's worked example is a finite-field Fermi-contact calculation on PbX (`one_electron_operators.rst`), including the scaling factor `3/(4π g_e) = 1/8.3872954891254192`. This is the pattern Denis 2015 and Fleig 2017 used (both say their implementation is described in their Refs.). |
| `c_I` | `**PROPERTIES .SPIN-ROTATION` | linear response; open-shell Ω = 1 applicability UNVERIFIED (§1.3). Related: `.ROTG` (rotational g-tensor), `.NMR`. |
| `eQq₂` | **nothing predefined** | see §5. `.ESR` ("g-tensors and hyperfine coupling tensors, using first-order quasi-degenerate perturbation theory based on configuration interaction") is the closest existing machinery — QDPT over a CI space is exactly how you get off-diagonal-in-Ω elements — but the manual names only *magnetic* tensors, not the EFG. |
| open-shell reference | `.KRCI` (Kramers-restricted CI), `.KRMCSCF`, `.RELCCSD`, `.EXACC` | KRCI is Fleig's module and is what Denis 2015 used (GASCI). |

**Pitfalls the manual states or implies:**

- **Gaunt/Breit.** `hamiltonian.rst`, verbatim: "The Gaunt term, which contains
  the spin-other orbit interaction, **is implemented at the SCF level in
  DIRAC**." NL-eEDM 2022 (PRA 105, 052811) says the same in practice: "a rough
  estimate of the Breit contribution can be extracted by including the Gaunt
  interaction in the calculation. However, this contribution can presently be
  calculated only on the Hartree-Fock level"; they measured ≈ 0.2 % at DHF for
  the BaF ground state and dropped it from the budget. For Th (Z = 90) the Gaunt
  contribution to a core-region property is larger; Skripnikov & Titov 2015 warn
  that "for some properties even taking account of Breit interaction (mainly
  between valence and core electrons of Th) can be important". Fabbro/Saue 2026
  *do* include Gaunt for the ⁸⁷Sr EFG.
- **Finite nuclear model.** `hamiltonian.rst` states point-nucleus is the default
  only for `.NONREL`; the Gaussian finite-nucleus model is DIRAC's default for
  relativistic Hamiltonians. I did **not** locate the `.NUCMOD`-equivalent
  keyword in `molecule.rst`/`general.rst` this session — **UNVERIFIED**, check it
  before running. For a core-region property at Z = 90 the nuclear model is not
  optional.
- **Virtual space.** `.RELCCSD` defaults to correlating orbitals between −10 and
  +20 hartree; both Petrov 2018 (2000 a.u. cutoff in a related paper) and Ng 2022
  ("virtual orbitals higher than 100 hartree" frozen) reset this. Uncontracted
  Dyall basis sets generate very large virtual spaces; Yuan, Visscher & Gomes,
  JCP 156, 224108 (2022) = arXiv:2202.01146 show MP2 frozen natural orbitals cut
  the virtual space to "around half" with controlled error on EFGs specifically.
  That is the cost lever if the calculation does not fit.

### 2.2 CFOUR + X2CAMF — *the route that already exists for this molecule*

Ng et al., PRA 105, 022823 (2022), §II D, verbatim:

> We perform numerical differentiation of coupled-cluster singles and doubles
> augmented with a non-iterative triples correction [CCSD(T)] energies to obtain
> d_mf, A∥, and G∥. … These calculations treat relativistic effects using an
> exact two-component (X2C) Hamiltonian with atomic mean-field spin-orbit (AMF)
> integrals. **We use the CFOUR program package for all the electronic structure
> calculations presented here.** … Calculations of d_mf and A∥ use uncontracted
> ANO-RCC basis sets. … All CCSD(T) calculations freeze sixty-four core electrons
> and virtual orbitals higher than 100 hartree.

Result: `A∥(¹⁹F) = −21.5 MHz` vs measured `−20.1(1)`; `d_mf = 3.46 D` vs
`3.37(9)` (their Table I). The calculation was done by **Lan Cheng** (Department
of Chemistry, Johns Hopkins University), a co-author, with computational support
acknowledged from NSF PHY-2011794.

CFOUR computes EFG tensors — the CFOUR wiki has a dedicated module `xefgiso` for
transforming a computed EFG tensor between isotopomers
(`https://cfour.uni-mainz.de/cfour/index.php?n=Main.DeterminationElectric-fieldGradientsForOtherIsotopomers`,
retrieved 2026-09-05). The exact keyword (`PROP=FIRST_ORDER` / `PROP=1`, giving
the correlated density from analytic gradients) comes from a web-search summary
and the ASH interface docs, **not** from the CFOUR manual, which I could not
reach — **UNVERIFIED**. Whether the X2CAMF property machinery in CFOUR covers
the EFG (as opposed to the hyperfine and dipole operators Ng used) is likewise
**UNVERIFIED and is the single most valuable thing to ask Lan Cheng.**

### 2.3 Skripnikov's two-step GRECP + restoration

Skripnikov & Titov 2015 §"Two-step approach", and Skripnikov, JCP (2016)
"Combined 4-component and relativistic pseudopotential study…"
(`docs/lit/skripnikov2016-combined-4c.pdf`), verbatim from the latter:

> The two-step method allows one to consider high-order correlation effects and
> large basis sets with rather modest requirements to computer resources in
> comparison to four-component approaches. However, some uncertainty remains due
> to the impossibility to consider the full version of the GRECP operator in the
> currently available codes and neglect of the inner-core correlation effects.

Their combined scheme, verbatim, is a four-way decomposition:

> (i) the main correlation contributions within the 52-electron four-component
> Dirac-Coulomb coupled cluster with single, double, and noniterative triple
> cluster amplitudes [CCSD(T)] theory; (ii) the inner-core correlation
> contributions; (iii) correction on inclusion of the Gaunt interaction; (iv) the
> contribution of high-order correlation effects up to … CCSDT(Q) for the valence
> electrons within the two-component … two-step approach

and the codes, verbatim: "4-component Dirac-Coulomb(-Gaunt) Hartree-Fock
calculations were performed within the **dirac12** code. Scalar-relativistic
coupled cluster … correlation calculations were performed within the **cfour**
code. All 4-component coupled cluster calculations as well as scalar-relativistic
CC and CI calculations with the treatment of the high-order cluster amplitudes
and excitations were performed within the **mrcc** code. The nonvariational
restoration code developed by us … and interfaced to these program packages was
used to restore the 4-component electronic structure near the Th nucleus."

For ThF⁺ specifically (Skripnikov & Titov 2015, Results): 1s–4f inner-core Th
electrons excluded via GRECP; main calculation **38e-2c-CCSD(T)** with MBas
`(30,20,10,11,4,6,5)/[30,8,10,4,4,2,1]` on Th and aug-cc-pVQZ-minus-2g on F;
CCSDT(Q)−CCSD(T) correction with 20 outer-core electrons frozen (Th 5s²5p⁶5d¹⁰,
F 1s²) in the CBasSO ANO basis; basis-enlargement corrections; a vibrational
contribution at R = 3.75 a.u. Their Table II decomposition for `A∥`:

| step | `A∥` (μ_Th/μ_N · MHz) | `G∥` |
|---|---|---|
| 38e-2c-CCSD | −4214 | 0.039 |
| 38e-2c-CCSD(T) | −4164 | 0.033 |
| correlation correction | +13 | +0.001 |
| basis set correction | −14 | — |
| vibrational contribution | +2 | — |
| **FINAL(ThF⁺)** | **−4163** | **0.034** |

`[derived]` Note how small every correction past CCSD(T) is: triples move `A∥` by
1.2 %, and everything after that by < 0.4 %. **The property is converged at
CCSD(T)**, which is the strongest argument that a single CCSD(T)-level run is
enough for our purposes. Denis 2015 independently reports the same insensitivity:
"The hyperfine interaction constant A∥ is insensitive to these higher
excitations, allowing triple excitations to the virtual space changes the value
by only 0.2 %" — while E_eff in the same calculation "exhibit[s] a strong
dependence on higher excitations. The inclusion of triple excitations yields a
drop of 25 %". **A∥ is an easy property; E_eff is the hard one.** That is good
news for us and is the reason to be optimistic about cost.

### 2.4 Fleig's KR-CI in DIRAC, and the Groningen FSCC line

Denis et al. 2015 (ThF⁺) used DIRAC12/DIRAC11 with the exact-two-component
Hamiltonian of Iliaš and Saue, in both atomic-mean-field (`amf`) and
molecular-mean-field (`mmf`) variants; Th basis sets uncontracted Dyall
`cv2z [26s23p17d13f1g]`, `cv3z [33s29p20d15f5g1h]`, `cv4z
[37s34p26d23f9g5h1i]`; F aug-cc-pVnZ; wavefunctions from GASCI/MRCI and
intermediate-Hamiltonian Fock-space CC (IHFSCC); geometry R = 1.981 Å.

This is the natural home for `eQq₂` (§5) because KRCI carries both Ω components
in one CI space.

The Groningen/NL-eEDM line (Borschevsky, Hao, Denis, Haase, Chamorro et al.,
PRA 105, 052811 (2022)) is the group that publishes **uncertainty budgets** —
basis-set, correlation, virtual-space cutoff, higher-order relativistic, summed
in quadrature — which is what a number entering `heff` actually needs.

### 2.5 EXP-T / relativistic Fock-space CC

NL-eEDM 2022 cites Oleynichenko et al. as the only recent example of accurate
molecular HFS by FSCC (for KCs). Zaitsevskii, Oleynichenko et al. use EXP-T for
RaF/RaCl (arXiv:2012.09818 and related). This is a live option for a
multireference ³Δ₁ but I found no ThF⁺ application, so it would be new work.

---

## 3. Cost, feasibility, and the collaboration route

### 3.1 Cost — what I know and what I do not

**I do not have core-hour or memory numbers.** None of Petrov 2018, Skripnikov &
Titov 2015, Denis 2015, Ng 2022 or Fabbro/Saue 2026 report timings, node counts,
or memory. Anyone who tells you "X core-hours" for ThF⁺ without a benchmark is
guessing. What I can give is the structure of the cost.

`[derived]` **Basis size.** Denis 2015 states the uncontracted Dyall cv3z set for
Th as `[33s29p20d15f5g1h]`. In spherical harmonics that is
33·1 + 29·3 + 20·5 + 15·7 + 5·9 + 1·11 = **381 large-component functions on Th
alone**; cv4z `[37s34p26d23f9g5h1i]` gives
37 + 102 + 130 + 161 + 81 + 55 + 13 = **579**. In a four-component calculation
kinetic balance roughly doubles this for the small component `[inference]`. So
a 4c CCSD(T) on ThF⁺ at cv3z is a ~1000-spinor problem with ~98 electrons
available to correlate.

`[inference]` **Scaling.** CCSD(T) is O(N⁷) in basis size and steeply increasing
in correlated-electron count. Petrov 2018 correlated all 80 electrons of HfF⁺ at
uncontracted Dyall CVTZ in DIRAC12 in 2018 — so an 80-electron 4c-CCSD(T) at TZ
was already a publishable-scale, not heroic, calculation eight years ago. ThF⁺
has 98 electrons `[derived]`: Th(90) + F(9) − 1. Two mitigations are already
standard practice: freeze the deep core (Petrov froze Hf 1s–3d; Ng froze 64 core
electrons; Skripnikov removes Th 1s–4f entirely via GRECP), and truncate the
virtual space (Ng: > 100 E_h; Yuan/Visscher/Gomes MP2FNO: half).

**Honest bracket, marked as inference, not measurement:**

- X2CAMF-CCSD(T) with ~34 correlated electrons and ANO-RCC-unc, the Ng 2022
  recipe: this ran on ordinary academic resources for a PRA paper. Plausibly
  **hundreds to a few thousand core-hours**, i.e. a large workstation
  (32–64 cores, 256 GB) over days, or a small cluster allocation.
- Full-electron 4c-DC-CCSD(T) in the Petrov 2018 style with QZ basis-set
  corrections and Gaunt: **cluster**, not workstation. Tens of thousands of
  core-hours is my guess and I would not defend it.
- The Fabbro/Saue 1 %-accuracy protocol (iterative triples, distributed tensor
  library, vibrational corrections): a dedicated methods project, not a
  side-calculation.

**Realistic uncertainty on the deliverable:**

| quantity | expected uncertainty | basis for the claim |
|---|---|---|
| `A∥(Th)` magnitude | **5–10 %** | Ng 2022 got 7 % on `A∥(¹⁹F)` in this molecule; Skripnikov quotes 7 % for his ThF⁺ properties ([HAM] §8) |
| `A∥(Th)` sign | **unambiguous**, if `A∥(¹⁹F)` is computed alongside (§4d) | — |
| `eQq₀(Th)` electronic factor `q₀` | **10–20 %** at CCSD(T)/CVTZ; ~1 % achievable but only with the full Fabbro/Saue-class protocol | `[inference]` from the 4 % triples contribution in Petrov 2018 and the 0.8 % achieved by Fabbro/Saue with far more machinery |
| `eQq₀(Th)` as an MHz number | **inherits the Q(²²⁹Th) dispute**: 3.11(2) e·b (Porsev 2021) vs 2.95(7) e·b (Li 2021) — a 5 % spread larger than the quoted error bars | [DIGEST] §1.3; Li et al. PRA 104, 062808 (2021) abstract |
| `c_I` | unknown; no precedent for an actinide diatomic | — |
| `eQq₂` | order of magnitude only, by any route currently available | §5 |

**The Q dispute is the reason to demand `q₀` in atomic units** rather than
`eQq₀` in MHz. `heff` can then rescale when the nuclear moment settles, exactly
as [DIGEST] §2.2 does for `A∥/g_N`.

### 3.2 Who has the pipeline

Affiliations as of the cited publications; verify before writing.

| group | what they have already published on this exact problem | what to ask them for |
|---|---|---|
| **Lan Cheng** (Chemistry, Johns Hopkins) — with **JILA (Ng, Cornell, Ye)** | The ThF⁺ X³Δ₁ `A∥(¹⁹F)`, `d_mf`, `G∥` calculation itself: X2CAMF-CCSD(T)/ANO-RCC-unc in CFOUR, Ng et al. PRA 105, 022823 (2022) §II D | **The shortest path.** He has the converged ThF⁺ reference, the basis, the frozen-core recipe, and an existing JILA collaboration. Ask: can the same run report the EFG at Th and `A∥(²²⁹Th)`? |
| **Skripnikov & Petrov** (PNPI Gatchina / St Petersburg State U) | `A∥(²²⁹Th)`, `W_M`, `G∥` for ThF⁺ (PRA 91, 042504); `eQq₀` and `eQq₂` for HfF⁺ (PRA 98, 042502); the ²³²ThF⁺ E-field-dependent g-factor paper (arXiv:2503.02840) | **The only group that has ever computed a ³Δ₁ `eQq₀`, and the only one that has ever written down `eQq₂`.** They also already have the ThF⁺ effective-Hamiltonian basis (Petrov 2025 Eq. 3 lists the interacting states including ³Π₀±). Ask for `eQq₀`, `eQq₂`, and the sign convention behind their −4163. |
| **Fleig** (LCPQ/IRSAMC, Toulouse) with **Nayak** (BARC), **Knecht**, **Gomes** (Lille), **Jensen** (Odense) | The other ThF⁺ paper (NJP 17, 043005) — DIRAC KRCI/GASCI/IHFSCC, and Fleig's HfF⁺ PRA 96, 040502 | The natural home for **`eQq₂`**, because KRCI carries both Ω = ±1 components in one wavefunction (§5). Also the other side of the `A∥` sign disagreement. |
| **Borschevsky / Hao / Chamorro** (Van Swinderen Institute, Groningen; NL-eEDM) | FSCC hyperfine with full uncertainty budgets, PRA 105, 052811 (2022) | If you want a **number with a defensible error bar** rather than a number. |
| **Saue / Fabbro / Brandejs** (Toulouse; DIRAC developers) | The current best molecular-EFG protocol, arXiv:2602.00855 | If 1 %-level `eQq₀` ever matters. |
| **Zaitsevskii / Oleynichenko** (EXP-T) | FS-RCC for RaF/RaCl | Alternative multireference route; no ThF⁺ precedent. |

### 3.3 What a request must specify

A collaborator can waste months on the wrong convention. Specify, in writing:

1. **State:** ThF⁺ **X ³Δ₁**, v = 0, the Ω = 1 component. (Note this is the
   *ground* state — Denis 2015 and Skripnikov 2015 predate/argue the assignment;
   Ng 2022 and later JILA work treat it as X.)
2. **Geometry:** R(Th–F) = **3.75 a₀ = 1.984 Å** (Skripnikov & Titov 2015,
   38e-2c-CCSD(T) equilibrium, and the value they used for all properties). Note
   Denis 2015 used **1.981 Å** (the CCSDT(Q) equilibrium of the ¹Σ⁺ state from
   Barker 2012). The 0.15 % difference is below the property uncertainty
   `[inference]` but the number should be stated, not assumed.
3. **Isotope and nuclear moments:** ²²⁹Th, I = 5/2, and **ask for the electronic
   factors, not the isotope products**: `q₀` in a.u. (E_h/e a₀²) and
   `A∥/g_N` in MHz. If they insist on MHz, state which μ and Q they used
   (recommend μ = 0.366(6) μ_N, Porsev 2021; and give *both* Q = 3.11(2) and
   Q = 2.95(7)).
4. **Axis convention:** see §4e — the whole point.
5. **Also compute `A∥(¹⁹F)`** in the same wavefunction, as the convention anchor
   and accuracy check against the measured −20.1(1) MHz.
6. **Nuclear model:** finite (Gaussian) nucleus for Th, stated explicitly.
7. **Deliverable format:** the decomposition table, not just the final number —
   Skripnikov & Titov's Table II layout (CCSD → CCSD(T) → correlation correction
   → basis correction → vibrational) is exactly right, because it lets us read
   off the convergence and set our own error bar.

---

## 4. Recommended minimal plan

### 4a. The one calculation

**A single X2CAMF-CCSD(T) run on ThF⁺ X³Δ₁ at R = 3.75 a₀ that reports
`A∥(²²⁹Th)`, `A∥(¹⁹F)`, the EFG `q₀` at Th in a.u., and `G∥`, in one stated axis
convention.**

Why this and not something bigger:

- It closes **G4 (the `A∥` sign) for free** via the measured ¹⁹F anchor (§4d).
- It closes the **electronic half of G2** (`eQq₀`), replacing [DIGEST] §4.3's
  ±40 % HfF⁺-anchored estimate with a ~15 % number.
- `A∥` is demonstrably converged at CCSD(T) (§2.3: triples move it 1.2 %,
  everything beyond 0.4 %) — so the expensive high-order corrections that E_eff
  needs are not needed here.
- It reuses the exact recipe that already reproduced ThF⁺'s measured `A∥(¹⁹F)`
  to 7 % (Ng 2022), so the accuracy is anchored, not asserted.
- `G∥` is measured (|G∥| = 0.0476, [HAM] §2.8) and computed twice (0.034
  Skripnikov 2015; 0.047 fitted by Petrov 2025), so it is a third free check on
  the wavefunction.

If DIRAC rather than CFOUR: `**HAMILTONIAN` Dirac-Coulomb (or `.X2Cmmf`),
`**WAVE FUNCTION .RELCCSD`, `**PROPERTIES .EFG .NQCC` with `.RDCCDM` to use the
CC density, plus a user-defined `.OPERATOR` (or finite-field FC-style
perturbation) for `A∥`. Uncontracted Dyall cv3z on Th + cc-pVTZ on F for the main
number, cv4z/aug-cc-pVQZ for the basis correction, Th 1s–3d (or deeper) frozen,
Gaunt at DHF as a correction only.

### 4b. What to ask a collaborator

**First ask, to Lan Cheng (cc: Ng / Cornell if the JILA route is open):** you
already have the converged ThF⁺ X³Δ₁ X2CAMF-CCSD(T) wavefunction from PRA 105,
022823. Can that same calculation output (i) `A∥` for a ²²⁹Th nucleus, (ii) the
electric field gradient at the Th centre in atomic units, in a stated sign
convention, and (iii) `A∥(¹⁹F)` from the same run so the convention is anchored
on your published −21.5 MHz? This is a re-run of an existing input with extra
property requests, not a new project.

**Second ask, to Skripnikov / Petrov:** you computed `eQq₀` and estimated `eQq₂`
for ¹⁷⁷/¹⁷⁹HfF⁺ (PRA 98, 042502 §IV) and `A∥` for ²²⁹ThF⁺ (PRA 91, 042504
Table II). Would you do the same pair for ²²⁹ThF⁺? Also: what molecular-axis and
Ω sign convention underlies the −4163 in your Table II — we get the opposite sign
from Denis et al. NJP 17, 043005.

**Third ask, to Fleig (and/or Nayak, Knecht, Gomes):** the `eQq₂` question of §5
— is the off-diagonal ΔΩ = 2 element of the EFG operator between the two ³Δ₁
components accessible in a KRCI property calculation?

### 4c. What can be estimated meanwhile, and how good it is

| constant | interim value | how good | improvement available *now*, no collaborator |
|---|---|---|---|
| `A∥(²²⁹Th)` | \|A∥\| = 1.51(6) GHz, sign unknown ([DIGEST] §2.2) | magnitude good to ~7 % (two independent calculations agree to 2.2 %); **sign unknown** | none — needs §4a |
| `eQq₀(²²⁹Th)` | −2.6(1.0) GHz ([DIGEST] §4.3) | factor ~1.5 either way | none |
| `eQq₂(²²⁹Th)` | 200–400 MHz ([DIGEST] §4.4) | order of magnitude | **yes** — see below |
| `c_I(Th)` | unconstrained, 1 kHz – 1 MHz ([DIGEST] §4.6) | three decades | none; treat as a free parameter |
| `c_I(¹⁹F)` | ~20 kHz, factor 3 ([HAM] §2.6) | factor 3 | none |

**The one free improvement: `eQq₂`.** [DIGEST] §4.4 carries a guessed bracket
`⟨1/r³⟩_6d,Th / 4.86 a.u. ≈ 0.6–1.0`. That bracket is removable by **one atomic
Dirac–Hartree–Fock calculation of ⟨1/r³⟩ for the Th⁺ 6d orbital** — the exact
analogue of what Petrov did for Hf⁺ (his 4.86 a.u.). That is minutes of compute
in any relativistic atomic code and turns a factor-1.7 guess into a number.

Two further points that firm up [DIGEST] §4.4 and that I can now assert:

1. **The sign of `w` is settled and positive.** Petrov 2018 Eq. (25) gives
   `w = G∥ + 0.002319`. Skripnikov & Titov 2015 Table II computes
   `G∥(ThF⁺) = +0.034`, and Petrov et al. 2025 (arXiv:2503.02840, local copy)
   state verbatim: "We also used G∥ = 0.047 to reproduce the g = 0.0149 value for
   the g-factor of the F = 3/2 state. Our value is slightly different from
   G∥ = 0.048 estimated in Ref. [11]". So `G∥ > 0`, hence `w = 0.049 > 0`
   `[derived]`, hence `eQq₂ > 0` within Petrov's model, exactly as [DIGEST] §4.4
   assumed. Good — this was a real risk (a weight `w` cannot be negative, so a
   negative `G∥` would have broken the model, not just flipped the sign).
2. **The admixed state has an identified analogue.** Petrov 2018 attributes the
   HfF⁺ `eQq₂` to a Π state of leading configuration |5s5dπ|; Petrov et al. 2025
   list ³Π₀∓ states in the ThF⁺ basis with "the same configuration (∼7s¹6d¹ …) as
   ³Δ₁". So the ThF⁺ analogue is not hypothetical. `[inference]` Fleig 2017
   attributes the HfF⁺ `G∥` deviation instead to a ¹Π₁ state (6s→6p on Hf) with
   `G∥(¹Π₁) = +1` — *the two papers name different admixed Π states*. Since
   Petrov's Eq. (24)/(25) calibration assumes one and the same state drives both
   `G∥` and `eQq₂`, this is a genuine soft spot in the estimate. Flag it; it
   caps the accuracy of the `w`-from-`G∥` route regardless of how well
   `⟨1/r³⟩` is computed.

### 4d. Fixing the `A∥` sign — the free-rider trick

[DIGEST] §5 G4 records that Skripnikov prints −4163 and Denis prints +1833 from
the same defining equation, agreeing to 2.2 % in magnitude, with neither stating
its axis convention.

`[derived]` **A sharpening that changes the diagnosis.** `A∥` as both papers
define it — `(μ/(IΩ))⟨(α × r/r³)_ζ⟩` — is **invariant** under the simultaneous
flip n̂ → −n̂, Ω → −Ω: the ζ projection changes sign and so does the Ω in the
denominator. So *if both papers applied their conventions self-consistently, the
sign cannot be a convention difference at all* — it is either a genuine
disagreement between the two wavefunctions, a typo, or one paper referring ζ and
Ω to opposite senses of the axis. [DIGEST] §2.2 already suspected the last; the
point here is that "different but equally valid conventions" is *not* an
available explanation. Something is actually wrong in one of the two papers.

`[derived]` **`G∥` cannot arbitrate.** `G∥ = (1/Ω)⟨L̂_e·n̂ − g_S Ŝ_e·n̂⟩` is
*also* invariant under n̂ → −n̂, Ω → −Ω, so the fact that Skripnikov's +0.034
matches the experimental +0.047 in sign tells us nothing about his ζ direction.
Do not use it as a check; it is a false friend.

**What does arbitrate: `A∥(¹⁹F)`.** It is measured, `−20.1(1) MHz`
([HAM] §2.4, Ng 2022 Table I), and it has already been reproduced as `−21.5 MHz`
by X2CAMF-CCSD(T) (Ng 2022). Any calculation that reports `A∥(¹⁹F)` and
`A∥(²²⁹Th)` **from the same wavefunction in the same convention** has its
convention validated by the F value and transfers that validation to the Th
value. This is the cheapest possible resolution of G4 and it is why §4a insists
on computing both. Note the F and Th constants need not have the same sign — in
HfF⁺ they do (`A∥(F) = −62.0`, `A∥(¹⁷⁷Hf) = −1429`, Petrov 2018 §III) — but the
*relative* sign is the physical output and it is convention-free.

### 4e. The input choices that must be fixed, verbatim, in any request

1. **n̂.** `heff` fixes `n̂ = F → Th` ([HAM] §7). Note that Petrov et al. 2025
   use the opposite: "D = −0.133 a.u. (with respect to center of nuclear mass;
   the molecular axis directed from Th to F)". State which one you want and
   restate it in every table caption.
2. **Ω.** Ω = +1 means the component with Λ = +2, Σ = −1 **with respect to that
   same n̂**. Every `1/Ω` in a property definition must use this Ω, not |Ω|.
3. **`A∥`:** report `A∥ = (μ_Th/(IΩ))⟨Ψ_Ω| Σ_i (α_i × r_i/r_i³)_ζ |Ψ_Ω⟩` with
   ζ the projection on that same n̂ and with the same Ω in numerator state and
   denominator. Report `A∥/g_N` as well.
4. **`q₀`:** report the raw EFG component in a.u. **and** say which sign
   convention: `q = φ^[2]_zz/e` (DIRAC's `.NQCC`) versus B&C's "q₀ is the
   negative of the electric field gradient" versus Petrov's Eq. (22). These
   differ by −1 relative to each other in ways that are documented as a known
   trap: Derevianko et al. arXiv:2601.07098 exists specifically because codes
   disagree here. [DIGEST] §3.4 already verified the B&C factor of −1 numerically
   against the Casimir function — that check must be re-run against whatever
   convention the calculation returns.
5. **Origin.** The EFG and hyperfine operators are centred on the **Th nucleus**
   and are origin-independent. The molecule-frame dipole is *not*; if `d_mf` is
   reported as a cross-check, give the origin (Petrov 2025: centre of nuclear
   mass). ThF⁺ is a charged species, so this matters for `d_mf` only.
6. **Nuclear moments** used, and preferably the moment-free electronic factors
   (§3.3 item 3).

---

## 5. Can `eQq₂` be computed directly?

**Short answer: yes in principle, no as a keyword.** It is a research task, not a
production run. Petrov's admixture estimate is not the *only* route, but it is
currently the only *published* route, and the fact that the group that wrote down
Eq. (23) chose to estimate rather than compute it is itself evidence about how
routine it is not.

### 5.1 Why no expectation-value method can give it

Petrov 2018 Eq. (23) is
`eQq₂ = 2√6 eQ ⟨³Δ₁| Σ_i √(2π/5) Y₂₂/r_1i³ |³Δ₋₁⟩` — bra Ω = +1, ket Ω = −1.

Every predefined EFG facility computes `⟨Ψ|Ô|Ψ⟩` for one state: DIRAC `.EFG` and
`.NQCC` are listed under "Expectation values" in the manual's own taxonomy, and
`.RDCCDM` feeds them a single CC density matrix. There is no density matrix that
represents a transition between two components. So the `.EFG`/`.NQCC`/CFOUR-`PROP`
route gives `q₀` and **cannot** give `q₂`, no matter how good the correlation
treatment.

There is a second, deeper obstruction. Petrov 2018 §IV: "eQq₂ has no nonzero
matrix elements within a main nonrelativistic term ³Δ₁." `[derived]` The reason:
|³Δ_{Ω=+1}⟩ has (Λ, Σ) = (+2, −1) and |³Δ_{Ω=−1}⟩ has (−2, +1). The quadrupole
operator is purely spatial, so it requires ΔΣ = 0, but here ΔΣ = 2. The element
therefore vanishes identically in any **scalar-relativistic or spin-free**
treatment. **A spin-free calculation of `eQq₂` returns exactly zero and is not a
bug.** Whatever route is used must carry spin–orbit coupling explicitly (4c, X2C
with AMF/mmf, or 2c with the SO term retained), and must resolve the SO-admixed
Π character that supplies the nonzero value.

### 5.2 Routes that could work

**(a) KRCI / KRMCSCF with an off-diagonal property — the principled route.**
DIRAC's `.KRCI` (Fleig's module, the one Denis 2015 used as GASCI) constructs a
Kramers-restricted CI wavefunction in which both Ω = ±1 components live in the
same CI space. Evaluating a one-electron operator *between* two CI vectors is
standard CI machinery (it is what transition dipole moments are). Denis 2015
already computed "electric dipole transition moments between different electronic
states" this way. The pieces needed are: (i) the two ³Δ₁ components as separate
CI roots or as a degenerate pair, (ii) the q = ±2 EFG component as a one-electron
operator, which DIRAC can build from user-defined `.OPERATOR` combinations
(`one_electron_operators.rst` documents the `.OPERATOR` syntax with operator-type
keywords, factors and `COMFACTOR`), and (iii) a driver that contracts them.
`[inference]` — I have **not** verified that DIRAC exposes a "transition property
between KRCI roots for an arbitrary user operator" path; that is precisely the
question for Fleig. Note also that `Y_{2,±2} ∝ (x ± iy)²/r²`, so the required
operator is the `xx − yy` and `xy` combination of the same EFG integrals `.EFG`
already computes — the integrals exist; the contraction is what is missing.

**(b) `**PROPERTIES .ESR` — the nearest existing machinery, probably wrong
tensor.** Verbatim: "Evaluate ESR parameters — g-tensors and hyperfine coupling
tensors — using first-order quasi-degenerate perturbation theory based on
configuration interaction." QDPT over a degenerate CI manifold is *exactly* the
formalism that produces off-diagonal-in-Ω matrix elements. But the manual names
only magnetic tensors. Whether the same QDPT driver can be pointed at an electric
rank-2 operator is **UNVERIFIED** and would be a code question, possibly a small
patch. Worth asking; the infrastructure is 90 % there.

**(c) Finite-field symmetry breaking — the sneaky route.** `[inference]`, mine,
not in any source. Within the degenerate 2-D space {Ω = +1, Ω = −1}, an added
perturbation λÔ_{q=±2} with off-diagonal element V splits the pair by 2|V|. So:
add `λ·Ô(q = ±2)` to the Hamiltonian via `**HAMILTONIAN .OPERATOR` (the manual's
Fermi-contact example shows exactly this pattern, with `FACTORS 0.000000001`),
run a correlated method that resolves both components, and read `eQq₂` off the
linear-in-λ splitting. Attractions: it works with CCSD(T)-quality correlation, and
the manual's own worked examples are finite-field. Obstructions, all real:
(i) the perturbation is not totally symmetric, so the point group must be lowered
and the run gets much more expensive; (ii) a single-reference CC on one component
of a broken degeneracy is delicate — this wants KRMCSCF/CASPT2 or a
state-averaged reference; (iii) the manual warns "Don't forget to decrease the
symmetry of your system" for far milder perturbations. I would propose this only
as a cross-check on (a).

**(d) A better version of Petrov's own estimate — cheap and worth doing now.**
Keep Eq. (24) but replace its two guessed inputs with computed ones:
`⟨1/r³⟩_6d` for Th⁺ from an atomic DHF calculation (§4c), and `w` from the
*computed* CI composition of the ³Δ₁ wavefunction rather than back-fitted from
`G∥`. Denis 2015's GASCI wavefunctions contain exactly that composition
information. This does not remove the model dependence — in particular the
Petrov/Fleig disagreement about *which* Π state is admixed (§4c.2) — but it
would take `eQq₂` from "order of magnitude" to maybe ±30 %.

### 5.3 Recommendation on `eQq₂`

Do **not** put a real `eQq₂` calculation on the critical path. Do:

1. Run the atomic DHF `⟨1/r³⟩_6d(Th⁺)` now (§4c) — minutes, removes the largest
   guessed factor.
2. Carry `eQq₂` in `heff` as a **free parameter** with the 200–400 MHz range as
   the prior, and make the ²²⁹ThF⁺ level structure explicitly a function of it —
   Petrov 2018 §V does exactly this for HfF⁺ (their Figs. 1–3 plot the T,P-odd
   shifts *as functions of* `A∥`, `eQq₀`, `eQq₂`). That is the right shape for a
   quantity this uncertain, and it makes the eventual measurement a determination
   of `eQq₂` rather than a test of it.
3. Ask Fleig whether route (a) is a week or a year.

Also outstanding from [DIGEST] §3.2: the **normalisation bridge** between
Petrov's Eq. (23) `eQq₂` and B&C Eq. (9.52) at q = ±2 is still UNVERIFIED. That
is pure algebra, needs no collaborator, and should be settled before any computed
`eQq₂` is substituted into `heff` — otherwise a correct ab initio number gets
inserted with the wrong prefactor.

---

## 6. What I could not verify this session

Listed so absence is not mistaken for evidence.

1. **`diracprogram.org` was down.** All DIRAC keyword text is from the master
   branch of the git repo's manual source, not from a released version's docs.
   Keyword availability in an installed release must be re-checked.
2. **The DIRAC finite-nucleus keyword.** I confirmed point-nucleus is the default
   only for `.NONREL`, but did not find the keyword that selects the Gaussian
   model. Check `molecule.rst` / the `.mol` file format.
3. **DIRAC `.SPIN-ROTATION` for an open-shell Ω = 1 state.** Existence of the
   keyword is verified; applicability to ³Δ₁ is not. Also unverified: whether
   the tensor it returns equals B&C's `c_I` in Eq. (8.7) for a state with
   electronic angular momentum.
4. **CFOUR's `PROP` keyword and whether X2CAMF property machinery covers the
   EFG.** The CFOUR manual pages were unreachable; the EFG capability is
   evidenced only by the existence of the `xefgiso` module page.
5. **Flygare, Chem. Rev. 74, 653 (1974).** Title and year from a web search; the
   paper was not opened.
6. **Petrov 2018's `483` prefactor** being a pure unit conversion (my 483/234.9647
   = 2.056 consistency check is suggestive, not proof).
7. **Whether DIRAC can evaluate a user-defined one-electron operator between two
   KRCI roots** — the crux of §5.2(a).
8. **Current affiliations** of everyone in §3.2.
9. **No compute-cost measurement of any kind.** Every core-hour statement in §3.1
   is labelled inference and should be treated as such.
