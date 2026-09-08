# Audit of the ²²⁹ThF⁺ and ²²⁷ThF⁺ nuclear/hyperfine estimates

**Scope.** Source audit of the numerical estimates in `heff.params.thf_v2`
at worktree HEAD `5db8df7dce01707731752cd379f11acf1d8100c1`. This report
distinguishes a usable central estimate from a supported uncertainty bound. It
does not change numerical code or model metadata. Sources were checked through
2026-09-08; the search was targeted to the cited nuclear data, the two 2015
ThF⁺ electronic-structure calculations, the recent ThF⁺ calculation already in
the repository, and recent primary ²²⁹Th/²²⁷Th papers. It is not a claim that no
other calculation exists.

## Verdicts

| Quantity | Current value | Verdict | Supported action |
|---|---:|---|---|
| ²²⁹Th `g_N_Th` | `0.1464(24)` | Compatible, but superseded | Use `0.1460(12)`, derived from `mu = 0.365(3) mu_N` and `I = 5/2`. |
| ²²⁹ThF⁺ `A_par_Th` | `-1510(60) MHz` | Central value and negative sign remain usable; `60 MHz` is not a hard source-supported error | A source-informed update is about `-1.50(11) GHz`. If the central value is retained, use an uncertainty scale of about `110 MHz` and say that it is an ab-initio envelope, not a calibrated 1-sigma interval. |
| ²²⁷Th spin | tentative `I = 1/2` | Correctly marked tentative | Retain the tentative status. Do not infer a pure spherical `s_1/2` orbital from it. |
| ²²⁷Th `g_N_Th` | `-3.826` Schmidt | Arithmetic is correct under the stated free-neutron assumption; it is not an empirical estimate for this deformed actinide | Preserve it only as the prior user-selected stress test. A 2024 direct HFBCS alternative is `g_N = -0.1720`, but it has no quoted uncertainty. |
| ²²⁷ThF⁺ `A_par_Th` | `+39821 MHz` Schmidt | Same conditional verdict; no uncertainty bound is supported | Preserve only as the Schmidt stress test. The 2024 HFBCS alternative gives `+1790 MHz` using the same molecular electronic factor, still as a placeholder with no interval. |
| `E_eff` | `35.0(2.45) GV/cm` | Usable adopted ab-initio value | Keep. The central value follows the 35.2/37.3-GV/cm calculations and the error is 7% of the adopted value; document that it is an adoption, not a literal quoted datum. |
| `W_TP` | `50.0(3.5) kHz` | Source-supported | Keep. Skripnikov--Titov give 50 kHz and state a 7% theoretical uncertainty; Denis et al. give 48.35 kHz. |

## ²²⁹Th nuclear moment and nuclear g factor

Zitzer et al. measured the trapped-ion hyperfine structure of ²²⁹Th³⁺ and,
using atomic electronic factors, extracted in their Table I (printed p. 3)

```
mu(²²⁹Th ground state) = 0.365(3) mu_N
I = 5/2
```

Their stated uncertainty combines experimental and atomic-theory contributions.
The repository's definition is the usual nuclear factor

```
g_N = (mu/mu_N) / I,
```

so the update is

```
g_N = 0.365 / 2.5 = 0.1460,
sigma(g_N) = 0.003 / 2.5 = 0.0012.
```

Thus `0.1464(24)` is statistically compatible, but `0.1460(12)` is the better
current input. The definition is independently printed in Porsev and Safronova's
2026 atomic calculation, Eq./text preceding Table II. That paper reports the
reduced *atomic* constant `A_t(7s) = A/g = 40563(364) MHz`; it does not determine
a new nuclear moment and cannot replace the molecular ThF⁺ electronic factor.

Sources:

- G. Zitzer et al., [Phys. Rev. A **111**, L050802 (2025)](https://doi.org/10.1103/PhysRevA.111.L050802), Table I and printed pp. 3--4.
- S. G. Porsev and M. S. Safronova, [arXiv:2608.13700v1 (2026)](https://arxiv.org/abs/2608.13700), Sec. IV and Table II.

## ²²⁹ThF⁺ parallel hyperfine constant

### Rescaling and units

Skripnikov and Titov's Table II labels the `A_parallel` column in
`(mu_Th/mu_N) MHz` and gives the final electronic coefficient `-4163`.
Consequently this source scales with the magnetic moment itself:

```
A_S = -4163 * 0.365 = -1519.495 MHz
sigma_mu(A_S) = 4163 * 0.003 = 12.489 MHz.
```

Denis et al. give `A_parallel = +1833 MHz` after explicitly adopting
`mu/mu_N = 0.45`, again at `I = 5/2`. Holding their electronic calculation fixed,

```
A_D = +1833 * (0.365 / 0.45) = +1486.767 MHz
sigma_mu(A_D) = 12.220 MHz.
```

The two conventions use opposite molecular-axis directions, but a consistent
axis reversal changes both the body-fixed component and Omega, leaving
`A_parallel` invariant. The negative Skripnikov value is therefore not converted
into Denis's positive value by an axis change. The repository's complete sign
audit, including the HfF⁺ cross-check, supports the negative branch; Denis's
printed positive sign remains an unexplained outlier. See
[the sign-convention audit](../../lit/lookup-apar-th-sign-convention.md).

### What uncertainty is supported

Skripnikov and Titov state immediately before Table II that the theoretical
uncertainties of `E_eff`, `W_TP`, and `A_parallel` are within 7%. For the updated
Skripnikov value this is 106.4 MHz, already larger than the current 60-MHz field.
Combining that scale with the nuclear-moment uncertainty gives 107.1 MHz. A
simple negative-sign synthesis of the two magnitudes gives

```
mean magnitude = (1519.495 + 1486.767) / 2 = 1503.131 MHz
half-spread = 16.364 MHz
quadrature(7%, moment error, half-spread) = 107.2 MHz.
```

Rounded appropriately, `A_parallel = -1.50(11) GHz` is a defensible
source-informed estimate. The `11` should not be presented as a statistically
calibrated 1-sigma error: the primary source says “within 7%,” and Denis et al.
do not quote a comparably explicit error model. The current `-1510 MHz` central
value differs by only 7 MHz and need not move if preserving continuity matters;
the current `60 MHz` can describe an earlier spread estimate, but not a hard
uncertainty bound.

Sources:

- L. V. Skripnikov and A. V. Titov, [Phys. Rev. A **91**, 042504 (2015)](https://doi.org/10.1103/PhysRevA.91.042504), Sec. V and Table II (printed pp. 7--8).
- M. Denis et al., [New J. Phys. **17**, 043005 (2015)](https://doi.org/10.1088/1367-2630/17/4/043005), Secs. 2 and 4, Tables 6--7.

## ²²⁷Th: attribution, spin, and moment models

### What is established

Thorium has even `Z = 90`; ²²⁷Th has `N = 227 - 90 = 137`, so its odd nucleon
is a neutron. The live IAEA/ENSDF API extraction on 2026-09-08 reports the
ground state as parenthesized `(1/2+)`, with `(5/2+)` at `9.3(3) keV` and
`(3/2+)` at `24.38(3) keV`; both magnetic-dipole and quadrupole fields are blank.
The underlying evaluation has an ENSDF publication cutoff of 2016-01-15, so the
live endpoint is current access to an older evaluation, not a new 2026 spin
measurement.

Kovalík et al.'s 2021 conversion-electron measurement finds mixed M1+E2 character
for the 9.2-keV transition and says its rotational/Coriolis calculation prefers a
`1/2+, 3/2+, 3/2+` sequence over the adopted `1/2+, 5/2+, 3/2+` sequence. It
therefore reinforces that the low-level assignments and band interpretation are
not settled.

Sources:

- [IAEA LiveChart ground-state row for ²²⁷Th](https://nds.iaea.org/relnsd/v0/data?fields=ground_states&nuclides=227th) and [level rows](https://nds.iaea.org/relnsd/v0/data?fields=levels&nuclides=227th), extracted 2026-09-08.
- A. Kovalík et al., [Phys. Lett. B **820**, 136593 (2021)](https://doi.org/10.1016/j.physletb.2021.136593), abstract and Secs. 4--5.

### Why the Schmidt value is only a stress test

The Schmidt arithmetic is internally correct. The 2022 CODATA free-neutron
values are `mu_n/mu_N = -1.91304276(45)` and `g_n = -3.82608552(90)`; rounding
to `-1.913` and `-3.826` gives

```
mu_Schmidt = -1.913 mu_N
g_N = mu/(I mu_N) = -3.826
A_parallel = (-10408 MHz) * (-3.826) = +39821 MHz.
```

The invalid step is the inference that a tentative nuclear `I^pi = 1/2+`
bandhead forces the odd neutron into a pure spherical `s_1/2` state. In a
deformed actinide, `K = 1/2` is a projection quantum number. The intrinsic state
contains multiple Nilsson/spherical components, and Coriolis mixing and the
decoupling parameter affect the laboratory moment. Kovalík et al. explicitly
find strong Coriolis mixing and a decoupling parameter whose sign is sensitive
to mixing. Therefore the Schmidt value is useful for exercising a large-
hyperfine numerical regime, but it is not a sourced estimate of the actual
²²⁷Th moment and supports no error bar.

Source for the free-neutron constants: [2022 CODATA recommended values, NIST](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf), constants table.

### A direct 2024 theory alternative

The earlier repository claim that no measured **or estimated** ²²⁷Th moment
exists is now false. Minkov et al. performed blocked Skyrme-SIII HFBCS
calculations for neighboring actinides. Their reflection-unconstrained,
octupole-deformed ²²⁷Th solution makes the `1/2<+>` state the ground state, in
agreement with the adopted ground-state assignment, and Table IV (printed p. 10)
gives

```
mu(²²⁷Th, 1/2<+>) = -0.0860 mu_N
g_N = -0.0860 / (1/2) = -0.1720
A_parallel = (-10408 MHz) * (-0.1720) = +1790.176 MHz.
```

This is a more relevant central sensitivity point than a free neutron because
it includes deformation, pairing, core polarization, and spin-gyromagnetic
quenching. It is not a precision prediction. The paper quotes no uncertainty;
it says Coriolis mixing and collective coupling remain for future work. Its own
²²⁹Th calibration is mixed: the octupole solutions predict ground-state moments
`0.4181--0.4254 mu_N`, above the 2025 `0.365(3) mu_N`, while the
reflection-symmetric solution predicts `0.7366 mu_N`. Accordingly, the
`+1.790 GHz` value should remain `status="placeholder"`, `uncertainty=None`, and
its source/model assumptions should be carried in the note. It is a recommended
alternative to the prior Schmidt stress test, not an automatic replacement of
that user-selected test value.

Source: N. Minkov et al., [Phys. Rev. C **110**, 034327 (2024)](https://doi.org/10.1103/PhysRevC.110.034327), Tables II--IV and Secs. III--IV; [open manuscript](https://arxiv.org/abs/2408.11010).

## `E_eff` and `W_TP` spot check

Skripnikov and Titov's final Table II values are `E_eff = 37.3 GV/cm` and
`W_TP = 50 kHz`, with the same stated 7% theoretical-uncertainty scale used
above. Denis et al. report `E_eff = -35.2 GV/cm` in their signed axis convention
and `W_PT = 48.35 kHz`. The effective-Hamiltonian registry uses the adopted
magnitude `E_eff = 35.0 GV/cm`, consistent with the latter central value and
within 7% of the former, and `W_TP = 50 kHz` directly from Skripnikov--Titov.
The stored uncertainties `2.45 GV/cm` and `3.5 kHz` are exactly 7% of the adopted
central values. Both are reasonable ab-initio inputs as long as the `E_eff`
metadata says that 35.0 is an adopted value rather than a verbatim
Skripnikov--Titov result.

The recent ThF⁺ paper checked in this audit, Petrov and Skripnikov
[arXiv:2503.02840](https://arxiv.org/abs/2503.02840), concerns electric-field-
dependent molecular g factors for ²³²ThF⁺. It supplies no updated odd-Th nuclear
moment or `A_parallel(Th)`. The 2026 Porsev--Safronova paper is atomic Th³⁺, as
described above. This is a bounded statement about the recent papers read, not
a global absence claim.

## Concrete metadata corrections implied by the evidence

1. Cite Zitzer et al. 2025 for `mu(²²⁹Th) = 0.365(3) mu_N`; define `g_N` explicitly
   as `(mu/mu_N)/I` and update it to `0.1460(12)`.
2. Describe ²²⁹ThF⁺ `A_par_Th` as an ab-initio estimate with the negative sign
   selected by the sign audit. Replace a hard `60 MHz` interpretation with the
   roughly `110 MHz` source-informed theoretical scale, or remove the numeric
   uncertainty and state the 7% envelope in the note.
3. Remove every maintained claim that no ²²⁷Th moment estimate exists. The
   narrower verified statement is: the 2026-09-08 IAEA/ENSDF extraction has no
   evaluated magnetic moment, while Minkov et al. 2024 give a model value.
4. Describe the ²²⁷Th Schmidt values as a prior user-selected free-neutron stress
   test. Remove the assertion that `(1/2+)` forces a pure spherical `s_1/2`
   orbital. Record `mu = -0.0860 mu_N`, `g_N = -0.1720`, and
   `A_parallel = +1.790 GHz` as the direct-HFBCS alternative, with no uncertainty
   and no silent change to the selected stress-test default.
5. Keep `E_eff = 35.0(2.45) GV/cm` and `W_TP = 50.0(3.5) kHz`; preserve their
   ab-initio status and the distinction between an adopted central value and a
   literal source value.

The unresolved experimental question is the ²²⁷Th ground-state magnetic moment;
the unresolved molecular question is the unexplained sign disagreement in the
two 2015 ²²⁹ThF⁺ `A_parallel` calculations. Neither uncertainty is repaired by
assigning a precise error bar to a placeholder.
