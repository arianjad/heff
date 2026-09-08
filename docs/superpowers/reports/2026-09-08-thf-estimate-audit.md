# ThF+ parameter estimate audit

Status: complete; scientific limits and final regression result recorded below. Audit date: 2026-09-08. Starting implementation: `5db8df7`.

The missing 232Th TOML provenance is repaired. **The audit does not validate all
existing defaults as precision predictions.** It separates source-supported
inputs from transferred estimates and uncalibrated sensitivity placeholders.
All central values, units, and physics matrix elements are unchanged.

## Findings and applied fixes

| Input | Finding | Applied change |
|---|---|---|
| 232Th TOML records | Native `Param` metadata was missing | Preserve full value, unit, status, uncertainty, source, isotope, convention, and note |
| B0 | Division of the reported interval by four propagated its error incorrectly by a factor ten | Uncertainty corrected from 0.001 to 0.010 MHz; status `derived`; central value unchanged |
| Shared odd-isotope constants | 232Th measurement errors do not quantify isotope-transfer errors | Label B0, D0, omega_ef, A_par, d_mf, G_par as transferred estimates; retain original errors in notes, leave target uncertainty unquantified |
| F spin rotation | CsF scaling gives 19.94 kHz, but establishes neither the ThF sign nor a factor-three bound | Retain 20 kHz as a sensitivity estimate; remove the unsupported bound claim |
| 229Th magnetic hyperfine | Approximately -1.5 GHz is a usable theory scale; the old 60 MHz spread omits the cited 7% theoretical uncertainty | Retain -1510 MHz and sign selection; remove the incomplete numeric uncertainty and state the roughly 106 MHz theory scale in its note |
| 229Th nuclear g | 0.1464(24) is compatible with the newer 0.1460(12) extraction | Retain the documented older input; cite the newer value in the note |
| 229Th quadrupoles | Hf-to-Th electronic transfer and signed normalization are unvalidated for both q=0 and q=2 | Keep -2600 and +300 MHz as `placeholder`; remove unsupported +/-1000 and +/-100 MHz intervals |
| 227Th moment/hyperfine | Free-neutron Schmidt arithmetic is conditional on an unestablished spherical model | Preserve the user-selected stress test; correct the orbital inference and record the direct 2024 nuclear-theory alternative |

Th spin rotation held at zero remains missing input, not evidence of a small
interaction. The adopted E_eff and W_TP values remain reasonable ab-initio
inputs with their documented theory-error scales. Two-photon placeholders were
not converted into physical predictions by this audit.

## The newly located 227Th theory alternative

[Minkov et al., PRC 110, 034327 (2024), Table IV](https://arxiv.org/abs/2408.11010)
predicts a magnetic moment of -0.0860 nuclear magnetons for its octupole-deformed
1/2 ground solution. Parent inspection of the rendered table confirms the row.
For I=1/2 and the existing molecular electronic factor:

```text
g_N = -0.0860 / 0.5 = -0.1720
A_parallel = (-10408 MHz) * (-0.1720) = +1790.176 MHz
```

This is about +1.79 GHz, versus the retained Schmidt default of +39.821 GHz.
The paper gives no calibrated uncertainty and leaves Coriolis/collective mixing
for future work. It supplies a more relevant nuclear-model sensitivity point,
not a measured moment or a precision prediction. The previous claim that no
published moment estimate exists is corrected in the lookup, digest, open
questions, Hamiltonian note, and parameter metadata. A deformed nuclear I=1/2
does not imply a pure spherical s1/2 neutron.

## Shared-input source checks

The canonical Zotero Brown and Carrington PDF (storage CKZKCGXY, 1045 pages)
was overviewed before targeted text and rendered-page reads. Table 8.12,
printed p. 481 / PDF p. 513, gives CsF F-spin rotation 15.1 kHz and
B=0.183782 cm^-1. Scaling by the current ThF rotational constant gives
19.936385 kHz. Printed p. 421 / PDF p. 453 discusses opposite-sign first- and
second-order contributions. These sources support the analogy's arithmetic,
not its accuracy or sign for ThF+.

B&C Eq. (7.199), printed p. 345 / PDF p. 377, gives the leading Dunham mass
scaling. Using mass numbers solely to estimate the size of the effect,
B proportional to reduced-mass^-1 gives shifts of +7.213712 MHz (229) and
+12.128782 MHz (227) from 232. D proportional to reduced-mass^-2 gives
+0.007733 and +0.013006 kHz. These diagnostic values were not substituted for
an isotope fit; vibrational and Born-Oppenheimer corrections remain relevant.

[Ng et al. 2022, Sec. II B, p. 2](https://arxiv.org/abs/2202.01346)
reports the J=2-to-J=1 interval 29.09733(4) GHz. Division by four gives
7274.3325 MHz with uncertainty 0.010 MHz. However, the package includes an
independent centrifugal term: its rotation-plus-centrifugal gap is
4B-32D=29097.205296 MHz, 124.704 kHz below the reported interval. B0 is therefore
an interval-derived approximation, not a joint fit of all included terms.
Adding 8D alone would repair that two-term interval but would not constitute
a complete spectroscopic refit. No such refit was silently performed.

## Independent review and validation

Full source analyses:
- [Nuclear moments and magnetic hyperfine](2026-09-08-thf-nuclear-estimate-audit.md).
- [Electric quadrupole transfer and normalization](2026-09-08-thf-quadrupole-estimate-audit.md).

Parent source checks included Minkov Table IV, Petrov 2018 Eqs. (22)-(25),
B&C Eqs. (4.30)-(4.32) and (9.52)-(9.53), and the shared-input pages above.
The report review corrected an equation-number attribution and clarified that
a negative inferred admixture weight invalidates the single-channel model;
it is not a physical negative probability or a predicted quadrupole sign.
No unresolved conversion factor or phase was selected.

Regression-first checks reproduced the missing metadata, incorrect B0 error,
and unsupported uncertainty fields before their fixes. Intermediate checks:
32 model/isotope tests; 52 focused tests; 380 passed and 2 existing skips before
the final hyperfine-uncertainty metadata correction.
Final full regression: **381 passed, 2 existing skips in 67.54 seconds**.
`git diff --check` also passed.
A direct comparison against `5db8df7` confirms all native 232/229/227 parameter
central values and units are unchanged. Regression tests establish software
consistency, not the physical truth of a placeholder.

Local checkpoints preserve the fixes and source reports. No push, merge,
pull request, installation, or other external write was performed.
