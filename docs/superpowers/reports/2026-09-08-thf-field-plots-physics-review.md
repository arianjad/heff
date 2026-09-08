# Independent physics review: ThF+ isotope field plots

**Review date:** 2026-09-08  
**Reviewed revision:** `6213473` plus the uncommitted plotting/results work  
**Scope:** `scripts/plot_thf_isotopes.py` and `results/thf-fields-2026-09-08/`; no package matrix elements or native defaults were changed by this plotting work.

## Verdict

The calculation is suitable as an **exploratory effective-model visualization**. Each portable odd-isotope figure/PDF page now states that shared molecular constants are unscaled transfers from 232ThF+, that the unknown Th spin-rotation term is omitted, and that no physical uncertainty band is shown. The overview page carries the same scope boundary and also identifies the omitted 229Th quadrupole interaction.

The field units, symmetry blocks, state counts, character tracking, and revised sampled J-cutoff comparison are internally consistent. The numerical checks do not establish physical accuracy. In particular, 229Th quadrupole constants and both odd-isotope Th spin-rotation constants remain unknown here; the 227Th nuclear moment is a model result without a calibrated uncertainty; and shared molecular constants have not been isotope-scaled or refit.

## Source procedure

I used the existing estimate/nuclear/quadrupole reports only as maps. I checked the load-bearing claims against the local primary PDFs or primary paper text. PDF metadata and tables of contents were inspected first for the long Brown and Carrington and Jadbabaie (2025) sources; targeted pages were then extracted. The Brown and Carrington PDF has only page-number bookmarks rather than a semantic TOC, so its chapter map came from the local `brown-carrington` skill and the printed contents. Petrov et al. (2018) is a five-page paper; its metadata and relevant pp. 2–3 were read directly. APS/arXiv primary versions were used for Zitzer et al. (2025), Skripnikov and Titov (2015), and Minkov et al. (2024).

## Claim-by-claim review

[claim 1] The 229Th plotting values `g_N=+0.1460` and `A_parallel(Th)=-1519.495 MHz` follow from the cited inputs.
  verdict: OK
  source: Zitzer et al., *Phys. Rev. A* **111**, L050802 (2025), Table I, printed pp. 2–3; Skripnikov and Titov, *Phys. Rev. A* **91**, 042504 (2015), Sec. V and Table II, printed pp. 7–8.
  reason: Zitzer reports `mu=0.365(3) mu_N` for `I=5/2`, hence `g_N=mu/I=0.1460(12)`. Skripnikov and Titov report the ThF+ factor `A_parallel=-4163 (mu_Th/mu_N) MHz`; multiplication gives `-1519.495 MHz`. Their quoted roughly 7% electronic-structure scale is about 106 MHz, much larger than the 12.5 MHz moment-only uncertainty, so the exported extra digits are arithmetic traceability rather than physical precision.

[claim 2] The 229Th baseline with `eQq0=eQq2=0` is an acceptable quadrupole-omitted reference, but zero is not a physical estimate.
  verdict: OK
  source: Brown and Carrington, Sec. 9.4.3, Eqs. (9.52)–(9.53), printed pp. 604–605; Petrov et al., *Phys. Rev. A* **98**, 042502 (2018), Sec. II, Eqs. (19), (22)–(25), printed pp. 2–3.
  reason: The quadrupole interaction is an independent rank-2 hyperfine term and can mix `Delta J=0,+/-1,+/-2`. No validated 229ThF+ values or convention conversion were supplied. The script and README explicitly call zero an omission and separate the legacy coefficients as an unvalidated sensitivity case, which is the defensible interpretation.

[claim 3] Holding `c_I(Th)=0` supplies a reasonable numerical baseline.
  verdict: OK
  source: No loaded ThF+ source supplies an odd-Th spin-rotation constant; `results/thf-fields-2026-09-08/README.md`, “Adopted parameters and limits.”
  reason: Zero is usable as an operator switch-off for a missing input, not as a physical estimate, and no uncertainty envelope follows from it. That required context now appears in the README, parameter export, odd-isotope figure footers, and overview footer.

[claim 4] The 227Th plotting values `g_N=-0.1720` and `A_parallel(Th)=+1790.176 MHz` implement the cited deformed-nucleus model.
  verdict: OK
  source: Minkov et al., arXiv:2408.11010 version of *Phys. Rev. C* **110**, 034327 (2024), Table IV, PDF p. 10, and Sec. V, PDF pp. 11–12.
  reason: Their reflection-asymmetric solution lists a `1/2+` ground state with `mu=-0.0860 mu_N`, giving `g_N=-0.1720`; applying the repository's existing electronic factor `-10408 MHz` gives `+1790.176 MHz`. The paper states that Coriolis mixing and collective coupling are needed for a more realistic quantitative description, so this is correctly labeled a model estimate with unresolved uncertainty.

[claim 5] The copied `B`, `D`, dipole, omega doubling, F hyperfine, and electronic-g inputs are unscaled transfers rather than a precision isotope prediction.
  verdict: OK
  source: Brown and Carrington, Sec. 7.5.4, Eqs. (7.198)–(7.199), printed p. 345.
  reason: Even the leading Dunham coefficients scale with reduced mass, `Y_kl ~= mu^{-(k+2l)/2} U_kl`, with further isotope-dependent corrections. Other transferred constants can also have isotope or nuclear dependence. The README and regenerated figure/PDF footers now explicitly call these unscaled transfers, which correctly limits the claim.

[claim 6] The requested field ranges and plotted units are implemented consistently.
  verdict: OK
  source: Ng et al., *Phys. Rev. A* **105**, 022823 (2022), Table I, printed p. 5; Brown and Carrington, Sec. 5.5.6, Eqs. (5.182)–(5.186), printed pp. 174–175, and Sec. 9.4.4, Eqs. (9.54)–(9.55), printed pp. 605–606.
  reason: The script scans `0..10000 V/cm` and `0..100 G`, forms the longitudinal rank-1 Stark and Zeeman matrices in MHz, and converts only display axes to kV/cm/GHz. With Ng's center-of-mass dipole `3.37(9) D`, `dE/h=3.37*0.5034118*10000=16965 MHz` at the endpoint, matching the implementation's scale. The dipole origin is explicitly identified, which matters for an ion.

[claim 7] The plotted set contains every magnetic sublevel correlated with zero-field `J=1,2,3`.
  verdict: OK
  source: Brown and Carrington, Sec. 5.5.5, Eqs. (5.172), (5.174)–(5.176), printed pp. 173–174; `scripts/plot_thf_isotopes.py`, lines 154–184 and 312–319; `results/thf-fields-2026-09-08/validation.json`.
  reason: The coupled basis is partitioned by conserved `m_F`; using `m_F=1/2,3/2,...` for 232Th and `m_F=0,1,...` for the integer-F odd isotopes reaches every allowed F manifold. The code asserts equality between selected dominant-J parents and the primitive `J<=3` dimension and restores negative-m degeneracy in the exported count. The resulting totals are 60, 360, and 120 states.

[claim 8] At `B=0`, drawing only nonnegative `m_F` Stark curves does not omit distinct energies, and the generated negative-`m_F` Zeeman spectra obey the stated time-reversal identity.
  verdict: OK
  source: Jadbabaie, *Measuring Fundamental Symmetry Violation in Polyatomic Molecules* (2025), App. A.2.2, Eqs. (A.16)–(A.20), printed pp. 276–277; Petrov et al. (2018), Sec. II, printed p. 2; `scripts/plot_thf_isotopes.py`, lines 156–164 and 187–217.
  reason: Time reversal maps angular-momentum projections to their negatives, while an electric field is T-even and a magnetic field changes sign. The script additionally constructs negative-m matrices directly at `(E,B)=(0,100),(10000,0),(10000,100)` and finds spectral residuals of order `1e-9 MHz`. Parity commutes with the zero-field/Zeeman Hamiltonians and anticommutes with the Stark operator in the implemented basis, so parity blocking is applied only where appropriate.

[claim 9] Stark branch labels preserve zero-field character with controlled numerical tracking.
  verdict: OK
  source: `scripts/plot_thf_isotopes.py`, lines 90–121 and 169–186; per-case JSON block diagnostics.
  reason: Each point contains exact eigenvalues; the ordering is assigned by successive squared-overlap Hungarian matching, with adaptive bisection until plotted-state step overlaps are at least 0.90. The observed minimum is 0.9006. This supports character-continuation labels, not physical adiabatic evolution; the README states this and notes that only the eigenspace is unique at exact degeneracy.

[claim 10] The `J<=8` cutoff is numerically adequate for the displayed `J=1..3` branches at the sampled validation points.
  verdict: OK
  source: `scripts/plot_thf_isotopes.py`, lines 96–144 and 169–177; `results/thf-fields-2026-09-08/validation.json`.
  reason: Embedded `J<=7` eigenvectors are overlap-matched to actual tracked `J<=8` Stark snapshots at 23 fields spanning `0..10000 V/cm`, plus Zeeman checks at `B=0,1,10,100 G`. The maximum sampled shift is 0.0962 Hz for the supported baselines and 0.93 Hz for the unvalidated sensitivity case, far below the 1 kHz numerical target. This is a sampled cutoff diagnostic, not a bound between samples or a physical uncertainty; the README now says so.

[claim 11] Zeeman curves are energy ordered within conserved `(m_F, parity)` blocks and do not claim continuous character tracking across crossings.
  verdict: OK
  source: `scripts/plot_thf_isotopes.py`, lines 187–203 and 245–264; `results/thf-fields-2026-09-08/README.md`, “Reading the figures.”
  reason: The code deliberately energy-orders eigenvalues independently at every field within each conserved `(m_F, parity)` block. The figure annotation and README say this explicitly, so a same-color line through an exact crossing does not assert a character-tracked state identity.

## Acceptance boundary

I find no remaining load-bearing physics or numerical blocker for delivery as an exploratory plotting set. The 229Th quadrupole sensitivity page retains its **UNVALIDATED** label and is not described as a bound or uncertainty envelope. Quantitative isotope predictions remain contingent on measured or independently validated odd-isotope quadrupole/spin-rotation inputs and isotope-specific molecular constants.

summary: 11 total claims; 11 OK, 0 WRONG, 0 NEEDS-CONTEXT, 0 STOP
