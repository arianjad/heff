# ThF+ isotope level and field plots

Complete, 2026-09-08. [Open the five-page PDF](thf-isotopes-levels-stark-zeeman.pdf).

The figures contain the zero-field structure and shifts of every state
correlated with J=1–3: 60 states for 232Th, 360 for 229Th, and 120 for 227Th,
including magnetic degeneracy. Stark scans cover 0–10 kV/cm at B=0; Zeeman
scans cover 0–100 G at E=0. They are separate scans, not simultaneous fields.

## Figures

1. [Rotational overview](thf-level-overview.png): each isotope has its own energy zero.
2. [232Th detailed panels](thf-232-baseline.png): rows J=1,2,3; columns zero-field, Stark, Zeeman.
3. [229Th baseline panels](thf-229-baseline.png): electric quadrupole deliberately omitted.
4. [227Th theory panels](thf-227-baseline.png): the published deformed-nucleus moment replaces the free-neutron stress test in this plotting set.
5. [229Th quadrupole sensitivity case](thf-229-legacy-quadrupole.png): unvalidated coefficients, not a prediction or an uncertainty envelope.

## Adopted parameters and limits

These are separate exploratory parameter sets. Package defaults were not changed.
All values and provenance are exported in `parameters.csv` and the per-case JSON files.

| Input | 232Th | 229Th baseline | 227Th theory |
|---|---:|---:|---:|
| B0 (MHz) | 7274.3325 | same, unscaled transfer | same, unscaled transfer |
| D0 (kHz) | 3.897 | same, unscaled transfer | same, unscaled transfer |
| body-fixed dipole, center-of-mass origin (D) | 3.37 | same transfer | same transfer |
| A_parallel(Th) (MHz) | absent, I=0 | -1519.495 | +1790.176 |
| g_N(Th) | absent, I=0 | +0.1460 | -0.1720 |
| eQq0(Th), eQq2(Th) (MHz) | absent, I=0 | 0,0: explicitly omitted | absent, I=1/2 |

229Th uses mu=0.365(3) nuclear magnetons from
[Zitzer et al. 2025, Table I](https://doi.org/10.1103/PhysRevA.111.L050802),
with the negative electronic hyperfine coefficient -4163 MHz per nuclear
magneton from [Skripnikov and Titov 2015, Table II](https://arxiv.org/abs/1503.01001).
The cited 7% theoretical scale remains much larger than numerical errors.

227Th uses mu=-0.0860 nuclear magnetons from the octupole-deformed 1/2 ground
solution in [Minkov et al. 2024, Table IV](https://arxiv.org/abs/2408.11010).
The tentative I=1/2 assignment and nuclear-model uncertainty remain unresolved;
no calibrated moment error is quoted. A=(-10408 MHz)*g_N transfers the existing
molecular electronic factor. These figures use +1.79 GHz, not the +39.8 GHz
Schmidt stress-test default.

229Th has **no validated complete parameter set here**: its quadrupole terms
may substantially alter its spectrum. The baseline sets them to zero to show
what the supported magnetic-hyperfine model produces; zero is not their
physical estimate. The separate sensitivity figure uses -2600 and +300 MHz
as raw B&C-parameter scenarios. Neither signs, conversion, nor magnitudes are
validated. That single scenario is not a bound on the unknown quadrupole effect.

All three retain F spin rotation 20 kHz as an analogy-based estimate. Th spin
rotation is held at zero for missing input. Shared B, D, dipole, omega-doubling,
F hyperfine, and electronic g inputs are unscaled isotope transfers; the plots
are exploratory structure predictions, not precision spectroscopy fits.
See [the parameter audit](../../docs/superpowers/reports/2026-09-08-thf-estimate-audit.md).

## Reading the figures

J0 is the dominant J at zero field, not an exact conserved label at finite E.
The minimum zero-field parent-J weight is 0.9935 (229), 0.9992 (227), and
0.9999999 (232). Left-column energies subtract the rotational reference
E_rot/h=B0*J*(J+1)-D0*[J*(J+1)]^2, with D converted to MHz.
Blue/orange marks distinguish physical parity +1/-1; nearby doublets can be
unresolved on the GHz-scale thorium hyperfine axes. CSV files preserve the
individual numerical levels. F is the total angular momentum including both
nuclear spins, where present.

Stark shifts subtract each branch's own zero-field energy. Branches follow
stepwise eigenvector character using Hungarian assignment, with adaptive
bisection when a plotted state's squared overlap drops below 0.90. This is a
labeling procedure, not an assertion of dynamical adiabaticity. The minimum
accepted overlap is 0.9006. Positive/negative mF Stark partners are degenerate
at B=0 and drawn once. At exact degeneracies, only the eigenspace is unique.

Zeeman curves are energy ordered separately within conserved (signed mF, parity)
blocks, referenced to the same rank at B=0. Negative-mF curves are dashed.
The implementation obtains them from positive mF at negative B and independently
checks the time-reversal spectral identity using directly constructed negative-mF
matrices at three field endpoints, including simultaneous E and B.

## Verification and reproduction

The basis includes J=1–8. For each plotted mF block, embedded J<=7 eigenvectors
were overlap matched against the actual tracked J<=8 subspace at 23 electric
fields spanning the range and at B=0,1,10,100 G. Maximum energy differences:
0.0138 Hz (232), 0.0962 Hz (229 baseline), 0.0201 Hz (227); the unvalidated
quadrupole scenario remains below 0.93 Hz. These are sampled cutoff differences,
not certified bounds between samples or physical error bars. Time-reversal
spectral residuals are below 0.002 Hz. Parity commutation and Stark
anticommutation checks passed. See `validation.json` and per-block JSON records.

From the repository root, after installing `.[plot]`:

```shell
python scripts/plot_thf_isotopes.py
```

`--render-only` regenerates figures/tables from the saved NPZ and JSON data.
The calculation builds term matrices from the current source; no machine-local
checkpoint or disk matrix cache is required. Raw curve arrays and plotting
parameter sets are committed. The script uses existing native heff operators
and NumPy diagonalization; no package matrix element is replaced.
No GPU, installation, external write, or package-default update was needed.
