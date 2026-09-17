# Checkpoint: two-photon figures at two polarizability anisotropies (2026-09-17)

Task: produce the 232ThF+ X3Delta1 two-photon transition-matrix figure set for
a_xx - a_yy = 0 and a_xx - a_yy > 0, and a per-B-panel comparison figure.
Numerical deck unchanged from the 2026-09-15 approval (J_max = 7, J = 1..5,
E_z = 0, B panels 0.001/1/3/5/10 G, sweep 0.001..20 G at 201 points, six
polarization pairs). All alpha values are PLACEHOLDERS; the closure form is
valid only for detuning >> ~7 GHz, which is not the JILA regime (OPEN-23).

## Milestone 1 (done): script parametrized, default unchanged

`scripts/plot_thf_twophoton.py` gained `--case {default,xxyy0,xxyy_pos}` mapping
to (alphas, output directory):

| case | alphas | output |
|---|---|---|
| default | alpha_K2_dOm0 = 1, alpha_K2_dOm2 = 1 | results/thf-twophoton-2026-09-15 |
| xxyy0 | alpha_K2_dOm0 = 1 only | results/thf-twophoton-2026-09-17-anisotropy/xxyy0 |
| xxyy_pos | alpha_K2_dOm0 = 1, alpha_K2_dOm2 = 1 | results/thf-twophoton-2026-09-17-anisotropy/xxyy_pos |

Verified: `--help` reports `default=default`, and the `CASES["default"]` entry
carries the original alphas and the original output path, so a no-argument run
still reproduces the 2026-09-15 configuration. The default case was NOT run;
`results/thf-twophoton-2026-09-15/` is left untouched.

Probe (2026-09-17): omitting `alpha_K2_dOm2` and setting it to 0.0 give
bit-identical amplitudes (max abs difference 0.0 at B = 1 G, (sigma+, pi)).
The xxyy0 case OMITS the knob, so `TwoPhotonOperator._keys()` emits no
(2, 2, P) key at all and `channels.npz` carries only the five (2, 0, P)
matrices. That is the direct evidence the dOmega = 2 channel is absent.

## Milestone 2 (done): both cases run

`results/thf-twophoton-2026-09-17-anisotropy/{xxyy0,xxyy_pos}/`, 48 s and 39 s
wall clock. Each directory holds 5 heatmap grids, 5 matrix npz, 5 label CSV,
6 curve PNG, 6 sweep npz, channels.npz and its own README.

## Milestone 3 (done): comparison figures

`scripts/compare_thf_twophoton_anisotropy.py` writes `compare_B*.png`, six panels
each, background |M|^2 from xxyy_pos via `heff.plot_transition.heatmap`, overlay
marking cells the |dOmega| = 2 pathway opens, closes, enhances >2x or suppresses
<0.5x. Counts printed by the same code that draws the markers.

## Milestone 4 (done): sanity checks, tests, README

All numbers are in `results/thf-twophoton-2026-09-17-anisotropy/README.md`.
Headlines:

- xxyy_pos reproduces `thf-twophoton-2026-09-15/matrices_B*.npz` exactly, max abs
  amplitude difference 0.000e+00 at all five B panels. Control: xxyy0 against the
  same file differs by 0.4903892810042282 at B = 1 G, so the check can fail.
- The (2, 2, P) channel matrices have support only on |dOmega| = 2 basis pairs and
  the (2, 0, P) matrices have exactly zero there, so omitting the knob removes the
  whole pathway.
- Two corrections to the task's framing, both measured:
  1. "Dominant Omega" is undefined at E_z = 0. Every selected eigenvector carries
     population exactly 0.500000 on Omega = +1, being an exact parity eigenstate.
     The operator-level statement (channel support) replaces it.
  2. The 1e-12-RELATIVE threshold is below the float64 noise floor of these
     matrices and invents 244 "opened" and 282 "closed" cells that sit at
     ~1e-13 amplitude. At an absolute 1e-10 floor the B = 1 G support is
     identical between the two cases, 0 opened and 0 closed.
- The real effect of a_xx - a_yy is interference, not new lines: about 80% of
  populated cells move by more than 2x in |M|^2, split about evenly between
  enhancement and suppression. A handful of cells (up to 20 of ~900 per pair) do
  open at fields other than 1 G.
- Channel amplitudes add linearly, residual exactly 0.
- `pytest tests/test_transition.py tests/test_twophoton.py -q` -> `22 passed in 7.92s`.
- Two figures inspected directly (`compare_B0.001G.png`,
  `xxyy0/heatmaps_B1G.png`): axes labelled, J boundaries drawn, six panels
  populated, Delta m_F in panel titles, placeholder and closure caveats in the
  suptitle.

Two README claims were written wrong on the first pass and corrected after
measuring: the 1e-10 floor's margin (it is 1.7-1.9 decades below the smallest
element it keeps, not nine), and the sigma+/sigma- relation (the matrices are not
elementwise equal; they are related by m_F -> -m_F, exactly only at zero field,
with a residual linear in B).

`results/thf-twophoton-2026-09-15/` was not touched, and the default `--case`
was never run.
