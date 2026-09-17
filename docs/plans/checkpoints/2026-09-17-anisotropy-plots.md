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

## Milestone 2 (pending): both cases run

## Milestone 3 (pending): comparison figures

## Milestone 4 (pending): sanity checks, tests, README
