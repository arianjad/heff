# Greedy quantum-logic purification of 232ThF+ X3Delta1 from 4 K (2026-09-18)

Index over the four library variants. Model: `heff.purify` (sideband pi pulses on the `heff.graph` transition graph, binary phonon click, classical belief vector, greedy expected-information-gain policy). Deck: `scripts/purify_thf.py`. Common settings: J_max=10, manifold J<=8 (320 states, 0.9992 of thermal population, H(b0)=6.953 bits), bw=0.001 MHz (Lamb-Dicke eta=0.09, Omega_ref/2pi=2 kHz, from Pipi 2026 Eq. A9), readout eta=1, dark=0, target max-belief 0.99, budget 60 cycles, 100 trajectories, seed 0. PLACEHOLDER polarizabilities (K=2 only, no scalar K=0 shift). B_z=1 G throughout.

Libraries: `six` = all six polarization pairs; `pi` = Pipi's E1=pi, E2=sigma+/-. `_closed` keeps only closed two-level pulses (no open Zeeman or diagonal ladders; see `pulse_library(closed=True)`).

## Cycles to max belief >= 0.99 (mean +/- s.e. / median / budget hits)

| library | E=0 V/cm | E=60 V/cm |
|---|---|---|
| six, open | 16.7 +/- 1.1 / 14 / 0 | 25.3 +/- 2.0 / 18 / 14 |
| pi, open | 19.6 +/- 1.4 / 16 / 1 | 32.9 +/- 2.2 / 31 / 27 |
| pi, closed | 25.5 +/- 1.7 / 20 / 5 | 32.9 +/- 2.2 / 32 / 25 |
| six, closed | 26.7 +/- 1.9 / 22 / 10 | 33.4 +/- 2.1 / 31 / 24 |

Budget hits = trajectories stopped at 60 cycles; those are counted as 60, so every E=60 mean (14-27 % censored) and the closed E=0 means are lower bounds on the true mean. Floor: H(b0) = 7 cycles (one <= 1-quantum mode reset carries <= 1 bit).

## Reading

- E=0, open six-pair library is the best case at ~17 cycles, median 14, no budget hits. The histogram is bimodal: a lump at 4-7 cycles (near the entropy floor) and a long tail to 50.
- Closing the ladders costs ~10 cycles at E=0: the degenerate Zeeman ladder is the main multi-state resource. At E=60 the closed and open libraries are indistinguishable within s.e.
- E=60 roughly doubles the count. The Stark shift (J=1 doublet split d*E = 3.37 D x 60 V/cm = 102 MHz) spreads lines so a 1 kHz window addresses fewer states; the first-pulse gain drops from 0.32 to 0.14-0.26 bits.
- Comparable scale to Pipi 2026: 26.5 pulses to purity 0.98 for H3O+ at 20 K.

## Files

- `README_<variant>.md`: per-variant stats, first pulse, most-used pulses.
- `cycles_E{0,60}_B1_<variant>.png`: histograms, dashed line = H(b0) floor.
- `graph_J12_E{0,60}_B1.png`: transition graph, J<=2 blocks.
- `run_<variant>.log`: stdout of each run (provenance).

Known model gaps and next steps: `docs/plans/2026-09-18-thf-purify-handoff.md`.
