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

## Follow-up 2026-09-18 (evening): budget 200, B sweep, propagator model

All 100 trajectories, seed 0, budget 200 cycles (no censoring unless noted). `run()` now samples the joint (final state, phonon) table, so these are fresh realizations, not the budget-60 trajectories extended. Sibling directories: `../thf-purify-2026-09-18-budget200/`, `../thf-purify-2026-09-18-bsweep/B<G>/`, `../thf-purify-2026-09-18-prop/` (`B<G>/` for E=0 at other fields). Each holds `README_<variant>.md`, histograms and run logs as above; graph images are omitted where they would repeat this directory's.

### Budget 200, pi-pulse model, B=1 G (mean +/- s.e. / median / hit budget)

| library | E=0 V/cm | E=60 V/cm |
|---|---|---|
| six, open | 15.9 +/- 1.2 / 12 / 0 | 25.7 +/- 2.3 / 18 / 0 |
| pi, open | 18.4 +/- 1.5 / 15 / 0 | 49.7 +/- 4.7 / 34 / 3 |
| pi, closed | 25.7 +/- 2.3 / 17 / 0 | 50.1 +/- 4.4 / 42 / 1 |
| six, closed | 26.2 +/- 2.7 / 16 / 0 | 33.5 +/- 3.1 / 22 / 0 |

The budget-60 E=60 pi-library means (32.9) were censored by ~25 %; uncensored they are ~50. The six-pair E=60 numbers were not (25.7, 33.5 vs 25.3, 33.4).

### B sweep, pi-pulse model, E=0, six pairs

Zeeman rung anharmonicity |f(m+1->m+2) - f(m->m+1)| over the J<=8 manifold from heff at E=0: median / thermal-population-weighted mean, and the population fraction of rungs above eta*Omega = 180 Hz and above the 1 kHz window. "Open pulses" = open minus closed library size (Zeeman ladders plus the 2 diagonal light-shift drives over ~310 states).

| B (G) | anharm. Hz (median / wtd mean) | > 180 Hz / > 1 kHz | open pulses | open | closed |
|---|---|---|---|---|---|
| 1 | 0.18 / 3.6 | 0.00 / 0.00 | 58 | 15.9 +/- 1.2 | 26.2 +/- 2.7 |
| 3.6 | 2.3 / 46 | 0.13 / 0.00 | 90 | 15.8 +/- 1.1 | 27.4 +/- 2.4 |
| 10 | 18 / 360 | 0.39 / 0.13 | 146 | 22.0 +/- 1.8 | 29.9 +/- 3.0 |
| 30 | 170 / 3200 | 0.95 / 0.39 | 312 | 27.7 +/- 2.3 | 41.1 +/- 3.1 |
| 100 | 2000 / 36000 | 1.00 / 1.00 | 10 | 29.8 +/- 2.6 | 45.7 +/- 4.3 |

Anharmonicity scales as B^2 (second-order Zeeman through the 5.3 MHz Omega doublet and 15 MHz hyperfine gap at J=1). The ladders close, in the 1 kHz-window sense, only near 100 G; at Pipi's 3.6 G nothing changes. Raising B costs cycles in both libraries because resolving the Zeeman lines removes the degeneracy that lets one pulse address many states. At 100 G only 10 open pulses remain (8 ladders, 2 diagonal drives) and the 16-cycle open-closed gap is the diagonal light-shift drive, not the Zeeman ladder. Hypothesis (handoff gap 5) falsified. At E=60 the m levels are Stark-split by ~0.1 MHz and there is no Zeeman ladder at any B.

### Propagator model (`--model prop`, handoff gaps 1-2), budget 200

Per-pulse propagator on the addressed states x phonons 0..8 with detunings, ladders and displacement drives; same centers and per-pulse pi calibration (eta*Omega_ref/2pi = 180 Hz, t = 2.78 ms, square pulse), lines included when their peak off-resonant transfer exceeds 1e-3, phonon truncation converged (nph 8 vs 12: 1e-8 bits). No open/closed distinction: ladders are physical.

| library | E=0, B=1 G | E=60, B=1 G |
|---|---|---|
| six | 27.0 +/- 2.6 / 18 / 0 | 30.3 +/- 2.6 / 24 / 0 |
| pi | 38.5 +/- 4.1 / 25 / 1 | 50.4 +/- 5.0 / 34 / 3 |

six, E=0 vs B: 3.6 G 29.6 +/- 2.7 / 24; 30 G 31.4 +/- 2.8 / 22; 100 G 31.3 +/- 2.1 / 29.

Reading: the pi-model open-library 16 cycles at E=0 was an artifact of treating Zeeman ladders and light-shift drives as perfect pi pulses; with the real dynamics the six-pair library needs ~27 cycles at E=0 and ~30 at E=60 (a 3-cycle difference, within s.e.), and B does not help. Pipi's pi library costs ~40 % more at E=0 and ~65 % at E=60. The best first pulse under the propagator is a single rotational line at 43.6 GHz (0.26 bits), not the Zeeman ladder (0.32 bits in the pi model). Histograms stay bimodal: a lump at the 7-cycle floor and a tail beyond 100; the tail sets the mean (median 18 vs mean 27).
