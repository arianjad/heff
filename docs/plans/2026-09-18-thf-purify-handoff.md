# Handoff: ThF+ quantum-logic purification study (updated 2026-09-18 evening)

Status: budget-200 rerun, B sweep (gap 5) and propagator model (gaps 1-2) done and committed. Next decision is Arian's (see bottom).

## Goal (Arian's words)
Find the best pulse library for ML-optimized quantum-logic purification of ThF+ X3Delta1 from a 4 K thermal mixed state to a pure state: Raman sideband pulses coupling molecular transitions to an ion-crystal motional mode, binary phonon readout, recool, repeat; minimize cycles. Test E = 0 and E > 0.

## Model (heff.purify, committed a081508)
- Transition graph from `heff.graph.transition_graph` (nodes = eigenstates J<=8, edges = two-photon lines with strength S, frequency f, polarization pair). Thermal populations at 4 K: 320 states, H(b0) = 6.953 bits, 0.9992 of population.
- A `Pulse` = (pair, center f0, addressed subspace u, joint tables Mc/Mn: initial -> final with/without a phonon). Two libraries build those tables:
  - `pulse_library` (`--model pi`): top-hat window bw, every addressed line on resonance as a sideband pi pulse calibrated on the strongest line, P = sin^2(pi/2 sqrt(S/S_ref)); `closed=True` drops open Zeeman ladders and diagonal light-shift drives.
  - `propagator_library` (`--model prop`): exact propagator on (addressed states) x (phonons 0..nph) in the drive's rotating frame: detunings, ladders, diagonal drives as coherent displacement; lines enter if peak off-resonant transfer g^2/(g^2+D^2) > tol. Power set per pulse (pi on the reference line, t = 1/(2 eta Omega_ref) = 2.78 ms); square pulse; carrier and red sideband dropped (MHz detuned); molecular coherences dropped at readout. Gates: tests/test_graph.py (resonant = pi model, sin^2 law, detuned suppression, harmonic ladder climbs / anharmonic ladder closes, displacement click 1 - exp(-(pi/2)^2)).
- Binary click, ideal readout (eta = 1, dark = 0). Diagonal belief vector; greedy expected-information-gain policy; stop at max belief 0.99.
- Floor: >= H(b0) = 7 cycles on average (one <= 1-quantum mode reset carries <= 1 bit). No pulse in any library has expected gain > 1 bit.

## Operating point (agreed 2026-09-18, from Pipi et al. 2026 FNO paper + Arian)
- Lamb-Dicke eta = 0.09, Omega_ref/2pi = 2.000 kHz (Pipi Eq. A9) -> eta*Omega/2pi = 180 Hz, pi pulse 2.78 ms; grid/window bw = 0.001 MHz.
- All polarizations available (Arian): `six` (all six pairs) and `pi` (Pipi's E1 = pi, E2 = sigma+/-).
- E in {0, 60 V/cm}, B = 1 G (60 V/cm is the repo tutorial convention, not from Arian; Ng 2022 measured at 24 V/cm static, HfF+ eEDM ran 20-24 V/cm rotating).
- Placeholder polarizabilities (K = 2 only, dOmega 0 and 2 channels equal).

## Results (results/thf-purify-2026-09-18/README.md indexes the sibling dirs -budget200, -bsweep, -prop)
Cycles to max belief 0.99, 100 trajectories, budget 200, mean +/- s.e.:

| model / library | E=0 | E=60 |
|---|---|---|
| pi, six open | 15.9 +/- 1.2 | 25.7 +/- 2.3 |
| pi, six closed | 26.2 +/- 2.7 | 33.5 +/- 3.1 |
| pi, pi open / closed | 18.4 / 25.7 | 49.7 / 50.1 |
| propagator, six | 27.0 +/- 2.6 | 30.3 +/- 2.6 |
| propagator, pi | 38.5 +/- 4.1 | 50.4 +/- 5.0 |

- The budget-60 E=60 pi-library means were censored (~25 % at budget); uncensored ~50.
- B sweep (pi model, E=0, six): open 15.9 / 15.8 / 22.0 / 27.7 / 29.8 and closed 26.2 / 27.4 / 29.9 / 41.1 / 45.7 at B = 1 / 3.6 / 10 / 30 / 100 G. Rung anharmonicity from heff scales as B^2 (median 0.2 Hz at 1 G, 2 kHz at 100 G); ladders close only near 100 G; raising B costs cycles because it resolves the Zeeman degeneracy. The residual open-closed gap at 100 G is the 2 diagonal light-shift drives, not ladders. Gap-5 hypothesis falsified.
- Propagator: the pi-model open-library 16 cycles at E=0 was an artifact; honest dynamics gives ~27 (E=0) and ~30 (E=60) for six pairs, B-independent (3.6 G 29.6, 30 G 31.4, 100 G 31.3). Comparable to Pipi 2026: 26.5 pulses to purity 0.98 for H3O+ at 20 K. Histograms bimodal: lump at the 7-cycle floor, tail beyond 100; the tail sets the mean (median 18).

## Known model gaps (ranked, after this round)
1. Tail diagnosis: which true states generate the >50-cycle trajectories (weakly connected states? states reachable only through weak lines?). Cheap: log the true state per trajectory and histogram cycles by (J, F) of the true state. This targets library design directly.
2. Greedy is myopic and fully adaptive; Pipi plan an open-loop sequence along the no-click branch. A lookahead-2 or learned policy (torch in env `structure`; qutip, jax, gymnasium absent) could shorten the tail.
3. Power per pulse is unbounded (a weak reference line gets S_max/S_ref more power); a laser-power cap would drop or lengthen weak pulses. Fixed-power variant: normalize Omega to Smax, t = pi/(eta Omega_l) per pulse, cap t.
4. Ideal readout, no off-resonant scattering (Sinhal 2020: ~1000 QND cycles at 10 GHz detuning, ~10 at 100 MHz), no BBR between cycles, molecular coherences dropped at readout, square pulses (shaped pulses would suppress off-resonant excitation and favor the multi-line pulses).
5. Placeholder polarizabilities; no scalar K = 0 shift.

## Git state
- Committed: 9a8bc16, 06d8c3e, 2e9cc3c (earlier); b407e34 (results index), a081508 (propagator model + deck flags), plus this round's results commit. Nothing pushed. Arian confirms before any `git push`.

## Provenance of numbers
- Pipi 2026 text: job tmp `pipi2026.txt` (and `pipi2024.txt`); probes `probe_*.py`/`.json` and `lit_sideband.md` in `/Users/arianjadbabaie/.claude/jobs/260d4de9/tmp/` (job-scoped, not durable, Mac).
- Anharmonicity table: scratch `zeeman_anharm.py` (Windows session temp, not durable); numbers copied into results/thf-purify-2026-09-18/README.md.
- Verified on PDFs: Chou 2017 sideband probe pulses 1 ms (p. 9); Sinhal 2020 eta ~ 0.1, Omega0/2pi ~ 90 kHz refer to the Ca+ readout sideband, not the molecular coupling; scattering budget p. 4.
