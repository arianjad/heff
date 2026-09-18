# Handoff: ThF+ quantum-logic purification study (paused 2026-09-18)

Session paused mid-run at Arian's request. Everything below was verified in-session unless marked otherwise.

## Goal (Arian's words)
Find the best pulse library for ML-optimized quantum-logic purification of ThF+ X3Delta1 from a 4 K thermal mixed state to a pure state: Raman sideband pulses coupling molecular transitions to an ion-crystal motional mode, binary phonon readout, recool, repeat; minimize cycles. Test E = 0 and E > 0.

## Model (heff.purify, committed)
- Transition graph from `heff.graph.transition_graph` (nodes = eigenstates J<=8, edges = two-photon lines with strength S, frequency f, polarization pair). Thermal populations at 4 K: 320 states, H(b0) = 6.953 bits, 0.9992 of population.
- Pulse = (pair, center f0, top-hat window bw). Sideband pi pulse on the strongest addressed line; weaker lines transfer P = sin^2(pi/2 sqrt(S/S_ref)). Binary click, ideal readout (eta = 1, dark = 0). Diagonal belief vector; greedy expected-information-gain policy; stop at max belief 0.99, budget 60 cycles.
- `pulse_library(..., closed=True)` keeps only closed two-level pulses (no state both source and destination, no diagonal light-shift drive). Open ladders (Zeeman band |m,0> -> |m+1,1> -> |m+2,2>, and diagonal |u,0> -> |u,1> -> |u,2>) are outside the pi-pulse model.
- Floor: the motional mode is the entropy sink and one reset of a mode holding <= 1 quantum carries <= 1 bit, so >= H(b0) = 7 cycles on average. Verified in the classical model: no pulse in either library has expected gain > 1 bit; no many-to-one (non-unitary) pulses exist.

## Operating point (agreed 2026-09-18, from Pipi et al. 2026 FNO paper + Arian)
- Lamb-Dicke eta = 0.09 (Arian: "eta like Pipi, 0.09"), Omega_ref/2pi = 2.000 kHz (Pipi Eq. A9) -> eta*Omega/2pi = 180 Hz, pi pulse ~2.8 ms -> window bw = 0.001 MHz.
- All polarizations available (Arian). Two libraries run: `six` (all six pairs) and `pi` (Pipi's E1 = pi, E2 = sigma+/-).
- E in {0, 60 V/cm}, B = 1 G. 60 V/cm is the repo tutorial convention, not from Arian; Ng 2022 measured at 24 V/cm static, HfF+ eEDM ran 20-24 V/cm rotating.
- Placeholder polarizabilities (K = 2 only, dOmega 0 and 2 channels equal). Varying their ratio moved cycle means by amounts comparable to 20-trajectory noise, not the order of magnitude.

## Results so far (results/thf-purify-2026-09-18/, 100 trajectories, seed 0, NOT yet committed)
Means from run logs; s.e. and budget hits are in each README_<variant>.md.

| library | E=0 mean / median | E=60 mean / median |
|---|---|---|
| six, open | 16.7 / 14 | 25.3 / 18 |
| pi, open | 19.6 / 16 | 32.9 / 31 |
| pi, closed | 25.5 / 20 | 32.9 / 32 |
| six, closed | 26.7 / 22 | RUNNING when paused |

Files present: cycles_E{0,60}_B1_{variant}.png (7 of 8), graph_J12_E{0,60}_B1.png, README_{six,pi,pi_closed}.md, run_*.log. Missing until the last process finishes: README_six_closed.md, cycles_E60_B1_six_closed.png. Check with `pgrep -f purify_thf.py` and `tail results/thf-purify-2026-09-18/run_sixclosed.log` (look for a line starting `wrote`).

Reading: E=0 open six-pair library is the best case at 17 cycles; realistic resolution (1 kHz) roughly doubled the earlier 50 kHz estimate; closing the ladders costs ~10 cycles at E=0 because the degenerate Zeeman ladder is the main multi-state resource; E=60 spreads lines by the Stark shift (J=1 doublet split d*E = 3.37 D x 60 V/cm = 102 MHz, consistent with heff) so pulses address fewer states. Comparable scale to Pipi 2026: 26.5 pulses to purity 0.98 for H3O+ at 20 K.

## Known model gaps (ranked)
1. Open ladders are treated as independent pi pulses; a real drive climbs the phonon ladder (coherent displacement for diagonal drives). Pipi truncate to phonons 0,1 too but integrate the real time-dependent Hamiltonian with detunings. Fix: per-pulse propagator on the addressed states x phonons {0,1,2}, precomputed once per pulse into a transfer/click table so greedy stays cheap.
2. Window and Rabi rate decoupled; detuning inside the window ignored. Same fix.
3. Ideal readout, no off-resonant scattering (Sinhal 2020: ~1000 QND cycles at 10 GHz detuning, ~10 at 100 MHz), no BBR between cycles.
4. Placeholder polarizabilities; no scalar K = 0 shift.
5. Hypothesis to test: raising B until the Zeeman ladder is anharmonic by more than eta*Omega closes the ladders (Pipi use 3.6 G). B sweep in heff.
6. Greedy is myopic and fully adaptive; Pipi plan an open-loop sequence along the no-click branch. Not the same protocol.

## Git state
- Committed: 9a8bc16 (graph/purify/deck), 06d8c3e (architecture.md), 2e9cc3c (closed option, deck at Pipi operating point). Working tree: only `results/thf-purify-2026-09-18/` untracked.
- Nothing pushed. Arian confirms before any `git push`.

## Next steps when resumed
1. Confirm the six_closed run finished; view the 8 histograms; write `results/thf-purify-2026-09-18/README.md` as an index over the four README_<variant>.md files (table above plus s.e.); delete run_*.log or keep as provenance (ask); commit results by pathspec.
2. Add the docs/architecture.md sentence for the `closed` option (one line).
3. Decide with Arian: build the per-pulse propagator model (gap 1-2), or the B sweep (gap 5), or a learned policy (torch available in env `structure`; qutip, jax, gymnasium absent).

## Provenance of numbers
- Pipi 2026 text: job tmp `pipi2026.txt` (and `pipi2024.txt`); probes `probe_*.py`/`.json` and `lit_sideband.md` in `/Users/arianjadbabaie/.claude/jobs/260d4de9/tmp/` (job-scoped, not durable).
- Verified on PDFs: Chou 2017 sideband probe pulses 1 ms (p. 9); Sinhal 2020 eta ~ 0.1, Omega0/2pi ~ 90 kHz refer to the Ca+ readout sideband, not the molecular coupling; scattering budget p. 4.
