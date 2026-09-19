# Quantum-logic purification of thermal 232ThF+ X3Delta1: pulse-library study

Written 2026-09-19 from the 2026-09-18 sessions. Code: `heff.graph`, `heff.purify`, `scripts/purify_thf.py`,
`scripts/purify_summary_figs.py`. Data and figures: `results/thf-purify-2026-09-18/` (index README) with the
sibling directories `-budget200`, `-bsweep/B<G>`, `-prop`. Working notes and ranked gaps:
`docs/plans/2026-09-18-thf-purify-handoff.md`.

## Question

Which library of Raman sideband pulses purifies a 4 K thermal ThF+ X3Delta1 ion into one known quantum state in the
fewest cycles, where one cycle is: pulse coupling a molecular transition to an ion-crystal motional mode, binary
phonon readout, recool, repeat (Pipi et al. 2026 style)? Compare E_z = 0 and E_z > 0, and the full set of six
polarization pairs against Pipi's pi / sigma+- library.

## Model

- **Manifold.** Eigenstates of the heff 232ThF+ X3Delta1 Hamiltonian (J_max = 10 basis), J <= 8 kept: 320 states
  holding 0.9992 of the 4 K Boltzmann population, initial entropy H(b0) = 6.953 bits. Placeholder polarizabilities
  (K = 2 only, dOmega = 0 and 2 channels equal, no scalar K = 0 shift).
- **Lines.** Two-photon (rank <= 2) transition strengths S_ij = |M_ij|^2 and frequencies f_ij = E_j - E_i for each
  polarization pair, as a networkx graph (`heff.graph.transition_graph`). Six pairs: (s+,s-), (s-,s+), (s+,s+),
  (s+,pi), (s-,pi), (pi,pi). Lines below 1e-5 Smax dropped.
- **Pulse.** One polarization pair and a center frequency f0 on a 500 Hz grid, calibrated as a sideband pi pulse on
  the strongest line within +-500 Hz of f0 (S_ref). Power is set per pulse. Operating point (agreed 2026-09-18, from
  Pipi et al. 2026 Eq. A9 and Arian): Lamb-Dicke eta = 0.09, Omega_ref/2pi = 2 kHz, so eta*Omega/2pi = 180 Hz and
  t_pi = 2.78 ms.
- **Two transfer models**, both stored as joint tables Mc/Mn (initial -> final, with/without a phonon) over the
  addressed subspace:
  - *pi model* (`pulse_library`): top-hat 1 kHz window, every addressed line on resonance, one destination per
    initial state, P = sin^2(pi/2 sqrt(S/S_ref)). `closed=True` drops open Zeeman ladders and diagonal
    light-shift drives, which this model cannot describe.
  - *propagator model* (`propagator_library`): exact propagator of each pulse on (addressed states) x (phonons
    0..8) in the drive's rotating frame, with detunings, ladder climbing, and diagonal drives as coherent
    displacement. A line enters if its peak off-resonant transfer g^2/(g^2 + D^2) > 1e-3. Square pulse; carrier
    and red sideband (MHz detuned) dropped; molecular coherences dropped at readout. Phonon truncation converged
    (nph 8 vs 12: 1e-8 bits). Gates in `tests/test_graph.py`.
- **Readout and policy.** Ideal binary click (efficiency 1, dark 0). Diagonal belief vector over the 320 states;
  greedy expected-information-gain pulse choice; stop at max belief >= 0.99; budget 200 cycles.
- **Floor.** One reset of a mode holding <= 1 quantum carries <= 1 bit, so the mean cycle count is >= H(b0) = 7.
  No pulse in any library has expected gain > 1 bit (checked).

## Results

100 trajectories, seed 0, budget 200, B_z = 1 G; mean +- s.e. / median cycles to max belief 0.99.

| model / library | E_z = 0 | E_z = 60 V/cm |
|---|---|---|
| pi model, six pairs, open | 15.9 +- 1.2 / 12 | 25.7 +- 2.3 / 18 |
| pi model, six pairs, closed | 26.2 +- 2.7 / 16 | 33.5 +- 3.1 / 22 |
| pi model, pi library, open / closed | 18.4 / 25.7 | 49.7 / 50.1 |
| **propagator, six pairs** | **27.0 +- 2.6 / 18** | **30.3 +- 2.6 / 24** |
| propagator, pi library | 38.5 +- 4.1 / 25 | 50.4 +- 5.0 / 34 |

B sweep at E_z = 0, six pairs (B = 1 / 3.6 / 10 / 30 / 100 G): pi open 15.9 / 15.8 / 22.0 / 27.7 / 29.8;
pi closed 26.2 / 27.4 / 29.9 / 41.1 / 45.7; propagator 27.0 / 29.6 / - / 31.4 / 31.3.

An earlier budget-60 run censored ~25 % of the E = 60 pi-library trajectories and reported ~33; uncensored is ~50.
`run()` now samples the joint table, so the budget-200 trajectories are fresh realizations, not extensions.

## What we learned

1. **~27 cycles at E = 0 and ~30 at E = 60 with all six pairs** is the honest number. The pi model's 16 cycles at
   E = 0 came from treating Zeeman ladders and light-shift drives as perfect pi pulses; with real detuned dynamics
   the open and closed libraries collapse to the same answer. Comparable to Pipi 2026 (26.5 pulses to purity 0.98,
   H3O+ at 20 K).
2. **E_z barely matters once the dynamics are honest**: +3 cycles at 60 V/cm, within error. The pi model had said +10.
3. **Raising B does not help and mostly hurts.** Rung anharmonicity scales as B^2 (median 0.2 Hz at 1 G, 2 Hz at
   Pipi's 3.6 G, 2 kHz at 100 G), so Zeeman ladders close only near 100 G; but resolving the Zeeman lines removes
   the degeneracy that lets one pulse address many states. The residual open-vs-closed gap at 100 G is the two
   diagonal light-shift drives, not ladders. Hypothesis that a modest B closes the ladders: falsified.
4. **Pipi's pi / sigma+- library costs 40-65 % more** than all six pairs (38 vs 27 at E = 0, 50 vs 30 at E = 60).
5. **Every cycle distribution is bimodal**: a lump on the 7-cycle entropy floor and a tail past 100 cycles. The tail
   sets the mean (median 18 vs mean 27) and has not yet been attributed to specific initial states.

## Structure of the line graph (what the pulses have to work with)

Figures: `summary_strength_matrix.png` (S_ij heatmap, J <= 2), `summary_lines_3d.png` (all lines as (i, j, f_ij)),
`summary_lines_dJ0.png` and `summary_lines_dJ12.png` (the Delta J planes unfolded). Counts: 11796 lines at E = 0
(2616 with Delta J = 0, 2 x 2822 with |Delta J| = 1, 2 x 1768 with |Delta J| = 2), 22058 at E = 60; no |Delta J| >= 3.

- **Shelves.** Two E1 photons give Delta J = 0, +-1, +-2, at 2B(J+1) and 2B(2J+3) with B_rot = 7274.3 MHz.
- **Rows within a shelf, E = 0.** Parity: e/f alternates with J and two photons conserve parity, so Delta J = +1
  connects e(J) <-> f(J+1) only, at residual +-[Delta_ef(J) + Delta_ef(J+1)]/2, and Delta J = +2 connects e <-> e,
  f <-> f at +-[Delta_ef(J+2) - Delta_ef(J)]/2. Delta_ef = 5.29 J(J+1)/2 MHz (J(J+1) law verified in the
  eigenvalues). Each row splits by the four F -> F' combinations; hyperfine F = J-1/2 <-> J+1/2 falls from 15.0 MHz
  (J = 1) to 2.2 MHz (J = 8). Inside the Delta J = 0 plane only Zeeman rungs (21 kHz per rung at J = 1 falling to
  1.7 kHz at J = 8) and hyperfine lines exist; e <-> f is forbidden.
- **Streaks.** One initial state reaches a whole (J', F', e/f) multiplet whose members differ by kHz Zeeman shifts.
  Since g_F differs between J and J+1, the rotational rows resolve into single m -> m' lines at 1 kHz; the Delta J = 0
  Zeeman ladder is the only degenerate multi-line resource at E = 0.
- **E = 60 V/cm.** Parity mixing turns on the forbidden rows (line count doubles). J <= 3 is linear Stark
  (dE Omega m/J(J+1) = 51, 17, 8.5 MHz per unit m at J = 1, 2, 3): m-states fan out, rows dissolve into dense
  1-100 MHz clusters, the J = 1 block spans 104 MHz = dE. J >= 4 is quadratic Stark (doublet dominates): rung
  spacing linear in m, rows bend into arcs of ~0.7 MHz (Delta J = 0, J = 8) to ~10 MHz (Delta J = 1, J = 3-5).
  Formula-level row positions were checked against the figures by eye; the splittings are from the eigenvalues.

## Known gaps, ranked

1. Attribute the >50-cycle tail to initial states (log the true state per trajectory; histogram by J, F). Cheapest
   and most directly useful for library design.
2. Greedy is myopic and adaptive; Pipi plan an open-loop no-click sequence. Lookahead or a learned policy.
3. Per-pulse power is unbounded; a fixed-power variant would lengthen or drop weak-line pulses.
4. Ideal readout, no off-resonant scattering or BBR, coherences dropped at readout, square pulses.
5. Placeholder polarizabilities; no scalar K = 0 shift.

## Provenance

Pipi et al. 2026 numbers from the paper text (job-scoped notes, Mac). Anharmonicity table and structure splittings
from `heff` diagonalizations in-session (scratch scripts, numbers copied into `results/thf-purify-2026-09-18/README.md`
and this note). E_z = 60 V/cm is the repo tutorial convention, not a chosen operating point (Ng 2022 measured at
24 V/cm; HfF+ eEDM ran 20-24 V/cm rotating). Seed and trajectory count are per-run headers in each `README_<variant>.md`.
