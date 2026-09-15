# Transition-matrix pipeline — design

Agreed with Arian 2026-09-15 (this session). Status: approved for planning.

## 1. Goal

A generator that, given an operator, an initial manifold, a final manifold and a
field point (E_z, B_z), returns the complex transition matrix

    M_ij(E, B) = <i(E,B)| T |j(E,B)>

with labelled rows and columns, plus plots of |M_ij|^2 and of its change with B.
First test: two-photon (2 x E1) transitions between eigenstates of X 3Delta1,
v = 0, 232Th19F+, initial and final manifolds J = 1..5, for every polarization
combination of two photons, at a few B near a few gauss and as a B sweep.

The operator, the manifolds, the isotopologue and the fields are all user input.
Nothing in the core is ThF+-specific.

## 2. What already exists (verified this session)

- `heff/twophoton.py`: the rank-K effective polarizability geometry per channel
  (K, |dOmega|, P) between any two ket arrays, and the polarization dyad.
  Audited 2026-09-15 against B&C 5.141/5.142/5.172/5.174/5.186 and against the
  C2v AC-Stark tensor (Caldwell 2020 App. A): every claim confirmed, five
  cosmetic defects, see `docs/lit/2026-09-15-twophoton-operator-audit.md`.
- `heff/spectra._strengths_from_matrices`: the eigenvector sandwich
  `sum_k w_k V_a^H M_k V_b`, amplitudes summed over channels before squaring.
- `heff/track.py`: deterministic sign gauge (`apply_gauge('pinned')`) and
  adiabatic ordering along a sweep (`order_states`).
- `heff/engine.sweep`, `heff/model.load_model`, `heff/spec.block_by_mF`.
- Probe results (job tmp `probe2g.py`, `probe_pad.py`, 232ThF+ J <= 3):
  K = 0 channel is exactly the identity; every K = 2 channel obeys
  <a|T_P|b> = (-1)^P <b|T_-P|a> to 0.0; the 232 basis padded to the two-spin
  dtype at I_Th = 0 equals the native one-spin basis ket-for-ket, its rank-2
  geometry equals an independent one-spin B&C formula elementwise and signed to
  6e-17 on all 11 channel matrices, and its Hamiltonian equals the native one
  elementwise to 1e-14.

## 3. Decisions

### 3.1 Operator = channel matrices + weights

An operator is any object exposing

    channels(kets_a, kets_b, ctx) -> {key: (n_a, n_b) real array}   # field-free, built once
    weights(**pol)                -> {key: complex}                 # per polarization choice

`transition.py` owns the sandwich; `spectra` (E1) and `twophoton` (rank-2
polarizability) supply operators. Nothing else in the pipeline knows what T is.
This is the Molecule-Structure `Calculate_TDM_evecs` pattern, generalized from
"one key per p" to arbitrary channel keys, with the coherent sum kept.

### 3.2 Two-photon operator scope

- Transitions use K = 2 only: `alpha_K2_dOm0` (T^2_0, ~ 2 a_zz - a_xx - a_yy,
  dOmega = 0) and `alpha_K2_dOm2` (T^2_{+-2}, ~ a_xx - a_yy, dOmega = +-2).
  Two general scalars; the Cartesian mapping is documentation, not code, since
  for a case (c) Omega = +-1 state the operator statement is the reduced element
  per q. Default alphas are 1.0 placeholders and every plot says so.
- K = 0 is the identity (probe) and is excluded from transition matrices; it
  stays available as a light-shift knob only.
- K = 1 is out of scope (vanishes in common-denominator closure; survives for a
  single resolved intermediate, OPEN-21). Not built.
- Photon reading: both photons absorbed, Delta m_F = p1 + p2. `dyad_weights`
  gets `reading="ladder"` as the default; `"raman"` (conjugate second vector)
  stays as an explicit option. Tests, notebook and [HAM] S9.5 follow the new
  default; stale sentences are deleted, not annotated.
- Polarizations: named (`sigma+`, `sigma-`, `pi`, `x`, `y`) or any Cartesian
  Jones vector, all relative to lab z, where E and B point. Panels are
  labelled by the pair and its Delta m_F. Because the K = 1 part is absent the
  operator is symmetric in the two photons, so the default grid is the six
  unordered pairs ++ , -- , 00 , +- , +0 , -0.

### 3.3 232 padding lives in the package

`twophoton` requires the two-spin ket dtype (it reads `F1` and `ctx.spins[0]`).
A package-level helper pads a one-spin problem to the two-spin dtype at
I_Th = 0. Gate: padded geometry == independent one-spin B&C formula, elementwise
and signed, for every (K, dOmega, P); padded H == native H elementwise. The
analytic reason (the I_Th = 0 recoupling line is exactly +1) goes in the
docstring.

### 3.4 Diagonalization per signed-m_F block, deterministic gauge

Fields along z conserve m_F, so each block is diagonalized separately (no
cross-m_F degeneracy at B = 0; m_F exact). Eigenvector phases from `eigh` are
arbitrary whether or not blocks are merged. For a single element M_ij every
channel shares the same two eigenvectors, so the coherent channel sum carries
one overall phase and |M_ij|^2 is gauge invariant. Relative phases between
different elements are fixed by a deterministic gauge: `track.apply_gauge
('pinned')` (dominant Condon-Shortley basis component real positive) applied
identically in every block, plus `track.order_states` along B. The complex
gauge-fixed M is stored; |M|^2 is what is plotted by default.

Gates: (i) random per-eigenvector phases leave |M|^2 unchanged; (ii) blockwise
and full-basis diagonalization give the same |M|^2.

### 3.5 Manifolds

A manifold is a filter on zero-field dominant labels: J, F1, F, e/f, m_F
(any subset). Initial and final are independent filters over the same
problem. Labels at field come from adiabatic tracking from the zero-field
assignment. Row/column order: by (J, F1, F, e/f, m_F).

### 3.6 Field dependence

`M` is recomputed per field point from cached channel matrices and the
eigenvectors at that point (the generator/eigensystem split Arian asked for).
Deliverables: (a) the full labelled matrix at a few B values; (b) |M_ij|^2 vs
B curves for selected elements, tracked; (c) B grid starts at a small nonzero
value (1e-3 G default) as a convention, not a necessity, since blockwise
diagonalization never sees the m_F degeneracy.

### 3.7 Plots

One heatmap per polarization pair, rows = initial manifold, columns = final
manifold, hierarchical ticks: J blocks as major ticks always; F, e/f, m_F minor
labels when the block is small enough to read. Zero entries drawn as blank.
Output to `results/thf-twophoton-<date>/` with a README naming the model,
J_max, field points, alphas (placeholder) and the closure caveat.

## 4. Gates (each names its unique failure mode; PASS and FAIL both reachable)

1. Delta m_F of every nonzero element equals p1 + p2 for the pair
   (catches a reading/conjugation slip).
2. At E = 0 the matrix connects e->e and f->f only (catches a parity phase).
3. sum_j |M_ij|^2 over a complete final manifold is independent of the initial
   m_F for pure sigma pairs, in the eigenbasis at field (catches a lost
   channel or a wrong Wigner-Eckart projection).
4. Random eigenvector phases leave |M|^2 unchanged (catches a gauge leak).
5. Blockwise == full-basis |M|^2 (catches a cross-block bookkeeping error).
6. x-polarized M from one Jones vector == phased sigma+ / sigma- amplitude sum,
   and != the incoherent sum somewhere (catches squaring before summing).
7. Zero-field |M|^2 reproduces `two_photon_line_strengths` (regression).
8. Padded 232 geometry == one-spin B&C formula, elementwise and signed; padded
   H == native H (catches a padding phase).
9. Along B, tracked |M_ij|^2 is continuous and label permutations are absent
   (catches a tracking failure at an avoided crossing).

## 5. Files

- `heff/twophoton.py`: `reading` argument; named polarizations; docstring sweep.
- `heff/transition.py` (new): operator protocol, manifold filter, padding
  helper, `transition_matrix`, `sweep_transition_matrix`.
- `heff/plot_transition.py` (new): heatmap grid, vs-B curves.
- `scripts/plot_thf_twophoton.py` (new): the 232 J = 1..5 run.
- `tests/test_transition.py` (new), edits to `tests/test_twophoton.py`.
- Docs: `docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5, the
  isotopologue notebook builder S11-13, `docs/open-questions.md`,
  `docs/models.md` if knobs change. Audit defects 1-4 fixed in passing.

## 6. Left out, with the trigger that adds it

- Resolved intermediate-state sum and K = 1: when the ThF+ intermediate ladder
  (energies, 0+/0- labels, dipoles) is settled. Every plot carries the caveat.
- Non-collinear E and B, or light not referenced to z: when an experiment
  geometry needs it (the operator already takes arbitrary Jones vectors).
- Interactive B slider: small multiples cover (b); add if the grid is not enough.
- Physical alpha values: none exist for ThF+ (digest gap 1).
