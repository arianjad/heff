# M_F-resolved level plots per J manifold — ThF+ X3Delta1 isotopologues

Status: in progress, 2026-09-09. Requested by Arian.

## Goal

For each ThF+ isotopologue (232, 229, 227) and each rotational manifold J = 1, 2, 3,
plot every M_F sublevel as a static level diagram, at four field configurations:

| config | E_z (V/cm) | B_z (G) |
|---|---:|---:|
| field free | 0 | 0 |
| magnetic only | 0 | 100 |
| electric only | 100 | 0 |
| both | 100 | 100 |

These are **static diagonalizations at fixed field**, not sweeps. The existing
`scripts/plot_thf_isotopes.py` produces field *scans* (Stark at B=0, Zeeman at E=0)
and never computes a combined E and B point; the combined config here is new.

## Physics facts established before implementation

Verified 2026-09-09 by running the package in the `structure` env.

- The basis is already blocked by signed M_F: `heff.spec.block_by_mF`, exact for
  collinear E || B || z. `plot_thf_isotopes.py:60` already slices one signed-M_F
  block at a time. No change to the Hamiltonian is required for this task.
- **Every (J, F, M_F) carries two states**, not one: X3Delta1 has Omega = +/-1, so
  each M_F has an e/f (Omega-doublet) parity pair. At E = 100 V/cm the doublet is
  Stark-split by roughly 170 MHz at J = 1, which is the largest structure in the
  figure. Parity must be encoded or the levels read as duplicates.
- Basis content at J_max = 3 (`problem.kets`):

  | iso | I(Th) | dim | F at J=1 | F at J=2 | F at J=3 | M_F range |
  |---|---|---:|---|---|---|---|
  | 232 | 0 | 60 | 1/2, 3/2 | 3/2, 5/2 | 5/2, 7/2 | -3.5 .. 3.5 |
  | 229 | 5/2 | 360 | 1,2,3,4 | 0..5 | 0..6 | -6 .. 6 |
  | 227 | 1/2 | 120 | 0,1,2 | 1,2,3 | 2,3,4 | -4 .. 4 |

- Per-J state counts for 232: J=1 -> 12, J=2 -> 20, J=3 -> 28 (sum 60). Confirmed
  that selecting parent-J <= 3 states out of a J_max = 8 basis reproduces exactly
  these counts, so the parent-J label is a usable selector.
- Above the Omega-doublet the J label is a correlation label, not a quantum number.
  At E = 100 V/cm the molecule is essentially polarized. Figures must say so.

## Decisions (per Arian, 2026-09-09)

- **Layout 3**: one axis per J manifold, three axes per figure (J = 1, 2, 3).
- **All M_F shown**, signed, -F .. +F.
- **Energy reference: per-J zero-field centroid** — subtract the mean zero-field
  energy of that J manifold. The same constant is used for all four configs so
  panels are comparable across configs.
- Bars are drawn at **true M_F**, with no sideways F-offset. Different F levels at
  the same M_F separate vertically on their own; the offset would make the x-axis
  misreport M_F. Colour encodes parity; F is annotated.
- Recomputation is authorized. Parameter sets are the exploratory ones already in
  `plot_thf_isotopes.py::parameters` (229 quadrupole omitted, 227 deformed-nucleus
  estimate), not the package defaults.

## Deliverables and write scopes

Two disjoint scopes, worked in parallel.

### Scope A — M_F tick formatting

Files: `scripts/_mf_ticks.py`, `tests/test_mf_ticks.py`.

Half-integer and integer M_F must render as physicists write them:

| value | rendered |
|---|---|
| 1.0 | `1` |
| -1.0 | `-1` |
| 0.0 | `0` |
| 0.5 | `1/2` |
| -0.5 | `-1/2` |
| 1.5 | `3/2` |
| -3.5 | `-7/2` |

Never `1.0`, never `0.5`. Public surface: `format_mF(value) -> str` and a
matplotlib tick formatter built on it.

### Scope B — data layer and figures

Files: `scripts/_thf_params.py`, `scripts/plot_thf_mf_resolved.py`,
`tests/test_mf_resolved.py`, output under `results/thf-mf-resolved-2026-09-09/`.

`parameters()` is extracted from `plot_thf_isotopes.py` into `_thf_params.py` and
imported by both scripts, so the exploratory parameter sets have one home.

## Required checks (TDD — each written and seen to fail first)

Physical invariants, not smoke tests. Each must be able to both pass and fail.

1. **Completeness.** Parent-J <= 3 selection from a J_max = 8 basis reproduces the
   full J_max = 3 state counts per (J, F) in the table above, for all three
   isotopologues. Fails if any M_F block is dropped or double-counted.
2. **Zero-field M_F degeneracy.** At E = B = 0, states sharing (J, F, parity) are
   degenerate across M_F. Tolerance justified from the numerics, not tuned to pass.
3. **Stark +/-M_F degeneracy.** At B = 0 and E = 100 V/cm, E(+M_F) = E(-M_F) exactly.
4. **Time reversal.** E(M_F, E_z, B_z) = E(-M_F, E_z, -B_z) for the combined config.
5. **Zeeman slope cross-check.** At E = 0, B = 100 G, the splitting within a (J, F)
   is linear in M_F with slope matching `heff.observe.g_factors` computed
   independently. This is the strongest check: it compares the plotted numbers
   against a separate package API rather than against themselves.
6. **Parity is exactly +/-1 at zero field**, and `parity_operator` commutes with the
   B term and anticommutes with the E term (as asserted in the existing script).
7. **Energy reference.** The per-J subtracted constant is identical across all four
   field configs, and equals the zero-field mean of that manifold.
8. **True M_F positions.** Drawn bar centres equal the record M_F exactly — the
   regression guard for the rejected offset layout.
9. **Basis convergence.** J_max = 6 versus J_max = 8 agree within a stated tolerance
   at the strongest config.
10. **Parameter provenance.** The exploratory values are actually applied — e.g.
    229 `g_N_Th` is 0.1460, not the package default. Guards a silent fallback.

## Boundaries

- Do not modify anything under `heff/` — the package excludes plotting by design
  (`heff/__init__.py` docstring). This work lives in `scripts/` and `tests/`.
- Do not modify `scripts/plot_thf_isotopes.py` beyond removing `parameters()` into
  `_thf_params.py` and importing it back. Its outputs must not change.
- Do not touch `results/thf-fields-2026-09-08/`. New output goes in its own dated
  directory.
- The existing suite (351 passed, 2 skipped) must stay green.
