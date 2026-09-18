# 2026-09-18: Regenerate two-photon results after Raman-default flip (24ccf28)

Regenerated all three result sets under the new default reading (Raman: photon 1
absorbed, photon 2 emitted, Delta m_F = p1 - p2) and rewrote the anisotropy README.

## Step 1: regeneration

Ran, in order:
```
conda run -n structure python scripts/plot_thf_twophoton.py
conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy0
conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy_pos
conda run -n structure python scripts/compare_thf_twophoton_anisotropy.py
```
Each of the three plot_thf_twophoton.py runs took roughly 10-15 minutes wall
clock (140x140 manifold, 5 B panels + 201-point sweep, 6 polarization pairs);
the compare script ran in under a minute. Regenerated in place: `channels.npz`,
`heatmaps_B*.png`, `labels_B*.csv`, `matrices_B*.npz`, `curves_*.png`,
`sweep_*.npz`, and each directory's `README.md`, in
`results/thf-twophoton-2026-09-15/`, `results/thf-twophoton-2026-09-17-anisotropy/xxyy0/`,
and `.../xxyy_pos/`, plus `compare_B*.png` in
`results/thf-twophoton-2026-09-17-anisotropy/`.

New pair `(sigma-, sigma+)` produced new files `curves_sigma-_sigma+.png` and
`sweep_sigma-_sigma+.npz` in all three case directories. The old pair
`(sigma-, sigma-)` is no longer in `PAIRS`, so its per-pair files were stale
(untouched by the regen run) and were `git rm`'d from all three directories:
- `results/thf-twophoton-2026-09-15/curves_sigma-_sigma-.png`
- `results/thf-twophoton-2026-09-15/sweep_sigma-_sigma-.npz`
- `results/thf-twophoton-2026-09-17-anisotropy/xxyy0/curves_sigma-_sigma-.png`
- `results/thf-twophoton-2026-09-17-anisotropy/xxyy0/sweep_sigma-_sigma-.npz`
- `results/thf-twophoton-2026-09-17-anisotropy/xxyy_pos/curves_sigma-_sigma-.png`
- `results/thf-twophoton-2026-09-17-anisotropy/xxyy_pos/sweep_sigma-_sigma-.npz`

Aggregate per-B files (`matrices_B*.npz`, `heatmaps_B*.png`, `labels_B*.csv`,
`channels.npz`) hold all six current pairs' keys and were simply overwritten in
place, not renamed.

Confirmed via `git status --short`: after the `git rm`, every remaining tracked
file under both result directories was either regenerated (`M`) or newly added
(`??` for the new `sigma-_sigma+` files, then staged), i.e. each directory holds
exactly the new run's outputs.

## Step 2: identity check against pre-regen ladder snapshot

Snapshot files (saved before regeneration, old ladder-reading run):
`ladder_2026-09-15_matrices_B1G.npz`, `ladder_xxyy0_matrices_B1G.npz` in
`/Users/arianjadbabaie/.claude/jobs/260d4de9/tmp/ladder_snapshot/`.

Checked `new["key"]` against `sign * old["oldkey"]` at B = 1 G, default case,
via `np.max(np.abs(...))`:

| new key | relation | old key | max abs diff |
|---|---|---|---|
| sigma+_sigma- | = -1 * | sigma+_sigma+ | 0.000e+00 |
| sigma-_sigma+ | = -1 * | sigma-_sigma- | 0.000e+00 |
| sigma+_sigma+ | = -1 * | sigma+_sigma- | 0.000e+00 |
| sigma+_pi | = +1 * | sigma+_pi | 0.000e+00 |
| sigma-_pi | = +1 * | sigma-_pi | 0.000e+00 |
| pi_pi | = +1 * | pi_pi | 0.000e+00 |

All six signs matched exactly what was stated in the task brief (no sign
surprises). All differences are exactly 0.0 (bit-identical), confirming
Raman(eps1, eps2) = -Ladder(eps1, eps2*) elementwise, not just in magnitude.

`xxyy_pos` vs the new default-run `matrices_B1G.npz`: every key (`sigma+_sigma-`,
`sigma-_sigma+`, `sigma+_sigma+`, `sigma+_pi`, `sigma-_pi`, `pi_pi`, `evals`,
`rows`, `cols`) matches with diff exactly 0 -- same alphas, same run.

Additionally verified directly (not requested by name, but needed for the
README's mirror-pair paragraph): Raman(sigma+, sigma+) equals Raman(sigma-,
sigma-) **exactly, same sign** (max abs diff 0.0 at B = 0.001/1/10 G, both
xxyy0 and xxyy_pos) -- this follows from the same identity plus the
eps1<->eps2 symmetry of the K=1-absent operator, both collapsing onto
-Ladder(sigma+, sigma-).

## Step 3: README rewrite

`results/thf-twophoton-2026-09-17-anisotropy/README.md` rewritten for the Raman
reading and the new pair set: reading statement, the pair list and Delta m_F
values, and sanity-check sections 5 (nonzero counts), 6 (support at 3 fields),
7 (linearity), 8 (interference ratios), and the comparison-figure count table
were all recomputed from the regenerated outputs (script:
`/Users/arianjadbabaie/.claude/jobs/260d4de9/tmp/readme_recompute.py`, plus a
follow-up script for the sigma+_sigma+/sigma-_sigma- identity). The
"(sigma+, sigma+) and (sigma-, sigma-) count rows coincide" paragraph was
replaced: under Raman those two are the exact same matrix (not merely matching
counts), so there is no separate (sigma-, sigma-) row; the m_F -> -m_F mirror
pairs are now (sigma+, sigma-) vs (sigma-, sigma+) and (sigma+, pi) vs
(sigma-, pi), which coincide in the comparison-figure counts but are not
elementwise equal (residual grows linearly with B, from the Zeeman term).
Placeholder-alpha and closure caveats kept verbatim.

## Step 4: PNG spot check

Read `results/thf-twophoton-2026-09-15/heatmaps_B1G.png` and
`results/thf-twophoton-2026-09-17-anisotropy/compare_B1G.png` directly. Panel
titles read `(sigma+, sigma-) Delta m_F = +2`, `(sigma+, sigma+) Delta m_F = +0`,
`(sigma-, pi) Delta m_F = -1`, and all six panels in both figures are populated.

## Step 5: tests

`conda run -n structure python -m pytest tests/test_transition.py
tests/test_twophoton.py -q` -> `22 passed in 8.31s`.

## Step 6: commit

Committed by explicit pathspec: `results/thf-twophoton-2026-09-15`,
`results/thf-twophoton-2026-09-17-anisotropy`, and this checkpoint file. No
push.
