# 232ThF+ two-photon transition matrices vs polarizability anisotropy, 2026-09-17

Two runs of the 2026-09-15 numerical deck at two settings of the rank-2 anisotropy
knob, plus one comparison figure per B panel.

| case | alphas passed | physical setting | directory |
|---|---|---|---|
| xxyy0 | `{"alpha_K2_dOm0": 1.0}` | a_xx = a_yy, so a_xx - a_yy = 0 | `xxyy0/` |
| xxyy_pos | `{"alpha_K2_dOm0": 1.0, "alpha_K2_dOm2": 1.0}` | a_xx - a_yy > 0 | `xxyy_pos/` |

ALL ALPHA VALUES ARE PLACEHOLDERS. No ThF+ polarizability value exists, so the
strengths here are geometry in units of alpha^2, not rates, and the ratio between
the two channels is a choice (both set to unit magnitude, same sign), not a
measurement. Closure caveat, unchanged from 2026-09-15: the effective two-photon
operator is the exact-closure form, valid only for detuning >> the intermediate
rotational structure ~7 GHz, which the 0.16-1.5 GHz JILA detunings do NOT satisfy
(`docs/open-questions.md` OPEN-23). K = 1 is absent (OPEN-21).

Reading: Raman (photon 1, eps1, absorbed; photon 2, eps2, emitted), the default
since 2026-09-18. Delta m_F = p1 - p2. `reading='ladder'` (both photons absorbed,
Delta m_F = p1 + p2) is available but not used here.

## Generating commands

```
conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy0
conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy_pos
conda run -n structure python scripts/compare_thf_twophoton_anisotropy.py
```

Numerical deck, identical in both cases and unchanged from the 2026-09-15 approval:
`padded_thf('232', J_max=7)`, params `thf_v2('232')`, J = 1..5 initial and final,
E_z = 0 V/cm, fixed-B panels at 0.001/1/3/5/10 G, B sweep 0.001..20 G at 201 points
for J=1 -> J=1,2, six polarization pairs: (sigma+, sigma-), (sigma-, sigma+),
(sigma+, sigma+), (sigma+, pi), (sigma-, pi), (pi, pi), reaching Delta m_F =
+2, -2, 0, +1, -1, 0. Each case directory carries its own README, the heatmap
grids, the complex amplitudes (`matrices_B*.npz`), the state labels
(`labels_B*.csv`), the B-sweep curves and arrays, and `channels.npz`.

## What differs physically

Both channels are rank K = 2. They differ in the molecule-frame projection
q = Omega' - Omega:

- **(K, |dOmega|) = (2, 0)**, knob `alpha_K2_dOm0` ~ (2 a_zz - a_xx - a_yy). Present
  in both cases. Connects states of the same signed Omega.
- **(K, |dOmega|) = (2, 2)**, knob `alpha_K2_dOm2` ~ (a_xx - a_yy). Present only in
  xxyy_pos. This is the Omega-doublet-flipping pathway: |q| = 2 requires K >= 2, and
  two E1 legs bound K <= 2, so it exists at K = 2 only. Physically it proceeds
  through an Omega = 0 intermediate state. Both q = +2 and q = -2 belong to this one
  channel and share the one scalar, because parity maps each into the other.

Setting a_xx = a_yy therefore removes the Omega-doublet-flipping pathway entirely
while leaving the dOmega = 0 pathway intact.

The knob is OMITTED in xxyy0 rather than set to zero. `TwoPhotonOperator._keys()`
selects channels by key presence, so omitting it means no (2, 2, P) channel matrix
is ever built and `xxyy0/channels.npz` carries only the five (2, 0, P) matrices.
A present-but-zero knob gives bit-identical amplitudes (probe 2026-09-17, max abs
difference 0.0 at B = 1 G, (sigma+, pi)) but leaves five dead keys in the file.
Neither the channel definitions nor this equivalence depend on the reading, since
the reading only rescales/relabels how eps1, eps2 combine into the dyad weights.

## Sanity checks

**1. The |dOmega| = 2 channel is absent from xxyy0.** `xxyy0/channels.npz` holds
five keys, `K2_dOm0_P{+0,+1,+2,-1,-2}`, and no `dOm2` key. `xxyy_pos/channels.npz`
holds ten. Unaffected by the reading (this is about the channel operator, not the
polarization contraction).

**2. The channel matrices carry exactly the support their labels claim** (from
`xxyy_pos/channels.npz`, against the basis-pair mask |Omega' - Omega| = 2):

| matrix | max abs on |dOmega|=2 pairs | max abs on all other pairs |
|---|---|---|
| K2_dOm0_P+0 | 0 | 0.389706 |
| K2_dOm0_P±1 | 0 | 0.388456 |
| K2_dOm0_P±2 | 0 | 0.554700 |
| K2_dOm2_P+0 | 0.504307 | 0 |
| K2_dOm2_P±1 | 0.424264 | 0 |
| K2_dOm2_P±2 | 0.489898 | 0 |

The |dOmega| = 2 basis pathway lives entirely in the (2, 2, P) matrices, which
xxyy0 never builds. That is the well-defined form of the "dOmega = 2 elements are
exactly zero in xxyy0" check.

**3. "Dominant Omega" is not a usable label at E_z = 0.** Every one of the 140
selected eigenvectors carries population exactly 0.500000 on Omega = +1 at
E_z = 0, B = 1 G, because the eigenstates are exact parity eigenstates and
therefore equal-weight superpositions of Omega = +1 and Omega = -1. So `argmax`
over the eigenvector picks arbitrarily between two equal amplitudes, and a check
phrased as "elements connecting states whose dominant Omega differs" has no
content at zero field. Check 2 is the statement that survives. A consequence:
both pathways connect the same pairs of eigenstates, and the dOmega = 2 channel
mostly rescales existing elements rather than opening new ones.

**4. Both pathways are parity-even.** Max abs amplitude over opposite-parity state
pairs, worst of the six polarizations at B = 1 G: 1.723e-12 for the dOmega = 0
pathway alone and 1.743e-12 for the dOmega = 2 pathway alone. Both are float64
noise, so neither channel drives an e <-> f two-photon line at zero field, as
gate V24 requires.

**5. Nonzero-element counts at B = 1 G.** The threshold matters, and the
1e-12-relative threshold named in the task is below the float64 noise floor of
these matrices, so it is reported first and then discarded. Recomputed for the
Raman polarization pairs (values shifted from the pre-2026-09-18 ladder-reading
numbers because the amplitude at fixed pair labels is a different matrix under
the new reading):

At |amp| > 1e-12 * max|amp| per matrix:

| pair | nonzero xxyy0 | nonzero xxyy_pos | apparently opened | apparently closed |
|---|---|---|---|---|
| sigma+_sigma- | 958 | 968 | 42 | 32 |
| sigma-_sigma+ | 958 | 968 | 42 | 32 |
| sigma+_sigma+ | 1120 | 1108 | 40 | 52 |
| sigma+_pi | 1108 | 1091 | 40 | 57 |
| sigma-_pi | 1108 | 1091 | 40 | 57 |
| pi_pi | 1120 | 1108 | 40 | 52 |
| total | 6372 | 6334 | 244 | 282 |

Every one of those 526 opened/closed cells is numerical dust (same conclusion as
2026-09-15, only the per-pair labels moved). The counts are an artifact of
scaling the threshold by a max that differs between the two cases.

At an absolute floor |amp| > 1e-10, which sits well above the noise cells and
below the smallest genuine element (see check 13 below), so the gap on either
side is wide and nothing lands near the cut:

| pair | nonzero xxyy0 | nonzero xxyy_pos | opened | closed |
|---|---|---|---|---|
| sigma+_sigma- | 852 | 852 | 0 | 0 |
| sigma-_sigma+ | 852 | 852 | 0 | 0 |
| sigma+_sigma+ | 988 | 988 | 0 | 0 |
| sigma+_pi | 952 | 952 | 0 | 0 |
| sigma-_pi | 952 | 952 | 0 | 0 |
| pi_pi | 988 | 988 | 0 | 0 |

At B = 1 G the two cases have identical support. Every matrix has 19600 cells, so
roughly 95% of them are zero in both cases by the shared |dF|, |dJ|, |dF1|,
|dm_F| <= 2 and even-parity selection rules.

**6. The support is not identical at every field.** Comparing each pathway alone
for the (sigma+, sigma+) and (pi, pi) pairs, the dOmega = 2 pathway is nonzero in
a few cells where the dOmega = 0 pathway is not, and never the reverse:

| field | pair | nonzero dOm0 only | nonzero dOm2 only | dOm2-and-not-dOm0 | dOm0-and-not-dOm2 |
|---|---|---|---|---|---|
| 0.001 G | sigma+_sigma+ | 924 | 932 | 8 | 0 |
| 0.001 G | pi_pi | 940 | 956 | 16 | 0 |
| 1 G | sigma+_sigma+ | 988 | 988 | 0 | 0 |
| 1 G | pi_pi | 988 | 988 | 0 | 0 |
| 10 G | sigma+_sigma+ | 1012 | 1028 | 16 | 0 |
| 10 G | pi_pi | 1028 | 1028 | 0 | 0 |

**7. Channel amplitudes add linearly, as the operator requires.** With
`amp_dOm2` the amplitude from an operator carrying only `alpha_K2_dOm2`, the
residual max|amp_xxyy_pos - (amp_xxyy0 + amp_dOm2)| is exactly 0 for all six
polarization pairs at B = 1 G. Max abs dOmega = 2 amplitude per pair:

| pair | max abs amp, dOmega = 2 pathway alone |
|---|---|
| sigma+_sigma- | 0.490389 |
| sigma-_sigma+ | 0.490389 |
| sigma+_sigma+ | 0.192300 |
| sigma+_pi | 0.300051 |
| sigma-_pi | 0.300051 |
| pi_pi | 0.384601 |

**8. The interference is large where the support is shared.** Ratio
|M_xxyy_pos|^2 / |M_xxyy0|^2 over cells above the 1e-10 floor, sampled at
(sigma+, sigma+) and (pi, pi):

| field | pair | shared cells | enhanced >2x | suppressed <0.5x | min ratio | max ratio | median ratio |
|---|---|---|---|---|---|---|---|
| 0.001 G | sigma+_sigma+ | 932 | 378 | 362 | 0.0333 | 11.94 | 1.983 |
| 0.001 G | pi_pi | 960 | 398 | 370 | 0.0333 | 11.94 | 1.983 |
| 1 G | sigma+_sigma+ | 988 | 402 | 394 | 0.0333 | 11.94 | 1.983 |
| 1 G | pi_pi | 988 | 402 | 394 | 0.0333 | 11.94 | 1.983 |
| 10 G | sigma+_sigma+ | 1028 | 430 | 406 | 0.0333 | 24.85 | 1.983 |
| 10 G | pi_pi | 1028 | 430 | 406 | 0.0333 | 24.85 | 1.983 |

About 80% of the populated cells move by more than a factor 2 in |M|^2, split
almost evenly between enhancement and suppression: the real effect of a_xx - a_yy
at these placeholder alphas is coherent interference that reshuffles the
strength distribution, not a new set of lines. The median ratio is 1.983 at
every field and every polarization pair, a structural coincidence of setting
both alphas to the same unit magnitude, not a physical prediction. The
(sigma+, sigma+)/(pi, pi) numbers coincide at B = 1 G and B = 10 G but diverge at
other fields (e.g. B = 3 G: 1008/414/402 vs 1028/430/406) -- that agreement is
itself incidental, not a symmetry (see check 11).

**9. xxyy_pos reproduces the 2026-09-15 run bit-for-bit.** Per B panel, over all
six polarization pairs:

| B | np.allclose, all 6 pairs | max abs amplitude difference | max eigenvalue difference |
|---|---|---|---|
| 0.001 G | True | 0.000e+00 | 0.000e+00 MHz |
| 1 G | True | 0.000e+00 | 0.000e+00 MHz |
| 3 G | True | 0.000e+00 | 0.000e+00 MHz |
| 5 G | True | 0.000e+00 | 0.000e+00 MHz |
| 10 G | True | 0.000e+00 | 0.000e+00 MHz |

The comparison can fail: xxyy0 against the same 2026-09-15 file differs by
0.4903892810042282 at B = 1 G. The eigensystem is identical between the two cases
(max eigenvalue difference 0.0, identical row and column index arrays), as it
must be, since the alphas enter only the transition operator and not the
Hamiltonian. This check is insensitive to the reading, since both files were
regenerated under the same (Raman) default.

**10. Tests.** `conda run -n structure python -m pytest tests/test_transition.py
tests/test_twophoton.py -q` reports `22 passed in 8.31s`.

**11. (sigma+, sigma+) and (sigma-, sigma-) under the Raman reading.** These are
no longer both in the default `PAIRS` list (only (sigma+, sigma+) is), so the
identity was checked directly: computing (sigma-, sigma-) on the fly and
comparing to the saved (sigma+, sigma+) matrix gives max|amp(sigma+,sigma+) -
amp(sigma-,sigma-)| = 0.000e+00 at B = 0.001, 1, 10 G, in both xxyy0 and
xxyy_pos. **They are the exact same matrix, same sign** (not merely equal up to
a sign flip): with sigma+* = sigma-, Raman(eps1, eps2) = -Ladder(eps1, eps2*)
gives Raman(sigma+, sigma+) = -Ladder(sigma+, sigma-) and, using the eps1<->eps2
symmetry from the absent K = 1 part, Raman(sigma-, sigma-) = -Ladder(sigma-,
sigma+) = -Ladder(sigma+, sigma-) as well -- the two collapse onto the same
ladder amplitude. This is why the 2026-09-18 PAIRS list only needs one of them.

The other four Raman pairs pair up under the physical m_F -> -m_F mirror instead:
(sigma+, sigma-) <-> (sigma-, sigma+) and (sigma+, pi) <-> (sigma-, pi). That
relation is exact only at zero field, since the Zeeman term is the only part of
the Hamiltonian odd under m_F -> -m_F. Residual
max | |amp(sigma+,X)|^2 - |amp(sigma-,X')|^2[mF->-mF] | at B = 1 G, xxyy_pos:
1.560e-03 for (sigma+,sigma-) vs (sigma-,sigma+) and 5.295e-04 for (sigma+,pi) vs
(sigma-,pi), against matrix maxima of order 0.5-0.7 -- small but nonzero, growing
with B as in the 2026-09-15 ladder-reading check.

**12. Floor margin (unaffected by the reading).** At B = 1 G: xxyy0 min|amp|
above the 1e-10 floor 8.074e-09, max|amp| 5.222e-01 (1.91 decades of margin);
xxyy_pos min 5.462e-09, max 7.354e-01 (1.74 decades of margin). Both sit
comfortably clear of the float64 noise floor identified in check 5.

## Comparison figures

`compare_B{0.001,1,3,5,10}G.png`, six panels each, one per polarization pair.
Background is |M|^2 for xxyy_pos drawn by `heff.plot_transition.heatmap`. Overlays,
against an absolute floor |amp| > 1e-10: red squares for cells the |dOmega| = 2
pathway opens, magenta plus signs for cells it closes, orange circles for cells it
enhances by more than 2x in |M|^2, blue crosses for cells it suppresses below 0.5x.
Opened and closed counts per panel, recomputed for the Raman pairs:

| B | pair | shared | opened | closed | enhanced >2x | suppressed <0.5x |
|---|---|---|---|---|---|---|
| 0.001 G | sigma+_sigma- | 794 | 16 | 6 | 302 | 314 |
| 0.001 G | sigma-_sigma+ | 794 | 16 | 6 | 302 | 314 |
| 0.001 G | sigma+_sigma+ | 924 | 8 | 0 | 370 | 362 |
| 0.001 G | sigma+_pi | 890 | 18 | 6 | 352 | 348 |
| 0.001 G | sigma-_pi | 890 | 18 | 6 | 352 | 348 |
| 0.001 G | pi_pi | 932 | 20 | 8 | 378 | 362 |
| 1 G | all six | 852-988 | 0 | 0 | 328-402 | 346-394 |
| 3 G | sigma+_sigma- | 864 | 14 | 4 | 334 | 352 |
| 3 G | sigma-_sigma+ | 864 | 14 | 4 | 334 | 352 |
| 3 G | sigma+_sigma+ | 996 | 8 | 4 | 406 | 398 |
| 3 G | sigma+_pi | 964 | 18 | 6 | 386 | 388 |
| 3 G | sigma-_pi | 964 | 18 | 6 | 386 | 388 |
| 3 G | pi_pi | 1004 | 20 | 4 | 410 | 402 |
| 5 G | sigma+_sigma- | 873 | 14 | 3 | 340 | 355 |
| 5 G | sigma-_sigma+ | 873 | 14 | 3 | 340 | 355 |
| 5 G | sigma+_sigma+ | 1004 | 16 | 4 | 410 | 402 |
| 5 G | sigma+_pi | 972 | 14 | 4 | 392 | 390 |
| 5 G | sigma-_pi | 972 | 14 | 4 | 392 | 390 |
| 5 G | pi_pi | 1012 | 16 | 0 | 414 | 406 |
| 10 G | sigma+_sigma- | 884 | 8 | 0 | 349 | 358 |
| 10 G | sigma-_sigma+ | 884 | 8 | 0 | 349 | 358 |
| 10 G | sigma+_sigma+ | 1012 | 16 | 0 | 414 | 406 |
| 10 G | sigma+_pi | 988 | 4 | 0 | 404 | 394 |
| 10 G | sigma-_pi | 988 | 4 | 0 | 404 | 394 |
| 10 G | pi_pi | 1028 | 0 | 0 | 430 | 406 |

The (sigma+, sigma-) and (sigma-, sigma+) count rows coincide at every field, as
do (sigma+, pi) and (sigma-, pi). (sigma+, sigma+) and (pi, pi) have no separate
mirror row: (sigma+, sigma+) mirrors onto (sigma-, sigma-), which check 11 shows
is the exact same matrix (not merely a matching count row), and (pi, pi) mirrors
onto itself trivially, pi being real. The mirrored-row matrices are NOT
elementwise equal (max abs amplitude difference 0.735 between (sigma+,sigma-)
and (sigma-,sigma+) at B = 1 G xxyy_pos): they are related by reversing the sign
of m_F on both the initial and the final state (check 11), a relation that is
exact only at zero field. The residual grows linearly with B, consistent with
the Zeeman term being the only m_F -> -m_F-odd piece of the Hamiltonian.

Gauge, unchanged: the dominant Condon-Shortley component of each eigenvector is
real positive, so relative phases between matrix elements are deterministic.
Raman reading: photon 1 (eps1) absorbed, photon 2 (eps2) emitted, so
Delta m_F = p1 - p2.
