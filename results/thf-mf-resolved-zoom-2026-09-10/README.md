# ThF+ M_F-resolved level plots — zoomed per-manifold grids

Calculations from 2026-09-14, using the body-frame Zeeman tensor
(G_zz, G_xx, G_yy) and the e-below-f Ω-doublet ordering at every J. The
zoomed companion set to
[../thf-mf-resolved-2026-09-09](../thf-mf-resolved-2026-09-09): same static
diagonalizations at the same four field configurations, one panel per
(J, F1) manifold instead of one axis per J. See
[the design note](../../docs/design/2026-09-09-mf-resolved-level-plots.md).

Each panel is referenced to its own (J, F1) zero-field centroid, printed as
its offset; differences of those offsets give the spacing between any two
manifolds. Panels in a row share one energy window, set by the widest
manifold in that row. Colour encodes zero-field parity (blue p = +1, orange
p = -1).

At B_z = 100 G, each F band carries a fitted g_F -- dE/dM_F of that band
divided by -mu_B B_z -- printed in the panel. For 232Th19F+, J = 1, F = 3/2
this reads g_F = -0.0149, the fit across both Ω-doublet parities; the
closed-form values from the tensor are g_F(e) = -0.015036 and
g_F(f) = -0.014762 (see
[the Hamiltonian summary](../../docs/thf-plus-heff-summary.md)), so the
parity-dependent Δg = g^u - g^l = +2.74e-4 is a real but sub-percent
splitting of that fitted slope, not visible at the printed precision.

12 PNG+PDF pairs, `thf-<iso>-zoom-E<E>-B<B>.{png,pdf}`. Parameter sets match
[../thf-mf-resolved-2026-09-09](../thf-mf-resolved-2026-09-09) and
[../thf-fields-2026-09-08](../thf-fields-2026-09-08).

Known cosmetic issue, not fixed: where two F bands are nearly degenerate
(F = 0 and F = 1 at the top of the 229Th J = 1 and J = 3 panels) their
inline labels overlap. Colour is already spent on parity, so separating them
needs a different channel.

## Reproduction

From the repository root, after installing `.[plot]`:

```shell
python scripts/plot_thf_mf_resolved.py
```

This writes both this set and
[../thf-mf-resolved-2026-09-09](../thf-mf-resolved-2026-09-09) in one run.
