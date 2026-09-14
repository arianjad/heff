# ThF+ M_F-resolved level plots

Calculations from 2026-09-14, using the body-frame Zeeman tensor
(G_zz, G_xx, G_yy) and the e-below-f Ω-doublet ordering at every J.
Static diagonalizations at fixed field, not the sweeps in
[../thf-fields-2026-09-08](../thf-fields-2026-09-08). See
[the design note](../../docs/design/2026-09-09-mf-resolved-level-plots.md) and
[the one-page Hamiltonian summary](../../docs/thf-plus-heff-summary.md).

For each isotopologue (232, 229, 227Th19F+) and each of four field
configurations, one figure with three axes (J = 1, 2, 3) shows every M_F
sublevel, drawn at its true M_F with no sideways offset. Colour encodes
zero-field parity (blue p = +1, orange p = -1); annotations give F.

| config | E_z (V/cm) | B_z (G) |
|---|---:|---:|
| field free | 0 | 0 |
| magnetic only | 0 | 100 |
| electric only | 100 | 0 |
| both | 100 | 100 |

12 PNG+PDF pairs, `thf-<iso>-E<E>-B<B>.{png,pdf}`.

Energies are referenced to the per-J zero-field centroid, the same constant
in all four configurations, so panels stay comparable across configs. At
E = 100 V/cm the molecule is essentially polarized: above the Ω-doublet, J is
a correlation label rather than a quantum number, and parity is the
zero-field label carried along the field ramp used to connect each state to
its finite-field energy.

Parameter sets are the exploratory ones in `scripts/_thf_params.py` (229
quadrupole omitted, 227 deformed-nucleus estimate) -- the same inputs as
`../thf-fields-2026-09-08`, not the package defaults.

## Reproduction

From the repository root, after installing `.[plot]`:

```shell
python scripts/plot_thf_mf_resolved.py
```

This also writes the zoomed companion set in
[../thf-mf-resolved-zoom-2026-09-10](../thf-mf-resolved-zoom-2026-09-10).
