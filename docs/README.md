# Documentation

Start with [getting started](getting-started.md) to install `heff` and calculate
a spectrum. Then work through the
[232ThF+ notebook](../notebooks/ThF_plus_X3Delta1_Tutorial.ipynb), where we connect
the basis and Hamiltonian terms to the level structure and field shifts.

| If you want to… | Read |
|---|---|
| Install, diagonalize, and plot your first spectrum | [Getting started](getting-started.md) |
| Choose an isotope/backend and understand its inputs | [Models, units, and conventions](models.md) |
| Learn how model files become matrices and spectra | [Architecture and API map](architecture.md) |
| Check which inputs and approximations limit a calculation | [Scientific limitations](open-questions.md) |
| Change code, refresh notebooks, or package a student copy | [Contributing](contributing.md) |

## Worked material

- [232ThF+ tutorial](https://github.com/arianjad/heff/blob/main/notebooks/ThF_plus_X3Delta1_Tutorial.ipynb): basis,
  Hamiltonian terms, parity, field shifts, observables, and E1 lines.
- [Isotopologues and two-photon tutorial](../notebooks/ThF_plus_Isotopologues.ipynb):
  two coupled nuclear spins and the scope of the closure operator.
- [Three-isotope field figures and data](../results/thf-fields-2026-09-08/README.md):
  J=1–3, fields to 10 kV/cm and 100 G, with plot-specific parameter choices.
- [Synthetic amide model](../examples/models/amide_synthetic.toml): the
  equivalent-proton backend exercised with illustrative coefficients.

## Physics and parameter provenance

[The ThF+ Hamiltonian reference](thf-plus-x3delta1-effective-hamiltonian.md)
contains the equations, phase conventions, source comparisons, and physical
approximations behind the case-(c) operators. The code cites its equation and
OPEN identifiers.

The [parameter audit](superpowers/reports/2026-09-08-thf-estimate-audit.md),
[nuclear-moment review](superpowers/reports/2026-09-08-thf-nuclear-estimate-audit.md),
and [quadrupole review](superpowers/reports/2026-09-08-thf-quadrupole-estimate-audit.md)
provide dated source evidence for current uncertainty labels. The
[amide source comparison](superpowers/reports/amide-source-comparison.md) records
the restricted port's matrix-element comparison.

Use the guides above for the API and supported models. The cited primary
publications contain the underlying measurements and calculations.
