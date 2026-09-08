# Documentation

Start with [the repository README](../README.md) and
[getting started](getting-started.md). The guides below describe the code that
is available now.

| If you want to… | Read |
|---|---|
| Install, diagonalize, and plot your first spectrum | [Getting started](getting-started.md) |
| Choose an isotope/backend and understand its inputs | [Models, units, and conventions](models.md) |
| Learn how model files become matrices and spectra | [Architecture and API map](architecture.md) |
| Understand limitations before making a physical claim | [Scientific limitations and current decisions](open-questions.md) |
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
approximations behind the implemented case-(c) operators. Equation and OPEN
identifiers are retained because the code cites them.

The [parameter audit](superpowers/reports/2026-09-08-thf-estimate-audit.md),
[nuclear-moment review](superpowers/reports/2026-09-08-thf-nuclear-estimate-audit.md),
and [quadrupole review](superpowers/reports/2026-09-08-thf-quadrupole-estimate-audit.md)
provide dated source evidence for current uncertainty labels. The
[amide source comparison](superpowers/reports/amide-source-comparison.md) records
the restricted port's matrix-element comparison.

The source reviews are dated evidence for parameter choices and uncertainty
labels. Use the guides above for the supported API and current model scope, and
consult the cited primary publications for the underlying measurements and
calculations.
