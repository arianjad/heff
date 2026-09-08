# heff — molecular effective Hamiltonians

`heff` builds and diagonalizes effective Hamiltonians for molecular rotation,
hyperfine structure, and external fields. Use it to learn how a model becomes
an energy spectrum, inspect individual interactions, and explore Stark and
Zeeman shifts with traceable parameter choices.

The central idea is

```text
H/h = sum_k c_k M_k       (frequencies in MHz)
```

The matrices `M_k` describe the operators in a chosen basis. The coefficients
`c_k` contain molecular parameters and field strengths. Changing a field or
parameter reuses the matrices; changing the basis or operator model rebuilds them.

## Start here

1. **Install** using Python 3.12 or newer. From the repository root, in your
   chosen Python environment:

   ```shell
   python -m pip install -e ".[notebooks,test]"
   ```

   [Environment setup for Windows, macOS, and Linux](docs/getting-started.md#1-install-in-a-python-environment).

2. **Calculate** the 232ThF+ spectrum:

   ```python
   import numpy as np
   from heff import load_model

   model = load_model("thf_plus")
   problem = model.problem(J_max=3)
   H = problem.hamiltonian(E_z=20.0, B_z=0.01)
   energies, states = np.linalg.eigh(H)

   print(H.shape)  # (60, 60)
   print(energies[:6] - energies[0])  # MHz above the lowest level
   ```

   Fields are **V/cm** and **gauss**. Columns of `states` are eigenvectors in
   `problem.kets`. The small cutoff is an example, not a convergence guarantee.

3. **Explore** [your first field sweep](docs/getting-started.md#3-plot-a-stark-sweep),
   then open [the 232ThF+ tutorial](notebooks/ThF_plus_X3Delta1_Tutorial.ipynb).
   Run `python -m jupyterlab` from this environment to work through it.

4. **Check assumptions** in [models and units](docs/models.md) before comparing
   to a measurement. Inspect `model.describe()` for parameter sources and status.

## What is available

| Model | What you can calculate | Physical scope |
|---|---|---|
| 232Th19F+ X 3Delta1 | Rotation, hyperfine/parity structure, Stark/Zeeman shifts, E1 transitions, selected observables | Effective Omega=±1 manifold; measured and adopted inputs with stated approximations |
| 229Th19F+ and 227Th19F+ | Two coupled nuclear spins, Th hyperfine, applicable quadrupole operators, field shifts | Exploratory isotope models; several inputs are transferred estimates or placeholders |
| Equivalent-proton amides | Asymmetric rotation, electron/nuclear spin interactions, exchange-filtered proton pair, optional metal spin | Restricted `amide_c2v` backend; the provided coefficients are synthetic |

Select an isotope with
`load_model("thf_plus", isotope="229Th19F")`, or load a custom TOML file with
`load_model("examples/models/amide_synthetic.toml")`. See
[the model guide](docs/models.md) for coupling order, supported terms, and units.
TOML describes inputs for existing backends; it does not generate missing physics.

The lower-level API also exposes E1 line strengths and rank-0/rank-2 effective
two-photon operators. Their validity conditions and parameter limitations are
explained in [the isotope tutorial](notebooks/ThF_plus_Isotopologues.ipynb).
Numerical diagonalization currently uses NumPy/SciPy on the CPU.

## Worked isotope figures

[Open the level, Stark, and Zeeman figure set](results/thf-fields-2026-09-08/thf-isotopes-levels-stark-zeeman.pdf):
J=1–3, electric fields to 10 kV/cm, and magnetic fields to 100 G.

![Zero-field rotational structure of three ThF+ isotopes](results/thf-fields-2026-09-08/thf-level-overview.png)

The [figure notes and data](results/thf-fields-2026-09-08/README.md) record the
actual inputs, state labels, and basis-convergence checks. These plotting
parameter sets are deliberately separate from the bundled defaults:

- The 229Th baseline omits unknown quadrupoles; an extra page shows an
  unvalidated sensitivity case.
- The 227Th plots use a deformed-nucleus theory estimate. The bundled 227Th
  model and advanced notebook retain the much larger Schmidt stress-test value.

## Before interpreting a result

Parameter status matters: a placeholder can produce a smooth, converged plot
without being a reliable molecular prediction. Odd-isotope constants are not
all independently measured or fitted. Increase the rotational cutoff for the
states and fields you study, and distinguish energy ordering from eigenstate
tracking at crossings. [Current scientific limitations](docs/open-questions.md)
and the [Hamiltonian reference](docs/thf-plus-x3delta1-effective-hamiltonian.md)
explain these boundaries and the source conventions.

## Tests and repository map

```shell
python -m pytest -q
```

The default suite checks assembly, symmetries, analytic limits, parameter
metadata, and model loading. Two literature comparisons are opt-in; commands
are in [contributing](docs/contributing.md). Tests do not certify unknown inputs.

| Directory | Purpose |
|---|---|
| `heff/` | Python implementation and bundled model data |
| `examples/` and `notebooks/` | Runnable models and teaching material |
| `docs/` | Current guides, physics references, and supporting source records |
| `results/` | Reproducible figures, parameter records, and numerical data |
| `tests/` and `scripts/` | Checks, plot reproduction, and local student packaging |

Use the [documentation index](docs/README.md) for the learning path and the
[architecture guide](docs/architecture.md) to navigate the implementation.
Development plans and handoffs are supporting records, not setup instructions.
To create a clean local copy to share with students, run
`python scripts/package_students.py`; see [what it includes](docs/contributing.md#build-a-local-wheel-or-student-archive).
