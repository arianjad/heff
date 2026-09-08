# Getting started

This guide takes you from a fresh checkout of `heff` to an energy spectrum and
a field sweep. You need basic Python, matrix diagonalization, and
angular-momentum notation.

## 1. Clone and install

Use Python 3.12 or newer. Clone the repository, enter its root directory, and
create a virtual environment:

```shell
git clone https://github.com/arianjad/heff.git
cd heff
python -m venv .venv
```

If your system provides Python as `python3`, use that command instead.

Activate it in PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Or in a macOS/Linux shell:

```shell
source .venv/bin/activate
```

Install the package, plotting/notebook tools, and tests:

```shell
python -m pip install -e ".[notebooks,test]"
python -m pytest -q
```

The editable install makes changes to local `heff/` files visible to Python.
Dependencies are declared in [pyproject.toml](../pyproject.toml). For calculations
without plots or notebooks, `python -m pip install -e .` is sufficient. For
scripts with plots, use `python -m pip install -e ".[plot]"`.

If PowerShell does not allow activation, use the environment's interpreter
directly, without changing your execution policy:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[notebooks,test]"
.\.venv\Scripts\python.exe -m pytest -q
```

An existing Conda environment with Python >=3.12 also works; activate it and run
the same `python -m pip` commands. Install this local checkout with `-e .`
so you use the code and model files described in this guide.

## 2. Calculate a spectrum

```python
import numpy as np
from heff import load_model

model = load_model("thf_plus")  # 232Th19F+, X 3Delta1
problem = model.problem(J_max=3)
H = problem.hamiltonian(E_z=20.0, B_z=0.01)
energies, vectors = np.linalg.eigh(H)

print(H.shape)  # (60, 60)
print(energies[:6] - energies[0])  # MHz above the lowest state
print(model.describe())  # basis, active terms, parameters, provenance
```

The matrix represents `H/h` in **MHz**. `E_z` is in **V/cm** and `B_z` in
**gauss**; 10 kV/cm is `E_z=10000`, not `10`. Columns of `vectors` are
eigenstates expressed in the ordered basis `problem.kets`. A common energy
offset has no effect on a transition frequency.

The basis contains all magnetic sublevels in the requested J range. The
high-level `Problem` assembles the full matrix even when the basis metadata
says `M="blocks"`; the low-level API provides explicit `block_by_mF` handling
for larger calculations. Collinear static fields preserve signed mF.

`J_max=3` is a small starting basis, not a convergence guarantee. Increase the
cutoff and compare the states you care about before reporting precise shifts.

## 3. Plot a Stark sweep

```python
import matplotlib.pyplot as plt

fields = np.linspace(0.0, 100.0, 101)  # V/cm
result = problem.sweep(E_z=fields, B_z=np.zeros_like(fields))

fig, ax = plt.subplots()
ax.plot(fields, result.evals[:, :12] - result.evals[0, 0])
ax.set(xlabel="Electric field (V/cm)", ylabel="Energy / h (MHz, common zero)")
fig.savefig("first_stark.png", dpi=180, bbox_inches="tight")
```

This displays the twelve lowest eigenvalues at each field, relative to one
common zero-field energy. `Problem.sweep` defaults to energy ordering: a line
index is an energy rank, not a tracked quantum-state identity. The low-level
`heff.sweep` supports explicit tracking and gauge choices; see
[architecture](architecture.md) and the tutorials for that workflow.

The ready-made [isotope figure set](../results/thf-fields-2026-09-08/README.md)
uses a larger basis, branch labels, and numerical convergence checks. It also
documents its parameter overrides; its odd-isotope plots are not identical
to loading the bundled defaults.

## 4. Explore a different model

```python
odd_model = load_model("thf_plus", isotope="229Th19F")
odd_problem = odd_model.problem(J_max=2)
print(odd_model.describe())

amide = load_model("examples/models/amide_synthetic.toml").problem()
amide_H = amide.hamiltonian(E_z=1.0, B_z=0.1)
```

The 229Th model contains uncalibrated quadrupole placeholders. The bundled
227Th model retains a Schmidt stress-test moment, whereas the saved field
plots use a deformed-nucleus theory estimate. The amide example contains
synthetic coefficients. Read [models and units](models.md) before interpreting
any of these as a spectrum of a real isotope.

To change a parameter without changing the bundled file, make a new problem:

```python
from dataclasses import replace

old = problem.params.params["c_I"]
without_spin_rotation = replace(
    problem,
    params=problem.params.with_(c_I=replace(
        old, value=0.0, status="held-fixed", uncertainty=None,
        source="Sensitivity comparison", note="Explicitly omitted interaction")),
)
```

This retains the unit and creates new metadata for the chosen assumption.
For a reusable model file, copy [thf_plus.toml](../heff/models/thf_plus.toml)
or the [amide example](../examples/models/amide_synthetic.toml), edit it, and
pass its path to `load_model`.

## Tutorials and saved results

Start Jupyter from the installed environment:

```shell
python -m jupyterlab
```

Open [ThF_plus_X3Delta1_Tutorial.ipynb](../notebooks/ThF_plus_X3Delta1_Tutorial.ipynb)
first. It explains the low-level basis, matrices, parity, observables, and E1
transitions. Then use
[ThF_plus_Isotopologues.ipynb](../notebooks/ThF_plus_Isotopologues.ipynb) for two
nuclear spins and effective two-photon operators. Select the kernel belonging
to your environment; `import sys; print(sys.executable)` identifies it.

To redraw the committed field figures from saved data:

```shell
python scripts/plot_thf_isotopes.py --render-only
```

Omit `--render-only` to recompute them. This overwrites the corresponding
generated figures/data; copy the result folder first if you are preserving a
separate experiment. The script needs no machine-specific launcher.

## When an example fails

| Symptom | Check |
|---|---|
| `No module named heff` | Install from this repository root using the same Python that runs the example. |
| Unknown model or missing TOML path | Use `heff.list_models()` for bundled names; custom paths are relative to the current directory. |
| Missing plotting or notebook module | Install the `plot` or `notebooks` extra in the active environment. |
| Slow odd-isotope calculation | Reduce the exploratory cutoff; use explicit mF blocks for large sweeps. Keep a separate convergence check for reported results. |
| Unexpected ordering at a crossing | Inspect `result.order` and state character; energy order does not follow an eigenvector through a crossing. |

Tests check code behavior and specified physical identities. Two comparisons to
published numbers are opt-in; the [contributing guide](contributing.md) gives
the commands. Passing tests does not establish the accuracy of an estimated
Hamiltonian parameter.
