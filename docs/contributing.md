# Working on heff

Use [getting started](getting-started.md) to install the local checkout. The
`test` extra supplies pytest; the `notebooks` extra supplies plotting and
notebook tools. Read [models](models.md) and [architecture](architecture.md)
before adding or changing a physical interaction.

## Changing a calculation

1. Identify the model, basis, units, and observable affected. Record the primary
   source and any phase or normalization conversion.
2. Make the change in the producing code or TOML file. Keep measured inputs,
   theoretical estimates, and sensitivity choices distinguishable in `Param`.
3. Run the relevant existing checks. For a new operator, add a meaningful
   analytic limit, symmetry check, or independent comparison.
4. Update the current guide and affected generated artifacts. Preserve earlier
   scientific runs as separate parameter/data records when inputs change.

A check is useful when it can distinguish a correct calculation from an
incorrect one. Hermiticity alone leaves many matrix-element errors undetected;
even correct matrix elements leave uncertainty in the molecular parameters.
Choose tolerances from the comparison being made, rather than adjusting them
or the constants to obtain a pass.

Keep comments short. Explain a physical assumption, cite a source, or give the
reason for a non-obvious choice. Development history belongs in Git or dated
research records; unresolved scientific questions belong with the relevant
model or convention.

## Test commands

Tests run when you call pytest, independently of ordinary calculations. Keep
checks of supported behavior, reproduced bugs, and independent physical
identities. Before adding a test, check whether an existing one covers the
same failure. Include invalid inputs where rejection matters, and preserve
signed matrix comparisons where a symmetry check would miss a sign error.

```shell
python -m pytest -q
```

Two literature-comparison tests are opt-in. In PowerShell:

```powershell
$env:HEFF_RUN_LITERATURE = "1"
python -m pytest tests/test_observe.py -q
Remove-Item Env:HEFF_RUN_LITERATURE
```

In a macOS/Linux shell:

```shell
HEFF_RUN_LITERATURE=1 python -m pytest tests/test_observe.py -q
```

These comparisons have stated convention and approximation assumptions in
their docstrings. Some measured quantities are also model inputs; reproducing
those numbers is a consistency check rather than independent validation.

## Refresh notebooks

Edit notebook prose in both the `.py` generator and the `.ipynb` notebook.
For prose-only edits, preserve the saved code and outputs. When code or
numerical inputs change, regenerate and execute the affected notebook:

```shell
python notebooks/build_tutorial.py
python notebooks/build_isotopologues.py
python -m jupyter nbconvert --to notebook --execute --inplace notebooks/ThF_plus_X3Delta1_Tutorial.ipynb
python -m jupyter nbconvert --to notebook --execute --inplace notebooks/ThF_plus_Isotopologues.ipynb
```

Regeneration replaces the notebook, so save any student edits first. Run these
commands in the environment where `heff` is installed, with a kernel that uses
the same environment. Inspect `sys.executable` if the kernel cannot import the
package. Saved outputs belong to the parameter sets stated in each notebook.

## Build a local wheel or student archive

```shell
python -m pip wheel --no-deps . --wheel-dir dist
python scripts/package_students.py
```

The wheel contains the Python package and bundled TOML model. The student
archive contains source, tests, examples, notebooks, guides, figures, data, and
supporting Markdown source notes. It omits Git metadata, caches, and build
products. The archive records its included-file list and source commit. Both
commands write local files under `dist/`.

## Sources and attribution

Cite the physical sources attached to the operators and parameter sets you use,
and record the heff revision and changed inputs alongside a result. Brown and
Carrington's *Rotational Spectroscopy of Diatomic Molecules* supplies the main
diatomic tensor conventions. Molecule-specific measurements and calculations
are cited in the [Hamiltonian reference](thf-plus-x3delta1-effective-hamiltonian.md).
The restricted amide implementation also records its C2V-Molecules provenance.
