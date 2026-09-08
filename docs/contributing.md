# Working on heff

Use [getting started](getting-started.md) to install the local checkout. The
`test` extra supplies pytest; the `notebooks` extra supplies plotting and
notebook tools. Read [models](models.md) and [architecture](architecture.md)
before adding or changing a physical interaction.

## Make one interpretable change

1. Identify the model, basis, units, and observable affected. Record the primary
   source and any phase or normalization conversion.
2. Make the change in the producing code or TOML file. Keep measured inputs,
   theoretical estimates, and sensitivity choices distinguishable in `Param`.
3. Run the relevant existing checks. For a new operator, add a meaningful
   analytic limit, symmetry check, or independent comparison.
4. Update the current guide and affected generated artifacts. Preserve earlier
   scientific runs as separate parameter/data records when inputs change.

Do not tune a constant or numerical tolerance merely to pass a check. A
Hermitian matrix does not establish correct matrix elements, and a correct
matrix does not establish accurate molecular parameters.

Comments and docstrings describe current behavior, physical assumptions, and
the reasons for non-obvious choices. Keep development history in Git or dated
research records; retain source citations and unresolved scientific questions.

## Test commands

Keep tests that exercise supported behavior, a reproduced bug, or an independent
physical identity. Check the existing coverage before adding another test.
Standalone demonstrations that deliberately corrupt a correct formula are
development diagnostics; they do not need permanent copies beside the direct
regression check. Retain negative cases that test actual validation behavior,
and do not replace signed matrix comparisons with symmetry checks alone.

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

The `.py` files are the editable notebook sources. Regeneration replaces the
corresponding notebook, so preserve any student edits before running it.

```shell
python notebooks/build_tutorial.py
python notebooks/build_isotopologues.py
python -m jupyter nbconvert --to notebook --execute --inplace notebooks/ThF_plus_X3Delta1_Tutorial.ipynb
python -m jupyter nbconvert --to notebook --execute --inplace notebooks/ThF_plus_Isotopologues.ipynb
```

Run these in the environment where `heff` is installed. The kernel must use
that environment; inspect its `sys.executable` if execution cannot import the
package. Saved outputs belong to the parameter sets stated in each notebook.

## Build a local wheel or student archive

```shell
python -m pip wheel --no-deps . --wheel-dir dist
python scripts/package_students.py
```

The wheel contains the Python package and bundled TOML model. The archive
contains the repository material students need: source, tests, examples,
notebooks, guides, figures/data, and supporting Markdown source notes. It omits
Git metadata, caches, machine-local coordination files, build products, and raw
literature/thesis extracts. Its included-file list and source commit are saved
inside the archive. Packaging uses local files and does not push or publish.

## Sources and attribution

Cite the physical sources attached to the operators and parameter sets you use,
and record the heff revision and changed inputs alongside a result. Brown and
Carrington's *Rotational Spectroscopy of Diatomic Molecules* supplies the main
diatomic tensor conventions. Molecule-specific measurements and calculations
are cited in the [Hamiltonian reference](thf-plus-x3delta1-effective-hamiltonian.md).
The restricted amide implementation also records its C2V-Molecules provenance.

The repository currently has no `LICENSE` file. The student archive preserves
that state and does not assign a new license to the code or cited source material.
