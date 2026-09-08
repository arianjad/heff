# Repository context

Read `README.md` and `docs/README.md` for the current user workflow. The
implementation supports the `case_c`, `case_c2`, and restricted `amide_c2v`
backends; `docs/models.md` describes their supported inputs and limitations.
Inspect live Git state before editing, and preserve unrelated user work.

## Physics and implementation

- Use the implemented basis, units, and convention definitions as the starting
  point. `docs/thf-plus-x3delta1-effective-hamiltonian.md` is the detailed
  diatomic reference; `docs/open-questions.md` states unresolved physical inputs.
- Source any new matrix element or convention conversion. Keep an unresolved
  sign or normalization explicit instead of choosing it silently.
- Preserve parameter provenance and distinguish measured inputs, theory,
  transfers, held-zero interactions, and placeholders. Bundled ThF defaults and
  the saved field-plot parameter sets intentionally differ.
- Keep numerical checks proportional to the change. Verify a relevant analytic
  limit, symmetry, or independent comparison for physical code changes; tests
  do not establish the accuracy of unknown molecular parameters.
- Preserve the `_EIGH` hook for the requested future PyTorch backend work;
  current execution uses NumPy/SciPy. Do not claim a GPU backend exists.

## Documentation and results

- Keep README and the current guides student-facing and consistent with the
  producing code. Use repository-relative paths and portable Python commands.
- Write code comments and docstrings as current contracts and rationale. Keep
  task history and debugging narratives in Git or dated evidence records;
  retain scientific citations, assumptions, and unresolved conventions.
- Edit notebook generators alongside notebook prose. Refresh outputs when
  numerical inputs or code change; preserve outputs for prose-only corrections.
- Preserve source-review evidence and original numerical runs. A changed
  parameter set is a new configuration, not a relabeling of an older result.
- Plans, briefs, dated reviews, and acceptance records describe their original
  scope. Their old task queues and approvals are not current instructions.
- Local changes and checkpoints do not authorize a push, publication, message,
  or change to another repository. Follow the user's current task scope.
