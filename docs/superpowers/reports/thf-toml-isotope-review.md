# Independent review: ThF+ TOML isotope integration

Reviewed commit: `84898cc` against baseline `96b4ad4` on
`feature/toml-model-foundation`.

## Verdict

Accepted within the requested Hamiltonian-integration scope. I found no
substantive defect in isotope selection, active parameter metadata, the native
two-spin adapter, or preservation of the existing low-level calculation.

This was a code-and-execution comparison, not a new physics or literature
audit. The implementation preserves the known uncertainty labels rather than
resolving them: in particular, the 229Th quadrupole estimates and OPEN-17
caveat remain estimates, while the 227Th magnetic moment and hyperfine inputs
remain explicitly marked placeholders.

## What I checked

- `select_manifold` applies overrides only to the selected isotope. Parameters
  and conventions merge by name; backend and the active term tuple replace the
  base values. It returns a new frozen manifold and leaves the parsed base
  definition unchanged. Unknown override fields and unknown manifold names are
  rejected.
- `load_model` selects the isotope before resolving the backend and parameters,
  so the case-c2 backend, term list, and parameter set are a coherent execution
  slice rather than a late mutation of the 232 model.
- The case-c2 adapter fixes the native coupling order as J + I_Th = F1 followed
  by F1 + I_F = F, reuses the existing case-c basis/electronic validation, and
  delegates enumeration, context construction, registry terms, and matrix
  elements to the existing low-level code. No `elements_c2` or other producing
  physics source changed in this commit.
- For both 229 and 227 at `J_max=2`, the high-level spec equals `thf_spec`, the
  ket arrays are identical, every selected term matrix is compared element by
  element to a fresh `REGISTRY_C2` build, and zero/finite-field Hamiltonians are
  compared to the native path. The tests compare each active parameter as a
  complete `Param` dataclass, covering value, unit, uncertainty, source,
  status, isotopologue tag, convention tag, and note.
- The 227 slice omits both quadrupole terms and both quadrupole parameters, as
  required for I(227Th) = 1/2. Its `A_par_Th` and `g_N_Th` retain placeholder
  status and the full existing caveat text.
- Selecting 229 does not mutate the default 232 execution slice. The bundled
  default remains 232Th19F on the existing case-c backend and base term list.

## Fresh verification

```text
conda run -n heff python -m pytest \
  tests/test_model_isotope_overrides.py \
  tests/test_backend_case_c2.py \
  tests/test_model_thf_isotopes.py -q

10 passed in 4.57s
```

An independent runtime probe at `J_max=2` also reported:

```text
229 spec_eq True terms 17 active_param_eq True hermitian True
227 spec_eq True terms 15 active_param_eq True hermitian True
```

## Deliberate boundary

The TOML execution slices serialize the parameters consumed by their selected
Hamiltonian terms. They do not expose the three low-level two-photon
polarizability placeholders, and the 229 slice does not expose provenance-only
`Q_Th`; those values remain available through `thf_v2` and are not consumed by
this Hamiltonian model. The exact native comparison therefore covers every
active Hamiltonian parameter, not every key carried by the broader low-level
`ParamSet`. This does not change the basis, matrices, Hamiltonian, or requested
232/229/227 model selection. If a later requirement is to route the separate
two-photon API or provenance-only quantities through model TOML, that is a
separate integration.

Override parameter validation errors currently identify the shared manifold
parameter path rather than the longer isotope-override source path. The value
is still validated and the computation is unaffected; this is diagnostic path
fidelity, not a correctness defect in the reviewed scope.

No load-bearing issue remains unchecked for this scoped verdict. I did not
repeat the full repository suite or wheel-import check; the focused test run
above directly exercises the changed selection, adapter, metadata, and native
matrix-comparison paths.
