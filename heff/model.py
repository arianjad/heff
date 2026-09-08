"""High-level composition of TOML model definitions and heff's numerical APIs."""
from dataclasses import dataclass
from functools import cached_property
from os import PathLike
from pathlib import Path

from .assemble import build_term_matrices, hamiltonian as _hamiltonian
from .backend_registry import get_backend
from .engine import sweep as _sweep
from .model_io import (list_bundled_models, read_bundled_model,
                       read_model_toml, resolve_param_set)
from .spec import enumerate_kets
from .terms import ctx_from, terms_for_case


def _with_source(source, message):
    return f"{source}: {message}"


def _case_insensitive_key(mapping, requested, *, kind, source):
    matches = [key for key in mapping if key.casefold() == str(requested).casefold()]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(_with_source(
        source, f"unknown {kind} {requested!r}; available: {tuple(mapping)}"))


def _bind_backend_path(path, manifold_id, isotope_id):
    scope, *parts = path
    if scope == "manifold":
        bound = f"manifolds.{manifold_id}"
    elif scope == "isotopologue":
        bound = f"isotopologues.{isotope_id}"
    else:
        return f"manifolds.{manifold_id}.backend"
    for part in parts:
        bound += f"[{part}]" if isinstance(part, int) else f".{part}"
    return bound


@dataclass(frozen=True)
class Problem:
    """One immutable, executable full-basis Hamiltonian problem."""

    backend: object
    spec: object
    params: object
    term_names: tuple[str, ...]

    @cached_property
    def kets(self):
        return enumerate_kets(self.spec)

    @cached_property
    def term_matrices(self):
        return build_term_matrices(
            self.kets, ctx_from(self.spec, self.params),
            case=self.backend.case, registry=self.backend.registry,
            term_names=self.term_names)

    def hamiltonian(self, **knobs):
        return _hamiltonian(self.term_matrices, self.params, knobs)

    def sweep(self, **knob_arrays):
        return _sweep(self.term_matrices, self.params, knob_arrays)


@dataclass(frozen=True)
class MoleculeModel:
    """A selected manifold and isotopologue from an immutable definition."""

    definition: object
    manifold: object
    isotopologue: object
    backend: object
    source: str

    def _problem(self, basis_overrides):
        basis = dict(self.manifold.basis)
        basis.update(basis_overrides)
        try:
            spec = self.backend.make_spec(
                basis, self.manifold.electronic, self.isotopologue.spins)
        except (KeyError, TypeError, ValueError) as exc:
            detail = str(exc)
            raw_paths = getattr(exc, "_heff_input_paths", ())
            paths = tuple(_bind_backend_path(
                path, self.manifold.id, self.isotopologue.id)
                for path in raw_paths)
            location = ", ".join(paths) if paths else (
                f"manifolds.{self.manifold.id}.backend")
            raise ValueError(
                _with_source(self.source, f"{location}: {detail}")) from exc

        try:
            params = resolve_param_set(self.manifold, self.backend)
        except ValueError as exc:
            raise ValueError(_with_source(self.source, str(exc))) from exc

        selected = self.manifold.metadata.get("terms", self.backend.default_terms)
        term_names = tuple(selected)
        try:
            terms = terms_for_case(
                self.backend.case, names=term_names, registry=self.backend.registry)
        except ValueError as exc:
            raise ValueError(_with_source(
                self.source, f"manifolds.{self.manifold.id}.terms: {exc}")) from exc

        available = (set(params.params) | set(self.backend.runtime_knobs)
                     | set(self.backend.optional_zero_parameters))
        for term in terms:
            for symbol in term.param:
                if symbol not in available:
                    path = f"manifolds.{self.manifold.id}.parameters.{symbol}"
                    raise ValueError(_with_source(
                        self.source,
                        f"{path}: missing parameter required by active term "
                        f"{term.name!r}"))
        return Problem(self.backend, spec, params, term_names)

    def problem(self, **basis_overrides):
        return self._problem(basis_overrides)

    def validate(self):
        """Validate executable structure without evaluating term matrices."""
        self._problem({})

    def describe(self):
        """Return the selected execution inputs in a compact readable form."""
        problem = self._problem({})
        basis = ", ".join(
            f"{key}={value!r}" for key, value in self.manifold.basis.items())
        terms = ", ".join(problem.term_names)
        return "\n".join((
            f"model: {self.definition.model_id}",
            f"manifold: {self.manifold.id}",
            f"isotopologue: {self.isotopologue.id}",
            f"backend: {self.backend.id}",
            f"basis: {basis}",
            f"terms: {terms}",
            "parameters:",
            problem.params.table(),
        ))


def list_models() -> tuple[str, ...]:
    """Return the bundled model IDs without initializing a backend."""
    return list_bundled_models()


def load_model(model_or_path, *, manifold=None, isotope=None) -> MoleculeModel:
    """Load a bundled model ID or TOML path and select its execution slice."""
    candidate = Path(model_or_path)
    is_path = (isinstance(model_or_path, PathLike) or candidate.exists()
               or candidate.suffix.casefold() == ".toml"
               or candidate.parent != Path("."))
    if is_path:
        source = str(candidate.resolve())
        try:
            definition = read_model_toml(candidate)
        except ValueError as exc:
            raise ValueError(_with_source(source, str(exc))) from exc
    else:
        bundled = list_bundled_models()
        model_id = _case_insensitive_key(
            {name: name for name in bundled}, model_or_path,
            kind="bundled model", source="bundled model catalogue")
        source = f"bundled model {model_id!r}"
        try:
            definition = read_bundled_model(model_id)
        except ValueError as exc:
            raise ValueError(_with_source(source, str(exc))) from exc

    manifold_id = _case_insensitive_key(
        definition.manifolds,
        definition.default_manifold if manifold is None else manifold,
        kind="manifold", source=source)
    isotope_id = _case_insensitive_key(
        definition.isotopologues,
        definition.default_isotopologue if isotope is None else isotope,
        kind="isotopologue", source=source)
    manifold_definition = definition.manifolds[manifold_id]

    from .backends import load_bundled_backends
    load_bundled_backends()
    try:
        backend = get_backend(manifold_definition.backend)
    except ValueError as exc:
        path = f"manifolds.{manifold_id}.backend"
        raise ValueError(_with_source(source, f"{path}: {exc}")) from exc

    return MoleculeModel(
        definition=definition,
        manifold=manifold_definition,
        isotopologue=definition.isotopologues[isotope_id],
        backend=backend,
        source=source,
    )
