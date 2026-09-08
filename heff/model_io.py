"""Lean TOML model definitions and backend-specific parameter resolution."""
from dataclasses import dataclass, field, fields, replace
from importlib import resources
from pathlib import Path
from types import MappingProxyType
from collections.abc import Mapping
import tomllib
import warnings

from .params import Param, ParamSet, STATUSES


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _empty_mapping():
    return MappingProxyType({})


@dataclass(frozen=True)
class LoadedParameter:
    """One unconverted parameter value, retained until a backend resolves it."""
    raw: object
    metadata: Mapping = field(default_factory=_empty_mapping)


@dataclass(frozen=True)
class IsotopologueDefinition:
    id: str
    spins: tuple[Mapping, ...]
    metadata: Mapping = field(default_factory=_empty_mapping)


@dataclass(frozen=True)
class ManifoldDefinition:
    id: str
    backend: str
    electronic: object
    basis: Mapping
    parameters: Mapping[str, LoadedParameter]
    conventions: Mapping = field(default_factory=_empty_mapping)
    convention_paths: Mapping = field(default_factory=_empty_mapping)
    metadata: Mapping = field(default_factory=_empty_mapping)


@dataclass(frozen=True)
class ModelDefinition:
    schema_version: int
    model_id: str
    default_manifold: str
    default_isotopologue: str
    isotopologues: Mapping[str, IsotopologueDefinition]
    manifolds: Mapping[str, ManifoldDefinition]
    metadata: Mapping = field(default_factory=_empty_mapping)


_TOP_LEVEL_KEYS = frozenset({
    "schema_version", "model_id", "default_manifold", "default_isotopologue",
    "isotopologues", "manifolds", "conventions",
})
_MANIFOLD_KEYS = frozenset({
    "backend", "electronic", "basis", "parameters", "conventions",
})


def _required(record, key, path):
    if key not in record:
        location = f"{path}.{key}" if path else key
        raise ValueError(f"missing required key {location}")
    return record[key]


def _mapping(value, path):
    if not isinstance(value, Mapping):
        raise ValueError(f"{path} must be a table")
    return value


def _loaded_parameter(raw):
    metadata = {}
    if isinstance(raw, Mapping):
        known = {"value", "unit", "uncertainty", "status", "source",
                 "isotopologue", "convention", "note"}
        metadata = {key: value for key, value in raw.items() if key not in known}
    return LoadedParameter(raw=_freeze(raw), metadata=_freeze(metadata))


def read_model_toml(path) -> ModelDefinition:
    """Read a model file without treating descriptive provenance as required."""
    path = Path(path)
    with path.open("rb") as stream:
        document = tomllib.load(stream)

    schema_version = _required(document, "schema_version", "")
    if schema_version != 1:
        raise ValueError(f"unsupported schema_version {schema_version!r}; expected 1")
    model_id = _required(document, "model_id", "")
    default_manifold = _required(document, "default_manifold", "")
    default_isotopologue = _required(document, "default_isotopologue", "")

    global_conventions = _mapping(document.get("conventions", {}), "conventions")
    global_convention_paths = {
        key: f"conventions.{key}" for key in global_conventions
    }
    isotopologue_records = _mapping(
        _required(document, "isotopologues", ""), "isotopologues")
    isotopologues = {}
    for isotope_id, raw_isotope in isotopologue_records.items():
        isotope_path = f"isotopologues.{isotope_id}"
        isotope = _mapping(raw_isotope, isotope_path)
        spins = _required(isotope, "spins", isotope_path)
        if not isinstance(spins, list):
            raise ValueError(f"{isotope_path}.spins must be an array of tables")
        isotopologues[isotope_id] = IsotopologueDefinition(
            id=isotope_id,
            spins=tuple(_freeze(_mapping(spin, f"{isotope_path}.spins")) for spin in spins),
            metadata=_freeze({key: value for key, value in isotope.items() if key != "spins"}),
        )

    manifold_records = _mapping(_required(document, "manifolds", ""), "manifolds")
    manifolds = {}
    for manifold_id, raw_manifold in manifold_records.items():
        manifold_path = f"manifolds.{manifold_id}"
        manifold = _mapping(raw_manifold, manifold_path)
        backend = _required(manifold, "backend", manifold_path)
        electronic = _required(manifold, "electronic", manifold_path)
        basis = _mapping(_required(manifold, "basis", manifold_path), f"{manifold_path}.basis")
        if backend == "case_c":
            for key in ("J_min", "J_max", "M", "frame"):
                _required(basis, key, f"{manifold_path}.basis")
        parameter_records = _mapping(
            _required(manifold, "parameters", manifold_path),
            f"{manifold_path}.parameters")
        manifold_conventions = _mapping(manifold.get("conventions", {}),
                                        f"{manifold_path}.conventions")
        conventions = dict(global_conventions) | dict(manifold_conventions)
        convention_paths = dict(global_convention_paths)
        convention_paths.update({
            key: f"{manifold_path}.conventions.{key}"
            for key in manifold_conventions
        })
        manifolds[manifold_id] = ManifoldDefinition(
            id=manifold_id,
            backend=backend,
            electronic=_freeze(electronic),
            basis=_freeze(basis),
            parameters=_freeze({symbol: _loaded_parameter(raw)
                                for symbol, raw in parameter_records.items()}),
            conventions=_freeze(conventions),
            convention_paths=_freeze(convention_paths),
            metadata=_freeze({key: value for key, value in manifold.items()
                              if key not in _MANIFOLD_KEYS}),
        )

    return ModelDefinition(
        schema_version=schema_version,
        model_id=model_id,
        default_manifold=default_manifold,
        default_isotopologue=default_isotopologue,
        isotopologues=_freeze(isotopologues),
        manifolds=_freeze(manifolds),
        metadata=_freeze({key: value for key, value in document.items()
                          if key not in _TOP_LEVEL_KEYS}),
    )


def list_bundled_models() -> tuple[str, ...]:
    """Return packaged model IDs without requiring backend initialization."""
    models = resources.files("heff").joinpath("models")
    if not models.is_dir():
        return ()
    return tuple(sorted(resource.name[:-5] for resource in models.iterdir()
                        if resource.name.endswith(".toml")))


def read_bundled_model(model_id) -> ModelDefinition:
    """Read one packaged TOML model by its stable identifier."""
    resource = resources.files("heff").joinpath("models", f"{model_id}.toml")
    if not resource.is_file():
        raise ValueError(f"unknown bundled model {model_id!r}; available: {list_bundled_models()}")
    with resources.as_file(resource) as path:
        return read_model_toml(path)


def _parameter(raw, *, path, canonical_unit):
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        return Param(float(raw), canonical_unit)
    if not isinstance(raw, Mapping) or "value" not in raw:
        raise ValueError(
            f"{path} must be a number or a table containing value")
    allowed = {"value", "unit", "uncertainty", "status", "source",
               "isotopologue", "convention", "note"}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        warnings.warn(f"{path}: preserving unknown metadata {unknown}")
    status = raw.get("status", "unspecified")
    if status not in STATUSES:
        warnings.warn(
            f"{path}: unknown optional status {status!r}; "
            "using 'unspecified'")
        status = "unspecified"
    try:
        param = Param(
            float(raw["value"]),
            unit=raw.get("unit", canonical_unit),
            uncertainty=raw.get("uncertainty"),
            status=status,
            source=raw.get("source"),
            isotopologue=raw.get("isotopologue"),
            convention=raw.get("convention"),
            note=raw.get("note"),
        )
        _ = param.canonical
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path}: {exc}") from exc
    return param


def resolve_param_set(manifold, backend) -> ParamSet:
    """Normalize one manifold's parameter records against its chosen backend."""
    if manifold.backend != backend.id:
        raise ValueError(
            f"manifolds.{manifold.id}.backend is {manifold.backend!r}, not {backend.id!r}")
    defaults = backend.default_conventions()
    recognized = {item.name for item in fields(defaults)}
    unknown = sorted(set(manifold.conventions) - recognized)
    if unknown:
        key = unknown[0]
        path = manifold.convention_paths.get(
            key, f"manifolds.{manifold.id}.conventions.{key}")
        raise ValueError(f"{path}: unknown convention key {key!r}")
    conventions = defaults
    for key, value in manifold.conventions.items():
        path = manifold.convention_paths.get(
            key, f"manifolds.{manifold.id}.conventions.{key}")
        try:
            conventions = replace(conventions, **{key: value})
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}: {exc}") from exc
    params = {
        symbol: _parameter(loaded.raw,
                           path=f"manifolds.{manifold.id}.parameters.{symbol}",
                           canonical_unit=backend.canonical_units.get(symbol, "MHz"))
        for symbol, loaded in manifold.parameters.items()
    }
    tagged = {
        symbol: parameter.convention
        for symbol, parameter in params.items()
        if parameter.convention
    }
    tags = set(tagged.values())
    if len(tags) > 1:
        paths = ", ".join(
            f"manifolds.{manifold.id}.parameters.{symbol}.convention"
            for symbol in tagged)
        raise ValueError(
            f"{paths}: conflicting convention tags in one ParamSet: "
            f"{sorted(tags)}")
    try:
        return ParamSet(params, conventions)
    except ValueError as exc:
        if "d_mf" in params:
            parameter_path = f"manifolds.{manifold.id}.parameters.d_mf"
            convention_path = manifold.convention_paths.get(
                "dipole_origin", "conventions.dipole_origin")
            raise ValueError(
                f"{parameter_path}, {convention_path}: {exc}") from exc
        raise ValueError(f"manifolds.{manifold.id}.parameters: {exc}") from exc
