"""Small registry of native Hamiltonian backends.

Backends describe an existing implementation; they do not own matrix elements.
"""
from dataclasses import dataclass
from typing import Callable, Mapping


@dataclass(frozen=True)
class Backend:
    id: str
    case: str
    registry: Mapping
    default_terms: tuple[str, ...]
    canonical_units: Mapping[str, str]
    runtime_knobs: frozenset[str]
    optional_zero_parameters: frozenset[str]
    make_spec: Callable
    default_conventions: Callable


_BACKENDS: dict[str, Backend] = {}


def register_backend(backend: Backend) -> Backend:
    if backend.id in _BACKENDS:
        raise ValueError(f"backend {backend.id!r} is already registered")
    _BACKENDS[backend.id] = backend
    return backend


def get_backend(backend_id: str) -> Backend:
    try:
        return _BACKENDS[backend_id]
    except KeyError as exc:
        raise ValueError(
            f"unknown backend {backend_id!r}; available: {sorted(_BACKENDS)}") from exc


def list_backends() -> tuple[str, ...]:
    return tuple(sorted(_BACKENDS))
