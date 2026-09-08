"""Bundled native Hamiltonian backend adapters."""
from ..backend_registry import get_backend, register_backend
from .case_c import CASE_C_BACKEND


def load_bundled_backends() -> None:
    """Register adapters shipped with heff exactly once per process."""
    try:
        get_backend(CASE_C_BACKEND.id)
    except ValueError:
        register_backend(CASE_C_BACKEND)
