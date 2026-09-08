"""Bundled native Hamiltonian backend adapters."""
from ..backend_registry import get_backend, register_backend
from .case_c import CASE_C_BACKEND


def load_bundled_backends() -> None:
    """Register adapters shipped with heff exactly once per process."""
    from .amide import AMIDE_BACKEND
    for backend in (CASE_C_BACKEND, AMIDE_BACKEND):
        try:
            get_backend(backend.id)
        except ValueError:
            register_backend(backend)
