import pytest

from heff.backend_registry import Backend, get_backend, list_backends
from heff.backends import load_bundled_backends


def test_case_c_backend_is_registered_by_stable_id():
    """Catches a missing/wrong bundled case-(c) registration."""
    load_bundled_backends()
    backend = get_backend("case_c")
    assert backend.id == "case_c"
    assert backend.case == "c"
    assert "rotation" in backend.default_terms
    assert backend.canonical_units["B0"] == "MHz"


def test_unknown_backend_names_the_available_ids():
    """Catches an unknown-backend error that omits usable registered IDs."""
    load_bundled_backends()
    with pytest.raises(ValueError, match="case_c"):
        get_backend("case_z")


def test_case_c_adapter_reproduces_thf_state_spec():
    """Catches incorrect mapping of the v1 TOML-shaped data onto StateSpec."""
    load_bundled_backends()
    backend = get_backend("case_c")
    spec = backend.make_spec(
        {"J_min": 1, "J_max": 4, "M": "blocks", "frame": "rotating"},
        {"label": "X3Delta1", "Omega": 1.0, "S": 1.0,
         "Lambda": 2.0, "T0": 0.0},
        ({"label": "19F", "I": 0.5, "couple_to": "J"},),
    )
    from heff.spec import thf_spec
    assert spec == thf_spec()


def test_case_c_adapter_rejects_an_unknown_basis_key():
    """Catches a misspelled basis input being silently ignored."""
    load_bundled_backends()
    backend = get_backend("case_c")
    with pytest.raises(ValueError, match="J_stop"):
        backend.make_spec(
            {"J_min": 1, "J_max": 4, "M": "blocks", "frame": "rotating",
             "J_stop": 5},
            {"label": "X3Delta1", "Omega": 1.0, "S": 1.0,
             "Lambda": 2.0, "T0": 0.0},
            ({"label": "19F", "I": 0.5, "couple_to": "J"},),
        )


def test_case_c_adapter_rejects_a_reversed_j_range():
    """Catches an invalid rotational truncation before StateSpec construction."""
    load_bundled_backends()
    backend = get_backend("case_c")
    with pytest.raises(ValueError, match="J_min"):
        backend.make_spec(
            {"J_min": 4, "J_max": 1, "M": "blocks", "frame": "rotating"},
            {"label": "X3Delta1", "Omega": 1.0, "S": 1.0,
             "Lambda": 2.0, "T0": 0.0},
            ({"label": "19F", "I": 0.5, "couple_to": "J"},),
        )
