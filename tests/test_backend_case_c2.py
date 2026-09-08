"""Native two-spin ThF+ case-(c) adapter coverage."""
import numpy as np
import pytest

from heff.assemble import build_term_matrices, hamiltonian
from heff.params import thf_v2
from heff.spec import enumerate_kets, thf_spec
from heff.terms import ctx_from

from heff.backends.case_c2 import CASE_C2_BACKEND


def _thf_inputs():
    return (
        {"J_min": 1, "J_max": 1, "M": "blocks", "frame": "rotating"},
        {"label": "X3Delta1", "Omega": 1.0, "S": 1.0,
         "Lambda": 2.0, "T0": 0.0},
        ({"label": "229Th", "I": 2.5, "couple_to": "J"},
         {"label": "19F", "I": 0.5, "couple_to": "F1"}),
    )


def test_case_c2_adapter_reproduces_the_229thf_low_level_basis_and_hamiltonian():
    """The TOML adapter is a composition layer over the native C2 kernels."""
    basis, electronic, spins = _thf_inputs()
    spec = CASE_C2_BACKEND.make_spec(basis, electronic, spins)
    low_level_spec = thf_spec("229", J_max=1)
    assert spec == low_level_spec

    params = thf_v2("229")
    kets = CASE_C2_BACKEND.enumerate_kets(spec)
    low_level_kets = enumerate_kets(low_level_spec)
    np.testing.assert_array_equal(kets, low_level_kets)

    matrices = build_term_matrices(
        kets, CASE_C2_BACKEND.make_context(spec, params),
        case=CASE_C2_BACKEND.case, registry=CASE_C2_BACKEND.registry)
    low_level_matrices = build_term_matrices(
        low_level_kets, ctx_from(low_level_spec, params),
        case="c2", registry=CASE_C2_BACKEND.registry)
    np.testing.assert_allclose(
        hamiltonian(matrices, params, {"E_z": 12.0, "B_z": 3.0}),
        hamiltonian(low_level_matrices, params, {"E_z": 12.0, "B_z": 3.0}),
        rtol=0.0, atol=0.0)


def test_case_c2_adapter_rejects_a_non_inner_first_spin_chain():
    """The outer fluorine spin must couple to the inner Th-coupled F1."""
    basis, electronic, spins = _thf_inputs()
    malformed = (spins[0], {**spins[1], "couple_to": "J"})
    with pytest.raises(ValueError, match="couple_to"):
        CASE_C2_BACKEND.make_spec(basis, electronic, malformed)


def test_case_c2_adapter_requires_exactly_two_spin_records():
    """A one-spin TOML record cannot select the two-spin native kernels."""
    basis, electronic, spins = _thf_inputs()
    with pytest.raises(ValueError, match="exactly two"):
        CASE_C2_BACKEND.make_spec(basis, electronic, spins[:1])
