"""Native two-spin ThF+ case-(c) adapter coverage."""
import pytest

from heff.backends.case_c2 import CASE_C2_BACKEND


def _thf_inputs():
    return (
        {"J_min": 1, "J_max": 1, "M": "blocks", "frame": "rotating"},
        {"label": "X3Delta1", "Omega": 1.0, "S": 1.0,
         "Lambda": 2.0, "T0": 0.0},
        ({"label": "229Th", "I": 2.5, "couple_to": "J"},
         {"label": "19F", "I": 0.5, "couple_to": "F1"}),
    )


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
