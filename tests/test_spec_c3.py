"""Three-spin case-(c) basis gates for |(((J I1) F1, I2) F2, I3) F, mF>."""

import numpy as np
import pytest

from heff.spec import KET_C2, KET_C3, ElecState, Spin, StateSpec, check_basis_invariants, enumerate_kets, thf_spec


def _three_spin_spec(*, I3=1.5):
    return StateSpec(
        case="c",
        electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
        I=0.5,
        J_range=(1, 2),
        spins=(
            Spin("A", 0.5, "J"),
            Spin("B", 1.0, "F1"),
            Spin("C", I3, "F2"),
        ),
    )


def test_ket_c3_dtype_is_named_float64_fields():
    assert KET_C3.names == ("J", "Om", "F1", "F2", "F", "mF")
    assert all(KET_C3[name] == np.dtype("f8") for name in KET_C3.names)


def test_c3_dimension_matches_the_unsymmetrized_tensor_product_at_nonzero_omega():
    spec = _three_spin_spec()
    kets = enumerate_kets(spec)
    expected = 2 * sum(2 * J + 1 for J in range(spec.J_range[0], spec.J_range[1] + 1))
    for spin in spec.spins:
        expected *= int(round(2 * spin.I + 1))
    assert len(kets) == expected


def test_c3_intermediate_and_final_triangles_hold():
    spec = _three_spin_spec()
    kets = enumerate_kets(spec)
    I1, I2, I3 = (spin.I for spin in spec.spins)
    J, F1, F2, F, mF = (kets[name] for name in ("J", "F1", "F2", "F", "mF"))
    assert np.all(F1 >= np.abs(J - I1)) and np.all(F1 <= J + I1)
    assert np.all(F2 >= np.abs(F1 - I2)) and np.all(F2 <= F1 + I2)
    assert np.all(F >= np.abs(F2 - I3)) and np.all(F <= F2 + I3)
    assert np.all(np.abs(mF) <= F)
    check_basis_invariants(kets, spec)


def test_c3_order_is_deterministic_and_keeps_omega_partners_adjacent():
    spec = _three_spin_spec()
    kets = enumerate_kets(spec)
    assert np.array_equal(kets, enumerate_kets(spec))
    assert np.array_equal(kets, np.sort(kets, order=["J", "F1", "F2", "F", "mF", "Om"]))
    for i in range(0, len(kets), 2):
        a, b = kets[i], kets[i + 1]
        assert (a["J"], a["F1"], a["F2"], a["F"], a["mF"]) == (
            b["J"], b["F1"], b["F2"], b["F"], b["mF"])
        assert {a["Om"], b["Om"]} == {-1.0, 1.0}


def test_zero_outer_spin_reduces_c3_to_the_existing_c2_basis():
    c2_spec = thf_spec("229", J_max=2)
    c3_spec = StateSpec(
        case=c2_spec.case,
        electronic=c2_spec.electronic,
        I=c2_spec.I,
        J_range=c2_spec.J_range,
        spins=(*c2_spec.spins, Spin("spectator", 0.0, "F2")),
    )
    c2, c3 = enumerate_kets(c2_spec), enumerate_kets(c3_spec)
    assert c2.dtype == KET_C2 and c3.dtype == KET_C3
    assert len(c3) == len(c2)
    assert np.array_equal(c3["F2"], c2["F"])
    for name in KET_C2.names:
        assert np.array_equal(c3[name], c2[name])


def test_invariants_reject_an_out_of_triangle_f2():
    spec = _three_spin_spec()
    bad = np.array([(1.0, 1.0, 0.5, 3.5, 2.0, 0.5)], dtype=KET_C3)
    with pytest.raises(ValueError, match="F2"):
        check_basis_invariants(bad, spec)


def test_more_than_three_spins_is_refused_clearly():
    spec = StateSpec(
        case="c",
        electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
        I=0.5,
        J_range=(1, 1),
        spins=(Spin("A", 0.5, "J"), Spin("B", 0.5, "F1"),
               Spin("C", 0.5, "F2"), Spin("D", 0.5, "F3")),
    )
    with pytest.raises(NotImplementedError, match="at most three"):
        enumerate_kets(spec)
