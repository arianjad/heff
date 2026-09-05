"""Gates V16a (basis half of the master reduction) and V17 (v2 enumerator invariants).

V16a uniquely catches: a two-spin enumerator that does not degenerate to the v1
basis at I_Th = 0 -- which would silently invalidate every v1 closed form as a
check on the v2 elements, because the two would no longer be comparing the same
states. FAIL is reachable and was demonstrated by mutation: shortening the F1
loop by one step (stop at `twoJ + twoI1 - 1`) makes the ket-for-ket comparison
in test_I_Th_zero_reduces_to_the_v1_basis fail.

The sort-key mutation originally proposed for V16a -- ordering the v2 kets by
('J','F','mF','Om') instead of ('J','F1','F','mF','Om') -- was PROBED and does
NOT break V16a: at I_Th = 0 the F1 column is identically J, so within a J the
two orders are the same permutation. It breaks test_omega_partners_are_adjacent
on the 229 basis instead, where F1 genuinely varies; that is that test's job,
and it is why the ordering guarantee needs its own gate rather than riding on
the reduction gate.

V17 uniquely catches: an off-by-one in the added F1 loop. The dimension closed
form (2I_Th+1)(2I_F+1) sum_J 2(2J+1) is DERIVED (the F1/F recoupling is a change
of basis at fixed (J, Omega, m_total), [SPEC-v2] S2.2), not transcribed from a
document, so this is a check and not a transcription. FAIL is reachable in
test_invariants_reject_an_out_of_triangle_F1 and was demonstrated by mutation:
deleting the F1 triangle raise from check_basis_invariants makes it fail, so the
F1 check is the only thing that catches that ket.
"""
import numpy as np
import pytest

from heff.spec import (KET_C, KET_C2, ElecState, Spin, StateSpec, block_by_mF,
                       check_basis_invariants, enumerate_kets, thf_spec)


def _closed_form_dim(I_Th, I_F, J_min, J_max):
    """dim = (2I_Th+1)(2I_F+1) sum_J 2(2J+1) -- derived, [SPEC-v2] S2.2."""
    n_spin = int(round((2 * I_Th + 1) * (2 * I_F + 1)))
    return n_spin * sum(2 * (2 * J + 1) for J in range(J_min, J_max + 1))


def test_ket_c2_dtype_is_named_float64_fields():
    assert KET_C2.names == ("J", "Om", "F1", "F", "mF")
    assert all(KET_C2[n] == np.dtype("f8") for n in KET_C2.names)


def test_thf_spec_with_no_argument_is_unchanged():
    """The v1 entry point still yields the v1 96-ket KET_C basis, element for element."""
    assert thf_spec() == thf_spec("232")
    kets = enumerate_kets(thf_spec())
    assert kets.dtype == KET_C
    assert len(kets) == 96
    assert np.array_equal(kets, enumerate_kets(thf_spec("232")))


def test_v2_dimension_matches_the_derived_closed_form():
    for iso, I_Th, J_max in (("229", 2.5, 4), ("227", 0.5, 4), ("229", 2.5, 6)):
        spec = thf_spec(iso, J_max=J_max)
        kets = enumerate_kets(spec)
        assert kets.dtype == KET_C2
        assert len(kets) == _closed_form_dim(I_Th, 0.5, spec.J_range[0], spec.J_range[1])


def test_I_Th_zero_reduces_to_the_v1_basis():
    """V16a: at I_Th = 0 the v2 basis IS the v1 basis, ket for ket."""
    spec = StateSpec(case="c",
                     electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
                     I=0.5, J_range=(1, 4),
                     spins=(Spin("232Th", 0.0, "J"), Spin("19F", 0.5, "F1")))
    kets = enumerate_kets(spec)
    v1 = enumerate_kets(thf_spec())
    assert len(kets) == len(v1)
    assert np.array_equal(kets["F1"], kets["J"])
    for name in KET_C.names:
        assert np.array_equal(kets[name], v1[name])


def test_F1_and_F_triangles_hold():
    spec = thf_spec("229")
    kets = enumerate_kets(spec)
    I_Th, I_F = spec.spins[0].I, spec.spins[1].I
    J, F1, F, mF = (kets[n] for n in ("J", "F1", "F", "mF"))
    assert np.all(F1 >= np.abs(J - I_Th)) and np.all(F1 <= J + I_Th)
    assert np.all(F >= np.abs(F1 - I_F)) and np.all(F <= F1 + I_F)
    assert np.all(np.abs(mF) <= F)
    assert np.all(np.abs(kets["Om"]) <= J)
    check_basis_invariants(kets, spec)


def test_every_value_is_a_multiple_of_one_half():
    kets = enumerate_kets(thf_spec("229"))
    for name in KET_C2.names:
        v = kets[name]
        assert np.array_equal(2.0 * v, np.round(2.0 * v))


def test_blocking_is_still_a_partition_in_the_v2_basis():
    kets = enumerate_kets(thf_spec("229"))
    blocking = block_by_mF(kets)
    allidx = np.concatenate([blocking.index[label] for label in blocking.labels])
    assert np.array_equal(np.sort(allidx), np.arange(len(kets)))
    seen = set()
    for label in blocking.labels:
        idx = set(blocking.index[label].tolist())
        assert not (idx & seen)
        seen |= idx
    labels = (blocking.labels[0], blocking.labels[-1])
    assert np.array_equal(blocking.merge(labels),
                          np.concatenate([blocking.index[l] for l in labels]))


def test_omega_partners_are_adjacent():
    """The pairwise-swap property parity_operator relies on."""
    kets = enumerate_kets(thf_spec("229"))
    assert len(kets) % 2 == 0
    for i in range(0, len(kets), 2):
        a, b = kets[i], kets[i + 1]
        assert (a["J"], a["F1"], a["F"], a["mF"]) == (b["J"], b["F1"], b["F"], b["mF"])
        assert {a["Om"], b["Om"]} == {-1.0, 1.0}


def test_invariants_reject_an_out_of_triangle_F1():
    """FAIL demonstration: F1 = J + I_Th + 1 is outside |J - I_Th| .. J + I_Th."""
    spec = thf_spec("229")
    J, I_Th = 1.0, spec.spins[0].I
    bad = np.array([(J, 1.0, J + I_Th + 1, J + I_Th + 1.5, 0.5)], dtype=KET_C2)
    with pytest.raises(ValueError, match="F1"):
        check_basis_invariants(bad, spec)


def test_a_spin_chain_that_is_not_a_chain_raises():
    with pytest.raises(ValueError, match="couple_to"):
        StateSpec(case="c", electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
                  I=0.5, J_range=(1, 4),
                  spins=(Spin("229Th", 2.5, "J"), Spin("19F", 0.5, "J")))


def test_227_spec_carries_the_tentative_one_half_assignment():
    spec = thf_spec("227")
    assert spec.spins[0].I == 0.5 and spec.spins[1].I == 0.5
    assert "ENSDF" in thf_spec.__doc__ and "9.3 keV" in thf_spec.__doc__


def test_an_unknown_isotopologue_is_refused():
    with pytest.raises(ValueError, match="isotopologue"):
        thf_spec("230")
