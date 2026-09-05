"""Gate A7 / [HAM] V11 -- basis enumerator invariants and blocking partition.

Uniquely catches: enumerator bugs, which are otherwise invisible until a
spectrum looks slightly odd -- double-counting Omega, admitting F outside
|J-I|..J+I, an m_F range off by one, or a blocking that is not a partition.
Every count below is COMPUTED from the enumerator; none is asserted from
docs/thf-plus-x3delta1-effective-hamiltonian.md S1.3, which is what makes this
a check rather than a transcription.

FAIL is reachable and demonstrated in test_invariants_reject_a_bad_basis.
"""
from dataclasses import replace

import numpy as np
import pytest

from heff.spec import (KET_C, ElecState, StateSpec, block_by_mF, blocks_for,
                       check_basis_invariants, enumerate_kets, thf_spec)


def test_ket_dtype_is_named_float64_fields():
    assert KET_C.names == ("J", "Om", "F", "mF")
    assert all(KET_C[n] == np.dtype("f8") for n in KET_C.names)


def test_every_value_is_an_exact_multiple_of_one_half():
    kets = enumerate_kets(thf_spec())
    for name in KET_C.names:
        v = kets[name]
        assert np.array_equal(2.0 * v, np.round(2.0 * v))


def test_basis_invariants_hold():
    spec = thf_spec()
    check_basis_invariants(enumerate_kets(spec), spec)


def test_no_duplicate_kets():
    kets = enumerate_kets(thf_spec())
    rows = {tuple(k) for k in kets}
    assert len(rows) == len(kets)


def test_total_dimension_is_ninety_six_computed_not_asserted():
    """96 = 2 Omega x sum over J of sum over F of (2F+1), recomputed here."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    expected = 0
    for J in range(spec.J_range[0], spec.J_range[1] + 1):
        for twoF in range(int(round(2 * abs(J - spec.I))), int(round(2 * (J + spec.I))) + 1, 2):
            expected += 2 * (twoF + 1)
    assert len(kets) == expected == 96


def test_signed_mF_block_dimensions():
    kets = enumerate_kets(thf_spec())
    blocking = block_by_mF(kets)
    dims = {abs(label): len(idx) for label, idx in blocking.index.items()}
    assert sorted(dims) == [0.5, 1.5, 2.5, 3.5, 4.5]
    assert [dims[m] for m in sorted(dims)] == [16, 14, 10, 6, 2]
    assert 2 * sum(dims[m] for m in sorted(dims)) == len(kets)


def test_blocking_is_a_partition():
    kets = enumerate_kets(thf_spec())
    blocking = block_by_mF(kets)
    allidx = np.concatenate([blocking.index[label] for label in blocking.labels])
    assert np.array_equal(np.sort(allidx), np.arange(len(kets)))
    seen = set()
    for label in blocking.labels:
        idx = set(blocking.index[label].tolist())
        assert not (idx & seen)
        seen |= idx


def test_merge_concatenates_in_label_order():
    kets = enumerate_kets(thf_spec())
    blocking = block_by_mF(kets)
    labels = (0.5, -0.5)
    merged = blocking.merge(labels)
    assert np.array_equal(merged, np.concatenate([blocking.index[0.5], blocking.index[-0.5]]))
    assert set(kets["mF"][merged].tolist()) == {0.5, -0.5}


def test_ket_order_is_deterministic_and_pairs_omega():
    kets = enumerate_kets(thf_spec())
    again = enumerate_kets(thf_spec())
    assert np.array_equal(kets, again)
    for i in range(0, len(kets), 2):
        a, b = kets[i], kets[i + 1]
        assert (a["J"], a["F"], a["mF"]) == (b["J"], b["F"], b["mF"])
        assert {a["Om"], b["Om"]} == {-1.0, 1.0}


def test_invariants_reject_a_bad_basis():
    """FAIL demonstration: a fabricated basis with F outside |J-I| .. J+I."""
    spec = thf_spec()
    bad = np.array([(1.0, 1.0, 2.5, 0.5)], dtype=KET_C)
    with pytest.raises(ValueError, match="triangle"):
        check_basis_invariants(bad, spec)

    worse = np.array([(1.0, 2.0, 1.5, 0.5)], dtype=KET_C)
    with pytest.raises(ValueError, match="Omega"):
        check_basis_invariants(worse, spec)

    thirds = np.array([(1.0, 1.0, 1.5, 1.0 / 3.0)], dtype=KET_C)
    with pytest.raises(ValueError, match="multiple of 0.5"):
        check_basis_invariants(thirds, spec)


def test_blocks_for_blocks_is_the_per_mF_partition():
    """M='blocks' (the default) is exactly block_by_mF -- the collinear fast path."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    got, want = blocks_for(spec, kets), block_by_mF(kets)
    assert got.kind == want.kind == "signed_mF"
    assert got.labels == want.labels
    assert all(np.array_equal(got.index[m], want.index[m]) for m in got.labels)


def test_blocks_for_all_is_one_block_holding_the_whole_basis():
    """M='all' is the basis a Delta-m_F != 0 term has to be built on."""
    spec = replace(thf_spec(), M="all")
    kets = enumerate_kets(spec)
    blocks = blocks_for(spec, kets)
    assert blocks.kind == "all" and blocks.labels == ("all",)
    assert np.array_equal(blocks.index["all"], np.arange(len(kets)))
    assert blocks.n_total == len(kets) == 96


def test_blocks_for_none_names_its_trigger():
    """M='none' is declared but not implemented; the message says what unlocks it."""
    spec = replace(thf_spec(), M="none")
    with pytest.raises(NotImplementedError, match="field-free spectra path"):
        blocks_for(spec, enumerate_kets(spec))


def test_case_other_than_c_is_refused():
    spec = StateSpec(case="a", electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
                     I=0.5, J_range=(1, 4))
    with pytest.raises(NotImplementedError):
        enumerate_kets(spec)
