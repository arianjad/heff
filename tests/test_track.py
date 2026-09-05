"""Gate A4 -- tracking invariants.

Uniquely catches tracking corruption: an assignment that is not a permutation
(Molecule-Structure's greedy np.argmax(overlap, axis=1) at Energy_Levels.py:1354
can silently duplicate a trace at a near-degeneracy), a sign gauge that is not
idempotent, and the np.sign(0) edge that zeroes an eigenvector
(Energy_Levels.py:312) -- C2V's _PIN_TIE rule is the fix and is lifted, not
re-derived.
"""
import numpy as np
import pytest

from heff import track


def test_assign_returns_a_permutation_even_at_a_near_degeneracy():
    """Adaptive and Hungarian strategies both resolve a near-tied overlap row."""
    overlap = np.array([[0.99, 0.98, 0.0],
                        [0.98, 0.99, 0.0],
                        [0.0, 0.0, 1.0]])
    for strategy in ("adaptive", "hungarian"):
        idx = track.assign(overlap, mode="overlap", strategy=strategy)
        assert sorted(idx.tolist()) == [0, 1, 2], strategy


def test_max_overlap_strategy_can_duplicate_which_is_why_adaptive_is_the_default():
    """FAIL demo: the greedy strategy is allowed to return a non-permutation."""
    overlap = np.array([[0.9, 0.8], [0.95, 0.1]])
    idx = track.assign(overlap, mode="overlap", strategy="max_overlap")
    assert sorted(idx.tolist()) == [0, 0]


def test_assign_indices_are_ints_not_vectors():
    """The carry across a match step must be an int index array, not eigenvectors."""
    overlap = np.eye(4)
    idx = track.assign(overlap, mode="overlap")
    assert idx.dtype.kind in "iu"


def test_pin_col_is_idempotent_and_handles_an_exact_zero_amplitude():
    """A second application of the pinned gauge changes nothing (idempotence)."""
    rng = np.random.default_rng(3)
    v = rng.standard_normal((6, 5))
    v[0, :] = 0.0
    v[0, 2] = -1.0
    once = track.apply_gauge(v, "pinned")
    twice = track.apply_gauge(once, "pinned")
    assert np.allclose(once, twice)
    assert np.all(track.pin_col(once) == 1)


def test_pin_col_tie_rule_is_deterministic():
    """_PIN_TIE = 1e-6: plain argmax(|.|) picks different components when two are
    near-tied with opposite signs; the tie rule (lowest-index near-max) is what
    makes the pin machine-reproducible (C2V matching/_utils.py:453)."""
    v = np.array([[0.5, -0.5 * (1 + 1e-9), 0.1]])
    assert track.pin_col(v)[0] == 1
    assert track._PIN_TIE == 1e-6


def test_apply_gauge_none_returns_an_untouched_copy():
    """gauge='none' is a no-op that still returns a NEW array, not an alias."""
    rng = np.random.default_rng(11)
    v = rng.standard_normal((4, 4))
    out = track.apply_gauge(v, "none")
    assert np.array_equal(out, v) and out is not v


def test_order_states_energy_is_the_identity_and_is_sorted():
    """order='energy' is the identity permutation because eigh already sorts."""
    rng = np.random.default_rng(5)
    w = np.sort(rng.standard_normal((6, 4)), axis=1)
    v = np.tile(np.eye(4), (6, 1, 1))
    perm = track.order_states(w, v, order="energy")
    assert np.array_equal(perm, np.tile(np.arange(4), (6, 1)))


def test_order_states_adiabatic_step_follows_a_clean_swap():
    """Two levels cross: energy order relabels them, adiabatic_step does not."""
    v = np.stack([np.eye(3), np.eye(3), np.eye(3)[:, [1, 0, 2]]])
    w = np.tile(np.array([0.0, 1.0, 2.0]), (3, 1))
    perm = track.order_states(w, v, order="adiabatic_step")
    assert np.array_equal(perm[0], [0, 1, 2])
    assert np.array_equal(perm[1], [0, 1, 2])
    assert np.array_equal(perm[2], [1, 0, 2])
    for i in range(3):
        assert sorted(perm[i].tolist()) == [0, 1, 2]


def test_order_states_adiabatic_zero_field_references_the_first_point():
    """adiabatic_zero_field matches every point against point 0, not its neighbor."""
    n = 4
    v = np.tile(np.eye(3), (n, 1, 1))
    v[-1] = np.eye(3)[:, [1, 0, 2]]
    w = np.tile(np.array([0.0, 1.0, 2.0]), (n, 1))
    perm = track.order_states(w, v, order="adiabatic_zero_field")
    assert np.array_equal(perm[0], np.arange(3))
    assert np.array_equal(perm[-1], np.array([1, 0, 2]))


def test_unknown_order_raises():
    """An unrecognised order string raises rather than silently defaulting."""
    with pytest.raises(ValueError, match="order"):
        track.order_states(np.zeros((2, 2)), np.zeros((2, 2, 2)), order="vibes")
