"""Gate A4 tracking invariants.

The tests require permutation-valued assignments near degeneracies, an
idempotent sign gauge, and a deterministic tie rule for zero amplitudes.
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
    """The greedy strategy may return a non-permutation."""
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
    """Plain argmax(|.|) picks different components when two are
    near-tied with opposite signs; the tie rule (lowest-index near-max) is what
    makes the pin machine-reproducible (C2V matching/_utils.py:453)."""
    v = np.array([[0.5, -0.5 * (1 + 1e-9), 0.1]])
    assert track.pin_col(v)[0] == 1


def test_apply_gauge_none_returns_an_untouched_copy():
    """gauge='none' is a no-op that still returns a NEW array, not an alias."""
    rng = np.random.default_rng(11)
    v = rng.standard_normal((4, 4))
    out = track.apply_gauge(v, "none")
    assert np.array_equal(out, v) and out is not v


def test_apply_gauge_global_is_one_flip_per_point_not_per_vector():
    """`global` applies one +/-1 factor to every state at each point.

    The factor uses the ground state's dominant component at point 0. This
    fixture has three points, two states, and four components; at point 1,
    component 0 vanishes while component 2 is negative. The gauge therefore
    flips the entire state pair using component 2.
    """
    v = np.zeros((3, 2, 4))
    v[0, 0] = [0.1, 0.2, 0.9, 0.1]     # point 0, state 0 (ground): dominant = index 2
    v[0, 1] = [0.9, 0.1, 0.1, 0.1]     # point 0, state 1
    v[1, 0] = [0.0, 0.2, -0.8, 0.1]    # point 1, state 0: component 0 is exactly zero,
    v[1, 1] = [0.7, 0.3, 0.2, 0.1]     # the shared point gauge flips this state too
    v[2, 0] = [0.1, 0.1, 0.9, 0.0]     # point 2, state 0: positive on component 2, no flip
    v[2, 1] = [0.6, 0.1, 0.1, 0.1]
    out = track.apply_gauge(v, "global")
    assert np.array_equal(out[0], v[0])            # point 0 already positive: untouched
    assert np.array_equal(out[1], -v[1])           # point 1 flipped WHOLE-SET
    assert np.array_equal(out[2], v[2])            # point 2 already positive: untouched
    assert out[1, 0, 2] > 0                        # the ground reference amplitude is now positive


def test_apply_gauge_global_requires_a_points_states_components_array():
    """A 2D array has no state axis to flip as a whole set; 'global' must say
    so rather than silently misinterpreting shape (unlike 'pinned'/'none',
    which are well-defined on any shape ending in the component axis)."""
    with pytest.raises(ValueError, match="global"):
        track.apply_gauge(np.eye(4), "global")


def test_apply_gauge_global_on_the_real_mF_three_halves_sweep():
    """On the real ThF+ |m_F|=3/2 block, `global` keeps the ground-state
    reference amplitude positive at every point. A whole-point sign flip is a
    valid eigensolution but violates this gauge contract."""
    from heff.assemble import build_term_matrices
    from heff.engine import sweep
    from heff.params import thf_v1
    from heff.spec import block_by_mF, enumerate_kets, thf_spec
    from heff.terms import ctx_from

    spec = thf_spec()
    kets = enumerate_kets(spec)
    idx = block_by_mF(kets).index[1.5]
    tm = build_term_matrices(kets[idx], ctx_from(spec, thf_v1()))
    E = np.linspace(0.0, 60.0, 5)
    res = sweep(tm, thf_v1(), {"E_z": E, "B_z": np.zeros_like(E)}, gauge="global")

    # res.evecs is (n_points, d, d) with eigenvector k in COLUMN k (engine.py
    # convention); reconstruct the same ground-reference index apply_gauge used.
    ground = res.evecs[:, :, 0]                          # (n_points, d): state 0 at every point
    ref_index = int(np.argmax(np.abs(ground[0])))
    ref_amp = ground[:, ref_index]
    assert np.all(ref_amp >= 0)                          # PASS: the real check

    corrupted = res.evecs.copy()
    corrupted[2] *= -1                                   # fabricate: flip point 2's WHOLE state set
    corrupted_amp = corrupted[:, ref_index, 0]
    assert corrupted_amp[2] < 0
    assert not np.all(corrupted_amp >= 0)                # negative control


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


def _rotation_fixture():
    """2x2 rotation eigenvectors [cosθ, sinθ], [-sinθ, cosθ] over 7 points from
    θ=0 to θ=π/2, crossing 45° BETWEEN grid points (40°, 50°) rather than on
    one -- exact 45° is a genuine argmax tie (cos45==sin45), not a bug, and
    landing on it would make the fixture's outcome depend on floating-point
    rounding rather than on the ordering policy under test."""
    theta = np.array([0, 12, 24, 40, 50, 66, 90]) * np.pi / 180
    v = np.stack([[[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]] for t in theta])
    w = np.zeros((7, 2))
    return w, v


def test_adiabatic_step_and_adiabatic_zero_field_disagree_past_a_crossing():
    """Gate A4: energy order alone cannot tell
    adiabatic_step and adiabatic_zero_field apart. On the rotation fixture,
    adiabatic_step only ever compares consecutive points (max step 24 deg <
    45 deg), so it tracks the physical state continuously and NEVER relabels.
    adiabatic_zero_field compares every point back to θ=0, so it DOES relabel
    once the cumulative rotation passes 45 deg. Both permutations are
    asserted explicitly, not just shown to differ."""
    w, v = _rotation_fixture()
    perm_step = track.order_states(w, v, order="adiabatic_step")
    perm_zero = track.order_states(w, v, order="adiabatic_zero_field")
    assert np.array_equal(perm_step, np.tile([0, 1], (7, 1)))
    assert np.array_equal(perm_zero, [[0, 1]] * 4 + [[1, 0]] * 3)
    assert not np.array_equal(perm_step, perm_zero)


def test_adiabatic_step_carries_the_accumulated_slot_through_a_held_swap():
    """Gate A4: in a [I, P, P] eigenvector stack, P is a
    0<->1 swap. adiabatic_step must reindex the PREVIOUS point by its own
    already-accumulated permutation (`evecs[i-1][:, perm[i-1]]`,
    heff/track.py) before overlapping against the next raw eigenvectors.
    Drop that reindex and point 2 overlaps raw P against raw P: P.T @ P == I
    for any orthogonal P, so the accumulated swap silently reverts to
    [0, 1, 2] instead of holding at [1, 0, 2] -- this fixture fails exactly
    that regression."""
    I = np.eye(3)
    P = np.eye(3)[:, [1, 0, 2]]
    v = np.stack([I, P, P])
    w = np.zeros((3, 3))
    perm = track.order_states(w, v, order="adiabatic_step")
    assert np.array_equal(perm, [[0, 1, 2], [1, 0, 2], [1, 0, 2]])


def test_adiabatic_zero_field_reference_names_which_point_is_zero_field():
    """`reference` selects the grid index that anchors adiabatic_zero_field.
    Anchoring at index 3 (θ=40 deg)
    instead of 0 on the same rotation fixture shifts the crossing so only the
    last point (θ=90 deg, |90-40|=50 deg > 45 deg) relabels."""
    w, v = _rotation_fixture()
    perm = track.order_states(w, v, order="adiabatic_zero_field", reference=3)
    assert np.array_equal(perm, [[0, 1]] * 6 + [[1, 0]])
