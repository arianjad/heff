"""Task 9 -- E1 spectra in the two-spin basis: `spectra.dipole_matrix`'s
`geometry=` keyword, F1 in `label_lines`, and the two-spin gates.

No new algebra: `geometry=None` keeps `heff.elements_c.dipole_geometry` (v1);
a v2 caller passes `functools.partial(heff.elements_c2.axial_geometry, k=1,
q=0.0)` (E1 within a fixed-Omega block has Delta Omega = 0).
"""
from functools import partial

import numpy as np
import pytest

from heff.assemble import build_term_matrices, hamiltonian
from heff.elements_c import dipole_geometry
from heff.elements_c2 import REGISTRY_C2, axial_geometry
from heff.params import thf_v1, thf_v2
from heff.spec import (Spin, StateSpec, block_by_mF, enumerate_kets,
                       thf_spec)
from heff.spectra import dipole_matrix, label_lines, line_strengths
from heff.terms import ctx_from

GEO = partial(axial_geometry, k=1, q=0.0)


def _v2_basis_at_I_Th_zero(J_max=4):
    """Same states as v1's `thf_spec(J_max=)`, KET_C2-shaped, I_Th = 0.

    Needed because `axial_geometry` reads `ket["F1"]` and `ctx.spins[0]`
    (heff/elements_c2.py `_spins`): the plain v1 `thf_spec()`/`ctx_from`
    (empty `ctx.spins`) has neither, so the reduction can only be exercised on
    this I_Th = 0 two-spin basis, per tests/test_elements_c2_reduction.py's
    `spec_with_I_Th_zero` (same construction, ket-for-ket identical to v1 --
    that file's own `test_the_two_bases_are_the_same_states_ket_for_ket`).
    """
    v1 = thf_spec(J_max=J_max)
    return StateSpec(case="c", electronic=v1.electronic, I=0.5,
                     J_range=v1.J_range, v=0, M="blocks", frame="rotating",
                     spins=(Spin(label="232Th", I=0.0, couple_to="J"),
                            Spin(label="19F", I=0.5, couple_to="F1")))


def test_dipole_matrix_default_is_byte_identical_to_v1():
    """`geometry=None` is exactly `elements_c.dipole_geometry` (atol=0), and
    `geometry=partial(axial_geometry, k=1, q=0.0)` agrees with it to 1e-12 on
    the I_Th = 0 two-spin basis (the analytic collapse, exercised through
    `spectra.dipole_matrix` rather than `axial_geometry` directly)."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    blocks = block_by_mF(kets)
    a, b = kets[blocks.index[0.5]], kets[blocks.index[1.5]]

    for p in (-1, 0, 1):
        got = dipole_matrix(a, b, ctx, p)
        want = np.zeros((len(a), len(b)))
        for i in range(len(a)):
            for j in range(len(b)):
                want[i, j] = dipole_geometry(a[i], b[j], ctx.I, p)
        assert np.array_equal(got, want)

    spec2 = _v2_basis_at_I_Th_zero()
    kets2 = enumerate_kets(spec2)
    ctx2 = ctx_from(spec2, thf_v2("232"))
    blocks2 = block_by_mF(kets2)
    a2, b2 = kets2[blocks2.index[0.5]], kets2[blocks2.index[1.5]]
    for p in (-1, 0, 1):
        default = dipole_matrix(a2, b2, ctx2, p)
        via_axial = dipole_matrix(a2, b2, ctx2, p, geometry=GEO)
        assert np.allclose(default, via_axial, atol=1e-12)
    # non-vacuous: p = -1 is the allowed direction (m_a = 0.5 = 1.5 + (-1))
    assert np.max(np.abs(dipole_matrix(a2, b2, ctx2, -1))) > 0.0


def test_B7_sum_rule_still_holds_in_the_two_spin_basis():
    """Gate B7, ported onto the 229ThF+ two-spin basis (J_max=2, I_Th = 5/2):
    total strength out of a J=1 basis state, summed over polarisation and
    over EVERY final state in the full basis, is exactly 1.0 (the closed
    form) and independent of which m_F sublevel of a given (J, F1, F) it is.

    Run at J_max = 3 so that TWO rungs are closed, not one: J = 1 needs J' in
    {0, 1, 2} (J' = 0 does not exist in this Omega = 1 molecule, so effectively
    {1, 2}) and J = 2 needs {1, 2, 3} -- all present. Only the top rung J = 3 is
    open, exactly as v1 excluded its own top rung (J = 4 in a J_max=4 basis);
    test_J2_is_truncated_and_does_not_hit_1_0 below keeps the J_max = 2 basis
    where J = 2 is the open rung, so the truncation boundary is still pinned.
    Looping J = 2 matters because every Delta J = +-1 path out of it exists in
    both directions, which J = 1 (no J = 0 partner) cannot test.

    Uniquely catches a recoupled dipole that lost its outer-spin (I_F) or
    Th-spectator (I_Th) 6j: either one breaks this exact identity, not merely
    a phase (test_spectra.py's B7b already isolates that case).
    """
    spec = thf_spec("229", J_max=3)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))

    for J, F1, F in ((1, 1.5, 1.0), (1, 1.5, 2.0), (1, 2.5, 2.0),
                     (2, 2.5, 2.0), (2, 3.5, 4.0), (2, 0.5, 1.0)):
        vals = []
        for mF in np.arange(-F, F + 0.5, 1.0):
            sel = ((kets["J"] == J) & (kets["F1"] == F1) & (kets["F"] == F)
                  & (kets["mF"] == mF) & (kets["Om"] == 1.0))
            i = int(np.flatnonzero(sel)[0])
            tot = 0.0
            for p in (-1, 0, 1):
                row = dipole_matrix(kets[i:i + 1], kets, ctx, p, geometry=GEO)
                tot += float(np.sum(np.abs(row) ** 2))
            vals.append(tot)
        assert np.allclose(vals, vals[0], rtol=1e-10), f"J={J} F1={F1} F={F}: {vals}"
        assert vals[0] == pytest.approx(1.0, rel=1e-10), (
            f"J={J} F1={F1} F={F}: sum rule = {vals[0]}, expected the closed "
            f"form 1.0")


def test_J2_is_truncated_and_does_not_hit_1_0():
    """Companion to the sum rule above: J = 2 needs the truncated J' = 3, so
    its total is the SAME fraction (7/15-style) for every m_F but not 1.0 --
    documents the truncation boundary rather than silently excluding it."""
    spec = thf_spec("229", J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    F1, F = 4.5, 4.0
    vals = []
    for mF in np.arange(-F, F + 0.5, 1.0):
        sel = ((kets["J"] == 2) & (kets["F1"] == F1) & (kets["F"] == F)
              & (kets["mF"] == mF) & (kets["Om"] == 1.0))
        i = int(np.flatnonzero(sel)[0])
        tot = 0.0
        for p in (-1, 0, 1):
            row = dipole_matrix(kets[i:i + 1], kets, ctx, p, geometry=GEO)
            tot += float(np.sum(np.abs(row) ** 2))
        vals.append(tot)
    assert np.allclose(vals, vals[0], rtol=1e-10)
    assert vals[0] < 0.99


def test_B7_fails_if_a_polarisation_is_dropped():
    """FAIL demo for B7 (v1's precedent test_B7_fails_if_a_polarisation_is_
    dropped, ported): every element is the real `axial_geometry` output, but
    the m-sum over polarisation is deliberately incomplete (p = -1 dropped),
    which breaks the m_F-independence the closed form otherwise guarantees."""
    spec = thf_spec("229", J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    F1, F = 1.5, 1.0
    vals = []
    for mF in np.arange(-F, F + 0.5, 1.0):
        sel = ((kets["J"] == 1) & (kets["F1"] == F1) & (kets["F"] == F)
              & (kets["mF"] == mF) & (kets["Om"] == 1.0))
        i = int(np.flatnonzero(sel)[0])
        tot = 0.0
        for p in (0, 1):                       # p = -1 dropped on purpose
            row = dipole_matrix(kets[i:i + 1], kets, ctx, p, geometry=GEO)
            tot += float(np.sum(np.abs(row) ** 2))
        vals.append(tot)
    assert not np.allclose(vals, vals[0], rtol=1e-10), (
        f"dropping p=-1 should have broken m_F-independence, got {vals}")


def test_E1_selection_rules_in_the_two_spin_basis():
    """Delta F1 = 0, +-1, Delta F = 0, +-1, Delta m_F = p, read off the
    COMPUTED matrix over the whole 229ThF+ J <= 2 basis (not declared): every
    non-zero element obeys all three, and every one of the seven (Delta F1,
    Delta F) classes this basis structurally reaches has at least one
    non-zero element (measured: (-1,-1) (-1,0) (0,-1) (0,0) (0,1) (1,0) (1,1);
    (-1,+1) and (+1,-1) never co-occur here -- 720/360/420/840/420/360/720
    non-zero elements respectively, so it is not vacuous)."""
    spec = thf_spec("229", J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    classes = set()
    for p in (-1, 0, 1):
        D = dipole_matrix(kets, kets, ctx, p, geometry=GEO)
        for i, j in np.argwhere(np.abs(D) > 1e-12):
            dF1 = float(kets["F1"][i] - kets["F1"][j])
            dF = float(kets["F"][i] - kets["F"][j])
            dmF = float(kets["mF"][i] - kets["mF"][j])
            assert dmF == float(p), (i, j, dmF, p)
            assert dF1 in (-1.0, 0.0, 1.0), (i, j, dF1)
            assert dF in (-1.0, 0.0, 1.0), (i, j, dF)
            classes.add((dF1, dF))
    assert len(classes) >= 7, sorted(classes)


def test_label_lines_reports_F1_and_its_purity():
    """Each eigenvector column gets its dominant F1 and the |amplitude|^2
    ("purity") of that dominant component -- the same argmax `label_lines`
    already uses for J and F, just read one field further (controller
    ruling), so the notebook can *show* F1 is a good (not exact -- [TH] S4.5
    predicts a further second-order Delta F1 mixing shift <= 25 kHz)
    quantum number rather than asserting it."""
    spec = thf_spec("229", J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    blocks = block_by_mF(kets)
    sub = kets[blocks.index[0.0]]
    tm = build_term_matrices(sub, ctx, case="c2", registry=REGISTRY_C2)
    w, v = np.linalg.eigh(hamiltonian(tm, thf_v2("229"), {"E_z": 0.0, "B_z": 0.0}))
    labels = label_lines(sub, v, ctx.S, rule=ctx.conventions.ef_rule, ell=0.0, s=0.0)

    assert len(labels) == len(sub)
    assert all("F1" in l and "purity" in l for l in labels)
    purities = [l["purity"] for l in labels]
    # well above a uniformly-mixed state (1/len(sub)); < 1 because the same
    # +-Omega doubling that already caps a v1 J/F label's dominant-component
    # weight near 1/2 caps this one too (label_lines never groups amplitude
    # across kets, by design -- see its docstring).
    assert all(5.0 / len(sub) < p < 1.0 for p in purities), purities


def test_label_lines_v1_dict_is_unchanged_without_F1():
    """v1 KET_C basis (no F1 field) -> label_lines' dict gains neither key."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    blocks = block_by_mF(kets)
    sub = kets[blocks.index[0.5]]
    tm = build_term_matrices(sub, ctx)
    w, v = np.linalg.eigh(hamiltonian(tm, thf_v1(), {"E_z": 0.0, "B_z": 0.0}))
    labels = label_lines(sub, v, ctx.S, rule=ctx.conventions.ef_rule, ell=0.0, s=0.0)
    assert all(set(l) == {"J", "F", "parity", "ef"} for l in labels)


def test_line_strengths_forwards_geometry_to_the_two_spin_basis():
    """`line_strengths` end to end with `geometry=` on a real two-spin block:
    same shape/positivity contract as v1, non-vacuous, and the default is not
    merely different on this basis but INAPPLICABLE to it: `geometry=None`
    keeps v1's `dipole_geometry`, which reads J, Omega, F and m_F but never
    F1, so on the two-spin basis it feeds an INTEGER F where its 6j expects
    J +- 1/2 and sympy rejects the triad. That is computed below (the raise is
    asserted) rather than claimed in prose -- and it is why `geometry=` had to
    be a keyword rather than a silently-tolerant default."""
    spec = thf_spec("229", J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    blocks = block_by_mF(kets)
    sub = kets[blocks.index[0.0]]
    tm = build_term_matrices(sub, ctx, case="c2", registry=REGISTRY_C2)
    w, v = np.linalg.eigh(hamiltonian(tm, thf_v2("229"), {"E_z": 0.0, "B_z": 0.0}))

    freqs, strengths = line_strengths(w, v, sub, w, v, sub, ctx, geometry=GEO)
    assert freqs.shape == strengths.shape == (len(w), len(w))
    assert np.all(strengths >= 0.0)
    assert np.max(strengths) > 0.0

    with pytest.raises(ValueError, match="triangle"):
        line_strengths(w, v, sub, w, v, sub, ctx)
