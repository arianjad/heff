"""Gate V16 -- the master reduction. With I_Th = 0 every v2 matrix element must
equal its v1 counterpart to 1e-12 on the same block.

Uniquely catches: ANY wrong spectator phase or 6j column order in any of the
seven re-coupled terms. Nothing else in v2 sees these at once: the symmetry
gates (parity, hermiticity, tracelessness) all survive a transposed 6j, and the
Th-side gates never exercise the outer-spin recoupler at all.


"""
import numpy as np
import pytest

from heff.assemble import build_term_matrices
from heff.elements_c import dipole_geometry
from heff.elements_c2 import REGISTRY_C2, axial_geometry
from heff.params import thf_v1, thf_v2
from heff.spec import Spin, StateSpec, block_by_mF, enumerate_kets, thf_spec
from heff.terms import ctx_from, terms_for_case

# v2 name -> v1 name. The five delegating terms keep their v1 names; the four
# 19F terms are suffixed _F because Task 5 adds the Th partner under the bare
# name (spec-v2 S2.3).
V2_TO_V1 = {
    "rotation": "rotation",
    "centrifugal": "centrifugal",
    "omega_doubling": "omega_doubling",
    "hyperfine_A_par_F": "hyperfine_A_par",
    "hyperfine_A_par_F_dJ1": "hyperfine_A_par_dJ1",
    "spin_rotation_cI_F": "spin_rotation_cI",
    "stark_z": "stark_z",
    "zeeman_Gpar": "zeeman_Gpar",
    "zeeman_nuclear_F": "zeeman_nuclear",
    "pt_odd_edm": "pt_odd_edm",
    "pt_odd_scalar_pseudoscalar": "pt_odd_scalar_pseudoscalar",
}

# The v2 terms that have NO v1 counterpart, stated EXPLICITLY so that a term
# added to REGISTRY_C2 without a line in V2_TO_V1 fails this file instead of
# slipping through a subset check. Task 5's six Th terms are the members; every
# one of them is identically zero at I_Th = 0 (asserted in
# tests/test_elements_c2_th.py::test_every_th_term_vanishes_at_I_Th_zero),
# which is why V16 is unaffected by their arrival.
EXPECTED_NO_V1_COUNTERPART = {
    "hyperfine_A_par_Th", "hyperfine_A_par_Th_dJ1", "spin_rotation_cI_Th",
    "zeeman_nuclear_Th", "quadrupole_eQq0_Th", "quadrupole_eQq2_Th",
}


def spec_with_I_Th_zero(J_max=4):
    """The v2 two-spin spec at I_Th = 0 -- the master gate's subject.

    Built through StateSpec directly and NOT through thf_spec('232'), which
    deliberately returns the one-spin v1 object (controller ruling R9): the
    point of this gate is to run the KET_C2 dtype and the two-spin enumerator
    against the v1 answer, so the Th spin has to be present-and-zero, not absent.
    """
    v1 = thf_spec(J_max=J_max)
    return StateSpec(case="c", electronic=v1.electronic, I=0.5,
                     J_range=v1.J_range, v=0, M="blocks", frame="rotating",
                     spins=(Spin(label="232Th", I=0.0, couple_to="J"),
                            Spin(label="19F", I=0.5, couple_to="F1")))


@pytest.fixture
def v1_setup():
    spec = thf_spec()
    return enumerate_kets(spec), ctx_from(spec, thf_v1())


@pytest.fixture
def v2_setup():
    spec = spec_with_I_Th_zero()
    return enumerate_kets(spec), ctx_from(spec, thf_v2("232"))


def _dense(t, kets, ctx):
    """Every (i, j) of one term, with NO selection-rule mask and no hermiticity
    check -- so a corrupted 6j shows up as a number, not as an exception from
    the assembler."""
    d = len(kets)
    M = np.zeros((d, d), dtype=float)
    for i in range(d):
        for j in range(d):
            M[i, j] = float(np.real(t.fn(kets[i], kets[j], ctx)))
    return M


def test_the_two_bases_are_the_same_states_ket_for_ket(v1_setup, v2_setup):
    """Precondition of the whole gate: comparing matrices only means something
    if row i of each is the same physical state ([SPEC-v2] S2.1, gate V16a)."""
    k1, _ = v1_setup
    k2, _ = v2_setup
    assert len(k1) == len(k2)
    for f in ("J", "Om", "F", "mF"):
        assert np.array_equal(k1[f], k2[f]), f
    assert np.array_equal(k2["F1"], k2["J"]), "at I_Th = 0 the F1 column must be J"


def test_the_two_param_sets_share_one_convention_stamp():
    """A convention divergence between thf_v1() and thf_v2('232') would break
    the reduction for a reason that has nothing to do with the recoupling, so
    it is separated out and reported on its own.

    `version` is excluded on purpose: it labels the parameter set ('thf-v1' vs
    'thf-v2'), not a physics convention, and no element reads it."""
    a = dict(thf_v1().conventions.stamp())
    b = dict(thf_v2("232").conventions.stamp())
    a.pop("version", None)
    b.pop("version", None)
    assert a == b


def test_V16_every_v2_term_reduces_to_its_v1_counterpart(v1_setup, v2_setup):
    """The master gate, term by term, on every signed-m_F block.

    It is also the gate that discriminates the B&C (5.173) ket-F1 vs bra-F1
    phase, which V19 cannot see: at I_Th = 0 the F1 column is J, so the two
    readings differ exactly on Delta J = +-1, and hyperfine_A_par_F_dJ1 is the
    only term here that has such elements -- the bra-F1 variant deviates from
    v1 hyperfine_A_par_dJ1 by 0.9682 in units of A_par^F (measured; [HAM] S9.3
    settles the phase by a decoupled-basis rebuild).
    """
    k1, c1 = v1_setup
    k2, c2 = v2_setup
    v2_names = {t.name for t in terms_for_case("c2", registry=REGISTRY_C2)}
    v1_names = {t.name for t in terms_for_case("c")}
    assert set(V2_TO_V1) <= v2_names, sorted(set(V2_TO_V1) - v2_names)
    assert v2_names - set(V2_TO_V1) == EXPECTED_NO_V1_COUNTERPART, (
        "a v2 term appeared with no v1 counterpart and no entry in "
        "EXPECTED_NO_V1_COUNTERPART; add it to one or the other so the master "
        f"gate's coverage stays explicit: {sorted(v2_names - set(V2_TO_V1))}")
    assert set(V2_TO_V1.values()) == v1_names, (
        "every v1 case-(c) term needs a v2 counterpart in this gate; missing "
        f"{sorted(v1_names - set(V2_TO_V1.values()))}")
    b1, b2 = block_by_mF(k1), block_by_mF(k2)
    assert b1.labels == b2.labels
    for label in b1.labels:
        sub1, sub2 = k1[b1.index[label]], k2[b2.index[label]]
        tm1 = build_term_matrices(sub1, c1, case="c")
        tm2 = build_term_matrices(sub2, c2, case="c2", registry=REGISTRY_C2)
        m1 = dict(zip(tm1.names, tm1.mats))
        m2 = dict(zip(tm2.names, tm2.mats))
        for v2_name, v1_name in V2_TO_V1.items():
            dev = np.max(np.abs(np.asarray(m2[v2_name]) - np.asarray(m1[v1_name])))
            assert dev < 1e-12, (
                f"m_F = {label}: v2 term {v2_name!r} deviates from v1 "
                f"{v1_name!r} by {dev:.3e}")


def test_axial_geometry_reduces_to_dipole_geometry_at_one_spin(v2_setup):
    """[HAM] S9.1's analytic collapse, element by element, independent of the
    term layer -- so a kernel bug is separated from a registration bug."""
    kets, ctx = v2_setup
    worst = 0.0
    for p in (-1, 0, 1):
        for i in range(len(kets)):
            for j in range(len(kets)):
                a = axial_geometry(kets[i], kets[j], ctx, k=1, q=0, p=p)
                b = dipole_geometry(kets[i], kets[j], ctx.I, p)
                worst = max(worst, abs(a - b))
    assert worst < 1e-12, f"max |axial_geometry - dipole_geometry| = {worst:.3e}"


def test_axial_geometry_refuses_q_greater_than_k(v2_setup):
    """Structural: a rank-k operator has no |q| > k component, so a caller that
    asks for one has a bug, not a zero."""
    kets, ctx = v2_setup
    with pytest.raises(ValueError):
        axial_geometry(kets[0], kets[0], ctx, k=1, q=2, p=0)
    with pytest.raises(ValueError):
        axial_geometry(kets[0], kets[0], ctx, k=0, q=-1, p=0)


def _19F_doublet_splittings(spec, J):
    """E(F = F1+1/2) - E(F = F1-1/2) of hyperfine_A_par_F, per F1, in units of
    A_par^F. Keyed by F1, at the given J."""
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229" if spec.spins[0].I == 2.5 else "232"))
    fn = REGISTRY_C2["hyperfine_A_par_F"].fn
    out = {}
    sel = (kets["J"] == J) & (kets["Om"] == 1.0)
    for F1 in sorted({float(v) for v in kets["F1"][sel]}):
        e = {}
        for F in (F1 - 0.5, F1 + 0.5):
            hit = np.flatnonzero(sel & (kets["F1"] == F1) & (kets["F"] == F)
                                 & (kets["mF"] == (0.0 if F % 1 == 0 else 0.5)))
            if len(hit) == 1:
                e[F] = float(fn(kets[hit[0]], kets[hit[0]], ctx))
        if len(e) == 2:
            out[F1] = e[F1 + 0.5] - e[F1 - 0.5]
    return out


def test_V19_the_19F_doublet_ordering_flips_between_F1_manifolds():
    """[TH] S4.5 / [HAM] S9.3 Table 2: at I_Th = 5/2, J = 1 the 19F doublet
    splittings over F1 = 3/2, 5/2, 7/2 are -0.4000, +0.1714, +0.5714 in units of
    A_par^F, against +0.7500 at I_Th = 0.

    THIS IS A SIGN TEST. The F1 = 3/2 entry INVERTS relative to the one-spin
    answer because [HAM] S9.3's closed-form bracket J(J+1)+F1(F1+1)-I_Th(I_Th+1)
    = 2 + 3.75 - 8.75 = -3 there, i.e. J is anti-aligned with F1.

    What it catches, per [HAM] S9.3's falsification block: a TRANSPOSED 6j
    COLUMN ORDER in either recoupler. Transposing the I_F 6j sends the three
    numbers to -0.1581, 0, 0; transposing the I_Th 6j sends them to 0, 0, 0 --
    and the I_Th one is the corruption the S9.1 I_Th = 0 collapse check is
    blind to by construction.

    What it does NOT catch: the B&C (5.173) ket-F1 vs bra-F1 phase. The same
    falsification block prints -0.4000, +0.1714, +0.5714 for the bra-F1 variant
    -- Table 2 reproduced identically, because every entry here is Delta F1 = 0
    and the two readings differ only on Delta F1 = +-1. That phase is pinned by
    [HAM] S9.3's decoupled-basis rebuild and, inside this file, by V16 through
    hyperfine_A_par_F_dJ1.
    """
    got = _19F_doublet_splittings(thf_spec("229", J_max=1), J=1)
    assert sorted(got) == [1.5, 2.5, 3.5]
    # [TH] S4.5 quotes four decimals, so the comparison against it is absolute
    # to half the last printed digit; [HAM] S9.3 Table 2 prints five, and that
    # one is compared at the brief's rel = 1e-4. Both are asserted -- 6/35 and
    # 4/7 round to 0.1714 and 0.5714, which rel = 1e-4 alone would reject.
    for F1, th, ham in ((1.5, -0.4000, -0.40000),
                        (2.5, +0.1714, +0.17143),
                        (3.5, +0.5714, +0.57143)):
        assert got[F1] == pytest.approx(th, abs=5e-5), f"[TH] S4.5 F1={F1}"
        assert got[F1] == pytest.approx(ham, rel=1e-4), f"[HAM] S9.3 F1={F1}"
    ref = _19F_doublet_splittings(spec_with_I_Th_zero(J_max=1), J=1)
    assert ref[1.0] == pytest.approx(+0.7500, rel=1e-4)
    assert got[1.5] < 0.0 < ref[1.0], "the F1 = 3/2 doublet ordering must invert"
