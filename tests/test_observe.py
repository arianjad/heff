"""Gate A3 and the ThF+ g-factor closed forms ([HAM] V5, V6, V7).

A3 is Hamiltonian-agnostic (random Hermitian matrices). V5/V6 are exact closed
forms in the PARAMETERS. V7 is a labelled, opt-in comparison to a published
model number.

No measured number appears in a default-on test: every gate that runs by
default checks a closed form or a symmetry, and the comparisons to Ng's
measured |g| and to his 32-level delta-g sit behind HEFF_RUN_LITERATURE=1,
each naming the convention block it assumes.
"""
import os
from dataclasses import replace

import numpy as np
import pytest
from _helpers import _gamma, _kappa

from heff import elements_c  # noqa: F401
from heff.assemble import build_term_matrices, hamiltonian
from heff.conventions import Conventions
from heff.observe import (expectation, g_factors, multi_curvature, offdiag,
                          pair_differential)
from heff.params import MU_B, MU_N, thf_v1
from heff.spec import ElecState, StateSpec, block_by_mF, enumerate_kets
from heff.terms import ctx_from

I_F = 0.5   # I(19F). _gamma/_kappa take it keyword-only with no default.


def _block(mF, J_range=(1, 4), pset=None):
    """One signed-m_F block of the v1 ThF+ spec, with the Ctx it was built from."""
    pset = thf_v1() if pset is None else pset
    spec = StateSpec(case="c",
                     electronic=(ElecState("X3Delta1", 1.0, 1.0, 2.0),),
                     I=I_F, J_range=J_range)
    kets = enumerate_kets(spec)
    idx = block_by_mF(kets).index[mF]
    ctx = ctx_from(spec, pset)
    return build_term_matrices(kets[idx], ctx), ctx


def _random_hermitian(rng, d):
    """A real symmetric (n, n) -- the Hamiltonian-agnostic fixture gate A3 needs."""
    A = rng.standard_normal((d, d))
    return A + A.T


def test_A3_first_derivatives_match_central_finite_difference():
    """Gate A3. Uniquely catches derivative-kernel sign and self-term errors --
    the self-term is excluded BY INDEX in the lifted kernel, and a version that
    excludes it by an energy-gap threshold gives subtly wrong answers near a
    small gap. FAIL is reachable: negate the vertex passed to multi_curvature
    and every first derivative flips sign against the finite difference."""
    rng = np.random.default_rng(2026)
    d = 6
    H0 = _random_hermitian(rng, d)
    verts = {"a": _random_hermitian(rng, d), "b": _random_hermitian(rng, d)}
    mc = multi_curvature(H0, verts)
    h = 1e-5
    for name, V in verts.items():
        up = np.linalg.eigvalsh(H0 + h * V)
        dn = np.linalg.eigvalsh(H0 - h * V)
        assert np.allclose(mc["d"][name], (up - dn) / (2 * h), atol=1e-6)


def test_A3_second_derivatives_match_central_finite_difference():
    """Gate A3. Uniquely catches a factor-of-two or a dropped conjugate in the
    sum-over-states kernel: d2 is the SECOND DERIVATIVE d2E/dx^2 itself, not the
    perturbation coefficient E2 = d2E/dx^2 / 2, so a kernel that returned the
    coefficient would land a clean factor 2 off the finite difference. FAIL is
    reachable by exactly that factor."""
    rng = np.random.default_rng(11)
    d = 5
    H0 = _random_hermitian(rng, d)
    V = _random_hermitian(rng, d)
    mc = multi_curvature(H0, {"a": V})
    h = 1e-3
    w0 = np.linalg.eigvalsh(H0)
    up = np.linalg.eigvalsh(H0 + h * V)
    dn = np.linalg.eigvalsh(H0 - h * V)
    fd = (up - 2 * w0 + dn) / h ** 2
    # d2 is the full second derivative, NOT the second-order perturbation
    # coefficient -- compare it to the finite difference with no factor of 2.
    assert np.allclose(mc["d2"][("a", "a")], fd, atol=1e-4)


def test_A3_role_swap_is_symmetric():
    """Gate A3. Uniquely catches a bra/ket asymmetry in the mixed second
    derivative: d2E/dx_a dx_b must not depend on which vertex was named first.
    Keys follow verts ORDER, so the swapped call is read at ("b", "a"). FAIL is
    reachable if curv() dropped the conjugate on one argument."""
    rng = np.random.default_rng(3)
    d = 5
    H0 = _random_hermitian(rng, d)
    A, B = _random_hermitian(rng, d), _random_hermitian(rng, d)
    ab = multi_curvature(H0, {"a": A, "b": B})["d2"][("a", "b")]
    ba = multi_curvature(H0, {"b": B, "a": A})["d2"][("b", "a")]
    assert np.allclose(ab, ba)


def test_A3_a_vertex_proportional_to_the_identity_has_zero_curvature():
    """Gate A3. Common-mode vertex: it shifts every eigenvalue equally, so its
    first derivative is exactly 1 and its second derivative exactly 0. Uniquely
    catches a self-term that leaked back into the sum-over-states denominator --
    that would divide by a zero gap and return inf/NaN here. FAIL is reachable
    by removing the dW[i, i] = inf line."""
    rng = np.random.default_rng(4)
    d = 5
    H0 = _random_hermitian(rng, d)
    mc = multi_curvature(H0, {"c": np.eye(d)})
    assert np.allclose(mc["d"]["c"], 1.0)
    assert np.allclose(mc["d2"][("c", "c")], 0.0)


def test_expectation_and_offdiag_match_the_explicit_quadratic_forms():
    """expectation is <k|O|k> and offdiag is the named-pair <i|O|j>."""
    rng = np.random.default_rng(6)
    d = 4
    H = _random_hermitian(rng, d)
    _, v = np.linalg.eigh(H)
    op = _random_hermitian(rng, d)
    ex = expectation(v[None, ...], op)
    for k in range(d):
        assert ex[0, k] == pytest.approx(v[:, k] @ op @ v[:, k])
    od = offdiag(v[None, ...], op, 0, 1)
    assert od[0] == pytest.approx(v[:, 0] @ op @ v[:, 1])


def test_offdiag_is_the_only_route_to_an_imaginary_hermitian_operator():
    """The PTV_shift trap (spec S3.6; Molecule-Structure Energy_Levels.py:481-483).

    For O = i x (real antisymmetric), O is Hermitian but purely imaginary, so
    <psi|O|psi> vanishes identically on the REAL eigenvectors of a real
    symmetric H -- the diagonal mode is structurally blind to such an operator
    and the named-pair off-diagonal element is the observable. Uniquely catches
    an implementation that answered a PT-odd question with expectation().
    Both outcomes are reachable and both are asserted here: the diagonal is
    exactly zero, the off-diagonal is not.
    """
    rng = np.random.default_rng(17)
    d = 4
    H = _random_hermitian(rng, d)
    _, v = np.linalg.eigh(H)
    A = rng.standard_normal((d, d))
    op = 1j * (A - A.T)
    assert np.allclose(op, op.conj().T)                       # Hermitian
    assert np.allclose(np.real(op), 0.0)                      # and purely imaginary
    assert np.allclose(expectation(v[None, ...], op), 0.0)
    assert abs(offdiag(v[None, ...], op, 0, 1)[0]) > 1e-3


def _g_by_JF(tm, ctx, pset, knobs):
    """g grouped by the (J, F) of each eigenstate's dominant basis ket."""
    res = g_factors(tm, pset, knobs, ctx=ctx)
    out = {}
    for k in range(len(res["g"])):
        dom = tm.kets[int(np.argmax(np.abs(res["V"][:, k])))]
        out.setdefault((float(dom["J"]), float(dom["F"])), []).append(float(res["g"][k]))
    return out


def test_V5_g_factor_closed_form_is_exact_in_a_single_J_basis():
    """[HAM] V5: g(J,F) = -G_par gamma_F + g_N (mu_N/mu_B) kappa_F.

    Exact only when the basis holds one J -- the B&C 9.51 hyperfine term mixes
    J by ~3e-4, and the Delta-J = +-1 Zeeman cross term that mixing opens up
    moves g at the 1e-5 level (see the next test). Ng's model has no Delta-J
    hyperfine, which is why his closed form is exact.

    Uniquely catches the Ng Eq. C.6 sign error of [HAM] S2.8: with Ng's printed
    minus the model gives 23.5 kHz/G instead of 20.85 kHz/G -- a 13 % error no
    other check sees. FAIL is reachable, and is
    test_V5_the_wrong_zeeman_sign_gives_the_wrong_g below.
    """
    pset = thf_v1()
    G, gN = pset.value("G_par"), pset.value("g_N")
    checked = 0
    for mF, J in ((1.5, 1), (0.5, 1), (2.5, 2), (1.5, 2)):
        tm, ctx = _block(mF, J_range=(J, J))
        for (Jd, Fd), gs in _g_by_JF(tm, ctx, pset, {"E_z": 0.0, "B_z": 0.0}).items():
            want = (-G * _gamma(Jd, Fd, I=I_F)
                    + gN * (MU_N / MU_B) * _kappa(Jd, Fd, I=I_F))
            for g in gs:
                assert g == pytest.approx(want, abs=1e-12), f"J={Jd} F={Fd}"
            checked += 1
    assert checked >= 5
    tm, ctx = _block(1.5, J_range=(1, 1))
    g = g_factors(tm, pset, {"E_z": 0.0, "B_z": 0.0}, ctx=ctx)["g"][0]
    assert g * MU_B == pytest.approx(-0.020853, abs=1e-6)   # -20.85 kHz/G, [HAM] S2.8


def test_V5_holds_in_the_full_basis_to_the_dJ1_hyperfine_level():
    """Same closed form in the J = 1-4 basis, where it is approximate.

    The residual is the cross term 2c <J|dH/dB|J+1> opened up by the B&C 9.51
    J-mixing amplitude c ~ 3e-4 -- of order 1e-5 in g. A tolerance tighter than
    that would be tuned; a tolerance looser than 1e-3 would stop catching the
    sign error, whose residual is 0.03 in g.
    """
    pset = thf_v1()
    G, gN = pset.value("G_par"), pset.value("g_N")
    tm, ctx = _block(1.5)
    for (Jd, Fd), gs in _g_by_JF(tm, ctx, pset, {"E_z": 0.0, "B_z": 0.0}).items():
        want = (-G * _gamma(Jd, Fd, I=I_F)
                + gN * (MU_N / MU_B) * _kappa(Jd, Fd, I=I_F))
        for g in gs:
            assert g == pytest.approx(want, abs=1e-4), f"J={Jd} F={Fd}"


def test_V5_the_wrong_zeeman_sign_gives_the_wrong_g():
    """The fabricated failure that makes V5 non-vacuous: with Ng's printed
    minus_Gpar operator the same block returns 23.52 kHz/G, not 20.85."""
    bad = replace(thf_v1(), conventions=Conventions(zeeman_sign="minus_Gpar"))
    tm, ctx = _block(1.5, pset=bad)
    g = g_factors(tm, bad, {"E_z": 0.0, "B_z": 0.0}, ctx=ctx)["g"][0]
    assert abs(g * MU_B) == pytest.approx(0.02352, abs=1e-4)   # 23.52 kHz/G, wrong


def test_V6_zeeman_is_even_in_omega_in_a_single_J_basis():
    """[HAM] V6: at E = 0 with only J = 1 in the basis (so no Delta-J = +-1
    Stark or Zeeman path), the two Omega-doublet components have identical
    g-factors -- exactly. Uniquely catches an odd-in-Omega contamination of the
    Zeeman operator, which would fake an eEDM signal in the four-way chop.
    FAIL is reachable: drop the ket['Om'] factor from zeeman_Gpar and the two
    g's differ.
    """
    tm, ctx = _block(1.5, J_range=(1, 1))
    res = g_factors(tm, thf_v1(), {"E_z": 0.0, "B_z": 0.0}, ctx=ctx)
    assert len(res["g"]) == 2
    assert res["g"][0] == pytest.approx(res["g"][1], abs=1e-14)
    dg, label = pair_differential(res["g"][1], res["g"][0], thf_v1().conventions)
    assert dg == pytest.approx(0.0, abs=1e-14)
    assert "g^u - g^l" in label


def test_pair_differential_honours_the_dg_def_convention():
    """Delta = g^u - g^l by default, delta = half that under dg_def='delta'."""
    conv = Conventions()
    assert pair_differential(0.3, 0.1, conv)[0] == pytest.approx(0.2)
    half = replace(conv, dg_def="delta")
    val, label = pair_differential(0.3, 0.1, half)
    assert val == pytest.approx(0.1) and "delta" in label


def test_pair_differential_keeps_the_sign_of_a_named_pair():
    """Swapping the partners flips the sign -- there is no fold to |.| here."""
    conv = Conventions()
    assert pair_differential(0.1, 0.3, conv)[0] == pytest.approx(-0.2)
    arr = pair_differential(np.array([0.3, -0.1]), np.array([0.1, 0.1]), conv)[0]
    assert np.allclose(arr, [0.2, -0.2])


def test_fully_polarised_stark_slope_is_gamma_F_m_F_d_mf():
    """[HAM] S2.7 / Leanhardt Eqs. 20-21: in the polarised limit the induced
    dipole saturates at gamma_F m_F d_mf, with gamma_(F=3/2)(J=1) = 1/3.

    Run in a J = 1-only basis so Stark J-mixing (which is 23 % at this field in
    the full basis) cannot contaminate the limit; the only correction left is
    the omega_ef one, ~5e-6 relative.
    """
    pset = thf_v1()
    tm, ctx = _block(1.5, J_range=(1, 1))
    res = g_factors(tm, pset, {"E_z": 1000.0, "B_z": 0.0}, ctx=ctx)
    want = 1.5 * (1.0 / 3.0) * pset.value("d_mf")
    assert _gamma(1, 1.5, I=I_F) == pytest.approx(1 / 3)
    assert sorted(abs(res["d_eff"]))[-1] == pytest.approx(want, rel=1e-4)


def _delta_g_over_g(pset, E_field, tm_ctx=None):
    """delta_g/g for the J = 1, F = 3/2, |m_F| = 3/2 Stark doublet.

    The pair is named by ENERGY ORDER, which is what "upper/lower Stark
    doublet" means, and the difference is taken signed through
    pair_differential -- no |.| fold anywhere.
    """
    tm, ctx = _block(1.5) if tm_ctx is None else tm_ctx
    res = g_factors(tm, pset, {"E_z": float(E_field), "B_z": 0.0}, ctx=ctx)
    order = np.argsort(res["W"])[:2]        # the J = 1 pair; J = 2 is 29 GHz away
    lower, upper = res["g"][int(order[0])], res["g"][int(order[1])]
    dg, _ = pair_differential(upper, lower,
                              replace(pset.conventions, dg_def="delta"))
    return float(dg / (0.5 * (upper + lower)))


def test_delta_g_scales_like_the_leanhardt_expression_within_a_factor_of_two():
    """Leanhardt Eq. 67: |delta_g/g| = 9 d_mf E / (40 B_e) = 0.00315 at 60 V/cm.

    LOOSE on purpose. Eq. 67 is a leading-order fully-polarised expression and
    [HAM] S2.9 records that the 96-state model gives -0.00223 at the same field
    -- a ~30 % gap that is physics, not a bug, so a tight tolerance here would
    be a tuned gate. What this check uniquely catches is a Stark Delta-J = +-1
    element that is wrong by ORDERS of magnitude or missing entirely: delta_g
    then collapses to ~0, which is the fabricated failure in the next test.
    """
    pset = thf_v1()
    got = _delta_g_over_g(pset, 60.0)
    want = 9 * pset.value("d_mf") * 60.0 / (40 * pset.value("B0"))
    assert got < 0.0, "delta_g/g is negative in the JILA data and in every model"
    assert 0.5 < abs(got) / want < 2.0


def test_delta_g_collapses_when_the_stark_J_mixing_is_removed():
    """The fabricated failure that makes the previous check non-vacuous."""
    pset = thf_v1()
    tm, ctx = _block(1.5)
    k = tm.names.index("stark_z")
    mats = np.array(tm.mats, copy=True)
    J = np.asarray(tm.kets["J"], dtype=float)
    mats[k][J[:, None] != J[None, :]] = 0.0     # kill Delta-J = +-1 Stark
    crippled = replace(tm, mats=tuple(mats))
    assert abs(_delta_g_over_g(pset, 60.0, tm_ctx=(crippled, ctx))) < 1e-4


@pytest.mark.skipif(os.environ.get("HEFF_RUN_LITERATURE") != "1",
                    reason="tier D: set HEFF_RUN_LITERATURE=1 to compare with published numbers")
def test_V7_delta_g_over_g_against_ng_thesis_opt_in():
    """TIER D, opt-in. Assumes the v1 conventions block (n_hat = F_to_Th,
    zeeman_sign = plus_Gpar, dg_def = delta for this quantity).

    At E = 60 V/cm the Omega = +-1, J = 1-4 model gives delta_g/g = -0.00223,
    identical to Ng thesis p.85's 32-level value ([HAM] S2.9, V7). The
    measurement is -0.00255(6); the 15 % gap is the KNOWN missing 3Delta2
    coupling and is a finding, not a failure. FAIL is reachable -- the same
    crippled-Stark block that collapses the loose gate collapses this one.
    """
    got = _delta_g_over_g(thf_v1(), 60.0)
    assert got == pytest.approx(-0.00223, rel=0.05)


@pytest.mark.skipif(os.environ.get("HEFF_RUN_LITERATURE") != "1",
                    reason="tier D: set HEFF_RUN_LITERATURE=1 to compare with published numbers")
def test_V7b_measured_g_and_omega_ef_are_reproduced_opt_in():
    """TIER D, opt-in: Ng 2022 Table I, under the v1 conventions block
    (zeeman_sign = plus_Gpar, n_hat = F_to_Th). |g_(F=3/2)| = 0.0149(3) and the
    J = 1 Omega-doublet splitting 5.29(5) MHz -- both are inputs to the fit, so
    this checks the round trip, not the physics. FAIL is reachable: the
    minus_Gpar block gives |g| = 0.0168.
    """
    pset = thf_v1()
    tm, ctx = _block(1.5)
    g = g_factors(tm, pset, {"E_z": 0.0, "B_z": 0.0}, ctx=ctx)["g"][0]
    assert abs(g) == pytest.approx(0.0149, abs=0.0003)
    w = np.linalg.eigvalsh(hamiltonian(tm, pset, {"E_z": 0.0, "B_z": 0.0}))
    assert (w[1] - w[0]) == pytest.approx(5.29, abs=0.05)
