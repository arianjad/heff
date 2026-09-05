"""Field and PT-odd case (c) elements against the closed forms of [HAM] S6
(V2, V5, V14) and the structural parity/Omega identities that V6 and V8 rest on.

Helpers live in tests/_helpers.py, so gamma_F and kappa_F are written once for
the whole suite. Nothing here compares against a stored spectrum.
"""
import numpy as np
import pytest

from _helpers import _elem, _find, _gamma, _kappa
from heff import elements_c
from heff.conventions import Conventions
from heff.params import MU_B, MU_N, Param, thf_v1
from heff.spec import enumerate_kets, thf_spec
from heff.terms import Ctx, ctx_from


@pytest.fixture
def basis():
    return enumerate_kets(thf_spec())


@pytest.fixture
def ctx():
    return ctx_from(thf_spec(), thf_v1())


def _reconvention(ctx, **kw):
    """The same Ctx with a different conventions block."""
    return Ctx(S=ctx.S, Lam=ctx.Lam, I=ctx.I, mu_B=ctx.mu_B, mu_N=ctx.mu_N,
               conventions=Conventions(**kw))


# ---------------------------------------------------------------------- Stark

def test_V2_stark_closed_form_for_every_diagonal_element(basis, ctx):
    """[HAM] V2: the Delta-J = 0, Delta-F = 0, p = 0 Stark element equals
    -Omega m_F gamma_F d_mf E, with
    gamma_F = [J(J+1)+F(F+1)-I(I+1)] / [2F(F+1)J(J+1)]  (Leanhardt Eqs. 20-21).

    Uniquely catches a wrong (-1)^(J'-Omega') phase, which leaves every
    MAGNITUDE right and corrupts only the (J, F) dependence, so no single-level
    check sees it. FAIL is reachable and was measured: dropping that phase
    leaves J = 1 and J = 3 exactly right and flips the sign at J = 2 and J = 4
    (Omega = +1 gives (-1)^(J-1)).

    [HAM] V2 also names "a wrong 6j column order"; that half is NOT reachable
    here and the wording is loose -- a column permutation is an exact symmetry
    of a 6j symbol, so it changes no number anywhere. A genuine 6j argument
    MISassignment collapses to the true symbol on the diagonal and needs the
    off-diagonal checks (hermiticity here, V7 later) to catch it.
    """
    d = thf_v1().value("d_mf")
    E = 1.0
    for i in range(len(basis)):
        J, Om, F, mF = (float(basis[n][i]) for n in ("J", "Om", "F", "mF"))
        got = d * E * _elem("stark_z", basis, ctx, i, i)
        assert got == pytest.approx(-Om * mF * _gamma(J, F, I=0.5) * d * E, abs=1e-12)


def test_stark_stretched_shift_at_the_jila_field(basis, ctx):
    """Cross-check against [HAM] S2.7, which checks it against Ng Table B.2's
    '2 pi x 50.9 MHz': J=1, F=3/2, m_F=3/2, Omega=+1 at 60 V/cm."""
    d = thf_v1().value("d_mf")
    i = _find(basis, 1, 1.0, 1.5, 1.5)
    assert d * 60.0 * _elem("stark_z", basis, ctx, i, i) == pytest.approx(-50.895, abs=1e-3)


def test_V14_stark_is_traceless_over_mF(basis, ctx):
    """[HAM] V14: sum over m_F of the diagonal Stark element vanishes at fixed
    (J, F, Omega).

    Uniquely catches a missing (-1)^(F-m_F) in the Wigner-Eckart 3j: that phase
    alternates with m_F, so dropping it turns the alternating sum into a
    non-zero one while every |element| stays right. FAIL is reachable and was
    measured: with the phase removed the J = 1, F = 3/2, Omega = +1 row sums to
    -2/3 and the J = 2, F = 5/2 row to -2/5, instead of 0.
    """
    for J in (1, 2, 3, 4):
        for F in (J - 0.5, J + 0.5):
            for Om in (-1.0, 1.0):
                idx = np.flatnonzero((basis["J"] == J) & (basis["F"] == F)
                                     & (basis["Om"] == Om))
                total = sum(_elem("stark_z", basis, ctx, i, i) for i in idx)
                assert total == pytest.approx(0.0, abs=1e-12)


def test_dipole_geometry_supports_p_plus_minus_one_for_spectra(basis, ctx):
    """heff.spectra (Task 10) needs p = +-1; the Stark term uses p = 0 only."""
    i = _find(basis, 1, 1.0, 1.5, 0.5)
    j = _find(basis, 1, 1.0, 1.5, 1.5)
    assert elements_c.dipole_geometry(basis[i], basis[j], ctx.I, 0) == 0.0
    assert abs(elements_c.dipole_geometry(basis[i], basis[j], ctx.I, -1)) > 0.0


# --------------------------------------------------------------------- Zeeman

def test_V5_zeeman_Gpar_diagonal_gives_minus_Gpar_gamma_F(basis, ctx):
    """[HAM] V5, G_par half. With E = -g mu_B B m_F the G_par term alone gives
    g = -G_par gamma_F, and g(1, 3/2) mu_B = -20.853 kHz/G once g_N is added --
    the number [HAM] S2.8 confirms numerically at G_par = 0.04756.

    Uniquely catches the Ng Eq. C.6 printed sign ([HAM] S2.8, OPEN-3): the
    minus sign gives g = +G_par gamma_F, i.e. 23.5 kHz/G rather than the
    measured-consistent 20.85 kHz/G. Both outcomes are reachable from the
    conventions block -- the 'minus_Gpar' branch is exercised below.
    """
    G, B = thf_v1().value("G_par"), 1.0
    for J in (1, 2, 3, 4):
        for F in (J - 0.5, J + 0.5):
            i = _find(basis, J, 1.0, F, F)
            shift = G * B * _elem("zeeman_Gpar", basis, ctx, i, i)
            g = -shift / (MU_B * B * F)
            assert g == pytest.approx(-G * _gamma(J, F, I=0.5), abs=1e-12)
    flipped = _reconvention(ctx, zeeman_sign="minus_Gpar")
    i = _find(basis, 1, 1.0, 1.5, 1.5)
    assert _elem("zeeman_Gpar", basis, flipped, i, i) == pytest.approx(
        -_elem("zeeman_Gpar", basis, ctx, i, i), abs=1e-12)


def test_V5_total_g_factor_closed_form(basis, ctx):
    """[HAM] V5 in full: both Zeeman terms together give
    g_F = -G_par gamma_F + g_N (mu_N/mu_B) kappa_F.

    Uniquely catches a RELATIVE sign or scale error between the two Zeeman
    terms, which the two single-term checks around it cannot see: they are the
    same 10 % of g_F that [HAM] S2.8 calls "the cleanest experimental handle on
    the nuclear contribution" (the ratio g(1,1/2)/g(1,3/2) is 2.19, not the
    G_par-only 2 of Leanhardt Eq. 24). FAIL is reachable and was measured:
    flipping the relative sign of the two terms moves the ratio to 1.830, and
    dropping the nuclear term moves it to exactly 2.

    Gated against the closed form at the parameter set's own G_par (0.04756),
    NOT against [HAM] S2.8's printed g_F table -- that table is a known
    erratum (it reproduces only at G_par = 0.048, not its own printed header
    value); see docs/open-questions.md "Erratum -- [HAM] S2.8 g_F table".
    """
    ps = thf_v1()
    G, gN, B = ps.value("G_par"), ps.value("g_N"), 1.0

    def g_of(J, F):
        i = _find(basis, J, 1.0, F, F)
        shift = B * (G * _elem("zeeman_Gpar", basis, ctx, i, i)
                     + gN * _elem("zeeman_nuclear", basis, ctx, i, i))
        return -shift / (MU_B * B * F)

    for J, F in ((1, 1.5), (1, 0.5), (2, 2.5), (4, 4.5)):
        closed = -G * _gamma(J, F, I=0.5) + gN * (MU_N / MU_B) * _kappa(J, F, I=0.5)
        assert g_of(J, F) == pytest.approx(closed, abs=1e-12)
    assert g_of(1, 1.5) * MU_B * 1e3 == pytest.approx(-20.853, abs=1e-3)
    assert g_of(1, 0.5) / g_of(1, 1.5) == pytest.approx(2.19, abs=5e-3)


def test_zeeman_nuclear_diagonal_gives_plus_gN_kappa_F(basis, ctx):
    """[HAM] S2.8: the nuclear term contributes g = +g_N (mu_N/mu_B) kappa_F,
    kappa_F = [F(F+1)-J(J+1)+I(I+1)] / [2F(F+1)]."""
    gN, B = thf_v1().value("g_N"), 1.0
    for J in (1, 2, 3, 4):
        for F in (J - 0.5, J + 0.5):
            i = _find(basis, J, 1.0, F, F)
            shift = gN * B * _elem("zeeman_nuclear", basis, ctx, i, i)
            g = -shift / (MU_B * B * F)
            assert g == pytest.approx(gN * (MU_N / MU_B) * _kappa(J, F, I=0.5), abs=1e-12)


def test_nuclear_zeeman_decouples_to_g_N_mu_N_per_gauss(basis, ctx):
    """[HAM] S2.8: g_N mu_N = 4.008 kHz/G. At fixed (J, Omega, m_F) the two
    F = J +- 1/2 kets span the m_I = +-1/2 space, so the eigenvalues of the
    nuclear-Zeeman block are exactly -/+ (1/2) g_N mu_N B at every J.

    The only check in this module that exercises the Delta-F = +-1 spectator
    element of B&C Eq. (5.175); the kappa_F check above is diagonal and blind
    to it. FAIL is reachable and was measured: dropping the off-diagonal element
    leaves the J = 2, m_F = 1/2 eigenvalues at +-0.1 g_N mu_N B and the J = 4,
    m_F = 1/2 pair at +-0.0556, instead of +-0.5 at every J.
    """
    gN, B = thf_v1().value("g_N"), 1.0
    for J in (1, 2, 3, 4):
        for mF in {0.5, J - 0.5}:
            idx = [_find(basis, J, 1.0, J - 0.5, mF), _find(basis, J, 1.0, J + 0.5, mF)]
            M = gN * B * np.array([[_elem("zeeman_nuclear", basis, ctx, i, j)
                                    for j in idx] for i in idx])
            assert np.linalg.eigvalsh(M) == pytest.approx(
                [-0.5 * gN * MU_N * B, 0.5 * gN * MU_N * B], abs=1e-12)
        s = _find(basis, J, 1.0, J + 0.5, J + 0.5)     # stretched: m_I = +1/2 only
        assert gN * B * _elem("zeeman_nuclear", basis, ctx, s, s) == pytest.approx(
            -0.5 * gN * MU_N * B, abs=1e-12)
    assert gN * MU_N * 1e3 == pytest.approx(4.008, abs=1e-3)


# ------------------------------------------------------- structure and PT-odd

def test_field_elements_are_hermitian(basis, ctx):
    """[HAM] V1 for the five field terms: a bra/ket swap in a spectator-theorem
    phase or a transposed 3j argument order shows up here and nowhere else."""
    for name in ("stark_z", "zeeman_Gpar", "zeeman_nuclear",
                 "pt_odd_edm", "pt_odd_scalar_pseudoscalar"):
        for i in range(0, len(basis), 7):
            for j in range(0, len(basis), 5):
                a = _elem(name, basis, ctx, i, j)
                b = _elem(name, basis, ctx, j, i)
                assert a == pytest.approx(b, abs=1e-12), f"{name} at ({i},{j})"


def test_stark_diagonal_block_is_odd_in_omega_and_zeeman_is_even(basis, ctx):
    """The structural fact behind [HAM] V6 and V8.

    Delta-J = 0: the Stark element is odd under Omega -> -Omega (it is
    -Omega m_F gamma_F d), while the G_par Zeeman element is EVEN (quadratic in
    n_hat: the explicit Omega and the geometry's Omega multiply out). That is
    exactly why the Zeeman gives g^u = g^l at leading order and the eEDM term,
    which is odd, does not.

    Uniquely catches an odd-in-Omega contamination of the Zeeman operator,
    which would fake an eEDM in the four-way chop; FAIL is reachable by
    dropping the explicit Omega factor from zeeman_Gpar, which makes it odd.
    """
    for J in (1, 2, 3):
        for F in (J - 0.5, J + 0.5):
            for mF in (0.5, F):
                p = _find(basis, J, 1.0, F, mF)
                m = _find(basis, J, -1.0, F, mF)
                assert _elem("stark_z", basis, ctx, p, p) == pytest.approx(
                    -_elem("stark_z", basis, ctx, m, m), abs=1e-12)
                assert _elem("zeeman_Gpar", basis, ctx, p, p) == pytest.approx(
                    _elem("zeeman_Gpar", basis, ctx, m, m), abs=1e-12)


def test_pt_odd_is_exactly_odd_in_omega_and_diagonal(basis, ctx):
    """[HAM] S2.12: H_eEDM = -d_e E_eff Omega/|Omega| -- +d_e E_eff on Omega = -1,
    -d_e E_eff on Omega = +1.

    Uniquely catches an eEDM term accidentally even in Omega, which would make
    it indistinguishable from a Zeeman background in the four-way chop. FAIL is
    reachable by dropping the sign(Omega): the Omega = -1 half then returns -1
    instead of +1.
    """
    for i in range(len(basis)):
        val = _elem("pt_odd_edm", basis, ctx, i, i)
        assert val == pytest.approx(-np.sign(basis["Om"][i]))
        assert _elem("pt_odd_scalar_pseudoscalar", basis, ctx, i, i) == pytest.approx(val)
    assert _elem("pt_odd_edm", basis, ctx, 0, 1) == 0.0     # the Omega partner


def test_leanhardt_half_convention_halves_the_pt_odd_element(basis, ctx):
    """[HAM] S2.12: Leanhardt Eq. 26 carries Omega/(2|Omega|); Ng removes the 1/2
    for consistency with the calculators' E_eff. The convention is a field, not a
    hand-entered factor."""
    half = _reconvention(ctx, edm_factor="leanhardt_half")
    i = 0
    assert _elem("pt_odd_edm", basis, half, i, i) == pytest.approx(
        0.5 * _elem("pt_odd_edm", basis, ctx, i, i))


def test_n_hat_convention_flips_the_omega_odd_terms_only(basis, ctx):
    """[HAM] OPEN-11: n_hat from Th to F flips the signed dipole and E_eff,
    hence every Stark and PT-odd element. The G_par Zeeman is quadratic in
    n_hat and does not move."""
    flipped = _reconvention(ctx, n_hat="Th_to_F")
    i = _find(basis, 1, 1.0, 1.5, 1.5)
    assert _elem("stark_z", basis, flipped, i, i) == pytest.approx(
        -_elem("stark_z", basis, ctx, i, i))
    assert _elem("pt_odd_edm", basis, flipped, i, i) == pytest.approx(
        -_elem("pt_odd_edm", basis, ctx, i, i))
    assert _elem("zeeman_Gpar", basis, flipped, i, i) == pytest.approx(
        _elem("zeeman_Gpar", basis, ctx, i, i))


def test_a_dipole_tagged_with_the_wrong_origin_is_refused():
    """[HAM] S2.7: d_mf is origin-dependent for an ION -- 3.37 D at the centre of
    nuclear mass is 2.74 D about the Th nucleus, a 20 % Stark error. A d_mf whose
    Param.convention disagrees with the conventions block's dipole_origin is a
    structural error and raises, rather than silently scaling every Stark
    element."""
    with pytest.raises(ValueError, match="dipole_origin"):
        thf_v1().with_(d_mf=Param(2.74, "D", convention="heavy_nucleus",
                                  status="ab-initio",
                                  source="Skripnikov & Titov 2015 Table II p.8"))
    ok = thf_v1().with_(d_mf=Param(3.46, "D", convention="center_of_mass",
                                   status="ab-initio", source="[HAM] S2.7"))
    assert ok.value("d_mf") == pytest.approx(3.46 * 0.5034118)
