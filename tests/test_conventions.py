"""Conventions module -- the phase contract, in one place.

Nothing else in the package writes a phase (spec S3.2). The load-bearing test
here is test_ef_rules_disagree_at_S_one: the thesis rule P = (-1)^(J-S-l) is the
S = 1/2 specialisation and gives the OPPOSITE e/f label to Brown 1975 at S = 1,
which would invert every e/f label in ThF+
(docs/thf-plus-x3delta1-effective-hamiltonian.md S2.3 TRAP, OPEN-2).
"""
import numpy as np
import pytest

from heff.conventions import (Conventions, ef_label, n_hat_sign,
                              parity_operator, parity_phase, superposition_parity)
from heff.spec import enumerate_kets, thf_spec


def test_defaults_are_the_agreed_v1_block():
    c = Conventions()
    assert (c.n_hat, c.ef_rule, c.zeeman_sign, c.dg_def, c.edm_factor, c.dipole_origin) == (
        "F_to_Th", "brown1975", "plus_Gpar", "Delta", "ng", "center_of_mass")


def test_unknown_convention_value_raises():
    with pytest.raises(ValueError, match="n_hat"):
        Conventions(n_hat="F_to_Xe")


def test_n_hat_sign_flips_for_the_petrov_convention():
    assert n_hat_sign(Conventions()) == 1.0
    assert n_hat_sign(Conventions(n_hat="Th_to_F")) == -1.0


def test_parity_phase_for_3delta_is_minus_one_to_the_J_minus_one():
    for J in (1, 2, 3, 4):
        assert parity_phase(J, S=1.0, ell=0.0, s=0.0) == pytest.approx((-1.0) ** (J - 1))


def test_parity_operator_squares_to_identity_and_flips_omega():
    kets = enumerate_kets(thf_spec())
    P = parity_operator(kets, S=1.0, ell=0.0, s=0.0)
    assert np.allclose(P @ P, np.eye(len(kets)))
    assert np.allclose(P, P.T)
    ev = np.linalg.eigvalsh(P)
    assert np.allclose(np.abs(ev), 1.0)
    assert sum(ev > 0) == sum(ev < 0) == len(kets) // 2
    # P must move amplitude from Omega = +1 to Omega = -1 and nowhere else
    for i in range(len(kets)):
        j = int(np.flatnonzero(np.abs(P[i]) > 0)[0])
        assert kets["Om"][j] == -kets["Om"][i]
        assert (kets["J"][j], kets["F"][j], kets["mF"][j]) == (
            kets["J"][i], kets["F"][i], kets["mF"][i])


def test_superposition_parity_matches_the_operator():
    kets = enumerate_kets(thf_spec())
    P = parity_operator(kets, S=1.0, ell=0.0, s=0.0)
    for i in range(0, len(kets), 2):
        J = kets["J"][i]
        for sym in (+1, -1):
            vec = np.zeros(len(kets))
            vec[i], vec[i + 1] = 1.0 / np.sqrt(2), sym / np.sqrt(2)
            assert np.allclose(
                P @ vec, superposition_parity(J, sym, S=1.0, ell=0.0, s=0.0) * vec)


def test_ef_rules_disagree_at_S_one():
    """OPEN-2: applying the thesis rule to a 3Delta inverts every e/f label."""
    for J in (1, 2, 3, 4):
        for parity in (+1, -1):
            bc = ef_label(J, parity, rule="brown1975", ell=0.0)
            th = ef_label(J, parity, rule="thesis_S_half", S=1.0, ell=0.0)
            assert bc != th


def test_ef_rules_agree_at_S_one_half():
    for twoJ in (1, 3, 5, 7):
        J = twoJ / 2.0
        for parity in (+1, -1):
            assert ef_label(J, parity, rule="brown1975", ell=0.0) == ef_label(
                J, parity, rule="thesis_S_half", S=0.5, ell=0.0)


def test_thesis_rule_without_S_raises():
    with pytest.raises(ValueError, match="S"):
        ef_label(1, +1, rule="thesis_S_half", ell=0.0)


def test_v2_fields_have_the_agreed_defaults():
    c = Conventions()
    assert (c.quadrupole_convention, c.eqq2_norm, c.two_photon_norm) == (
        "bc_q0_is_negative_efg", "bc_9p52_q2", "bc_5p142_reduced")


def test_parity_operator_works_on_the_two_spin_dtype():
    """[HAM]/spec-v2 S2.2: parity_operator must key on every dtype field (with
    Om negated), not a four-literal (J, Om, F, mF) key -- otherwise two F1
    branches that land on the same (J, Om, F, mF) collide and KeyError or
    silently mismap."""
    kets = enumerate_kets(thf_spec("229"))
    P = parity_operator(kets, S=1.0, ell=0.0, s=0.0)
    assert np.allclose(P @ P, np.eye(len(kets)))
    assert np.allclose(P, P.T)
    for i in range(len(kets)):
        j = int(np.flatnonzero(np.abs(P[i]) > 0)[0])
        assert kets["Om"][j] == -kets["Om"][i]
        assert (kets["J"][j], kets["F1"][j], kets["F"][j], kets["mF"][j]) == (
            kets["J"][i], kets["F1"][i], kets["F"][i], kets["mF"][i])
