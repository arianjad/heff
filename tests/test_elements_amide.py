"""Pinned numerical fixtures for the source-preserving amide kernels.

Fixture provenance: each literal was evaluated with the corresponding function
in ``C2V-Molecules/atm_core/physics.py`` at commit
``3b77b021ab2b3f9256f1af1b15ab3116eddb1068``.  The basis states came from its
``generate_states_and_blocks([0, 1, 2], [0, 1, 2], S=.5, I_N=1)`` m_F=1/2
block.  These tests retain literals rather than rebuilding expectations from
the port, and the committed suite has no dependency on C2V-Molecules.
"""

import pytest

from heff import elements_amide as amide


pytestmark = pytest.mark.filterwarnings(
    "ignore:bitcount function is deprecated:DeprecationWarning"
)


PARAMS = {
    "A": 10.0,
    "B": 4.0,
    "C": 3.0,
    "S": 0.5,
    "I_N": 1.0,
    "eps_xx": 0.2,
    "eps_yy": -0.1,
    "eps_zz": 0.3,
    "a_H": 1.2,
    "TH_aa": 0.4,
    "TH_bb": -0.15,
    "a_N": -0.7,
    "TN_aa": 0.5,
    "TN_bb": 0.2,
    "Q_aa": 0.8,
    "Q_bb": -0.25,
    "g_l": {
        0: {0: -0.12},
        1: {-1: 0.0, 0: 0.0, 1: 0.0},
        2: {-2: -0.04, -1: 0.0, 0: 0.07, 1: 0.0, 2: -0.04},
    },
}


def _matrix_args(q, q_prime):
    """Convert two 7-QN tuples to the source's alternating 14-QN arguments."""
    return tuple(value for pair in zip(q, q_prime) for value in pair)


# The chosen active cases jointly exercise off-diagonal K, N, J, F_N, and F.
ACTIVE_CASES = [
    (
        "rotational_hamiltonian_term",
        (2, -1, 1.5, 0.5, 1, 0.5, 0.5),
        (2, 1, 1.5, 0.5, 1, 0.5, 0.5),
        (PARAMS,),
        1.5,
    ),
    (
        "spin_rotation",
        (2, -2, 1.5, 0.5, 0, 0.5, 0.5),
        (1, 0, 1.5, 0.5, 0, 0.5, 0.5),
        (PARAMS,),
        -0.09185586535436917,
    ),
    (
        "spin_aniso_zeeman_hamiltonian_element",
        (2, -2, 1.5, 2.5, 0, 2.5, 0.5),
        (0, 0, 0.5, 1.5, 0, 1.5, 0.5),
        (PARAMS,),
        -0.010954451150103319,
    ),
    (
        "hydrogen_hyperfine",
        (1, -1, 0.5, 1.5, 1, 2.5, 0.5),
        (2, 1, 1.5, 2.5, 1, 2.5, 0.5),
        (PARAMS,),
        -0.014790199457749042,
    ),
    (
        "nitrogen_hyperfine",
        (2, -2, 1.5, 1.5, 0, 1.5, 0.5),
        (0, 0, 0.5, 1.5, 0, 1.5, 0.5),
        (PARAMS,),
        0.2053959590644373,
    ),
    (
        "quadrupole",
        (2, -2, 1.5, 0.5, 0, 0.5, 0.5),
        (0, 0, 0.5, 0.5, 0, 0.5, 0.5),
        (PARAMS,),
        -0.043301270189221946,
    ),
    (
        "stark_hamiltonian_element",
        (1, 0, 1.5, 2.5, 0, 2.5, 0.5),
        (2, 0, 2.5, 3.5, 0, 3.5, 0.5),
        (0.5, 1.0),
        0.4780914437337572,
    ),
    (
        "spin_iso_zeeman_hamiltonian_element",
        (2, -2, 1.5, 2.5, 0, 2.5, 0.5),
        (2, -2, 2.5, 3.5, 0, 3.5, 0.5),
        (0.5, 1.0),
        -0.47809144373375717,
    ),
    (
        "hydrogen_zeeman_hamiltonian_element",
        (1, -1, 1.5, 2.5, 1, 1.5, 0.5),
        (1, -1, 1.5, 2.5, 1, 2.5, 0.5),
        (),
        -0.7483314773547883,
    ),
    (
        "nitrogen_zeeman_hamiltonian_element",
        (2, -2, 2.5, 1.5, 0, 1.5, 0.5),
        (2, -2, 2.5, 2.5, 0, 2.5, 0.5),
        (1.0,),
        -0.7483314773547882,
    ),
]


@pytest.mark.parametrize("name,q,q_prime,extra,want", ACTIVE_CASES)
def test_amide_kernel_matches_pinned_source_fixture(name, q, q_prime, extra, want):
    """A sign, phase, recoupler, or parameter-semantics change breaks a fixture."""
    got = getattr(amide, name)(*_matrix_args(q, q_prime), *extra)
    assert got == pytest.approx(want, abs=3e-12)


@pytest.mark.parametrize("name,q,q_prime,extra,_want", ACTIVE_CASES)
def test_amide_kernel_rejects_delta_mf_for_lab_q_zero(name, q, q_prime, extra, _want):
    """Dropping a scalar delta or lab-q=0 3j selector makes this nonzero."""
    mismatched = (*q_prime[:-1], -0.5)
    got = getattr(amide, name)(*_matrix_args(q, mismatched), *extra)
    assert got == pytest.approx(0.0, abs=1e-15)


def test_edm_matches_pinned_off_diagonal_n_fixture_and_selectors():
    q = (1, -1, 1.5, 0.5, 1, 0.5, 0.5)
    q_prime = (2, -1, 1.5, 0.5, 1, 0.5, 0.5)
    assert amide.EDM(q, q_prime, 0.5) == pytest.approx(
        -0.4330127018922193, abs=3e-12
    )
    assert amide.EDM(q, (*q_prime[:-1], -0.5), 0.5) == 0.0


def test_delta_is_numeric_kronecker_delta():
    assert amide.delta(1.0, 1) == 1.0
    assert amide.delta(1.0, 1.5) == 0.0
