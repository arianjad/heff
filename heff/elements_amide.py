"""Low-level 14-QN matrix elements for the amide spectator basis.

This is a source-preserving port of the named functions in
``C:/Users/Arian/Code/C2V-Molecules/atm_core/physics.py`` at commit
``3b77b021ab2b3f9256f1af1b15ab3116eddb1068``.  The source comments cite
Sears (1984), Eqs. (21), (22), and (41), and Hirota, Eqs. (2.3.70),
(2.3.79)-(2.3.80), for the corresponding Zeeman and quadrupole conventions.
Only the Wigner backend import is changed here, to ``heff.wigner``.

The ``quadrupole`` expression is explicitly limited to nitrogen spin I_N=1:
its nuclear 3j normalization and 6j recoupler contain literal spin 1.  This
module does not generalize that expression or claim physics for another spin.
"""

import numpy as np

from .wigner import w3j as wigner_3j
from .wigner import w6j as wigner_6j
from .wigner import w9j as wigner_9j


def delta(x, y):
    return 1.0 if x == y else 0.0


def rotational_hamiltonian_term(N, N_prime, K, K_prime, J, J_prime,
                                F_N, F_N_prime, I_T, I_T_prime,
                                F, F_prime, m_F, m_F_prime, params):
    """Rotational Hamiltonian matrix element."""
    A = params['A']
    B = params['B']
    C = params['C']

    d = delta(m_F, m_F_prime) * delta(F_N, F_N_prime) * delta(J, J_prime) * delta(N, N_prime) * delta(F, F_prime) * delta(I_T, I_T_prime)
    c1 = (-1)**(N - K) * (-1) / np.sqrt(3) * (A + B + C) * wigner_3j(N, 0, N, -K, 0, K_prime) * N * (N + 1) * (2 * N + 1) * wigner_6j(N, N, 1, 0, 1, N)
    c2 = (-1)**(N - K) * 1 / np.sqrt(6) * (2 * A - B - C) * wigner_3j(N, 2, N, -K, 0, K_prime) * np.sqrt(5) * N * (N + 1) * (2 * N + 1) * wigner_6j(N, N, 1, 2, 1, N)
    c3 = (-1)**(N - K) * 1 / 2 * (B - C) * wigner_3j(N, 2, N, -K, -2, K_prime) * np.sqrt(5) * N * (N + 1) * (2 * N + 1) * wigner_6j(N, N, 1, 2, 1, N)
    c4 = (-1)**(N - K) * 1 / 2 * (B - C) * wigner_3j(N, 2, N, -K, 2, K_prime) * np.sqrt(5) * N * (N + 1) * (2 * N + 1) * wigner_6j(N, N, 1, 2, 1, N)
    return d * (c1 + c2 + c3 + c4)


def spin_rotation(N, N_prime, K, K_prime, J, J_prime,
                  F_N, F_N_prime, I_T, I_T_prime,
                  F, F_prime, m_F, m_F_prime, params):
    """Spin-rotation matrix element."""
    S = params['S']
    eps_xx = params['eps_xx']
    eps_yy = params['eps_yy']
    eps_zz = params['eps_zz']

    d = delta(m_F, m_F_prime) * delta(F_N, F_N_prime) * delta(F, F_prime) * delta(J, J_prime) * delta(I_T, I_T_prime)
    c1 = -np.sqrt(1 / 3) * (eps_xx + eps_yy + eps_zz) * np.sqrt(S * (S + 1) * (2 * S + 1)) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * (-1)**(J + N + S) * wigner_6j(N_prime, S, J, S, N, 1) * 1 / 2 * (np.sqrt(N * (N + 1) * (2 * N + 1)) * wigner_6j(1, 1, 0, N_prime, N, N) + np.sqrt(N_prime * (N_prime + 1) * (2 * N_prime + 1)) * wigner_6j(0, 1, 1, N_prime, N, N_prime)) * (-1)**(N - K) * wigner_3j(N, 0, N_prime, -K, 0, K_prime)
    c2 = np.sqrt(1 / 6) * (2 * eps_zz - eps_yy - eps_xx) * np.sqrt(5) * np.sqrt(S * (S + 1) * (2 * S + 1)) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * (-1)**(J + N + S) * wigner_6j(N_prime, S, J, S, N, 1) * 1 / 2 * (np.sqrt(N * (N + 1) * (2 * N + 1)) * wigner_6j(1, 1, 2, N_prime, N, N) + np.sqrt(N_prime * (N_prime + 1) * (2 * N_prime + 1)) * wigner_6j(2, 1, 1, N_prime, N, N_prime)) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 0, K_prime)
    c3 = 1 / 2 * (eps_xx - eps_yy) * np.sqrt(5) * np.sqrt(S * (S + 1) * (2 * S + 1)) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * (-1)**(J + N + S) * wigner_6j(N_prime, S, J, S, N, 1) * 1 / 2 * (np.sqrt(N * (N + 1) * (2 * N + 1)) * wigner_6j(1, 1, 2, N_prime, N, N) + np.sqrt(N_prime * (N_prime + 1) * (2 * N_prime + 1)) * wigner_6j(2, 1, 1, N_prime, N, N_prime)) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, -2, K_prime)
    c4 = 1 / 2 * (eps_xx - eps_yy) * np.sqrt(5) * np.sqrt(S * (S + 1) * (2 * S + 1)) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * (-1)**(J + N + S) * wigner_6j(N_prime, S, J, S, N, 1) * 1 / 2 * (np.sqrt(N * (N + 1) * (2 * N + 1)) * wigner_6j(1, 1, 2, N_prime, N, N) + np.sqrt(N_prime * (N_prime + 1) * (2 * N_prime + 1)) * wigner_6j(2, 1, 1, N_prime, N, N_prime)) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 2, K_prime)

    return d * (c1 + c2 + c3 + c4)


def spin_aniso_zeeman_hamiltonian_element(N, N_prime, K, K_prime, J, J_prime,
                                          F_N, F_N_prime, I_T, I_T_prime,
                                          F, F_prime, m_F, m_F_prime, params):
    """Anisotropic electron Zeeman matrix element."""
    S = params['S']
    I_N = params['I_N']
    g_l = params['g_l']

    d = delta(I_T, I_T_prime)
    c1 = (-1)**(F - m_F + F_prime + F_N + 1 + I_T + F_N_prime + J + 1 + I_N)
    c2 = np.sqrt((2 * F + 1) * (2 * F_prime + 1)) * np.sqrt((2 * F_N + 1) * (2 * F_N_prime + 1))
    c3 = wigner_3j(F, 1, F_prime, -m_F, 0, m_F_prime) * wigner_6j(F_N_prime, F_prime, I_T, F, F_N, 1) * wigner_6j(J_prime, F_N_prime, I_N, F_N, J, 1)
    rme = lambda _k: np.sqrt(3 * (2 * J + 1) * (2 * J_prime + 1)) * wigner_9j(N, N_prime, _k, S, S, 1, J, J_prime, 1) * sum([
        g_l[_k][_q] * (-1)**(N - K) * wigner_3j(N, _k, N_prime, -K, _q, K_prime) * np.sqrt((2 * N + 1) * (2 * N_prime + 1) * S * (S + 1) * (2 * S + 1)) for _q in np.arange(-_k, _k + 1,)
    ])
    c4 = sum([
        (-1)**_k * np.sqrt((1 / 3) * (2 * _k + 1)) * rme(_k) for _k in [0, 1, 2]
    ])
    return d * c1 * c2 * c3 * c4


def hydrogen_hyperfine(N, N_prime, K, K_prime, J, J_prime,
                       F_N, F_N_prime, I_T, I_T_prime,
                       F, F_prime, m_F, m_F_prime, params):
    """Hydrogen hyperfine matrix element."""
    S = params['S']
    I_N = params['I_N']
    a_H = params['a_H']
    TH_aa = params['TH_aa']
    TH_bb = params['TH_bb']

    d = delta(F, F_prime) * delta(m_F, m_F_prime) * delta(I_T, I_T_prime)
    c = (-1)**(F_N_prime + I_T + F) * np.sqrt((S * (S + 1) * (2 * S + 1)) * (I_T * (I_T + 1) * (I_T * 2 + 1)) * ((2 * J + 1) * (2 * J_prime + 1))) * wigner_6j(I_T, F_N_prime, F, F_N, I_T, 1) * np.sqrt((2 * F_N + 1) * (2 * F_N_prime + 1)) * (-1)**(F_N_prime + J + I_N + 1) * wigner_6j(J_prime, F_N_prime, I_N, F_N, J, 1)
    c1 = delta(N, N_prime) * delta(K, K_prime) * (-1)**(J + S + N + 1) * wigner_6j(S, J_prime, N, J, S, 1) * a_H
    c2 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 0, K_prime) * 0.5 * (TH_aa)
    c5 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 2, K_prime) * np.sqrt(1 / 24) * (TH_aa + 2 * TH_bb)
    c6 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, -2, K_prime) * np.sqrt(1 / 24) * (TH_aa + 2 * TH_bb)
    return d * c * (c1 - (c2 + c5 + c6))


def nitrogen_hyperfine(N, N_prime, K, K_prime, J, J_prime,
                       F_N, F_N_prime, I_T, I_T_prime,
                       F, F_prime, m_F, m_F_prime, params):
    """Nitrogen hyperfine matrix element."""
    S = params['S']
    I_N = params['I_N']
    a_N = params['a_N']
    TN_aa = params['TN_aa']
    TN_bb = params['TN_bb']

    d = delta(F, F_prime) * delta(m_F, m_F_prime) * delta(F_N, F_N_prime) * delta(I_T, I_T_prime)
    c = (-1)**(J_prime + I_N + F_N) * np.sqrt((S * (S + 1) * (2 * S + 1)) * (I_N * (I_N + 1) * (I_N * 2 + 1)) * ((2 * J + 1) * (2 * J_prime + 1))) * wigner_6j(I_N, J_prime, F_N, J, I_N, 1)
    c1 = delta(N, N_prime) * delta(K, K_prime) * (-1)**(J + S + N + 1) * wigner_6j(S, J_prime, N, J, S, 1) * a_N
    c2 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 0, K_prime) * 0.5 * (TN_aa)
    c5 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 2, K_prime) * np.sqrt(1 / 24) * (TN_aa + 2 * TN_bb)
    c6 = np.sqrt(30) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_9j(J, J_prime, 1, N, N_prime, 2, S, S, 1) * (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, -2, K_prime) * np.sqrt(1 / 24) * (TN_aa + 2 * TN_bb)
    return d * c * (c1 - (c2 + c5 + c6))


def quadrupole(N, N_prime, K, K_prime, J, J_prime,
               F_N, F_N_prime, I_T, I_T_prime,
               F, F_prime, m_F, m_F_prime, params):
    """Nitrogen quadrupole matrix element, limited to I_N=1 (14N)."""
    S = params['S']
    Q_aa = params['Q_aa']
    Q_bb = params['Q_bb']

    d = delta(F, F_prime) * delta(F_N, F_N_prime) * delta(m_F, m_F_prime) * delta(I_T, I_T_prime)
    c1 = 0.5 * (wigner_3j(1, 2, 1, -1, 0, 1))**(-1) * (-1)**(J_prime + 1 + F_N) * (-1)**(N + S + J_prime)
    c2 = wigner_6j(1, J, F_N, J_prime, 1, 2) * np.sqrt((2 * J_prime + 1) * (2 * J + 1)) * np.sqrt((2 * N_prime + 1) * (2 * N + 1)) * wigner_6j(N, J, S, J_prime, N_prime, 2)
    c3 = (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 0, K_prime) * Q_aa / 2
    c4 = ((-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, 2, K_prime) + (-1)**(N - K) * wigner_3j(N, 2, N_prime, -K, -2, K_prime)) * (2 * Q_bb + Q_aa) / (2 * np.sqrt(6))
    return d * c1 * c2 * (c3 + c4)


def stark_hamiltonian_element(N, N_prime, K, K_prime, J, J_prime,
                              F_N, F_N_prime, I_T, I_T_prime,
                              F, F_prime, m_F, m_F_prime, S, I_N):
    """Calculate the Stark effect Hamiltonian matrix element."""
    d = delta(I_T, I_T_prime)
    c1 = (-1)**(F - m_F) * wigner_3j(F, 1, F_prime, -m_F, 0, m_F_prime) * (-1)**(F_N + I_T + F_prime + 1) * np.sqrt((2 * F + 1) * (2 * F_prime + 1))
    c2 = wigner_6j(F_N_prime, F_prime, I_T, F, F_N, 1) * (-1)**(J + I_N + F_N_prime + 1) * np.sqrt((2 * F_N + 1) * (2 * F_N_prime + 1)) * wigner_6j(J_prime, F_N_prime, I_N, F_N, J, 1)
    c3 = (-1)**(N + S + J_prime + 1) * np.sqrt((2 * J + 1) * (2 * J_prime + 1)) * wigner_6j(N_prime, J_prime, S, J, N, 1)
    c4 = (-1)**(N - K) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * wigner_3j(N, 1, N_prime, -K, 0, K_prime)
    return d * c1 * c2 * c3 * c4


def spin_iso_zeeman_hamiltonian_element(N, N_prime, K, K_prime, J, J_prime,
                                        F_N, F_N_prime, I_T, I_T_prime,
                                        F, F_prime, m_F, m_F_prime, S, I_N):
    """Calculate the isotropic electron Zeeman matrix element."""
    d = delta(N, N_prime) * delta(K, K_prime) * delta(I_T, I_T_prime)
    c1 = (-1)**(F - m_F + F_prime + F_N + 1 + I_T + F_N_prime + J + 1 + I_N + J + N + S + 1)
    c2 = np.sqrt((2 * F + 1) * (2 * F_prime + 1)) * np.sqrt((2 * F_N + 1) * (2 * F_N_prime + 1)) * np.sqrt((2 * J + 1) * (2 * J_prime + 1)) * np.sqrt((2 * S + 1) * (S + 1) * S)
    c3 = wigner_3j(F, 1, F_prime, -m_F, 0, m_F_prime) * wigner_6j(F_N_prime, F_prime, I_T, F, F_N, 1) * wigner_6j(J_prime, F_N_prime, I_N, F_N, J, 1) * wigner_6j(S, J_prime, N, J, S, 1)
    return d * c1 * c2 * c3


def hydrogen_zeeman_hamiltonian_element(N, N_prime, K, K_prime, J, J_prime,
                                        F_N, F_N_prime, I_T, I_T_prime,
                                        F, F_prime, m_F, m_F_prime):
    """Calculate the hydrogen nuclear Zeeman matrix element."""
    d = delta(N, N_prime) * delta(K, K_prime) * delta(J, J_prime) * delta(F_N, F_N_prime) * delta(I_T, I_T_prime)
    c1 = (-1)**(F - m_F + F + F_N + 1 + I_T)
    c2 = np.sqrt((2 * F + 1) * (2 * F_prime + 1)) * np.sqrt(I_T * (I_T + 1) * (2 * I_T + 1))
    c3 = wigner_3j(F, 1, F_prime, -m_F, 0, m_F_prime) * wigner_6j(I_T, F_prime, F_N, F, I_T, 1)
    return d * c1 * c2 * c3


def nitrogen_zeeman_hamiltonian_element(N, N_prime, K, K_prime, J, J_prime,
                                        F_N, F_N_prime, I_T, I_T_prime,
                                        F, F_prime, m_F, m_F_prime, I_N):
    """Nitrogen nuclear Zeeman matrix element (lab-frame T^1_0(I_N))."""
    d = delta(N, N_prime) * delta(K, K_prime) * delta(J, J_prime) * delta(I_T, I_T_prime)
    c1 = (-1)**(F - m_F + F_prime + F_N + 1 + I_T + F_N + J + 1 + I_N)
    c2 = np.sqrt((2 * F + 1) * (2 * F_prime + 1)) * np.sqrt((2 * F_N + 1) * (2 * F_N_prime + 1)) * np.sqrt(I_N * (I_N + 1) * (2 * I_N + 1))
    c3 = wigner_3j(F, 1, F_prime, -m_F, 0, m_F_prime) * wigner_6j(F_N_prime, F_prime, I_T, F, F_N, 1) * wigner_6j(I_N, F_N_prime, J, F_N, I_N, 1)
    return d * c1 * c2 * c3


def EDM(q, q_prime, S):
    """EDM matrix element between two basis states."""
    N, K, J, F_N, I_T, F, m_F = q
    N_prime, K_prime, J_prime, F_N_prime, I_T_prime, F_prime, m_F_prime = q_prime
    d = delta(F, F_prime) * delta(m_F, m_F_prime) * delta(J, J_prime) * delta(F_N, F_N_prime) * delta(I_T, I_T_prime)
    c1 = (-1)**(N_prime + S + J) * wigner_6j(S, N_prime, J, N, S, 1)
    c2 = (-1)**(N - K) * np.sqrt((2 * N + 1) * (2 * N_prime + 1)) * wigner_3j(N, 1, N_prime, -K, 0, K_prime) * np.sqrt(S * (S + 1) * (2 * S + 1))

    return d * c1 * c2
