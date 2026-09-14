"""Numerical check behind audit §3 of
docs/lit/lookup-effective-zeeman-tensor.md. Run with
`conda run -n structure python docs/lit/2026-09-14-termC-second-order-model.py`.
This is documentation, not a test.

Does the constant body tensor B.G.J reproduce the 2nd-order Coriolis x Zeeman
effective operator on the 3Delta1 (Omega=+-1) manifold, off-diagonals included?

Model (Petrov 2014/2025 structure): electronic states X(+1), X(-1) and ONE
intermediate n with Omega_n = 0 or +-2, electronic energy T_n above X.
V = H_Cor + H_Z with
  H_Cor = -B_rot (J^b_+ J^e_- + J^b_- J^e_+)      (anomalous body ladders)
  H_Z   = mu_B B sum_q D^{1*}_{0q} G^b_q,  G = L + |g_S| S  (Petrov's G_perp element)
Second order, electronic denominators only (Petrov's approximation):
  Z_eff = P H_Z P - sum_n [P H_Cor Q_n H_Z P + P H_Z Q_n H_Cor P] / T_n
Fit (G_par-G_perp, G_perp, G_Delta) from three elements, check ALL elements.
"""
import numpy as np
from heff.wigner import w3j
from heff.elements_c import dipole_geometry, _ph

M = 1.0
JMAX = 5
B_ROT = 0.243
G_PAR = 0.047


def D1(Jp, Omp, J, Om, q, k=1, p=0):
    """<J' Om' M| D^{k*}_{pq} |J Om M> at M'=M (B&C 5.186 form)."""
    if abs(Omp - Om - q) > 1e-9 or abs(Omp) > Jp or abs(Om) > J:
        return 0.0
    return (_ph(Jp - M) * w3j(Jp, k, J, -M, p, M) * _ph(Jp - Omp)
            * np.sqrt((2 * J + 1) * (2 * Jp + 1)) * w3j(Jp, k, J, -Omp, q, Om))


def ladder(J, Om, dOm):
    """<J, Om+dOm | J^b_{-/+} | J, Om>: body ladder that shifts Om by dOm=+-1."""
    return np.sqrt(max(J * (J + 1) - Om * (Om + dOm), 0.0))


def run(Om_n, s, Gp, Dn, Tn):
    """Om_n: intermediate Omega (0 or 2). s: reflection sign (+-1) relating the
    Omega<0 side to the Omega>0 side. Gp, Dn: Petrov's G_perp and Delta for n."""
    # electronic states with Omega; intermediate has +-|Om_n| (one state if 0)
    el = [("X+", 1.0), ("X-", -1.0)]
    el += [("n0", 0.0)] if Om_n == 0 else [("n+", float(Om_n)), ("n-", -float(Om_n))]
    # <a|O_+|b> tables for O = J^e (times 2B => Delta) and O = G. Both share s.
    # Petrov: Delta = 2B <higher Om| J^e_+ |lower Om>, G_perp analogous.
    Jplus, Gplus = {}, {}
    if Om_n == 0:
        Jplus[("X+", "n0")] = Dn / (2 * B_ROT); Gplus[("X+", "n0")] = Gp   # 0 -> +1
        Jplus[("n0", "X-")] = s * Dn / (2 * B_ROT); Gplus[("n0", "X-")] = s * Gp  # -1 -> 0
    else:
        Jplus[("n+", "X+")] = Dn / (2 * B_ROT); Gplus[("n+", "X+")] = Gp   # +1 -> +2
        Jplus[("X-", "n-")] = s * Dn / (2 * B_ROT); Gplus[("X-", "n-")] = s * Gp  # -2 -> -1

    def Oplus(a, b, tab):  # <a|O_+|b>
        return tab.get((a, b), 0.0)

    def Ominus(a, b, tab):  # <a|O_-|b> = <b|O_+|a>* (real)
        return tab.get((b, a), 0.0)

    basis = [(n, Om, J) for (n, Om) in el for J in range(1, JMAX + 1) if J >= abs(Om)]
    if Om_n == 0:
        basis += [("n0", 0.0, 0)]
    idx = {b: i for i, b in enumerate(basis)}
    N = len(basis)
    Hc, Hz = np.zeros((N, N)), np.zeros((N, N))
    for (na, Oa, Ja) in basis:
        for (nb, Ob, Jb) in basis:
            i, j = idx[(na, Oa, Ja)], idx[(nb, Ob, Jb)]
            dq = Oa - Ob
            # Coriolis: same J; body ladder shifts Om by dq, electronic ladder the same
            if Ja == Jb and abs(dq) == 1:
                oe = Oplus(na, nb, Jplus) if dq == 1 else Ominus(na, nb, Jplus)
                Hc[i, j] = -B_ROT * ladder(Ja, Ob, dq) * oe
            # Zeeman: D^{1*}_{0,dq} times electronic spherical component
            if abs(dq) <= 1:
                if dq == 0:
                    ge = G_PAR * Oa if (na == nb and na.startswith("X")) else 0.0
                elif dq == 1:
                    ge = -Oplus(na, nb, Gplus) / np.sqrt(2)
                else:
                    ge = +Ominus(na, nb, Gplus) / np.sqrt(2)
                Hz[i, j] = D1(Ja, Oa, Jb, Ob, dq) * ge
    assert np.allclose(Hc, Hc.T) and np.allclose(Hz, Hz.T)
    P = np.array([b[0].startswith("X") for b in basis])
    Q = ~P
    Z2 = Hz.copy()
    Z2[np.ix_(P, P)] -= (Hc[np.ix_(P, Q)] @ Hz[np.ix_(Q, P)] + Hz[np.ix_(P, Q)] @ Hc[np.ix_(Q, P)]) / Tn
    Zeff = Z2[np.ix_(P, P)]
    xb = [b for b in basis if b[0].startswith("X")]

    # candidate operators on the X manifold
    A = np.zeros((len(xb), len(xb)))   # Omega * D^{1*}_{00}   (existing zeeman_Gpar)
    Bm = M * np.eye(len(xb))           # B.J = B M
    C = np.zeros((len(xb), len(xb)))   # ladder x D^{1*}_{0,-/+1}: Om +1 <-> -1
    K2 = np.zeros((len(xb), len(xb)))  # the spec's axial_geometry(k=2, q=+-2, p=0)
    for a, (na, Oa, Ja) in enumerate(xb):
        for b, (nb, Ob, Jb) in enumerate(xb):
            if Oa == Ob:
                A[a, b] = Oa * D1(Ja, Oa, Jb, Ob, 0)
            else:
                dq = Oa - Ob  # +-2
                # J^b ladder takes ket Om -> Om + dq/2 (stays in J=Jb), D^{1*} does the rest
                # C = D^{1*}_{0,-1} J^b_+  -  D^{1*}_{0,+1} J^b_-  (Hermitian: (D^{1*}_{0,-1})^+ = -D^{1*}_{0,+1})
                C[a, b] = -np.sign(dq) * ladder(Jb, Ob, dq / 2) * D1(Ja, Oa, Jb, Ob + dq / 2, dq / 2)
                K2[a, b] = D1(Ja, Oa, Jb, Ob, dq, k=2)
    # ordering check: D first then ladder must agree (they commute)
    C2 = np.zeros_like(C)
    for a, (na, Oa, Ja) in enumerate(xb):
        for b, (nb, Ob, Jb) in enumerate(xb):
            if Oa != Ob:
                dq = Oa - Ob
                C2[a, b] = -np.sign(dq) * D1(Ja, Oa - dq / 2, Jb, Ob, dq / 2) * ladder(Ja, Oa - dq / 2, dq / 2)
    order_gap = np.abs(C - C2).max()
    herm_gap = np.abs(C - C.T).max()
    # normalisation: <J, parity +-| C |J, parity +-> should be +-M (J-independent) if
    # G_Delta is to mean "g = ... +- G_Delta"
    cdiag = []
    for J in range(1, JMAX + 1):
        i, j = xb.index(("X+", 1.0, J)), xb.index(("X-", -1.0, J))
        cdiag.append(C[i, j] / M)

    def fit(ops):
        X = np.stack([o.ravel() for o in ops], axis=1)
        coef, *_ = np.linalg.lstsq(X, Zeff.ravel(), rcond=None)
        return coef, np.abs(X @ coef - Zeff.ravel()).max()

    (a, gperp, gdel), res = fit([A, Bm, C])
    _, res_k2 = fit([A, Bm, K2])
    _, res_noC = fit([A, Bm])
    # Petrov's closed form for this one-state model
    corr = Gp * Dn / (0.0 - Tn)
    out = dict(Om_n=Om_n, s=s, fit_Gpar_minus_Gperp=a, fit_Gperp=gperp, fit_GDelta=gdel,
               residual_ABC=res, residual_A_B_k2=res_k2, residual_AB_only=res_noC,
               C_ordering_gap=order_gap, C_hermiticity_gap=herm_gap,
               C_offdiag_over_M_by_J=[round(x, 9) for x in cdiag], petrov_corr=corr)
    # diagonal g in parity eigenstates for J=1..3 vs closed forms
    gs = {}
    for J in range(1, 4):
        i, j = xb.index(("X+", 1.0, J)), xb.index(("X-", -1.0, J))
        blk = Zeff[np.ix_([i, j], [i, j])]
        w = np.linalg.eigvalsh(blk)
        g = -w / M  # E = -g mu_B B M, mu_B=B=1
        if Om_n == 0:
            g_pet = (-G_PAR / (J * (J + 1)) + corr, -G_PAR / (J * (J + 1)))
        else:
            g_pet = (-G_PAR / (J * (J + 1)) + corr * (J + 2) * (J - 1) / (2 * J * (J + 1)),) * 2
        g_tensor = (-gperp - a / (J * (J + 1)) + gdel, -gperp - a / (J * (J + 1)) - gdel)
        gs[J] = dict(model=sorted(g), petrov=sorted(g_pet), tensor=sorted(g_tensor))
    out["g_by_J"] = gs
    # Field-free 2nd-order Omega-doubling from the same intermediate, same phases:
    E2 = -(Hc[np.ix_(P, Q)] @ Hc[np.ix_(Q, P)]) / Tn
    i, j = xb.index(("X+", 1.0, 1)), xb.index(("X-", -1.0, 1))
    doub = E2[i, j]                       # cm^-1, J=1 Omega-flip element
    zflip = Zeff[i, j]                    # per mu_B B, M=1
    blk = E2[np.ix_([i, j], [i, j])]
    w, v = np.linalg.eigh(blk)            # lower, upper parity level at J=1
    gl, gu = [-(v[:, k] @ Zeff[np.ix_([i, j], [i, j])] @ v[:, k]) / M for k in (0, 1)]
    out["J1_doubling_flip_elem_cm-1"] = doub
    out["J1_zeeman_flip_elem"] = zflip
    out["J1_split_MHz"] = (w[1] - w[0]) * 29979.2458
    out["J1_g_upper_minus_g_lower"] = gu - gl
    out["J1_ratio_zflip_over_doub"] = zflip / doub if doub else None
    return out


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    for Om_n, Gp, Dn, Tn in ((0, 0.48 * 1.088, 0.250 * 1.088, 314.0), (2, -2.41, -0.510, 1052.0)):
        for s in (+1, -1):
            r = run(Om_n, s, Gp, Dn, Tn)
            gs = r.pop("g_by_J")
            print({k: (round(v, 10) if isinstance(v, float) else v) for k, v in r.items()})
            for J, d in gs.items():
                print("   J=%d" % J, {k: [round(x, 9) for x in v] for k, v in d.items()})
    # m-dependence of the k=2 element: linear in M or not?
    vals = []
    for m in (1.0, 2.0):
        M = m
        vals.append(D1(2, -1.0, 2, 1.0, -2.0, k=2))
    print("k=2 element at M=1,2:", vals, " ratio", vals[1] / vals[0], "(linear would be 2)")
    print("k=2 has dJ=2 element J=1->3:", D1(3, -1.0, 1, 1.0, -2.0, k=2))
