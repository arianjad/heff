"""Rank-K effective two-photon transition geometry; see [HAM] S9.5 and B&C 5.141–5.142.

``REGISTRY_2G`` holds transition operators, never Hamiltonian terms. Registered
channels are (K, |ΔOmega|)=(0,0),(2,0),(2,2); K=1 has no exact-closure operator.
Geometry is dimensionless with ``two_photon_norm='bc_5p142_reduced'``. Closure
requires Δ much larger than intermediate rotational structure ([HAM] S9.5.5).
"""
import numpy as np

from .elements_c2 import axial_geometry
from .spectra import _strengths_from_matrices
from .terms import Rules, term
from .wigner import w3j

REGISTRY_2G = {}

_SQ2 = np.sqrt(2.0)


def _spherical(v):
    """Return Condon-Shortley spherical components ([HAM] S9.5.1(3))."""
    vx, vy, vz = (complex(c) for c in np.asarray(v).reshape(3))
    return {1: -(vx + 1j * vy) / _SQ2, 0: vz, -1: (vx - 1j * vy) / _SQ2}


def _leg(eps):
    """Return c[p]=(-1)^p eps[-p]; conjugate the Raman bra-side vector first."""
    s = _spherical(eps)
    return {p: (-1.0) ** p * s[-p] for p in (-1, 0, 1)}


def _cg(j1, m1, j2, m2, J, M):
    """Return the Clebsch-Gordan coefficient from B&C 5.141."""
    return (-1.0) ** (j1 - j2 + M) * np.sqrt(2 * J + 1) * w3j(j1, j2, J, m1, m2, -M)


def dyad_weights(eps1, eps2):
    """Return the Raman dyad weights for Cartesian Jones vectors ([HAM] S9.5.1).

    The bra-side vector is conjugated and first. K=1 weights are returned but
    have no registered exact-closure operator ([HAM] S9.5.3).
    """
    c_a = _leg(np.conj(np.asarray(eps2, dtype=complex)))
    c_b = _leg(np.asarray(eps1, dtype=complex))
    return {(K, P): sum(_cg(1, p_a, 1, P - p_a, K, P) * c_a[p_a] * c_b[P - p_a]
                        for p_a in (-1, 0, 1) if abs(P - p_a) <= 1)
            for K in (0, 1, 2) for P in range(-K, K + 1)}


def two_photon_geometry(bra, ket, ctx, *, K, P):
    """Return two-spin alpha^K geometry via ``axial_geometry`` ([HAM] S9.5.1).

    It is zero for |q|>K and requires unresolved intermediate hyperfine.
    """
    q = float(bra["Om"]) - float(ket["Om"])
    if abs(q) > float(K):
        return 0.0
    return axial_geometry(bra, ket, ctx, k=K, q=q, p=P)


def two_photon_matrix(kets_a, kets_b, ctx, *, K, dOmega, P):
    """Return one dimensionless channel matrix; ``dOmega=2`` includes both signs."""
    if dOmega not in (0, 2):
        raise ValueError(
            f"dOmega must be 0 or 2 (2 meaning |dOmega| = 2, both signs), got "
            f"{dOmega!r}; within X 3Delta1 q = Om' - Om is 0 or +-2 and |q| <= K "
            "([HAM] S9.5.2)")
    out = np.zeros((len(kets_a), len(kets_b)))
    for i in range(len(kets_a)):
        for j in range(len(kets_b)):
            if abs(float(kets_a["Om"][i]) - float(kets_b["Om"][j])) != dOmega:
                continue
            out[i, j] = two_photon_geometry(kets_a[i], kets_b[j], ctx, K=K, P=P)
    return out



_CITE = (
    "[HAM] S9.5, the rank-K effective two-photon operator. B&C Eqs. (5.141) "
    "and (5.142), PDF p.198 / book p.166: two rank-1 operators couple to "
    "K = 0, 1, 2 and nothing else, and with a COMMON detuning (5.142) read "
    "backwards makes the intermediate sum the reduced element of a single "
    "rank-K operator -- one scalar per (K, dOmega) channel times parameter-free "
    "geometry ([HAM] S9.5.1(1)). The geometry is [HAM] S9.5.1(2) = S9.1's "
    "two-spectator chain with k -> K: B&C (5.172) + (5.174) twice + (5.186), "
    "PDF pp.205-207 / book pp.173-175, i.e. elements_c2.axial_geometry at "
    "k = K -- no new algebra. Normalisation conventions.two_photon_norm = "
    "'bc_5p142_reduced': the RANK-K reduced element is unity per channel, "
    "<eta'||alpha^K||eta> = 1 for each (K, dOmega) -- not unit ONE-photon "
    "reduced elements, which is a different (and unmade) claim: the "
    "intermediate sum has already been closed away by (5.142). So a strength "
    "comes "
    "out in units of alpha^2. K = 1 is NOT registered ([HAM] S9.5.3, OPEN-21: "
    "identically zero in exact closure). TRANSITION operator in REGISTRY_2G, "
    "never summable into a Hamiltonian ([SPEC-v2] S3.2). Validity "
    "[HAM] S9.5.5: needs Delta >> the intermediate B_i ~ 7 GHz, which the "
    "0.16-1.5 GHz JILA detunings do not satisfy."
)

_K2 = (-2, -1, 0, 1, 2)


def _channel(K, dOmega):
    """Return the registered element with P fixed by Δm_F."""
    def fn(bra, ket, ctx):
        if abs(float(bra["Om"]) - float(ket["Om"])) != dOmega:
            return 0.0
        return two_photon_geometry(bra, ket, ctx, K=K,
                                   P=float(bra["mF"]) - float(ket["mF"]))
    fn.__name__ = f"two_photon_K{K}_dOm{dOmega}"
    return fn


two_photon_K0_dOm0 = term(
    name="two_photon_K0_dOm0", param=("alpha_K0_dOm0",), cases=("c2",),
    registry=REGISTRY_2G, hermitian=True, real=True,
    rules=Rules(dJ=(0,), dOm=(0.0,), dF1=(0,), dF=(0,), dmF=(0,)),
    cite="The SCALAR channel, K = 0, dOmega = 0: alpha^0 is the isotropic "
         "polarisability eps2*.eps1. At K = 0 every 3j and 6j collapses to a "
         "diagonal, so dJ = dOmega = dF1 = dF = dm_F = 0. " + _CITE
)(_channel(0, 0))

two_photon_K2_dOm0 = term(
    name="two_photon_K2_dOm0", param=("alpha_K2_dOm0",), cases=("c2",),
    # P=Δm_F gives odd-Δm_F antisymmetry by spherical-tensor reciprocity.
    registry=REGISTRY_2G, hermitian=False, real=True,
    rules=Rules(dJ=_K2, dOm=(0.0,), dF1=_K2, dF=_K2, dmF=_K2),
    cite="The rank-2 channel at dOmega = 0. Selection rules as data, "
         "[HAM] S9.5.4: |dF| <= K and |dm_F| = |P| <= K from B&C (5.172)'s "
         "triangle and projection (so NEVER +-3), |dJ| <= K from the "
         "molecule-frame 3j, |dF1| <= K from the I_Th spectator 6j, and parity "
         "EVEN -- at zero field it does not connect e to f ([2gamma] S3.3, "
         "gate V24). " + _CITE
)(_channel(2, 0))

two_photon_K2_dOm2 = term(
    name="two_photon_K2_dOm2", param=("alpha_K2_dOm2",), cases=("c2",),
    registry=REGISTRY_2G, hermitian=False, real=True,
    rules=Rules(dJ=_K2, dOm=(-2.0, 2.0), dF1=_K2, dF=_K2, dmF=_K2),
    cite="The |dOmega| = 2 channel, which exists at K = 2 ONLY: |q| = 2 needs "
         "K >= 2 and two E1 legs bound K <= 2 ([HAM] S9.5.2, verified there "
         "numerically -- the K = 0 and K = 1 molecule-frame 3j's are 0.000000 "
         "at Om' = -1, Om = +1). Physically it needs an Omega = 0 intermediate "
         "([2gamma] S3.3). Both q = +2 and q = -2 belong to this ONE channel: "
         "they share this alpha and parity maps each into the other, so a "
         "relative sign between them would fake an e <-> f two-photon line "
         "(gate V24). Same |dF|, |dJ|, |dF1|, |dm_F| <= K rules as the "
         "dOmega = 0 rank-2 channel. " + _CITE
)(_channel(2, 2))

#: The registered channels, (K, |dOmega|, alpha knob). [HAM] S9.5.2/S9.5.3.
CHANNELS = ((0, 0, "alpha_K0_dOm0"), (2, 0, "alpha_K2_dOm0"),
            (2, 2, "alpha_K2_dOm2"))


def two_photon_line_strengths(evals_a, evecs_a, kets_a, evals_b, evecs_b, kets_b,
                              ctx, *, eps1, eps2, alphas):
    """Line positions (MHz) and two-photon strengths between two separately
    diagonalised blocks. Mirrors heff.spectra.line_strengths.

    AMPLITUDES ARE SUMMED OVER EVERY (K, dOmega, P) CHANNEL AND THEN SQUARED --
    heff.spectra._strengths_from_matrices is reused, so the package implements
    the ordering in one place. Cossel's Fig. 6.18 (thesis p.224) reports a
    MEASURED cancellation of two sigma pathways "because of the signs of the
    Wigner 3j coefficients" ([2gamma] S3.0).

    alphas: {param name: complex} over the three registered knobs
    'alpha_K0_dOm0', 'alpha_K2_dOm0', 'alpha_K2_dOm2' -- the scalar per
    channel, in (MHz/(V/cm))^2/MHz. Missing keys are 0. The two field
    amplitudes are the caller's too: this function multiplies geometry by
    alpha and by the polarisation dyad, nothing else, so `strengths` comes out
    in units of alpha^2 (conventions.two_photon_norm = 'bc_5p142_reduced').

    eps1, eps2: the Cartesian Jones vectors of the ket-side and bra-side
    photons. RAMAN reading -- eps2 is conjugated inside `dyad_weights`; see its
    docstring for what that does to a "sigma+sigma+" label.
    """
    unknown = set(alphas) - {knob for _, _, knob in CHANNELS}
    if unknown:
        raise ValueError(
            f"unknown two-photon scalars {sorted(unknown)}; the registered "
            f"channels are {sorted(k for _, _, k in CHANNELS)} "
            "(K = 1 is not one of them -- [HAM] S9.5.3, OPEN-21)")
    w = dyad_weights(eps1, eps2)
    mats, weights = {}, {}
    for K, dOmega, knob in CHANNELS:
        alpha = complex(alphas.get(knob, 0.0))
        if abs(alpha) == 0.0:
            continue
        for P in range(-K, K + 1):
            c = alpha * w[(K, P)]
            if abs(c) == 0.0:
                continue
            mats[(K, dOmega, P)] = two_photon_matrix(kets_a, kets_b, ctx, K=K,
                                                     dOmega=dOmega, P=P)
            weights[(K, dOmega, P)] = c
    return _strengths_from_matrices(evals_a, evecs_a, evals_b, evecs_b, mats,
                                    weights=weights)
