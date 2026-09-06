"""The rank-K effective two-photon (2 x E1) operator within one electronic state.

EVERY form in this file is copied from docs/thf-plus-x3delta1-effective-
hamiltonian.md S9.5 -- the single source for the two-photon operator -- with its
Brown & Carrington equation number and PDF/book page carried into the decorator,
in the style heff/elements_c.py and heff/elements_c2.py use. Nothing here was
re-derived. If a check disagrees with one of these forms, report both with
citations -- do not adjust a sign (Arian's standing rule).

WHAT THE OPERATOR IS. Adiabatic elimination of a far-detuned intermediate
manifold ([HAM] S9.5, [2gamma] S3.1, Cossel PhD thesis (Colorado, 2014) Eqs.
(6.29), (6.34), pp. 213, 219) gives, between X-state levels,

    T_eff = sum_i (d.eps2*) |i><i| (d.eps1) / Delta_i .

Two rank-1 operators couple to K = 0, 1, 2 and nothing else (B&C Eq. (5.141),
PDF p.198 / book p.166). When the detuning is common, B&C Eq. (5.142) read
BACKWARDS says the intermediate sum IS the reduced element of a single rank-K
operator, so the whole process is ONE SCALAR PER (K, dOmega) CHANNEL times
parameter-free geometry ([HAM] S9.5.1(1)). That is what this module computes:
the geometry. The scalars are the placeholder alpha_K*_dOm* Params, and the
caller supplies them.

WHICH CHANNELS ARE REGISTERED: K in {0, 2} ONLY.
  * (K = 0, dOmega = 0), (K = 2, dOmega = 0), (K = 2, dOmega = +-2). Within
    X 3Delta1 both |Omega| are 1, so q = Om' - Om is 0 or +-2, and the
    molecule-frame 3j (J' K J; -Om' q Om) vanishes unless |q| <= K -- hence
    |dOmega| = 2 at K = 2 only ([HAM] S9.5.2).
  * K = 1 is NOT registered and there is no alpha_K1_* parameter. [HAM] S9.5.3
    (OPEN-21) rules that K = 1 is the antisymmetric part of the dyad, so in
    exact closure it is P_X [d_a, d_b] P_X / 2 = 0 identically -- the Cartesian
    components of d commute. It reappears at O(delta/Delta) in the resolved
    form, and at O(1) for a restricted intermediate manifold (a single
    intermediate electronic state -- the actual ThF+/HfF+ situation); either
    would need its own placeholder alpha with that order recorded, and neither
    is in v2's registered set.
  `dyad_weights` still returns the K = 1 weights: they are a property of the
  polarisation pair, not of the molecule (and the sin^2 theta law Cossel
  measured lives in them, [2gamma] S3.0). No operator in this module consumes
  them.

REGISTRY_2G IS A THIRD REGISTRY ([SPEC-v2] S3.2) -- not heff.terms.REGISTRY,
not heff.elements_c2.REGISTRY_C2 -- so a transition operator can never be
summed into an assembled Hamiltonian. The entries exist for their Rules, their
cite and gate A5; heff.assemble never sees them.

WHICH 6j THE RANK-K REDUCTION USES, AND WHERE IT LIVES. [HAM] S9.5.1(2) is
S9.1's two-spectator chain with k -> K: the reduction is the two spectator 6j's
{F1 F I_F; F' F1' K} and {J F1 I_Th; F1' J' K}, and BOTH are evaluated inside
heff.elements_c2.axial_geometry -- this module writes no recoupling algebra of
its own beyond the polarisation dyad. So the 6j patch site for a monkeypatched
FAIL demonstration (Task 8) is `heff.elements_c2.w6j` -- there is no
`heff.twophoton.w6j` to patch. The B&C (5.142) 6j {1 1 K; j' j j''} that a
RESOLVED intermediate sum would carry is exactly what closure removes; it
appears nowhere in this file, which is the whole content of the closure form.
w3j is used here directly, for the Clebsch-Gordan coefficients of the dyad.

UNITS. The geometry is dimensionless. conventions.two_photon_norm =
'bc_5p142_reduced': alpha^K_{dOmega} multiplies the geometry of B&C (5.142)
at <eta'||alpha^K||eta> == 1 per (K, dOmega) channel -- the closure relation
of [HAM] S9.5.1(1), sum_{j''} {1 1 K; j' j j''} <j||d||j''><j''||d||j'> =
(-1)^(K+j+j') (2K+1)^-1/2 <j||T^K(d,d)||j'>, alpha^K = T^K(d,d)/Delta, living
on the resolved side -- so an amplitude comes out in units of alpha and a
strength in units of alpha^2. The alpha Params carry (MHz/(V/cm))^2/MHz and
the two field amplitudes are the caller's, exactly as
heff.spectra.line_strengths leaves d_mf to the caller.

VALIDITY, in the sentence [HAM] S9.5.5 writes for the notebook: the closure form
requires a detuning large compared with the intermediate ROTATIONAL structure,
Delta >> B_i ~ 7 GHz for ThF+, whereas the JILA experiments run at 0.16-1.5 GHz
-- the right operator shape, the wrong limit for the current experiment.
"""
import numpy as np

from .elements_c2 import axial_geometry
from .spectra import _strengths_from_matrices
from .terms import Rules, term
from .wigner import w3j

REGISTRY_2G = {}

_SQ2 = np.sqrt(2.0)


def _spherical(v):
    """{p: v_p} for any complex Cartesian vector, [HAM] S9.5.1(3):

        v_{+1} = -(v_x + i v_y)/sqrt(2),  v_0 = v_z,  v_{-1} = +(v_x - i v_y)/sqrt(2)

    the Condon-Shortley phase fixed by [HAM] S2.
    """
    vx, vy, vz = (complex(c) for c in np.asarray(v).reshape(3))
    return {1: -(vx + 1j * vy) / _SQ2, 0: vz, -1: (vx - 1j * vy) / _SQ2}


def _leg(eps):
    """c[p] = (-1)^p eps_{-p}, the coefficients of one E1 leg. [HAM] S9.5.1(3).

    d.eps = sum_p (-1)^p eps_{-p} T^1_p(d), from A.B = sum_p (-1)^p A_p B_{-p}
    (a bilinear identity, valid for complex components).

    Pass the CONJUGATED Jones vector for the bra-side leg: `_leg(conj(eps2))`
    means conjugate the Cartesian vector first, then take spherical components.
    The anchor is the sigma+ Jones vector e_{+1} = -(x + iy)/sqrt(2), for which
    this returns c[+1] = +1 and nothing else -- absorbing a sigma+ photon raises
    m_F by one.
    """
    s = _spherical(eps)
    return {p: (-1.0) ** p * s[-p] for p in (-1, 0, 1)}


def _cg(j1, m1, j2, m2, J, M):
    """<j1 m1 j2 m2|J M> = (-1)^(j1-j2+M) sqrt(2J+1) (j1 j2 J; m1 m2 -M).

    B&C Eq. (5.141), PDF p.198 / book p.166, is exactly this coupling at
    k1 = k2 = 1; [HAM] S9.5.1(3) inverts it to get the dyad weights.
    """
    return (-1.0) ** (j1 - j2 + M) * np.sqrt(2 * J + 1) * w3j(j1, j2, J, m1, m2, -M)


def dyad_weights(eps1, eps2):
    """The coupled polarisation dyad {(K, P): complex} for K = 0, 1, 2.

    [HAM] S9.5.1(3), implemented DIRECTLY in the form printed there:

        w^K_P = sum_{p_a + p_b = P} <1 p_a 1 p_b|K P> c_a[p_a] c_b[p_b]

    with `eps_a = eps2*` in slot 1 (the bra-side leg) and `eps_b = eps1` in
    slot 2 (the ket-side leg), because T_eff = (d.eps2*)|i><i|(d.eps1)/Delta and
    B&C (5.142) puts the bra-side leg first. The contraction is then free of
    further phases: T_eff = sum_K sum_P w^K_P alpha^K_P, with Delta m_F = P.

    NEITHER THE CONJUGATION NOR THE SLOT ORDER IS COSMETIC. The equivalent
    dyad-notation form is w^K_P = (-1)^(K+P) (eps2* x eps1)^K_{-P}; the
    plausible-looking (-1)^P (eps1 x eps2)^K_{-P} is a different function of
    the two Jones vectors ([HAM] S9.5.1(3) measures the difference at 1.3e+01
    on random eps) and it swaps which physical polarisation pair reaches
    Delta m_F = +-2. That form is given in S9.5 only for recognition and is not
    what this function computes.

    THE READING IS RAMAN (photon 2 emitted, hence conjugated). Under it,
    eps1 = eps2 = sigma+ gives legs (p_a, p_b) = (-1, +1) and so Delta m_F = 0,
    while eps1 = sigma+ with eps2 = sigma- gives (+1, +1) and Delta m_F = +2 --
    conjugating eps2 flips the sign of its helicity label, so the beam-pair ->
    Delta m_F map is INVERTED relative to the ladder reading in which both
    photons are absorbed. The reachable SET under sigma+- only is
    {0, +-2} either way ([HAM] S9.5.4). A caller labelling a polarisation pair
    must say which reading it means; this one is Raman.

    K = 0 is the scalar eps2*.eps1; K = 1 is the antisymmetric part
    (prop. eps2* x eps1), non-zero only for non-parallel or elliptical
    polarisations; K = 2 is the symmetric traceless part. The K = 1 weights are
    returned because they are a property of the polarisation pair, but NO
    registered operator consumes them -- there is no K = 1 channel
    ([HAM] S9.5.3, OPEN-21).

    eps1, eps2 are Cartesian Jones 3-vectors, complex allowed, in the frame
    whose z is the quantisation axis.
    """
    c_a = _leg(np.conj(np.asarray(eps2, dtype=complex)))
    c_b = _leg(np.asarray(eps1, dtype=complex))
    return {(K, P): sum(_cg(1, p_a, 1, P - p_a, K, P) * c_a[p_a] * c_b[P - p_a]
                        for p_a in (-1, 0, 1) if abs(P - p_a) <= 1)
            for K in (0, 1, 2) for P in range(-K, K + 1)}


def two_photon_geometry(bra, ket, ctx, *, K, P):
    """<bra| alpha^K_P |ket> in the two-spin basis, in units of <eta'||alpha^K||eta>.

    [HAM] S9.5.1(2), the spectator reduction: alpha^K acts on the electronic-
    rotational part only, so S9.1's two-spectator chain applies UNCHANGED with
    k -> K and p -> P. That chain is heff.elements_c2.axial_geometry, and this
    function is that call at k = K, q = Om_bra - Om_ket, p = P -- no new
    recoupling algebra, which is the point of the closure form.

    Returns 0.0 when |q| > K: the molecule-frame 3j (J' K J; -Om' q Om) vanishes
    identically there, so the channel does not exist ([HAM] S9.5.2). This is
    the ONE place the |q| <= K closure is applied; `axial_geometry` itself
    RAISES on |q| > k, since for it such a call is a caller bug.

    Validity ([HAM] S9.5.1(2)): the spectator reduction is exact only when the
    intermediate hyperfine structure is unresolved -- if Delta_i depends on the
    intermediate F', the F' sum cannot be factored out and the reduction fails.
    """
    q = float(bra["Om"]) - float(ket["Om"])
    if abs(q) > float(K):
        return 0.0
    return axial_geometry(bra, ket, ctx, k=K, q=q, p=P)


def two_photon_matrix(kets_a, kets_b, ctx, *, K, dOmega, P):
    """<a| alpha^K_P |b> over one channel, shape (len(a), len(b)).

    Mirrors heff.spectra.dipole_matrix: dimensionless geometry, indexed purely
    by POSITION in the arrays the caller passes.

    `dOmega` selects the channel and is 0 or 2, where 2 means |dOmega| = 2 --
    BOTH q = +2 and q = -2. They are one channel, not two: they share one
    scalar (alpha_K2_dOm2, exactly as elements_c2's eQq2_Th covers both signs),
    and parity maps each into the other, so separating them would make
    [T, parity] = 0 false by construction (gate V24).
    """
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


# --------------------------------------------------------- the registry

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
    "'bc_5p142_reduced': unit one-photon reduced elements, so a strength comes "
    "out in units of alpha^2. K = 1 is NOT registered ([HAM] S9.5.3, OPEN-21: "
    "identically zero in exact closure). TRANSITION operator in REGISTRY_2G, "
    "never summable into a Hamiltonian ([SPEC-v2] S3.2). Validity "
    "[HAM] S9.5.5: needs Delta >> the intermediate B_i ~ 7 GHz, which the "
    "0.16-1.5 GHz JILA detunings do not satisfy."
)

_K2 = (-2, -1, 0, 1, 2)


def _channel(K, dOmega):
    """The registered element for one (K, |dOmega|) channel.

    A registered term is called as fn(bra, ket, ctx) (heff.terms), so the lab
    component is the one the kets themselves force, P = Delta m_F -- the sum
    over P at fixed weight. That is the right object for gate A5, which asks
    whether the FORMULA is non-zero inside its declared Rules and zero outside;
    a physical amplitude weights each P by the dyad, which is
    two_photon_line_strengths' job.
    """
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
    # hermitian=False: fn evaluates at P = Delta m_F, and by reciprocity
    # (M_P(a,b) = (-1)^P M_{-P}(b,a), test_rank_K_sum_rule_and_reciprocity)
    # that object is antisymmetric at odd Delta m_F -- measured max|M-Mt| =
    # 0.632 on the 229ThF+ J_max=2 basis (fix round 1, finding 1).
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
    # hermitian=False: same reciprocity argument as K2_dOm0 -- measured
    # max|M-Mt| = 0.741 on the same basis (fix round 1, finding 1).
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
    heff.spectra._strengths_from_matrices is reused unchanged, so the ordering
    is encoded in exactly one place in the package (gate B8 there, gate V28
    here). This is not a formality: Cossel's Fig. 6.18 (thesis p.224) is a
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
