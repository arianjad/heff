"""Gates A1, A2, B1, and [HAM] V8 (parity) and V9 (Kramers).

A1 and A2 run on RANDOM matrices with random coefficients: they test the
assembler's indexing and the fast path's agreement with the reference path, and
they must not depend on any physics. V8 and V9 run on the real ThF+ blocks.
"""
import numpy as np
import pytest

from heff import elements_c  # noqa: F401  (registers the terms)
from heff.assemble import (TermMatrices, active, build_term_matrices, hamiltonian,
                           hamiltonian_batch, sweep_coefficients, vertex)
from heff.conventions import parity_operator
from heff.params import Param, ParamSet, thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import ctx_from


@pytest.fixture
def block():
    spec = thf_spec()
    kets = enumerate_kets(spec)
    idx = block_by_mF(kets).index[1.5]
    return kets[idx], ctx_from(spec, thf_v1())


def _random_tm(rng, n_terms=5, d=7):
    mats = []
    for _ in range(n_terms):
        A = rng.standard_normal((d, d))
        mats.append(A + A.T)
    names = tuple(f"t{k}" for k in range(n_terms))
    return TermMatrices(names=names,
                        params=tuple((f"p{k}",) for k in range(n_terms)),
                        mats=np.array(mats), kets=np.zeros(d), manifest={})


def _pset_for(tm, values):
    from heff.conventions import Conventions
    return ParamSet({p[0]: Param.of(v) for p, v in zip(tm.params, values)}, Conventions())


def test_A1_assembly_is_exactly_the_weighted_sum():
    """Gate A1: H(c) == sum_k c_k M_k on random M_k and random c.

    Uniquely catches assembler indexing/ordering drift -- a term matrix paired
    with the wrong coefficient, or a transposed stack. Generalises
    Molecule-Structure's GATE A (export_crossing_subspace.py:323).
    """
    rng = np.random.default_rng(20260905)
    tm = _random_tm(rng)
    vals = rng.standard_normal(len(tm.names))
    pset = _pset_for(tm, vals)
    ref = sum(v * M for v, M in zip(vals, tm.mats))
    assert np.allclose(hamiltonian(tm, pset, {}), ref, atol=1e-13)


def test_A1_fails_if_the_stack_is_permuted():
    """FAIL demo for A1: permute the matrices, keep the coefficients."""
    rng = np.random.default_rng(7)
    tm = _random_tm(rng)
    vals = rng.standard_normal(len(tm.names))
    pset = _pset_for(tm, vals)
    from dataclasses import replace
    scrambled = replace(tm, mats=tm.mats[::-1])
    ref = sum(v * M for v, M in zip(vals, tm.mats))
    assert not np.allclose(hamiltonian(scrambled, pset, {}), ref)


def test_A2_batched_resum_equals_the_per_point_loop():
    """Gate A2: the tensordot fast path == the reference loop, random matrices.

    Uniquely catches the sweep fast path diverging from the reference path
    (an axis swap in tensordot, a chunk boundary dropping a set).
    """
    rng = np.random.default_rng(99)
    tm = _random_tm(rng, n_terms=4, d=6)
    c = rng.standard_normal((37, 4))
    fast = hamiltonian_batch(tm, c)
    assert fast.shape == (37, 6, 6)
    for n in range(37):
        ref = sum(c[n, k] * tm.mats[k] for k in range(4))
        assert np.allclose(fast[n], ref, atol=1e-13)


def test_sweep_coefficients_broadcasts_knob_arrays(block):
    """c[n_sets, n_terms] broadcasts a field sweep without re-reading records."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    pset = thf_v1()
    E = np.linspace(0.0, 60.0, 5)
    c = sweep_coefficients(tm, pset, {"E_z": E, "B_z": np.zeros(5)})
    assert c.shape == (5, len(tm.names))
    k = tm.names.index("stark_z")
    assert np.allclose(c[:, k], pset.value("d_mf") * E)
    j = tm.names.index("rotation")
    assert np.allclose(c[:, j], pset.value("B0"))


def test_B1_every_assembled_block_is_hermitian_and_real(block):
    """Gate B1. Uniquely catches a transposed 3j/6j argument order or a bra/ket
    swap in a spectator-theorem phase ([HAM] V1)."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    pset = thf_v1().with_(d_e=Param(1e-29, ""), k_TP=Param(1e-9, ""))
    for label, idx in block_by_mF(kets).index.items():
        tm = build_term_matrices(kets[idx], ctx)
        H = hamiltonian(tm, pset, {"E_z": 12.0, "B_z": 3.0})
        assert np.allclose(H, H.conj().T, atol=1e-12, rtol=0), f"block {label}"
        assert not np.iscomplexobj(H) or np.max(np.abs(H.imag)) == 0.0


def test_B1_check_is_not_loose_against_diagonals_at_the_1e5_MHz_scale():
    """Companion FAIL demo for B1 (fix round 1, finding 2): with numpy's
    DEFAULT rtol=1e-5, a block whose diagonals run to ~1e5 MHz (J up to 4)
    gets an effective tolerance of ~1 MHz, which swallows a real antisymmetric
    perturbation of only 1e-6 MHz on the largest off-diagonal pair. B1's own
    check must use rtol=0 so it actually FAILS on that perturbation.
    """
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    idx = block_by_mF(kets).index[3.5]
    tm = build_term_matrices(kets[idx], ctx)
    H = hamiltonian(tm, thf_v1(), {"E_z": 12.0, "B_z": 3.0})
    assert np.allclose(H, H.conj().T, atol=1e-12, rtol=0)

    off = np.abs(np.triu(H, k=1))
    i, j = np.unravel_index(np.argmax(off), off.shape)
    assert off[i, j] > 0.0, "need a genuinely nonzero off-diagonal pair to perturb"
    H_pert = H.copy()
    H_pert[i, j] += 1e-6  # MHz, antisymmetric: only one side of the pair moves

    # The bug: default rtol=1e-5 is loose enough to hide this at ~1e5 MHz diagonals.
    assert np.allclose(H_pert, H_pert.conj().T, atol=1e-12)
    # The fix: rtol=0 makes B1's own check correctly reject the perturbed matrix.
    assert not np.allclose(H_pert, H_pert.conj().T, atol=1e-12, rtol=0)


def test_a_term_declared_hermitian_that_is_not_raises(block):
    """Structural error -> raise (spec S3.3). Not a hash, not a snapshot."""
    from heff.terms import Rules, term

    kets, ctx = block
    reg = {}

    @term(name="not_hermitian", param=("B0",), cases=("c",),
          rules=Rules(dJ=(-1, 0, 1), dOm=(0.0,), dF=(-1, 0, 1), dmF=(0,)),
          hermitian=True, real=True, cite="toy", registry=reg)
    def bad(bra, ket, ctx):
        return 1.0 if bra["J"] > ket["J"] else 0.0

    with pytest.raises(ValueError, match="hermitian"):
        build_term_matrices(kets, ctx, registry=reg)


def test_inactive_terms_are_reported_not_hidden(block):
    """A term absent from the parameter set assembles with coefficient zero,
    which is indistinguishable from 'deliberately off' unless it is reported
    (spec S3.3 (ii); the thesis's own H_K = 0 deperturbation study is a
    legitimate use, so this must be visible, not forbidden)."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    names = active(tm, thf_v1(), {"E_z": 0.0, "B_z": 0.0})
    assert "rotation" in names and "omega_doubling" in names
    assert "stark_z" not in names and "pt_odd_edm" not in names
    assert "zeeman_Gpar" not in names
    names_on = active(tm, thf_v1(), {"E_z": 10.0, "B_z": 1.0})
    assert "stark_z" in names_on and "zeeman_Gpar" in names_on


def test_vertex_is_the_exact_derivative_of_H(block):
    """dH/dB_z is exact because H is linear in every knob: it is the sum over
    the terms containing B_z of (product of the other knobs) x M_k, checked
    here against that exact identity (not a finite difference -- see
    task-7-report.md ruling 1: a difference quotient at 1e5 MHz diagonals sits
    below the noise floor a 1e-9 atol would require)."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    pset = thf_v1()
    knobs = {"E_z": 7.0, "B_z": 0.0}
    V = vertex(tm, pset, knobs, "B_z")
    k1, k2 = tm.names.index("zeeman_Gpar"), tm.names.index("zeeman_nuclear")
    assert np.allclose(V, pset.value("G_par") * tm.mats[k1]
                       + pset.value("g_N") * tm.mats[k2], atol=1e-13)


def test_V8_parity_commutes_at_zero_field_and_not_with_the_stark_term():
    """[HAM] V8 / gate B2: at E = 0 with the PT-odd block off, H commutes with
    P = sigma_xz R_y(pi), P^2 = 1, and the block splits into two equal halves.

    Uniquely catches a parity-odd term (Stark, eEDM) leaking into the field-free
    Hamiltonian. FAIL is reachable in the same test: with E != 0 the commutator
    is large, which is what makes the E = 0 result non-vacuous.
    """
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    idx = block_by_mF(kets).index[0.5]
    sub = kets[idx]
    tm = build_term_matrices(sub, ctx)
    P = parity_operator(sub, S=ctx.S, ell=0.0, s=0.0)
    H0 = hamiltonian(tm, thf_v1(), {"E_z": 0.0, "B_z": 1.0})
    assert np.allclose(P @ P, np.eye(len(sub)))
    assert np.max(np.abs(H0 @ P - P @ H0)) < 1e-9
    ev = np.linalg.eigvalsh(P)
    assert sum(ev > 0) == sum(ev < 0) == len(sub) // 2
    Hs = hamiltonian(tm, thf_v1(), {"E_z": 10.0, "B_z": 1.0})
    assert np.max(np.abs(Hs @ P - P @ Hs)) > 1.0
    He = hamiltonian(tm, thf_v1().with_(d_e=Param(1e-20, "")), {"E_z": 0.0, "B_z": 0.0})
    assert np.max(np.abs(He @ P - P @ He)) > 0.0


def test_V9_kramers_degeneracy_at_zero_B():
    """[HAM] V9: at B = 0 (any E, no eEDM), E(F, m_F, Om) = E(F, -m_F, -Om), so
    the +m_F and -m_F block spectra coincide (Leanhardt Eq. 34, exact).

    Uniquely catches a sign error in an m_F-odd term. This is the degeneracy the
    whole experiment measures against. Non-vacuous: at B != 0 the two spectra
    differ, asserted below.
    """
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    blocks = block_by_mF(kets)
    pset = thf_v1()
    for m in (0.5, 1.5, 2.5, 3.5, 4.5):
        up = build_term_matrices(kets[blocks.index[m]], ctx)
        dn = build_term_matrices(kets[blocks.index[-m]], ctx)
        wa = np.linalg.eigvalsh(hamiltonian(up, pset, {"E_z": 25.0, "B_z": 0.0}))
        wb = np.linalg.eigvalsh(hamiltonian(dn, pset, {"E_z": 25.0, "B_z": 0.0}))
        assert np.allclose(wa, wb, atol=1e-9), f"|m_F| = {m}"
    up = build_term_matrices(kets[blocks.index[1.5]], ctx)
    dn = build_term_matrices(kets[blocks.index[-1.5]], ctx)
    wa = np.linalg.eigvalsh(hamiltonian(up, pset, {"E_z": 25.0, "B_z": 1.0}))
    wb = np.linalg.eigvalsh(hamiltonian(dn, pset, {"E_z": 25.0, "B_z": 1.0}))
    # This block spans J = 1..4, so diagonals run up to ~1.45e5 MHz. np.allclose's
    # default rtol=1e-5 alone gives a tolerance of ~1.5 MHz there -- large enough
    # to swallow the real ~0.01-0.1 MHz Zeeman asymmetry and make this assertion
    # vacuously pass regardless of B. Compare the absolute difference directly
    # (task-7-report.md: discovered during GREEN, not a physics or sign issue).
    assert np.max(np.abs(wa - wb)) > 1e-4


def test_manifest_records_terms_conventions_and_citations(block):
    """The manifest is a record (dimension, cites, conventions, cache key), not
    a gate -- rebuilding twice yields the same key."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    m = tm.manifest
    assert m["dimension"] == len(kets)
    assert set(m["cites"]) == set(tm.names)
    assert m["conventions"]["n_hat"] == "F_to_Th"
    assert m["wigner_backend"] == "sympy+lru_cache"
    assert isinstance(m["wigner_version"], str) and m["wigner_version"]
    assert isinstance(m["spec_hash"], str) and len(m["spec_hash"]) == 64  # sha256 hex
    # the cache key is a record, designed now, gating nothing
    import json
    key = json.loads(m["key"])
    assert key["dimension"] == len(kets) and "rotation" in key["terms"]
    rebuilt = build_term_matrices(kets, ctx).manifest
    assert rebuilt["key"] == m["key"]
    assert rebuilt["spec_hash"] == m["spec_hash"]


def test_cache_key_is_block_aware_not_just_case_dimension_terms_conventions():
    """Fix round 1, finding 1: the +m_F = 3/2 and -m_F = 3/2 blocks share case,
    dimension, term names and conventions, so the OLD key collided even though
    the two blocks' matrices differ (opposite-sign m_F-odd terms). The key
    must fold in a canonicalised description of the block's own kets."""
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    blocks = block_by_mF(kets)
    up = build_term_matrices(kets[blocks.index[1.5]], ctx)
    dn = build_term_matrices(kets[blocks.index[-1.5]], ctx)
    # same case/dimension/terms/conventions -- the OLD key's entire input
    assert up.manifest["dimension"] == dn.manifest["dimension"]
    assert up.manifest["terms"] == dn.manifest["terms"]
    assert up.manifest["conventions"] == dn.manifest["conventions"]
    # ... but the matrices differ (an m_F-odd term, e.g. Zeeman, flips sign)
    assert not np.allclose(up.mats[up.names.index("zeeman_Gpar")],
                           dn.mats[dn.names.index("zeeman_Gpar")])
    # so the key (and the dedicated spec_hash) must differ too
    assert up.manifest["key"] != dn.manifest["key"]
    assert up.manifest["spec_hash"] != dn.manifest["spec_hash"]


def test_unknown_knob_symbol_raises_naming_it_and_the_known_ones(block):
    """Fix round 1, finding 3: a typo'd knob symbol used to be silently
    dropped (`pset.value` / `_knob` fall through to a default), so
    `sweep_coefficients(tm, pset, {"NOPE": ...})` returned a valid-looking
    array with no warning. It must raise, naming the bad symbol and listing
    the known ones, in every entry point that accepts a knobs mapping."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    pset = thf_v1()

    with pytest.raises(ValueError, match="NOPE") as exc:
        sweep_coefficients(tm, pset, {"E_z": np.linspace(0.0, 60.0, 4),
                                      "NOPE": np.zeros(4)})
    assert "B0" in str(exc.value) or "E_z" in str(exc.value)  # known symbols listed

    with pytest.raises(ValueError, match="NOPE"):
        hamiltonian(tm, pset, {"E_z": 10.0, "NOPE": 1.0})

    with pytest.raises(ValueError, match="NOPE"):
        active(tm, pset, {"NOPE": 1.0})

    with pytest.raises(ValueError, match="NOPE"):
        vertex(tm, pset, {"NOPE": 1.0}, "B_z")

    with pytest.raises(ValueError, match="unknown knob symbol"):
        vertex(tm, pset, {}, "NOPE")


def test_sweep_coefficients_wrong_length_array_raises_a_clear_error(block):
    """Fix round 1, finding 3: two knob arrays with incompatible (non-
    broadcastable) lengths must raise a clear error naming the symbols and
    shapes involved, not an opaque bare numpy broadcast message."""
    kets, ctx = block
    tm = build_term_matrices(kets, ctx)
    pset = thf_v1()
    with pytest.raises(ValueError) as exc:
        sweep_coefficients(tm, pset, {"E_z": np.zeros(5), "B_z": np.zeros(3)})
    msg = str(exc.value)
    assert "E_z" in msg and "B_z" in msg
    assert "(5,)" in msg and "(3,)" in msg


def test_complex_term_promotes_H_only_when_its_coefficient_is_nonzero():
    """Fix round 1, finding 4: build_term_matrices used to compute one
    np.result_type across ALL terms and cast the whole stack, so a single
    complex term promoted every assembled H regardless of whether that term
    was even active. Each term now keeps its own dtype at build; hamiltonian
    (and hamiltonian_batch) promote only when an ACTIVE term is complex."""
    from heff.terms import Rules, term

    reg = {}
    d = 3

    @term(name="real_diag", param=("A",), cases=("c",), rules=Rules(),
          hermitian=True, real=True, cite="toy", registry=reg)
    def _real(bra, ket, ctx):
        return 0.0

    @term(name="fake_complex_hermitian", param=("C",), cases=("c",), rules=Rules(),
          hermitian=True, real=False, cite="toy", registry=reg)
    def _cplx(bra, ket, ctx):
        return 0.0

    # Build TermMatrices by hand (spec/rules plumbing isn't needed for this
    # dtype check): a real diagonal term and a complex Hermitian term.
    from heff.spec import KET_C
    kets_arr = np.zeros(d, dtype=KET_C)
    ctx = ctx_from(thf_spec(), thf_v1())
    tm = build_term_matrices(kets_arr, ctx, registry=reg)
    assert tm.mats[tm.names.index("real_diag")].dtype == np.float64
    assert tm.mats[tm.names.index("fake_complex_hermitian")].dtype == np.complex128

    pset = ParamSet({"A": Param.of(2.0), "C": Param.of(0.0)}, thf_v1().conventions)
    H_off = hamiltonian(tm, pset, {})
    assert not np.iscomplexobj(H_off)
    assert H_off.dtype == np.float64

    pset_on = pset.with_(C=Param.of(3.0))
    H_on = hamiltonian(tm, pset_on, {})
    assert np.iscomplexobj(H_on)
    assert H_on.dtype == np.complex128

    # Column order follows tm.names (sorted term names), not declaration order.
    i_A, i_C = tm.names.index("real_diag"), tm.names.index("fake_complex_hermitian")
    c_batch = np.zeros((2, len(tm.names)))
    c_batch[:, i_A] = 2.0  # A active in both sets, C inactive in both
    Hb_off = hamiltonian_batch(tm, c_batch)
    assert Hb_off.dtype == np.float64

    c_batch_on = c_batch.copy()
    c_batch_on[1, i_C] = 3.0  # second set turns C on
    Hb_on = hamiltonian_batch(tm, c_batch_on)
    assert Hb_on.dtype == np.complex128
