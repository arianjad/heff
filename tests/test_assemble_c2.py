"""Task 6: assembly and sweeps in the v2 (two-nuclear-spin) basis.

Gates A5, B1, [HAM] V8/V9, and V23 (the J-truncation report), on 229ThF+ and
227ThF+ through `build_term_matrices(..., case='c2', registry=REGISTRY_C2)` --
the v1 gates in tests/test_terms.py and tests/test_assemble.py, transplanted.
No new physics: this file wires the already-built v2 registry (Tasks 2-5)
through the already-supported assemble/engine path.

Controller rulings carried here (2026-09-05):
  R3  V23 compares J_max = 1 vs J_max = 4 (not 2 vs 4 -- dropping the whole
      Delta-J = +-1 coupling is the failure mode), thresholds 229 > 1 MHz,
      232 < 10 kHz ([HAM] S2.5 puts the 19F second-order shift at 1.6-2.6 kHz).
  R7  V23 gains a 227ThF+ row (same > 1 MHz threshold; its A_par_Th
      placeholder makes the shift far larger). j_convergence is reported for
      229 AND 227.
  R14 conventions.a_par_th_sign does not exist; thf_v2's own
      `a_par_th_sign=` keyword selects the signed A_par_Th Param.
  Task 5 review minor 4: test_A5_holds_for_every_c2_term and the hermiticity
      gate run at J_max = 3 (a J_max = 2 basis has no Delta-J = 2 pair, so a
      spurious eQq0_Th/eQq2_Th element there would be invisible).
"""
import numpy as np
import pytest

from heff.assemble import build_term_matrices, hamiltonian
from heff.conventions import parity_operator
from heff.elements_c2 import REGISTRY_C2, j_convergence
from heff.engine import sweep
from heff.params import Param, thf_v2
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import check_selection_rules, ctx_from, terms_for_case


def test_A5_holds_for_every_c2_term():
    """Gate A5 on the whole REGISTRY_C2, J_max = 3, 229ThF+ (I_Th = 5/2): every
    term is non-zero somewhere inside its declared rules and exactly zero
    outside them.

    229's ctx is the one that exercises every term including the two Th
    quadrupole terms (I_Th = 0.5 for 227 makes them structurally, correctly,
    zero everywhere -- [HAM] S9.4.1's I < 1 guard -- which would fail this
    same check for the wrong reason). Uniquely catches a dead operator, e.g.
    a rules/formula sign disagreement in eQq2_Th ([HAM] S9.4.4).
    """
    spec = thf_spec("229", J_max=3)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    for t in terms_for_case("c2", registry=REGISTRY_C2):
        report = check_selection_rules(t, kets, ctx)
        assert report["n_nonzero_inside"] > 0, t.name
        assert report["n_nonzero_outside"] == 0, t.name


def test_B1_every_v2_block_is_hermitian():
    """Gate B1 at J_max = 3, for 229ThF+ and 227ThF+, at a representative
    (E_z, B_z) point with the PT-odd block turned on so its Hermiticity is
    exercised too."""
    for iso in ("229", "227"):
        spec = thf_spec(iso, J_max=3)
        kets = enumerate_kets(spec)
        ctx = ctx_from(spec, thf_v2(iso))
        pset = thf_v2(iso).with_(d_e=Param(1e-29, ""), k_TP=Param(1e-9, ""))
        for label, idx in block_by_mF(kets).index.items():
            tm = build_term_matrices(kets[idx], ctx, case="c2", registry=REGISTRY_C2)
            H = hamiltonian(tm, pset, {"E_z": 12.0, "B_z": 3.0})
            assert np.allclose(H, H.conj().T, atol=1e-12, rtol=0), f"{iso} block {label}"
            assert not np.iscomplexobj(H) or np.max(np.abs(H.imag)) == 0.0


def test_V8_parity_commutes_at_zero_field_in_the_v2_basis():
    """[HAM] V8: at E = B = 0 with the PT-odd block off, H commutes with
    P = sigma_xz R_y(pi), including with eQq2_Th on (it is parity-EVEN, same
    as omega_doubling -- [HAM] S9.4.4). [H, P] != 0 once stark_z is on.

    Uniquely catches a parity-odd term leaking into the field-free
    Hamiltonian, and specifically a sign error in eQq2_Th that would make it
    parity-odd.
    """
    for iso in ("229", "227"):
        spec = thf_spec(iso)
        kets = enumerate_kets(spec)
        ctx = ctx_from(spec, thf_v2(iso))
        idx = block_by_mF(kets).index[0.0]
        sub = kets[idx]
        tm = build_term_matrices(sub, ctx, case="c2", registry=REGISTRY_C2)
        pset = thf_v2(iso)  # eQq0_Th/eQq2_Th at their non-zero defaults for 229
        P = parity_operator(sub, S=ctx.S, ell=0.0, s=0.0)
        H0 = hamiltonian(tm, pset, {"E_z": 0.0, "B_z": 0.0})
        assert np.allclose(P @ P, np.eye(len(sub)))
        # v1's V8 (tests/test_assemble.py) also asserts P's own spectrum is an
        # even split, +1 and -1 in equal numbers -- the statement that the block
        # really does carry both parities. Carried over here: an m_F block of
        # the two-spin basis is built from +-Omega pairs just as v1's was, so a
        # pairing bug that halved the -1 eigenvalues would leave [H, P] = 0
        # while quietly making the parity label meaningless.
        ev = np.linalg.eigvalsh(P)
        assert sum(ev > 0) == sum(ev < 0) == len(sub) // 2, iso
        assert np.max(np.abs(H0 @ P - P @ H0)) < 1e-9, iso
        Hs = hamiltonian(tm, pset, {"E_z": 10.0, "B_z": 0.0})
        assert np.max(np.abs(Hs @ P - P @ Hs)) > 1.0, iso


def test_V9_kramers_degeneracy_at_zero_B_in_the_v2_basis():
    """[HAM] V9: at B = 0 (any E, PT-odd off), E(F, m_F, Om) = E(F, -m_F, -Om),
    so the +m_F and -m_F block spectra coincide. Non-vacuous: they differ once
    B != 0."""
    for iso in ("229", "227"):
        spec = thf_spec(iso)
        kets = enumerate_kets(spec)
        ctx = ctx_from(spec, thf_v2(iso))
        pset = thf_v2(iso)
        blocks = block_by_mF(kets)
        for m in (m for m in blocks.labels if m > 0):
            up = build_term_matrices(kets[blocks.index[m]], ctx, case="c2",
                                     registry=REGISTRY_C2)
            dn = build_term_matrices(kets[blocks.index[-m]], ctx, case="c2",
                                     registry=REGISTRY_C2)
            wa = np.linalg.eigvalsh(hamiltonian(up, pset, {"E_z": 25.0, "B_z": 0.0}))
            wb = np.linalg.eigvalsh(hamiltonian(dn, pset, {"E_z": 25.0, "B_z": 0.0}))
            assert np.allclose(wa, wb, atol=1e-9), f"{iso} m_F = {m}"
        m = next(m for m in blocks.labels if m > 0)
        up = build_term_matrices(kets[blocks.index[m]], ctx, case="c2", registry=REGISTRY_C2)
        dn = build_term_matrices(kets[blocks.index[-m]], ctx, case="c2", registry=REGISTRY_C2)
        wa = np.linalg.eigvalsh(hamiltonian(up, pset, {"E_z": 25.0, "B_z": 1.0}))
        wb = np.linalg.eigvalsh(hamiltonian(dn, pset, {"E_z": 25.0, "B_z": 1.0}))
        assert np.max(np.abs(wa - wb)) > 1e-4, iso


def test_J_truncation_is_reported_and_the_report_is_non_vacuous():
    """V23. [HAM] S2.5 / [TH] S4.2: dropping Delta-J = +-1 Th hyperfine by
    truncating at J_max = 1 shifts the J = 1 levels by tens to thousands of
    MHz for 229/227ThF+ (real Th nuclear spin) but only ~kHz for 232ThF+ (no
    Th spin, 19F-only second-order shift) -- five decades of separation.

    Both outcomes are reachable in one test: 232ThF+ is the negative control,
    so this gate cannot pass by always reporting "large". n_levels is the
    full dimension of the J = 1, m_F = 0 (or 0.5 for 232) sub-block, so the
    comparison is exactly the lowest J = 1 levels (232ThF+ has no J = 0
    state either, so its lowest states are also its J = 1 manifold).
    """
    dims_J1 = {}
    for iso, mF in (("229", 0.0), ("227", 0.0), ("232", 0.5)):
        spec = thf_spec(iso, J_max=1)
        kets = enumerate_kets(spec)
        dims_J1[iso] = len(block_by_mF(kets).index[mF])

    shift_229 = j_convergence("229", J_maxes=(1, 4), mF=0.0,
                              n_levels=dims_J1["229"])["by_J_max"][1]["max_shift_MHz"]
    shift_227 = j_convergence("227", J_maxes=(1, 4), mF=0.0,
                              n_levels=dims_J1["227"])["by_J_max"][1]["max_shift_MHz"]

    # 232ThF+ has spins=() (thf_spec returns the v1 object): j_convergence
    # falls back to the v1 registry there, per its own docstring.
    shift_232 = j_convergence("232", J_maxes=(1, 4), mF=0.5,
                              n_levels=dims_J1["232"])["by_J_max"][1]["max_shift_MHz"]

    assert shift_229 > 1.0, shift_229          # positive control: MHz-scale
    assert shift_227 > 1.0, shift_227          # R7: also MHz-scale (larger)
    assert shift_232 < 1e-2, shift_232         # negative control: < 10 kHz = 1e-2 MHz


def test_j_convergence_reports_229_and_227():
    """R7: the reporting helper itself, at its own defaults, for both
    isotopologues that carry Th hyperfine. Reports; never raises on a
    magnitude."""
    for iso in ("229", "227"):
        report = j_convergence(iso)
        assert report["isotopologue"] == iso
        assert set(report["by_J_max"]) == {2, 4, 6}
        assert report["field_free"] is True
        for J_max, row in report["by_J_max"].items():
            assert row["dimension"] > 0
            assert row["max_shift_MHz"] >= 0.0
        # the largest J_max is its own reference: zero self-shift
        assert report["by_J_max"][report["reference_J_max"]]["max_shift_MHz"] == 0.0


def test_j_convergence_default_mF_differs_by_isotopologue():
    """229/227ThF+ have integer F (inner spin half-integer + 19F 1/2 ->
    F1 half-integer -> F integer), so m_F = 0 exists; 232ThF+'s v1 F = J +-
    1/2 is half-integer, so m_F = 0.5 is the smallest that exists there. No
    single default is a level of all three (checked numerically), so the
    default is per-isotopologue, not a shared literal."""
    assert j_convergence("229")["mF"] == 0.0
    assert j_convergence("227")["mF"] == 0.0
    assert j_convergence("232")["mF"] == 0.5


def test_j_convergence_rejects_an_mF_the_block_does_not_hold():
    with pytest.raises(ValueError, match="mF"):
        j_convergence("229", mF=0.25)


def test_a_sweep_in_the_v2_basis_returns_eigenvectors_and_the_active_term_list():
    """engine.sweep unchanged: returns eigenvectors, and SweepResult.manifest
    carries the v2 conventions stamp (version = 'thf-v2')."""
    spec = thf_spec("229")
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    idx = block_by_mF(kets).index[0.0]
    tm = build_term_matrices(kets[idx], ctx, case="c2", registry=REGISTRY_C2)
    pset = thf_v2("229")
    E = np.linspace(0.0, 30.0, 4)
    res = sweep(tm, pset, {"E_z": E, "B_z": np.zeros(4)})
    assert res.evals.shape == (4, len(idx))
    assert res.evecs.shape == (4, len(idx), len(idx))
    assert "rotation" in res.active_terms and "stark_z" in res.active_terms
    assert res.manifest["conventions"]["version"] == "thf-v2"
    assert res.manifest["case"] == "c2"
