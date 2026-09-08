"""Amide basis enumeration: sequential nuclear-spin coupling with NH2 exchange."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from heff.amide_basis import AMIDE_KET, AmideBasisSpec, enumerate_amide_kets
from heff.spec import block_by_mF


def _spec(*, I_M=0.0, vibronic_sign=1, M="all"):
    """Small mixed-K fixture: K=0 at N=0 and K=+-1 at N=1."""
    return AmideBasisSpec(
        S=0.5,
        I_N=1.0,
        i_H=0.5,
        I_M=I_M,
        N_range=(0, 1),
        K_values=(0, 1),
        vibronic_sign=vibronic_sign,
        M=M,
    )


def _coupled_values(a, b):
    """Hand-counting helper using doubled angular momenta, independent of the API."""
    two_a, two_b = int(round(2 * a)), int(round(2 * b))
    return [two / 2.0 for two in range(abs(two_a - two_b), two_a + two_b + 1, 2)]


def _legacy_rows_without_metal(spec):
    """Independent seven-column legacy construction for the I_M=0 reduction."""
    rows = []
    for K in (-1, 0, 1):
        for N in range(spec.N_range[0], spec.N_range[1] + 1):
            if N < abs(K):
                continue
            for I_T in _coupled_values(spec.i_H, spec.i_H):
                if (spec.vibronic_sign * (-1) ** K * (-1) ** int(2 * spec.i_H - I_T)
                        != (-1) ** int(2 * spec.i_H)):
                    continue
                for J in _coupled_values(N, spec.S):
                    for F_N in _coupled_values(J, spec.I_N):
                        for F in _coupled_values(F_N, I_T):
                            for two_mF in range(-int(2 * F), int(2 * F) + 1, 2):
                                rows.append((N, K, J, F_N, I_T, F, two_mF / 2.0))
    return rows


def _rows(kets, names):
    return sorted(tuple(float(row[name]) for name in names) for row in kets)


def test_enumerator_returns_named_float64_full_m_basis():
    """Catches a positional dtype or an m-free result that blocks cannot partition."""
    kets = enumerate_amide_kets(_spec())
    assert kets.dtype == AMIDE_KET
    assert AMIDE_KET.names == ("N", "K", "J", "F_N", "I_T", "F_core", "F", "mF")
    assert all(AMIDE_KET[name] == np.dtype("f8") for name in AMIDE_KET.names)
    assert len(kets) > 0


def test_spec_is_immutable_after_normalizing_user_sequences():
    """Catches a mutable basis definition whose inputs can change after validation."""
    spec = AmideBasisSpec(0.5, 1.0, 0.5, 0.0, [0, 1], [0, -1], 1)
    assert spec.N_range == (0, 1) and spec.K_values == (0, 1)
    with pytest.raises(FrozenInstanceError):
        spec.S = 1.0


def test_vibronic_exchange_filter_selects_opposite_singlet_triplet_k_parities():
    """Catches a missing vibronic or rotational exchange-sign factor."""
    plus, minus = enumerate_amide_kets(_spec(vibronic_sign=1)), enumerate_amide_kets(
        _spec(vibronic_sign=-1))
    for kets, expected in ((plus, {0.0: {0.0}, -1.0: {1.0}, 1.0: {1.0}}),
                            (minus, {0.0: {1.0}, -1.0: {0.0}, 1.0: {0.0}})):
        got = {K: set(kets["I_T"][kets["K"] == K]) for K in expected}
        assert got == expected


@pytest.mark.parametrize("I_M", (0.0, 1.0, 2.0))
def test_outer_metal_spin_has_the_tensor_product_dimension(I_M):
    """Catches skipping an outer F branch or adding it at the wrong point in the chain."""
    base = len(enumerate_amide_kets(_spec(I_M=0.0)))
    assert len(enumerate_amide_kets(_spec(I_M=I_M))) == base * int(2 * I_M + 1)


def test_zero_metal_spin_reduces_to_the_hand_constructed_legacy_basis():
    """Catches an I_M=0 branch that is only dimensionally, not ket-wise, equivalent."""
    spec = _spec(I_M=0.0)
    got = enumerate_amide_kets(spec)
    projected = _rows(got, ("N", "K", "J", "F_N", "I_T", "F_core", "mF"))
    assert np.array_equal(got["F_core"], got["F"])
    assert projected == sorted(_legacy_rows_without_metal(spec))


def test_all_and_blocks_modes_enumerate_the_same_full_basis():
    """Catches M='blocks' silently omitting mF components instead of leaving partitioning to callers."""
    assert np.array_equal(enumerate_amide_kets(_spec(M="all")),
                          enumerate_amide_kets(_spec(M="blocks")))


def test_blocks_mode_output_is_accepted_by_the_existing_mf_partitioner():
    """Catches a field-name mismatch that leaves M='blocks' unusable by the shared blocker."""
    kets = enumerate_amide_kets(_spec(M="blocks"))
    blocking = block_by_mF(kets)
    merged = np.concatenate([blocking.index[label] for label in blocking.labels])
    assert np.array_equal(np.sort(merged), np.arange(len(kets)))


def test_every_row_satisfies_the_coupling_triangles_and_is_unique():
    """Catches off-by-one angular-momentum endpoints and duplicate +/-K expansion."""
    kets = enumerate_amide_kets(_spec(I_M=1.0))
    assert np.all(kets["J"] >= np.abs(kets["N"] - 0.5))
    assert np.all(kets["J"] <= kets["N"] + 0.5)
    assert np.all(kets["F_N"] >= np.abs(kets["J"] - 1.0))
    assert np.all(kets["F_N"] <= kets["J"] + 1.0)
    assert np.all(kets["F_core"] >= np.abs(kets["F_N"] - kets["I_T"]))
    assert np.all(kets["F_core"] <= kets["F_N"] + kets["I_T"])
    assert np.all(kets["F"] >= np.abs(kets["F_core"] - 1.0))
    assert np.all(kets["F"] <= kets["F_core"] + 1.0)
    assert np.all(np.abs(kets["mF"]) <= kets["F"])
    assert np.all(kets["N"] >= np.abs(kets["K"]))
    assert len({tuple(row) for row in kets}) == len(kets)


@pytest.mark.parametrize(
    "kwargs",
    (
        {"S": 0.25},
        {"I_N": -0.5},
        {"i_H": float("nan")},
        {"I_M": float("inf")},
        {"N_range": (0, 1.5)},
        {"N_range": (2, 1)},
        {"K_values": (0, 0.5)},
        {"vibronic_sign": 0},
        {"M": "none"},
    ),
)
def test_spec_rejects_nonphysical_or_unsupported_quantum_number_inputs(kwargs):
    """Catches fractional/nonfinite spins and unsupported truncation or M inputs early."""
    values = dict(S=0.5, I_N=1.0, i_H=0.5, I_M=0.0, N_range=(0, 1),
                  K_values=(0, 1), vibronic_sign=1)
    values.update(kwargs)
    with pytest.raises(ValueError):
        AmideBasisSpec(**values)
