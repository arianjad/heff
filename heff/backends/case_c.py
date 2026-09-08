"""Adapter for heff's native Hund's case-(c) implementation."""
from collections.abc import Mapping
from operator import index as integer_index

from .. import elements_c as _elements_c  # populate the native term registry
from ..backend_registry import Backend
from ..conventions import Conventions
from ..spec import ElecState, Spin, StateSpec
from ..terms import REGISTRY, terms_for_case


CANONICAL_UNITS = {
    "B0": "MHz", "D0": "MHz", "omega_ef": "MHz", "A_par": "MHz",
    "c_I": "MHz", "d_mf": "MHz/(V/cm)", "G_par": "", "g_N": "",
    "E_eff": "MHz/(e cm)", "W_TP": "MHz", "d_e": "", "k_TP": "",
}

_BASIS_KEYS = frozenset({"J_min", "J_max", "M", "frame"})


class _InputValidationError(ValueError):
    """A backend input error with paths that the model layer can bind."""

    def __init__(self, message, *paths):
        super().__init__(message)
        self._heff_input_paths = tuple(paths)


def _required(record, key, path):
    if key not in record:
        raise _InputValidationError(
            f"missing case-c input key {key!r}", path + (key,))
    return record[key]


def _one_electronic_record(electronic):
    if isinstance(electronic, Mapping):
        return electronic
    try:
        records = tuple(electronic)
    except TypeError as exc:
        raise _InputValidationError(
            "case-c backend requires exactly one electronic record",
            ("manifold", "electronic")) from exc
    if len(records) != 1:
        raise _InputValidationError(
            "case-c backend requires exactly one electronic record",
            ("manifold", "electronic"))
    return records[0]


def make_spec(basis, electronic, spins):
    """Translate v1 TOML-shaped records into the existing case-(c) StateSpec."""
    unknown_basis_keys = set(basis) - _BASIS_KEYS
    if unknown_basis_keys:
        key = sorted(unknown_basis_keys)[0]
        raise _InputValidationError(
            f"unknown case-c basis key {key!r}",
            ("manifold", "basis", key))
    missing_basis_keys = _BASIS_KEYS - set(basis)
    if missing_basis_keys:
        key = sorted(missing_basis_keys)[0]
        raise _InputValidationError(
            f"missing case-c basis key {key!r}",
            ("manifold", "basis", key))
    range_paths = (("manifold", "basis", "J_min"),
                   ("manifold", "basis", "J_max"))
    for key in ("J_min", "J_max"):
        value = basis[key]
        try:
            if isinstance(value, bool):
                raise TypeError
            integer_index(value)
        except TypeError as exc:
            raise _InputValidationError(
                f"{key} must be a non-boolean integer usable as a range endpoint",
                ("manifold", "basis", key)) from exc
    reversed_range = basis["J_min"] > basis["J_max"]
    if reversed_range:
        raise _InputValidationError(
            "J_min must be less than or equal to J_max", *range_paths)

    electronic_record = _one_electronic_record(electronic)
    if not isinstance(electronic_record, Mapping):
        raise _InputValidationError(
            "case-c backend requires exactly one electronic record",
            ("manifold", "electronic"))
    electronic_path = ("manifold", "electronic")
    elec_state = ElecState(label=_required(electronic_record, "label", electronic_path),
                           Omega=_required(electronic_record, "Omega", electronic_path),
                           S=_required(electronic_record, "S", electronic_path),
                           Lam=_required(electronic_record, "Lambda", electronic_path),
                           T0=electronic_record.get("T0", 0.0))

    if len(spins) > 1:
        raise _InputValidationError(
            "high-level case-c foundation supports at most one spin",
            ("isotopologue", "spins", 1))

    spin_chain = []
    for index, record in enumerate(spins):
        spin_path = ("isotopologue", "spins", index)
        if not isinstance(record, Mapping):
            raise _InputValidationError(
                "case-c spins must be records", spin_path)
        spin = Spin(label=_required(record, "label", spin_path),
                    I=_required(record, "I", spin_path),
                    couple_to=_required(record, "couple_to", spin_path))
        expected_parent = "J" if index == 0 else f"F{index}"
        if spin.couple_to != expected_parent:
            raise _InputValidationError(
                f"spins must form a chain coupled inner-first: spins[{index}] "
                f"({spin.label!r}) must have couple_to == {expected_parent!r}, "
                f"got {spin.couple_to!r}",
                spin_path + ("couple_to",))
        spin_chain.append(spin)
    spin_chain = tuple(spin_chain)

    # The native v1 representation carries one nuclear spin through I.
    state_spins = ()
    I = spin_chain[-1].I if spin_chain else 0.0
    try:
        return StateSpec(case="c", electronic=(elec_state,), I=I,
                         J_range=(basis["J_min"], basis["J_max"]), M=basis["M"],
                         frame=basis["frame"], spins=state_spins)
    except ValueError as exc:
        # Range and spin-chain invariants were checked above; the remaining
        # StateSpec validation reachable here is the M mode.
        raise _InputValidationError(
            str(exc), ("manifold", "basis", "M")) from exc


CASE_C_BACKEND = Backend(
    id="case_c",
    case="c",
    registry=REGISTRY,
    default_terms=tuple(term.name for term in terms_for_case("c")),
    canonical_units=CANONICAL_UNITS,
    runtime_knobs=frozenset({"E_z", "B_z"}),
    optional_zero_parameters=frozenset({"d_e", "k_TP"}),
    make_spec=make_spec,
    default_conventions=Conventions,
)
