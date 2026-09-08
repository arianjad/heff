"""Adapter for heff's native Hund's case-(c) implementation."""
from collections.abc import Mapping

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


def _one_electronic_record(electronic):
    if isinstance(electronic, Mapping):
        return electronic
    records = tuple(electronic)
    if len(records) != 1:
        raise ValueError("case-c backend requires exactly one electronic record")
    return records[0]


def make_spec(basis, electronic, spins):
    """Translate v1 TOML-shaped records into the existing case-(c) StateSpec."""
    unknown_basis_keys = set(basis) - _BASIS_KEYS
    if unknown_basis_keys:
        raise ValueError(f"unknown case-c basis key {sorted(unknown_basis_keys)[0]!r}")
    missing_basis_keys = _BASIS_KEYS - set(basis)
    if missing_basis_keys:
        raise ValueError(f"missing case-c basis key {sorted(missing_basis_keys)[0]!r}")
    if basis["J_min"] > basis["J_max"]:
        raise ValueError("J_min must be less than or equal to J_max")

    electronic_record = _one_electronic_record(electronic)
    if not isinstance(electronic_record, Mapping):
        raise ValueError("case-c backend requires exactly one electronic record")
    elec_state = ElecState(label=electronic_record["label"],
                           Omega=electronic_record["Omega"],
                           S=electronic_record["S"],
                           Lam=electronic_record["Lambda"],
                           T0=electronic_record.get("T0", 0.0))

    spin_chain = tuple(Spin(label=record["label"], I=record["I"],
                            couple_to=record["couple_to"]) for record in spins)
    for index, spin in enumerate(spin_chain):
        expected_parent = "J" if index == 0 else f"F{index}"
        if spin.couple_to != expected_parent:
            raise ValueError(
                f"spins must form a chain coupled inner-first: spins[{index}] "
                f"({spin.label!r}) must have couple_to == {expected_parent!r}, "
                f"got {spin.couple_to!r}")

    # The native v1 representation carries exactly one nuclear spin through I.
    # Multiple spins use the existing explicit, inner-first StateSpec chain.
    state_spins = () if len(spin_chain) == 1 else spin_chain
    I = spin_chain[-1].I if spin_chain else 0.0
    return StateSpec(case="c", electronic=(elec_state,), I=I,
                     J_range=(basis["J_min"], basis["J_max"]), M=basis["M"],
                     frame=basis["frame"], spins=state_spins)


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
