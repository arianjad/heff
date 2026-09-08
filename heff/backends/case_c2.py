"""Adapter for ThF+ in the native two-spin Hund's case-(c) basis.

The input spin order is physical and fixed: ``spins[0]`` is the inner thorium
spin coupled to ``J`` and ``spins[1]`` is the outer fluorine spin coupled to
``F1``.  The native C2 kernels use those positions as Th and F, respectively;
labels document the isotope but are not used to infer the coupling order.
"""
from collections.abc import Mapping
from dataclasses import replace

from .. import elements_c2 as _elements_c2  # populate REGISTRY_C2
from ..conventions import Conventions
from ..elements_c2 import REGISTRY_C2
from ..spec import Spin, enumerate_kets
from ..terms import ctx_from, terms_for_case
from ..backend_registry import Backend
from . import case_c


CANONICAL_UNITS = {
    **case_c.CANONICAL_UNITS,
    "A_par_Th": "MHz",
    "c_I_Th": "MHz",
    "eQq0_Th": "MHz",
    "eQq2_Th": "MHz",
    "g_N_Th": "",
}


def _two_spin_records(spins):
    try:
        records = tuple(spins)
    except TypeError as exc:
        raise case_c._InputValidationError(
            "case-c2 backend requires exactly two spin records",
            ("isotopologue", "spins")) from exc
    if len(records) != 2:
        raise case_c._InputValidationError(
            "case-c2 backend requires exactly two spin records",
            ("isotopologue", "spins"))

    coupled = []
    for index, couple_to in enumerate(("J", "F1")):
        record = records[index]
        spin_path = ("isotopologue", "spins", index)
        if not isinstance(record, Mapping):
            raise case_c._InputValidationError(
                "case-c2 spins must be records", spin_path)
        spin = Spin(
            label=case_c._required(record, "label", spin_path),
            I=case_c._required(record, "I", spin_path),
            couple_to=case_c._required(record, "couple_to", spin_path),
        )
        if spin.couple_to != couple_to:
            raise case_c._InputValidationError(
                "spins must form a chain coupled inner-first: spins[{}] "
                "({!r}) must have couple_to == {!r}, got {!r}".format(
                    index, spin.label, couple_to, spin.couple_to),
                spin_path + ("couple_to",))
        coupled.append(spin)
    return tuple(coupled)


def make_spec(basis, electronic, spins):
    """Translate a ThF+ two-spin TOML record into the native C2 StateSpec."""
    # case_c owns the shared electronic and basis validation.  Its temporary
    # empty spin list lets this adapter apply the distinct two-spin contract.
    common = case_c.make_spec(basis, electronic, ())
    coupled = _two_spin_records(spins)
    return replace(common, I=coupled[1].I, spins=coupled)


def _thf_v2_conventions():
    return Conventions(version="thf-v2")


CASE_C2_BACKEND = Backend(
    id="case_c2",
    case="c2",
    registry=REGISTRY_C2,
    default_terms=tuple(
        term.name for term in terms_for_case("c2", registry=REGISTRY_C2)),
    canonical_units=CANONICAL_UNITS,
    runtime_knobs=frozenset({"E_z", "B_z"}),
    optional_zero_parameters=frozenset({"d_e", "k_TP"}),
    make_spec=make_spec,
    default_conventions=_thf_v2_conventions,
    enumerate_kets=enumerate_kets,
    make_context=ctx_from,
)
