"""heff -- effective Hamiltonians for molecules, term-matrix first.

Import stays light on purpose (gate A8): numpy only. sympy is imported inside
heff.wigner's kernels, scipy inside heff.track's assignment call, matplotlib
never (plotting lives in notebooks, spec S3.8).
"""

__version__ = "0.1.0"


def __getattr__(name):
    if name in {"load_model", "list_models", "MoleculeModel", "Problem"}:
        from . import model
        return getattr(model, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

from . import (assemble, conventions, elements_c, elements_c2, engine, formalism,
              observe, params, spec, spectra, terms, track, twophoton,
              wigner)  # noqa: F401
from .assemble import (TermMatrices, active, build_term_matrices, coefficients,
                       hamiltonian, hamiltonian_batch, sweep_coefficients, vertex)
from .conventions import (Conventions, ef_label, n_hat_sign,
                          parity_operator, parity_phase, superposition_parity)
from .elements_c import dipole_geometry
from .elements_c2 import (REGISTRY_C2, axial_geometry, j_convergence,
                          outer_spin_scalar)
from .engine import SweepResult, eigh_batch, sweep
from .formalism import convert_formalism
from .observe import (expectation, g_factors, multi_curvature, offdiag,
                      pair_differential)
from .params import MU_B, MU_N, Param, ParamSet, thf_v1, thf_v2
from .spec import (KET_C, KET_C2, Blocking, ElecState, Spin, StateSpec,
                   block_by_mF, blocks_for, check_basis_invariants,
                   enumerate_kets, thf_spec)
from .spectra import dipole_matrix, label_lines, line_strengths
from .terms import Ctx, Rules, Term, ctx_from, term, terms_for_case
from .twophoton import (REGISTRY_2G, dyad_weights, two_photon_geometry,
                        two_photon_line_strengths, two_photon_matrix)

__all__ = ["wigner", "spec", "KET_C", "KET_C2", "Blocking", "ElecState",
           "Spin", "StateSpec",
           "block_by_mF", "blocks_for", "check_basis_invariants",
           "enumerate_kets", "thf_spec",
           "conventions", "formalism", "params", "Conventions", "ef_label",
           "n_hat_sign", "parity_operator", "parity_phase",
           "superposition_parity",
           "convert_formalism", "dipole_geometry", "MU_B", "MU_N",
           "Param", "ParamSet", "thf_v1", "thf_v2", "terms", "Ctx", "Rules", "Term",
           "ctx_from", "term", "terms_for_case", "elements_c", "elements_c2",
           "REGISTRY_C2", "axial_geometry", "outer_spin_scalar", "j_convergence",
           "assemble",
           "TermMatrices", "build_term_matrices", "coefficients", "active",
           "hamiltonian", "hamiltonian_batch", "sweep_coefficients", "vertex",
           "engine", "SweepResult", "eigh_batch", "sweep",
           "track", "observe", "expectation", "g_factors",
           "multi_curvature", "offdiag", "pair_differential", "spectra",
           "dipole_matrix", "label_lines", "line_strengths",
           "twophoton", "REGISTRY_2G", "dyad_weights", "two_photon_geometry",
           "two_photon_matrix", "two_photon_line_strengths",
           "load_model", "list_models", "MoleculeModel", "Problem"]
