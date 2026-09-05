"""heff -- effective Hamiltonians for molecules, term-matrix first.

Import stays light on purpose (gate A8): numpy only. sympy is imported inside
heff.wigner's kernels, scipy inside heff.track's assignment call, matplotlib
never (plotting lives in notebooks, spec S3.8).
"""

__version__ = "0.1.0"

from . import (assemble, conventions, elements_c, engine, formalism, observe,
              params, spec, spectra, terms, track, wigner)  # noqa: F401
from .assemble import (TermMatrices, active, build_term_matrices, coefficients,
                       hamiltonian, hamiltonian_batch, sweep_coefficients, vertex)
from .conventions import (Conventions, ef_label, n_hat_sign, parity_operator,
                          parity_phase, superposition_parity)
from .elements_c import dipole_geometry
from .engine import SweepResult, eigh_batch, sweep
from .formalism import convert_formalism
from .observe import (expectation, g_factors, multi_curvature, offdiag,
                      pair_differential)
from .params import MU_B, MU_N, Param, ParamSet, thf_v1
from .spec import (KET_C, Blocking, ElecState, StateSpec, block_by_mF,
                   blocks_for, check_basis_invariants, enumerate_kets, thf_spec)
from .spectra import dipole_matrix, label_lines, line_strengths
from .terms import Ctx, Rules, Term, ctx_from, term, terms_for_case

__all__ = ["wigner", "spec", "KET_C", "Blocking", "ElecState", "StateSpec",
           "block_by_mF", "blocks_for", "check_basis_invariants",
           "enumerate_kets", "thf_spec",
           "conventions", "formalism", "params", "Conventions", "ef_label",
           "n_hat_sign", "parity_operator", "parity_phase", "superposition_parity",
           "convert_formalism", "dipole_geometry", "MU_B", "MU_N",
           "Param", "ParamSet", "thf_v1", "terms", "Ctx", "Rules", "Term",
           "ctx_from", "term", "terms_for_case", "elements_c", "assemble",
           "TermMatrices", "build_term_matrices", "coefficients", "active",
           "hamiltonian", "hamiltonian_batch", "sweep_coefficients", "vertex",
           "engine", "SweepResult", "eigh_batch", "sweep",
           "track", "observe", "expectation", "g_factors",
           "multi_curvature", "offdiag", "pair_differential", "spectra",
           "dipole_matrix", "label_lines", "line_strengths"]
