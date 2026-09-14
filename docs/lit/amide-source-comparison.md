# Amide source comparison

**Result.** The new `amide_c2v` registry reduces to the live C2V amide
Hamiltonian at `I_M=0` for the requested small full basis.  The maximum
elementwise difference is `5.329070518200751e-15 MHz`, with no matrix
symmetrization applied in either comparison path.

## Inputs and method

- Source: `C2V-Molecules/atm_core/physics.py`
  at `3b77b021ab2b3f9256f1af1b15ab3116eddb1068` (verified during this run;
  `git status --short -- atm_core/physics.py` was empty).
- Basis: `S=1/2`, `I_N=1`, `i_H=1/2`, `I_M=0`, `N=0..1`, requested
  `|K|={0,1}`, `vibronic_sign=+1`, and every `mF`.  Each representation has
  132 kets.  Source mF-blocks were concatenated and reordered by the seven
  quantum-number values to the named heff order, mapping `F_core -> F`.
- Constants: `tests/test_elements_amide.py::PARAMS`, including
  `g_l[0][0]=-0.12`, `g_l[2][0]=0.07`, and
  `g_l[2][+-2]=-0.04`; added `d_0=0.8`, `g_s=1.7`, `g_H=2.3`,
  `g_N=0.403761`, `E_z=17.0 V/cm`, and `B_z=0.19 G`.
- The probe AST-extracted the live C2V Wigner wrappers, basis generator, and
  the ten native Hamiltonian kernels.  It then evaluated every matrix entry
  directly.  The heff side evaluated every entry of all 23 registered terms
  directly before comparing the public registry-masked assembly.  The two
  metal terms were retained with nonzero `a_M=0.31` and `g_M=0.27`; their raw
  matrices are identically zero at `I_M=0`.

## Measurements

| Quantity | Value |
| --- | ---: |
| registered heff terms | 23 |
| basis content equal after value reorder | `True` |
| max source raw non-Hermiticity | `0.0` |
| max heff raw non-Hermiticity | `0.0` |
| max raw heff vs raw source difference | `5.329070518200751e-15 MHz` |
| max public assembled heff vs raw source difference | `5.329070518200751e-15 MHz` |
| max unmasked heff term vs registry-masked term difference | `0.0` |
| max `hyperfine_M_contact` / `zeeman_metal` element | `0.0` / `0.0` |

No term or total matrix was replaced by `(H + H.T)/2`; Hermiticity is an
observed raw property here.  This comparison establishes the `I_M=0`
reduction only.  It does not validate nonzero-metal recoupling, molecular
constants, or an end-to-end species model.

The original comparison probe is archived privately with the source checkout.
This report records that comparison; the package's maintained amide checks are
in `tests/test_elements_amide.py`.
