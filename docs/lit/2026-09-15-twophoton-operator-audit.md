# Citation-grade audit: `heff/twophoton.py` (effective 2xE1 operator)

Reviewer: spectroscopy-reviewer subagent. Read-only; nothing in `/Users/arianjadbabaie/Code/heff` was edited.

Sources read (targeted, citation-grade):

- Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules*, `/Users/arianjadbabaie/Zotero/storage/CKZKCGXY/Brown and Carrington, Rotational Spectroscopy of Diatomic Molecules.pdf`
  - PDF pp. 191-194 / book pp. 159-162: Eqs. (5.108)-(5.118), Table 5.2, the scalar product (5.111)/(5.112)
  - PDF pp. 197-199 / book pp. 165-167: Eqs. (5.136)-(5.142)
  - PDF pp. 205-207 / book pp. 173-175: Eqs. (5.172), (5.174), (5.186)
  - PDF p. 283 / book p. 251: Eq. (6.234) and the `(-1)^{J-S+s}` parity phase
- `heff` working tree: `heff/twophoton.py`, `heff/elements_c2.py` (`axial_geometry`, `_spins`), `heff/elements_c.py` (`_ph`), `heff/wigner.py`, `heff/conventions.py`, `heff/spectra.py` (`_strengths_from_matrices`), `tests/test_twophoton.py`, `docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5
- `/Users/arianjadbabaie/Code/C2V-Molecules/atm_core/physics.py` (AC-Stark block) and `atm_core/polarizability.py` (Caldwell 2020 App. A pins)

Test status: `conda run -n Structure python -m pytest tests/test_twophoton.py tests/test_twophoton_closure.py -q` -> **14 passed in 8.49 s**. The verdicts below do not rest on that; independent checks I ran myself are marked.

---

## 1. Rank content - CONFIRMED

B&C (5.141), PDF p. 198 / book p. 166, verbatim:

```
T^K_p(A1, B1) = (-1)^{k1-k2+p} (2K+1)^{1/2} sum_{p1 p2} (k1 k2 K; p1 p2 -p) T^{k1}_{p1}(A1) T^{k2}_{p2}(B1)
```

With `k1 = k2 = 1` the triangle admits `K = 0, 1, 2` and nothing else. The registered set `{0, 2}` plus the omission of `K = 1` is therefore a statement about `K = 1` alone.

**S9.5.3's argument is sound.** `<1 p1 1 p2|1 P>` is antisymmetric under `p1 <-> p2` (3j column-exchange symmetry at `j1 = j2 = 1`, `J = 1`: the exchange phase is `(-1)^{2j-J} = -1`). So the `K = 1` projection of the dyad is `alpha^1 ~ (1/2)(d_a P d_b - d_b P d_a)`. When `P` is the projector onto the *complete* opposite-parity space it acts as the identity on `d|g>`, leaving `P_X (1/2)[d_a, d_b] P_X = 0` because the Cartesian components of `d = -e sum_i r_i` commute. The doc states the scope of that premise explicitly and shows `K = 1` returns at first order in `delta/Delta`, and at **O(1) for a restricted (single-intermediate) manifold** - which is the actual ThF+/HfF+ situation. That caveat is correct and is the one that matters for the planned pipeline.

**K = 0 is exactly the identity.** Measured, not asserted. On 229ThF+ at `J_max = 2` (192 kets, `I_Th = 5/2`) the `K = 0, dOmega = 0` channel matrix has off-diagonal max `0.0` and **every** diagonal entry `1.000000`, uniform across `J = 1, 2`, both signs of `Omega`, and every `(J, F1, F)` group. The team lead's 232ThF+ probe (`I_Th = 0`) reproduces this in the spinless limit. So within one electronic state with unit reduced element per channel, the `K = 0` channel is a **pure, state-independent light shift**: no transitions, and no differential shift either. All the spectroscopic content of the registered operator sits in `K = 2`.

## 2. Photon reading - CONFIRMED

**Spherical components.** B&C Table 5.2, PDF p. 193 / book p. 161: `T^1_0(V) = V_Z`, `T^1_{+-1}(V) = -+(1/sqrt2)(V_X +- i V_Y)`. That is exactly `_spherical`: `{1: -(vx + i vy)/sqrt2, 0: vz, -1: (vx - i vy)/sqrt2}`.

**The scalar product.** B&C (5.111), same page: `T^k(R).T^k(S) = sum_p (-1)^p T^k_p(R) T^k_{-p}(S)`, and (5.112) confirms it equals the Cartesian `A.B` at `k = 1`. It is a bilinear identity, so it is valid for complex components - the doc says so and is right. Hence

```
d.eps = sum_p (-1)^p T^1_p(d) eps_{-p}  ==  sum_p c[p] T^1_p(d),    c[p] = (-1)^p eps_{-p}
```

which is `_leg` verbatim. **CONFIRMED.**

**The sigma+ Jones vector.** `(-1, -i, 0)/sqrt2` **is** the standard Condon-Shortley `e_{+1} = -(xhat + i yhat)/sqrt2`. I evaluated `_leg` on it directly: `c[+1] = 0.9999999999999999`, `c[0] = c[-1] = 0`. So absorbing a sigma+ photon raises `m_F` by one, and the anchor the whole convention rests on holds.

**The ladder reading.** Change one call in `dyad_weights`:

```python
c_a = _leg(np.asarray(eps2, dtype=complex))     # was: _leg(np.conj(np.asarray(eps2, dtype=complex)))
```

That is the whole change. Physically it is right: the positive-frequency (absorbing) part of `E` carries `eps`, not `eps*`, so two absorbed photons give `(d.eps2)(d.eps1)`. Verified by feeding a pre-conjugated sigma+ as `eps2` (which cancels the internal conjugate): the reachable set becomes `Delta m_F = {+2}`, as the ladder reading requires.

## 3. Coupling order and contraction - CONFIRMED

**`_cg`.** `(-1)^{j1-j2+M} sqrt(2J+1) w3j(j1, j2, J, m1, m2, -M)` is the standard CG-to-3j relation and is identical to the bracket in B&C (5.141). **CONFIRMED.**

**Slot assignment.** B&C (5.141) puts slot 1 on the leftmost operator factor, and (5.142) preserves that ordering in the reduced-element chain: `<j||T^{k1}(A1)||j''><j''||T^{k2}(B1)||j'>`, bra-side leg first. The leftmost factor of `T_eff = (d.eps2*)|i><i|(d.eps1)` is the `eps2*` leg. The code assigns `c_a = _leg(conj(eps2))` to slot 1 via `_cg(1, p_a, 1, P - p_a, K, P)`. **Consistent.**

**Contraction sign.** Inverting (5.141) by CG orthogonality gives `T^1_{p_a} T^1_{p_b} = sum_K <1 p_a 1 p_b|K P> T^K_P` with `P = p_a + p_b`, so

```
T_eff = sum_{K,P} w^K_P alpha^K_P ,     w^K_P = sum_{p_a + p_b = P} <1 p_a 1 p_b|K P> c_a[p_a] c_b[p_b]
```

with **no residual phase**. B&C (5.172) makes `<a|alpha^K_P|b>` nonzero only when `-m_a + P + m_b = 0`, i.e. `Delta m_F = +P`, and `axial_geometry` enforces exactly that (`if abs(m2 - m1 - p) > 1e-9: return 0.0`). So multiplying `geometry(K, P = m_bra - m_ket)` by `w[(K, P)]` **is the right pairing**. `w[(K, -P)]` with a `(-1)^P` would be wrong. `test_dyad_weights_reconstruct_the_ordered_polarization_product` exercises the inversion identity numerically and is a genuine slot-reversal detector.

**The amplitude path is correct too.** `two_photon_line_strengths` builds one matrix per fixed `(K, dOmega, P)` and hands the whole set to `heff.spectra._strengths_from_matrices`, which sums `c * (conj(evecs_a).T @ D @ evecs_b)` over channels **before** `abs(...)**2`. Amplitudes are summed then squared, as required for the arbitrary-polarization-pair use case.

## 4. Spectator chain in `axial_geometry` - CONFIRMED, line by line

Basis `|((J I_Th) F1, I_F) F, m_F>`; B&C's unprimed slot is the bra throughout. Code variables: `J1/G1/F1/m1` = ket `(J, F1, F, m_F)`, `J2/G2/F2/m2` = bra.

**Line 1** = B&C (5.172), PDF p. 205 / book p. 173: `<j,m|T^k_p(A)|j',m'> = (-1)^{j-m} (j k j'; -m p m') <j||T^k(A)||j'>`. Code: `_ph(F2 - m2) * w3j(F2, k, F1, -m2, p, m1)`. **Matches.**

**Lines 2 and 3** = B&C (5.174), same page: `<j1,j2,j12||T^{k1}(A1)||j1',j2',j12'> = delta_{j2 j2'} (-1)^{j12'+j1+k1+j2} [(2j12+1)(2j12'+1)]^{1/2} {j1' j12' j2; j12 j1 k1} <j1||T^{k1}(A1)||j1'>`.

- Line 2, outer layer, `(j1, j2, j12) = (F1, I_F, F)`: phase `(-1)^{F_ket + F1_bra + K + I_F}` -> `_ph(F1 + G2 + k + I_F)`; factor `sqrt((2F2+1)(2F1+1))`; 6j `{G1 F1 I_F; F2 G2 K}` -> `w6j(G1, F1, I_F, F2, G2, k)`. **Matches argument-for-argument.**
- Line 3, inner layer, `(j1, j2, j12) = (J, I_Th, F1)`: phase `(-1)^{F1_ket + J_bra + K + I_Th}` -> `_ph(G1 + J2 + k + I_Th)`; factor `sqrt((2G2+1)(2G1+1))`; 6j `{J1 G1 I_Th; G2 J2 K}` -> `w6j(J1, G1, I_Th, G2, J2, k)`. **Matches.**

**Line 4** = B&C (5.186), PDF p. 207 / book p. 175: `<J,Omega||D^{(k)}_{.q}(w)*||J',Omega'> = (-1)^{J-Omega} (J k J'; -Omega q Omega') [(2J+1)(2J'+1)]^{1/2}`. Code: `_ph(J2 - Om2) * sqrt((2J2+1)(2J1+1)) * w3j(J2, k, J1, -Om2, q, Om1)`. **Matches**, and the 3j projection forces `q = Omega_bra - Omega_ket`, which is what `axial_geometry` sets.

`_ph` raises on a non-integer exponent; each of the four phase arguments above is integer-valued for this basis, so the guard is live rather than decorative.

**Reciprocity.** I swept `<a|alpha^K_P|b> - (-1)^P <b|alpha^K_{-P}|a>` over strided pairs across all three registered channels: residual **exactly 0.0**. The team lead's 232ThF+ probe gives the same, `P = -2..+2`, both `|dOmega|` channels. The relation follows from (5.172) plus the 3j identity `(j' K j; -m' -P m) = (j K j'; -m P m')` (column reversal times overall projection-sign reversal, each contributing `(-1)^{j+K+j'}`) together with real reduced elements obeying `<j'||T^K||j> = (-1)^{j-j'} <j||T^K||j'>`.

Consequence for the declaration: since the registered `fn(bra, ket, ctx)` evaluates at `P = Delta m_F`, reciprocity gives `M(b,a) = (-1)^{Delta m_F} M(a,b)` - symmetric at even `Delta m_F`, antisymmetric at odd. `K = 2` reaches `|Delta m_F| = 1`, so its matrix genuinely is not symmetric. **`hermitian=False, real=True` is correct, and the "odd-Delta m_F antisymmetry" comment is the right explanation, not a mask for a phase error.** (`K = 0` is `Delta m_F = 0` only, so `hermitian=True` is fine - and it is the identity anyway.)

## 5. The `Delta Omega = +-2` channel - CONFIRMED, and forced by the convention

`heff/conventions.py:parity_phase` is `(-1)^{J - S - ell + s}`, pinned to B&C (6.234), PDF p. 283 / book p. 251, whose text names "the `(-1)^{J-S+s}` phase factor" in

```
|eta, Lambda_s; J, M; +> = (1/sqrt2){ |eta, Lambda_s; S, Sigma; J, Omega, M> + (-1)^p |eta, -Lambda_s; S, -Sigma; J, -Omega, M> }
```

`parity_operator` applies that phase while negating `Om` and leaving every other dtype field alone. **The phase carries no `Omega` dependence.** That is the load-bearing fact: parity maps the `q = +2` component onto the `q = -2` component with equal weight, so a single shared `alpha_K2_dOm2` with no relative sign is what parity conservation requires. A relative sign would fake an e <-> f two-photon line.

**Negative control I built** (the one the docstring promises but the repo does not contain): sign-flip only the `q = +2` component and re-measure the commutator with `parity_operator` on the 229ThF+ `J_max = 2` basis.

| `P` | `max\|T\|` | shipped `max\|[T,Par]\|` | `q=+2`-flipped `max\|[F,Par]\|` |
|---|---|---|---|
| -2 | 0.4899 | 0.000e+00 | 9.798e-01 |
| -1 | 0.3703 | 0.000e+00 | 7.407e-01 |
| 0 | 0.3513 | 0.000e+00 | 7.026e-01 |
| +1 | 0.3703 | 0.000e+00 | 7.407e-01 |
| +2 | 0.4899 | 0.000e+00 | 9.798e-01 |

Independent corroboration from a separately derived codebase: `C2V-Molecules/atm_core/physics.py:_alpha2_body` sets `a2[+2] = a2[-2] = (alpha_b - alpha_c)/2` - the same single scalar for both signs, reached by a different derivation (Caldwell 2020 App. A) for a different molecule class.

## 6. What `alpha_K_dOm` multiplies under `bc_5p142_reduced` - CONFIRMED

`alpha_K_dOm` multiplies the product of the four geometry lines above, and **is** the body-frame reduced element over the electronic-vibrational coordinate:

```
alpha_K_dOmega  ==  <eta' || alpha^K_q || eta>  ==  (1/Delta) <eta' || T^K_q(d, d) || eta> ,    q = Omega' - Omega
```

with the rotational structure supplied by B&C (5.186)'s D-matrix reduced element, not by `alpha`. Units: (dipole)^2 / energy, i.e. `(MHz/(V/cm))^2 / MHz` in package units, which is what the `two_photon_line_strengths` docstring states. A strength therefore comes out in units of `alpha^2`.

**Factors the geometry silently absorbs**, all of them square-root degeneracy factors and none of them `sqrt(2K+1)`:

| factor | where | source |
|---|---|---|
| `sqrt((2F'+1)(2F+1))` | line 2 | B&C (5.174) |
| `sqrt((2F1'+1)(2F1+1))` | line 3 | B&C (5.174) |
| `sqrt((2J'+1)(2J+1))` | line 4 | B&C (5.186) |

**`sqrt(2K+1)` appears nowhere in `axial_geometry`**, and that is correct given the definition above: the `(2K+1)^{1/2}` of B&C (5.142) lives on the *closure relation* that defines `T^K(d,d)`, not on the spectator geometry.

**Flag for the resolved-sum bridge.** A user relating `alpha_K_dOm` to a resolved sum over Denis 2015 dipoles must supply (5.142)'s full prefactor themselves:

```
<eta,j||T^K(d,d)||eta',j'> = (2K+1)^{1/2} (-1)^{K+j+j'} sum_{eta'' j''} {1 1 K; j' j j''} <eta,j||T^1(d)||eta'',j''> <eta'',j''||T^1(d)||eta',j'>
```

The `(-1)^{K+j+j'}` is **`j`-dependent**, so the map from a resolved sum onto one `j`-independent scalar `alpha` is constant only in the closure limit. That is the same limit S9.5.5 already flags as unmet at JILA detunings (`Delta ~ 0.16-1.5 GHz` against an intermediate `B_i ~ 7 GHz`).

## 7. Other findings

- **`K = 0` is exactly the identity** (item 1). Within X 3Delta1 it carries no spectroscopic information, and in a strength pipeline it produces a spurious zero-frequency self-line for every state. Neither the term's `cite` string nor doc S9.5 records this.
- **`dyad_weights` returns `K = 1` weights that no registered operator consumes.** Declared (OPEN-21), but a general transition pipeline will silently drop that channel - and per S9.5.3 it survives at O(1) for the single-intermediate case that ThF+/HfF+ actually is.
- **The Raman reading is hardwired**, while doc S9.5.1(3) states that "any polarization-labelled result must state which reading it uses."
- **A dangling test reference** and a **normalization sentence that contradicts the module** - see defects 1 and 2.
- No unit errors found: everything in the two-photon path is dimensionless geometry, and the only dimensioned quantity is the caller-supplied `alpha`, documented in MHz-anchored units.

---

## C2v / Caldwell cross-check

### (a) Contraction convention: equivalent, verified numerically

`C2V-Molecules/atm_core/physics.py` header states the AC-Stark operator as

```
H = -(E0^2/4) sum_{k=0,1,2} sum_p (-1)^p T^k_p(alpha) T^k_{-p}(eps, eps*)
```

which is B&C (5.111)'s scalar product `T^k(alpha) . T^k(eps, eps*)`. heff instead writes a plain `sum_{K,P} w^K_P alpha^K_P` with no `(-1)^P`. **These are the same object**, because the `(-1)^P` and the index reversal are already inside `w`:

```
w^K_P = (-1)^{K+P} (eps2* (x) eps1)^K_{-P}
```

I verified this over 60 random **complex** polarization pairs: max deviation `8.882e-16`. (Analytically: substituting `c[p] = (-1)^p eps_{-p}` pulls out `(-1)^P`, and `<1,-q_a 1,-q_b|K,-P> = (-1)^{2-K} <1 q_a 1 q_b|K P>` supplies the `(-1)^K`.)

### (b) `z_K` normalization: heff carries none

Caldwell (A.10b), via `atm_core/polarizability.py`: `(A_pm)^K_P = (1/z_K) T^K_P(d, R_pm d)` with `z_2 = sqrt(3/2)`. heff's `alpha^K ~ T^K(d,d)/Delta` is the **bare** coupled tensor with no `z_K`. So:

```
alpha_K2_dOm*  (heff)  =  sqrt(3/2)  x  a Caldwell-convention rank-2 polarizability
```

Convenient for the user: `C2V`'s `physics._alpha2_body` already folds `z_2` back out (its own docstring: "returns the BARE spherical tensor T^2_q (z_2 folded out), so `a2[0] = sqrt(3/2) * alpha2_q0`"). So **`C2V`'s `a2` dict and heff's `alpha_K2_*` are in the same convention**, up to the `1/Delta` that heff carries and the AC-Stark case does not:

| heff | C2V `_alpha2_body` |
|---|---|
| `alpha_K2_dOm0` | `a2[0] = (2 alpha_a - alpha_b - alpha_c)/sqrt(6)` |
| `alpha_K2_dOm2` | `a2[+2] = a2[-2] = (alpha_b - alpha_c)/2` |
| `alpha_K0_dOm0` | `alpha_bar = (alpha_a + alpha_b + alpha_c)/3`, no factor |

Caldwell's `alpha_2 = (2/3)(alpha_par - alpha_perp)` convention is the `z_2`-divided one; do not mix it with the bare form without the `sqrt(3/2)`.

### (c) Does `dyad_weights` reduce to Caldwell's tensor at `eps1 = eps2`? Yes, exactly

The AC-Stark case is degenerate (`eps` and `eps*`, same beam), i.e. the Raman reading with `eps1 = eps2 = eps`. Two `(-1)^K` factors cancel there - one from heff's conjugated slot order `(eps2*, eps1)`, one from swapping `(eps*, eps)` into Caldwell's `(eps, eps*)` (slot exchange on a coupled rank-K tensor of two rank-1s gives `(-1)^{1+1-K} = (-1)^K`). The result:

```
eps1 = eps2 = eps   =>   w^K_P = (-1)^P (eps (x) eps*)^K_{-P}
```

so `T_eff = sum_K T^K(alpha) . T^K(eps, eps*)`, which is **exactly** Caldwell's bracket, i.e. `-4/E0^2` times the C2v AC-Stark Hamiltonian. Verified over the same 60 random complex vectors: max deviation `8.951e-16`.

Sanity limit, lab-`Z` linear polarization, which the C2v header calls out ("only lab `p = 0` survives and the vector (`k = 1`) term vanishes"):

```
w^0_0 = -0.577350   (= -1/sqrt3 = <1 0 1 0|0 0>)
w^2_0 = +0.816497   (= sqrt(2/3) = <1 0 1 0|2 0>)
w^1_P = 0           for all P
```

All three match. The two codebases agree.

---

## Defects

| # | file:line | What | Fix |
|---|---|---|---|
| 1 | `tests/test_twophoton.py:68` | The V24 docstring cites `test_..._fails_if_the_q_plus_2_component_changes_sign` as the negative control. **No such test exists** anywhere (grepped `tests/` and `heff/`); the only hit is this mention of it. The gate's own falsifiability claim is unbacked. | Add the test - my hand version gives `max\|[F,Par]\|` of 0.70-0.98 against 0.000e+00 shipped - or delete the sentence. |
| 2 | `tests/test_twophoton.py:148` | Docstring reads "`conventions.two_photon_norm = 'bc_5p142_reduced'`, unit one-photon reduced elements". `heff/twophoton.py`'s `_CITE` explicitly rejects that reading: "not unit ONE-photon reduced elements, which is a different (and unmade) claim". A direct doc-vs-code contradiction about the normalization the whole channel scale rests on. | Replace with "unit rank-K reduced element per channel". |
| 3 | `heff/twophoton.py:118` (K=0 `cite`) and `docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5 | Neither records that the `K = 0` channel matrix is **exactly the identity**. The consequence for a strength pipeline (a zero-frequency self-line on every state, and no differential shift) is undocumented. | One line in the `cite` string plus a caller caveat in `two_photon_line_strengths`. |
| 4 | `heff/twophoton.py:37` (`def dyad_weights`) | The Raman conjugation is hardwired, yet doc S9.5.1(3) requires every polarization-labelled result to name its reading. For a general pipeline over arbitrary polarization pairs the caller cannot select the ladder reading without editing the module. | `def dyad_weights(eps1, eps2, *, raman=True)`; conjugate `eps2` only when true. One-line body change. |
| 5 | `heff/twophoton.py`, `CHANNELS` / OPEN-21 | `K = 1` is registered nowhere, but S9.5.3 establishes it survives at **O(1)** for a single intermediate electronic state - the actual ThF+/HfF+ configuration. Declared, not hidden, but it caps what the planned general transition pipeline can claim. | Escalate to Arian before building the pipeline on the closure operator alone. |

None of the five invalidates the operator: the geometry, the phases, the slot order, the contraction pairing, and the `q = +-2` sign are all correct as shipped.

---

**summary: 9 claims; 9 CONFIRMED, 0 WRONG, 0 CANNOT DETERMINE; 5 defects, none invalidating the operator.**
