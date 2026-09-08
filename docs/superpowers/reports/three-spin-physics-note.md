# Physics note: equivalent-proton amides and outer metal spin

Current verdict: the restricted `amide_c2v` implementation is accepted; see the
final bounded implementation review below. The initial generic derivation that
follows is supporting algebra, not the target molecule backend.

## Supporting sequential case-(c) derivation

Use one unsymmetrized, ordered coupling chain

\[
  K_0\equiv J,\qquad K_r=K_{r-1}+I_r\ (r=1,\ldots,n),\qquad
  K_n\equiv F,\qquad n\leq 3,
\]

with basis

\[
 |J,\Omega;K_1,\ldots,K_n,m_F\rangle .
\]

This is a source-grounded extension of the current v1/v2 case-(c) basis.  It
does not require a target species.  The basis spans the full labelled-spin
product space and is also a useful intermediate representation before a
separate permutation-symmetry projection.  It does **not** by itself supply the
permutation-inversion restrictions of equivalent nuclei or the rotational
Hamiltonian of an asymmetric top.  In particular, adding three-spin case-(c)
support does not establish an SrNH2 model.

For production, recurse the spectator 6j factors.  Keep a direct product-basis
Clebsch-Gordan (CG) transform as an independent test oracle.  The recursive
form is only O(n) Wigner factors per candidate element, exposes exact sparsity,
and reduces identically to the existing v1/v2 formulas.  A full product-basis
construction scales with `(2J+1) product_i (2I_i+1)` at each J and hides those
selection rules; it is better used to falsify phases than to assemble every
production matrix.

The equations below use primes for the bra, unlike Brown & Carrington (B&C) on
the cited pages, and `[x]=2x+1`.  The spin list order is physical basis data:
changing it defines a recoupled basis and is not a harmless relabelling.

## Three reusable recouplers

### 1. Tensor on the rotor/electronic part

For a lab rank-k component `p` whose molecule-frame component is `q`, define

\[
 L^{(k)}_p=(-1)^{K_n'-m_F'}
 \begin{pmatrix}K_n'&k&K_n\\-m_F'&p&m_F\end{pmatrix},
\]

and, for each coupled spectator spin,

\[
 R_r^{(k)}=(-1)^{K_r+K_{r-1}'+k+I_r}
 \sqrt{[K_r'][K_r]}
 \left\{\begin{matrix}
 K_{r-1}&K_r&I_r\\K_r'&K_{r-1}'&k
 \end{matrix}\right\}.
\]

The axial rotor reduced element is

\[
 D_{kq}=(-1)^{J'-\Omega'}\sqrt{[J'][J]}
 \begin{pmatrix}J'&k&J\\-\Omega'&q&\Omega\end{pmatrix},
 \qquad q=\Omega'-\Omega .
\]

Thus an axial rotor tensor has geometry

\[
 \langle\mathrm{bra}|T^{(k)}_p|\mathrm{ket}\rangle
 =L^{(k)}_p\left(\prod_{r=1}^{n}R_r^{(k)}\right)D_{kq}.
\]

This is B&C Eq. (5.172) followed by Eq. (5.174) once per spectator and
Eq. (5.186) at the rotor.  It is exactly S9.1 with one further `R_r` when
`n=3`.  It supplies Stark/electronic-Zeeman geometry at `k=1,q=0`, and any
other declared axial tensor rank without new recoupling algebra.

### 2. Lab tensor on nucleus i

At the layer where `I_i` is added, B&C Eq. (5.175), operator on the second
constituent, gives

\[
 N_i^{(k)}=(-1)^{K_i'+K_{i-1}+I_i+k}\sqrt{[K_i'][K_i]}
 \left\{\begin{matrix}
 I_i&K_i'&K_{i-1}\\K_i&I_i&k
 \end{matrix}\right\}
 \langle I_i\Vert T^{(k)}(I_i)\Vert I_i\rangle .
\]

The full lab element is

\[
 L^{(k)}_p\left(\prod_{r=i+1}^{n}R_r^{(k)}\right)N_i^{(k)},
\]

times deltas on `J, Omega, K_1,...,K_(i-1)`.  For nuclear Zeeman,
`k=1`, `p=0`,
`<I_i||T1(I_i)||I_i>=sqrt(I_i(I_i+1)[I_i])`, and the Hamiltonian multiplier is
`-g_N_i mu_N B_z`, matching the current v1 sign.  B&C Eq. (8.3), printed
p.378/PDF p.410, explicitly writes the nuclear Zeeman Hamiltonian as one term
per nucleus (with optional shielding factors); the minimal current contract
continues to omit those shielding corrections.

Selection rule: the inner labels before `i` are diagonal; `K_i,...,K_n` may
change subject to rank-k triangles; `m_F'=m_F+p`.

### 3. Scalar product of a rotor tensor with nucleus i

First lift the rotor reduced element through the spins inside `I_i`:

\[
 \langle K_{i-1}'\Vert T^{(k)}_\mathrm{rot}\Vert K_{i-1}\rangle
 =\left(\prod_{r=1}^{i-1}R_r^{(k)}\right)
 \langle J'\Vert T^{(k)}_\mathrm{rot}\Vert J\rangle .
\]

Then B&C Eq. (5.173) gives the scalar at layer `i`:

\[
 \begin{aligned}
 S_i^{(k)}={}&(-1)^{K_{i-1}+K_i+I_i}
 \left\{\begin{matrix}
 I_i&K_{i-1}&K_i\\K_{i-1}'&I_i&k
 \end{matrix}\right\}\\
 &\times\langle K_{i-1}'\Vert T^{(k)}_\mathrm{rot}\Vert K_{i-1}\rangle
 \langle I_i\Vert T^{(k)}(I_i)\Vert I_i\rangle .
 \end{aligned}
\]

It is multiplied by deltas on `K_i,...,K_n,m_F`.  These outer deltas are the
content of B&C Eq. (5.176): the scalar acts wholly inside `K_i`.  Inner totals
`K_1,...,K_(i-1)` can be off diagonal.  This rule reproduces the current
distinction between an inner-spin scalar (a v1 substitution) and an outer-spin
scalar (the S9.3 recoupling); it remains valid for a third spin at any position.

## Operator and parameter contract

Use positional parameter names tied to the declared coupling order, rather
than chemical labels:

- `A_par_i`: axial magnetic hyperfine.  The new backend may expose one unified
  `hyperfine_A_par_i` term with `Delta J=0,+-1`; the existing v1/v2 split terms
  remain untouched.  The unified operator is
  `(A_par_i/Omega) T1(n_hat).T1(I_i)`, with signed `Omega` and `q=0`.
  The projection theorem makes its `Delta J=0` block exactly
  `A_par_i T1(J).T1(I_i)/[J(J+1)]`, i.e. B&C Eq. (9.50); its off-diagonal block
  is Eq. (9.51).  This normalization is not defined at `Omega=0` and must not
  be silently reused for a Sigma/case-(b) model.
- `c_I_i`: `c_I_i T1(J).T1(I_i)`, using the scalar recoupler at `k=1` and the
  rotor angular-momentum reduced element.  B&C Eqs. (8.7), printed p.378/PDF
  p.410, and (8.20), printed p.382/PDF p.414.
- `g_N_i`: the lab-nucleus recoupler above, multiplied by `-g_N_i mu_N B_z`.
- `eQq0_i`: B&C Eq. (9.53), using the scalar recoupler at `k=2`, terminal
  axial `D_(2,0)`, and the nuclear normalization

  \[
    \frac{1}{4\begin{pmatrix}I_i&2&I_i\\-I_i&0&I_i\end{pmatrix}}.
  \]

  The plus `1/4` is for B&C's stated convention that `q0` is the negative of
  the electric-field gradient.  Return zero for `I_i<1` *before* evaluating
  the denominator: such a nucleus has no rank-2 quadrupole moment and the
  printed formula otherwise presents `0/0`.
- `eQq2_i`: the same rank-2 scalar geometry with `q=Omega'-Omega=+-2`, using
  the repository's existing B&C-internal definition
  `eq2Q=-2 eQ <T2_(+-2)(grad E)>`.  B&C do not define this named constant, and
  the current bridge to Petrov-style `eQq2` remains unresolved.  Do not accept
  a differently normalized literature value without an explicit converter.

The `q=+-2` term is specific to a basis containing `Omega=+-1`.  A truly
generic case-(c) manifold can require other molecule-frame components (for
example `q=+-1` between `Omega=+-1/2`); those need separately named constants
or a general tensor-component input.  `eQq2_i` must not be presented as that
general interface.

Rotation, centrifugal distortion, Omega doubling, and PT-odd diagonal terms
remain diagonal in every nuclear-coupling label.  Stark and electronic axial
Zeeman terms use recoupler 1.  No pairwise nuclear spin-spin term is implied by
the per-nucleus list; adding one requires its own two-spin tensor operator.

## Required validation

1. **Exact legacy reductions:** at `n=1`, compare every element with
   `heff.elements_c`; at `n=2`, map `(K1,K2)=(F1,F)` and compare every term with
   `heff.elements_c2`, including the `Delta J=+-1` hyperfine blocks and both
   nuclear Zeeman operators.  Compare elements, not only spectra.
2. **Independent CG oracle:** let `U` have product-basis rows and coupled-basis
   columns, with entries

   \[
   U=\prod_{r=1}^{n}
   \langle K_{r-1}M_{r-1},I_r m_r|K_rM_r\rangle,
   \quad M_r=M_{r-1}+m_r,
   \]

   then build primitive product-basis operators directly and compare the
   coupled result `O_coupled=U^dagger O_product U` element by element.  Use
   explicit ladder operators
   for `J.I_i`, direct spherical sums for `n_hat.I_i` and rank-2 quadrupole,
   and direct `I_(i,z)` for nuclear Zeeman.  The oracle must not call the 6j
   recouplers it tests.
3. Exercise all three spin positions, off-diagonal intermediate labels, all
   `p=-1,0,+1` for rotor rank 1, and `Delta J=+-1` axial hyperfine.  Include at
   least one synthetic `I_i>=1` at an outer position for quadrupole; this is an
   algebra fixture, not a species claim.
4. Verify representative off-diagonal elements that distinguish a bra/ket
   phase swap in Eq. (5.175) and a valid-triad 6j column permutation.
   Hermiticity can remain exact under the wrong phase, as S9's two-spin
   nuclear-Zeeman and outer-spin scalar checks already demonstrate.  CG
   unitarity alone also cannot validate operator phases.
5. Verify the zero-spin collapse: inserting `I_r=0` collapses its spectator
   factor to one and recovers the next-shorter chain; all direct interactions
   with that spin vanish.

## Source ledger and unresolved physics

The local Zotero attachment inspected was John M. Brown and Alan Carrington,
*Rotational Spectroscopy of Diatomic Molecules* (Cambridge, 2003), SHA-256
`a5405bf27439474fbebc403511746093badbc1525d132f88f0d3c7f9c0282cc4`.
The load-bearing equations were visually checked on printed p.173/PDF p.205
(5.172)-(5.176), printed pp.174-175/PDF pp.206-207 (5.179), (5.186), printed
p.378/PDF p.410 (8.3), (8.7), printed p.382/PDF p.414 (8.20), and printed
pp.604-605/PDF pp.636-637 (9.50)-(9.53).  The appendix forms (5.174) and
(5.175), cross-checked against direct CG sums, are the source contract; do not
substitute the known faulty earlier main-text Eqs. (5.136)/(5.138).

Still outside this result: target-specific constants; equivalent-nucleus
exchange symmetry and nuclear-spin statistical weights; asymmetric-top or
case-(b) rotational structure; non-axial hyperfine tensors; pairwise nuclear
spin-spin interactions; and conversion of non-B&C `eQq2` conventions.  These
are separate physical model choices, not reasons to hold the common
unsymmetrized sequential basis and recouplers.

## Target clarification: metal monoamides

The requested `226RaNH2`, `87SrNH2`, or `43CaNH2`-type target changes the
model class.  It is a planar C2v asymmetric top in a Hund-case-(b) basis, not
the axial case-(c) model above.  The live prior art is
`C2V-Molecules/atm_core/physics.py`, whose state generator
uses

\[
 |N,K;J(N,S),F_N(J,I_N),I_T,F_{\rm core}(F_N,I_T),m\rangle,
 \qquad I_T=I_{H1}+I_{H2}.
\]

In the source this is the seven-column tuple
`[N,K,J,F_N,I_T,F,mF]`; `F` there is renamed `F_core` above only to leave the
final symbol available after adding a metal spin.  The current coupling and
exchange filter are implemented at `physics.py:55-148`, not inferred from the
generic case-(c) work.

For two protons (`i_H=1/2`) in the metal-amide `X 2A1` vibronic ground state,
the live keep rule is

\[
 (+1)(-1)^K(-1)^{1-I_T}=-1.
\]

Consequently even `K=K_a` retains the antisymmetric proton singlet `I_T=0`,
and odd `K` retains the symmetric triplet `I_T=1`.  The parity doublet within
a fixed nonzero `K_a` sector is still present: its space-fixed inversion
parity is `(-1)^(K_c)`, so exchange symmetry fixes the parity of `K_a`, not a
single `K_c` component.  The repository's manuscript states these points for
the `X 2A1` metal amides at `2026-06-09_ASYM.tex:132-143,173,624-632`.

Hirota independently establishes the underlying equivalent-nucleus rule:
two spin-1/2 nuclei have symmetric `I_T=1` and antisymmetric `I_T=0` spin
functions, which must be paired with the rovibronic exchange character; their
sum-spin hyperfine matrix elements equal a single-nucleus form with averaged
site constants.  The difference operator `I_-=I_1-I_2` connects rotational
levels of different `K`.  This was visually checked in Hirota (1985), section
2.3.10, printed p.42/PDF p.52.  The bare-NH2 example on that page has a
different vibronic/axis assignment; its listed `(K_a,K_c)` pairings must not
be copied over the metal-amide `X 2A1` rule.  The live C2V comments that say
only "Hirota p.52" are citing the PDF page; the unambiguous citation is
printed p.42/PDF p.52.

The current hydrogen hyperfine kernel at `physics.py:270-290` is diagonal in
`I_T`.  The exchange-filtered basis can contain both singlet/even-`K` and
triplet/odd-`K` sectors when both `K` parities are requested, but that legacy
kernel cannot mix them.  This is the existing approximation that drops the
`I_-`/different-`K` interaction.  Adding a distinguishable metal spin neither
changes nor invalidates the proton exchange filter.

### What "three spins" means here

With the production default `14N`, an odd metal isotope gives four physical
spin-bearing nuclei but three coupled spin groups:

1. nitrogen `I_N`;
2. the symmetry-adapted proton pair `I_T=0 or 1`;
3. the metal `I_M`.

NNDC/ENSDF gives `87Sr` ground-state `I_M=9/2` and `43Ca` ground-state
`I_M=7/2`; `226Ra` is even-even with a `0+` ground state and hence `I_M=0`.
Thus `226Ra14NH2` has only the two nonzero coupled groups already represented
by the current C2V basis.  `87Sr14NH2` and `43Ca14NH2` require the third group.
The isotope values were checked against NNDC ENSDF records for 87Sr, 43Ca,
and 226Ra, rather than inferred from mass-number parity.

The least-disruptive exact extension is

\[
 |((((N S)J I_N)F_N I_T)F_{\rm core} I_M)F,m_F\rangle .
\]

Appending `I_M` outside the existing core preserves every proton-exchange
restriction and gives a one-to-one zero-spin collapse `I_M=0 => F=F_core`.
It is an exact representation when every allowed `F` and every off-diagonal
intermediate block required by the operators are retained.  It does not
assert that `F_core` is a good quantum number.  Hirota section 2.3.10,
printed p.41/PDF p.51, recommends coupling the nucleus with the larger
hyperfine interaction nearest `J` for useful approximate labels and warns
that comparable interactions require off-diagonal intermediate-coupling
matrix elements.  If metal hyperfine is dominant, coupling the metal first is
therefore a reasonable later recoupling for sparse matrices or spectroscopic
labels, but it is not required for a correct full diagonalization.  Outer
metal coupling is the minimal compatibility choice with the live C2V code.

### Minimal backend contract and its limit

The TOML representation needs a named C2v-amide backend and an explicit
equivalent-pair declaration (`site spin=1/2`, `equivalent_count=2`, and the
vibronic exchange sign).  Three ordinary labelled `spins` records alone are
not enough, because `I_T` is a symmetry-selected variable taking 0 or 1, not
a fixed third single-nucleus spin.  A minimal fixed coupling contract is
`N+S=J`, `J+I_N=F_N`, `F_N+I_T=F_core`, `F_core+I_M=F`.

For a first restricted odd-metal backend, it is physically consistent to port
the existing rotation, electron spin-rotation, N/H hyperfine, quadrupole,
Stark, electron/nuclear Zeeman, and EDM kernels unchanged on the core; lift
core lab tensors through the outer metal spectator; and add explicitly only
the supported metal terms.  An isotropic `a_M S.I_M` contact term and
`-g_M mu_N B I_(M,z)` are well-defined with user-supplied constants.  This is
not yet a complete precision Hamiltonian for `87SrNH2` or `43CaNH2`: both
metal nuclei have `I_M>=1`, so metal electric quadrupole, anisotropic metal
hyperfine, shielding/nuclear spin-rotation, or nuclear spin-spin terms may be
needed if supported by data.  They should remain refused term IDs until their
tensor conventions and constants are supplied, rather than being silently
zeroed or synthesized.

The available molecule-specific spectroscopy does not fill that gap.
Thompsen, Sheridan, and Ziurys, *Chemical Physics Letters* 330, 373-382
(2000), DOI `10.1016/S0009-2614(00)01113-1`, reports that the 87Sr analog was
probably present but too weak to study and that the measured high-N SrNH2
hyperfine structure was unresolved (article p.378/PDF p.6 of the inspected
Zotero attachment).  Brewster and Ziurys, *J. Chem. Phys.* 113, 3141-3149
(2000), confirms planar `X 2A1` CaNH2, the `J=N+S` case-(b) labeling, and
unresolved 14N/1H hyperfine at high N (article pp.3141-3143/PDF pp.1-3), but
does not supply a 43Ca metal-hyperfine fit.  Target constants therefore remain
inputs, not defaults justified by these papers.

For the outer lift, recovering a core reduced matrix element from an existing
`p=0` matrix is valid under B&C Eq. (5.172): divide the core matrix element by
`(-1)^(F_core'-m) (F_core' k F_core; -m 0 m)`, then apply Eq. (5.174) with
spectator `I_M` and the final lab 3j.  For rank 1, choosing
`m=min(F_core',F_core)` gives a nonzero denominator for every allowed pair
except the forbidden `0 <-> 0` pair; return zero there and for
`|Delta F_core|>1`.  Metal `S.I_M` follows by extracting the core reduced
matrix element of `S` from the unit isotropic-electron-Zeeman kernel and using
Eq. (5.173).  This reuse is sound only if the source kernel is exactly the
dimensionless `T^1_0(S)` geometry, with its Zeeman coefficient and sign kept
outside the extraction.

Decisive checks for this target backend are: exact matrix and spectrum
reduction to the live seven-label C2V implementation at `I_M=0`; the analytic
two-spin spectrum of `a_M S.I_M` with all N/H degrees as spectators in a small
`N=K=0` fixture; and the complete set of metal-Zeeman projections.  An
independent CG product-basis comparison remains the phase-sensitive check for
representative nonzero-`I_M` lifted core tensors.

### Bounded implementation review (2026-09-08)

Physics verdict: the implemented `amide_c2v` basis and restricted operator set
match the contract above.  No recoupling sign or normalization defect was
found.

- `heff/elements_amide.py` was compared by normalized Python AST against the
  live `C2V-Molecules/atm_core/physics.py` at commit
  `3b77b021ab2b3f9256f1af1b15ab3116eddb1068`.  All eleven ported functions
  (rotation, spin rotation, anisotropic electron Zeeman, H/N hyperfine, N
  quadrupole, Stark, isotropic electron Zeeman, H/N nuclear Zeeman, and EDM)
  had identical executable bodies; only their surrounding documentation and
  Wigner imports differ.
- The outer active-core lift in `heff/backends/amide.py` agrees with corrected
  B&C Eq. (5.174); the direct metal `I_M,z` element agrees with Eq. (5.175);
  and the metal scalar `S.I_M` agrees with Eq. (5.173).  A separate exact-CG
  scratch calculation compared 62 active-core rank-1 sectors, 153 metal-`I_z`
  sectors, and 23 scalar-dot-product sectors and found zero symbolic/numerical
  discrepancy.  The repository test additionally compares signed coupled
  matrices against an uncoupled-spin construction, including all seven lifted
  core rank-1 terms in an `N=1`, `|K|=1`, `I_T=1`, `I_M=1/2` sector.
- The initially invalid attempt to treat the nested C2V `g_l` tensor as one
  scalar TOML coefficient was removed.  The backend now exposes the three
  real C2v spherical components `g_l_00`, `g_l_20`, and `g_l_22`, with the last
  populating equal `q=+2` and `q=-2` entries.  This reproduces a diagonal real
  principal-axis tensor; the source conversion is
  `g_l_00=-(g_xx+g_yy+g_zz)/sqrt(3)`,
  `g_l_20=(2g_zz-g_xx-g_yy)/sqrt(6)`, and
  `g_l_22=(g_xx-g_yy)/2`, where source body `z` is the amide `a` axis.
- The focused amide basis, source-fixture, and model tests passed `49/49`
  before the final expanded spectator check; the final model subset passed
  `9/9` after it.  The coordinating run reports the full suite at `365 passed,
  2 skipped` before that test-only addition.

The acceptance is deliberately limited to the implemented model.  It
preserves the source approximation that omits `I_-`-driven different-`K`
singlet/triplet mixing, and it does not claim metal anisotropic hyperfine or
metal quadrupole support.  Since no fitted odd-metal constants were found in
the inspected SrNH2/CaNH2 papers, a TOML file becomes a molecule-specific
prediction only when the user supplies convention-matched constants.  The
outer coupling order is exact for the complete enumerated basis but should not
be advertised as a hyperfine hierarchy or as proof that `F_core` is a good
quantum number.
