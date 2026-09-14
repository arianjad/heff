# 229ThF+ electric-quadrupole estimate audit

Status: complete, 2026-09-08. The numerical values audited were present at
checkout `5db8df7dce01707731752cd379f11acf1d8100c1`. This report also checks the
subsequent metadata-only qualification of those values in the live worktree.

## Bottom line

Neither `eQq0_Th = -2600 +/- 1000 MHz` nor
`eQq2_Th = +300 +/- 100 MHz` is a convention-validated quantitative estimate
for the Brown-and-Carrington (B&C) Hamiltonian implemented by `heff`. Their
magnitudes are useful only as order-of-magnitude sensitivity placeholders.
The quoted `+/-1000` and `+/-100 MHz` ranges are unsupported as calibrated
uncertainties and should not define a likelihood, confidence interval, or
Monte-Carlo prior.

The live parameter metadata now implements that ruling: both values remain
numerically unchanged, have `status="placeholder"`, carry no numerical
uncertainty, and warn that the signed Petrov-to-B&C conversion is unresolved.
No Hamiltonian kernel or quadrupole number was changed by this audit.

## Sources inspected

Long references were overviewed before targeted page reads. Earlier project
digests and reports were treated only as pointers, not evidence.

- Brown and Carrington, *Rotational Spectroscopy of Diatomic Molecules*,
  printed pp. 132, 145, and 604-605, Eqs. (4.30)-(4.31), (5.33), and
  (9.52)-(9.53). Local
  Zotero attachment: `CKZKCGXY/Brown and Carrington, Rotational Spectroscopy
  of Diatomic Molecules.pdf`.
- Petrov et al., Phys. Rev. A **98**, 042502 (2018), printed pp. 2-4,
  Eqs. (19), (22)-(25), and Sec. IV. Local Zotero attachment:
  `2CWMNAD7/Petrov et al. - 2018 - Evaluation of CP violation in HfF +.pdf`.
- Skripnikov and Titov, Phys. Rev. A **91**, 042504 (2015), printed
  pp. 7-9 and Table II (arXiv:1503.01001).
- Porsev et al., Phys. Rev. Lett. **127**, 253001 (2021), printed pp. 1 and 4,
  especially Table IV (arXiv:2107.14723).
- Chakraborty and Sahoo, Phys. Rev. A **113**, L020801 (2026),
  *Comprehensive Assessment of Th3+ Properties for Nuclear Clock and
  Fundamental Physics Applications*, printed p. 3, Table IV
  (arXiv:2511.15346).
- Ng et al., Phys. Rev. A **105**, 022823 (2022), printed p. 5 and Table I
  (arXiv:2202.01346).

## Operative convention comparison

B&C define

```text
H_Q = -e T^2(grad E) . T^2(Q)                         (4.30)
T^2(grad E) = -(1/(4 pi epsilon_0))
              sum_i (e_i/R_i^3) C^2(theta_i,phi_i)   (4.31)
C^2_q = sqrt(4 pi/5) Y_2q                            (5.33, l=2)
```

and their case-(a) matrix element has the electronic coefficient
`-(1/2)eQ <T^2_q(grad E)>` [B&C, printed p. 604, Eq. (9.52)]. Its `q=0`
specialization instead carries `eq0Q/4`, and B&C state that `q0` is the
negative of the electric-field gradient [printed p. 605, Eq. (9.53) and the
sentence immediately following it]. Comparing Eqs. (9.52) and (9.53) gives

```text
eq0Q_B&C = -2 eQ <T^2_0(grad E)> .
```

Petrov et al. print

```text
eQq0 = 2 eQ <3Delta_1 | sum_i sqrt(2 pi/5) Y_20/r_i^3
                         | 3Delta_1>                  (22)
eQq2 = 2 sqrt(6) eQ <3Delta_1 |
                  sum_i sqrt(2 pi/5) Y_22/r_i^3
                         | 3Delta_-1>                 (23)
```

[Petrov 2018, printed p. 3]. Under the load-bearing assumption that Petrov and
B&C use the same molecular axis and electronic-state phase, literal use of
Petrov's printed `sqrt(2 pi/5) Y_2q` gives
`eq0Q_B&C = -sqrt(2) eQq0_Petrov`. If the printed factor was intended to be
the Racah-normalized `C^2_q = sqrt(4 pi/5)Y_2q`, the result is instead
`eq0Q_B&C = -eQq0_Petrov`. The ambiguity is load-bearing for the sign and
magnitude of `eQq0_Th`, even though the repository open item was originally
named for `q=2`.

For `q=2`, Petrov's Eq. (19) prints the same `q` on the nuclear and electronic
rank-2 tensors rather than the standard scalar-product `q,-q` pairing. Together
with the `sqrt(2 pi/5)` issue, this leaves two candidate conversions already
derived in the Hamiltonian note:

```text
eq2Q_B&C = -eQq2_Petrov/sqrt(3)   [literal printed normalization]
eq2Q_B&C = -eQq2_Petrov/sqrt(6)   [if the intended tensor is C^2_q].
```

The paper does not choose between them. One relative statement *is* fixed by
Petrov's Eqs. (22)-(23): if its `eQq0` is inserted directly as the B&C `q=0`
parameter, the correspondingly normalized `q=2` parameter is
`eQq2_Petrov/sqrt(6)`. The current raw ratio `300/(-2600) = -0.1154` would then
be `-0.0471`, smaller by `sqrt(6)`. Thus the two raw placeholder values cannot
simultaneously be interpreted as a convention-consistent B&C parameter pair.

## Claim verdicts

[claim 1] The HfF+ values are a real quantitative anchor, and nuclear-Q
scaling within that molecule is supported.

- **verdict: OK**
- **classification: justified quantitative anchor, limited to HfF+ and the
  nuclear-moment factor**
- **source:** Petrov 2018, printed p. 3, Eqs. (22)-(23) and the paragraph below
  them.
- **reason:** Petrov reports `eQq0=-2100 MHz` and `eQq2=110 MHz` for
  177HfF+, and `-2400 MHz` and `125 MHz` for 179HfF+. The paper explicitly
  says the isotope ratios correspond to `Q(179Hf)/Q(177Hf)`, using 3.793 b and
  3.365 b. It does not extend this scaling across HfF+ and ThF+ electronic
  wavefunctions.

For 229Th, Porsev 2021 recommends `Q=3.11(2) eb` from four atomic-HFS
extractions [printed p. 4, Table IV], giving the arithmetic factor
`3.11/3.365 = 0.9242`. A newer Th3+ calculation recommends `2.91(3) b`
[Chakraborty and Sahoo 2026, printed p. 3, Table IV], 6.4% lower, and tabulates
the spread among earlier determinations. The nuclear-Q input is consequently
method-dependent at a level larger than Porsev's internal quoted uncertainty,
although this disagreement is still much smaller than the molecular transfer
uncertainties below.

[claim 2] `eQq0_Th = -2600 +/- 1000 MHz` is a signed quantitative B&C
parameter estimate.

- **verdict: STOP**
- **classification: order-of-magnitude sensitivity only; signed value and
  uncertainty unsupported**
- **source:** B&C printed pp. 132, 145, and 604-605, Eqs. (4.30)-(4.31),
  (5.33), and (9.52)-(9.53);
  Petrov 2018 printed p. 3, Eq. (22); Skripnikov and Titov 2015 printed
  pp. 7-9 and Table II.
- **reason:** Scaling Petrov's `-2100 MHz` by 0.9242 and the assumed project
  range `R_el=1...1.65` produces `-1.94...-3.20 GHz`, so `-2.6 GHz` is inside
  that constructed range. No inspected source establishes that the ThF+/HfF+
  electric-field-gradient ratio lies between 1 and 1.65. The endpoint 1 is an
  assumption; the other endpoint is inferred from a magnetic hyperfine ratio,
  which probes a different electronic operator. Skripnikov and Titov explicitly
  observe that ThF/ThO trends in `A_parallel` behave inconsistently with the
  electron-EDM effective field [p. 7], and explain that the leading HFS matrix
  element is 7s-7s whereas the effective-field matrix element is 7s-7p [p. 9].
  This is direct evidence that core-sensitive operator ratios need not track
  across molecules; it does not supply an EFG ratio. Independently, if the
  molecular-axis and electronic-phase conventions are identified as assumed
  above, the Petrov-to-B&C conversion reverses the `q=0` sign under either
  candidate normalization. The printed sources do not settle the full
  convention identification or which magnitude factor applies, so this audit
  does not select or apply a numerical correction.

[claim 3] `eQq2_Th = +300 +/- 100 MHz` is a signed quantitative B&C parameter
estimate.

- **verdict: STOP**
- **classification: order-of-magnitude sensitivity only; signed value and
  uncertainty unsupported**
- **source:** Petrov 2018, printed p. 3, Eqs. (23)-(25); Ng et al. 2022,
  printed p. 5 and Table I; B&C printed p. 604, Eq. (9.52).
- **reason:** Petrov's HfF+ model is
  `eQq2 = 483 w Q <1/r^3>_5d MHz` [Eq. (24)], with
  `<1/r^3>_5d=4.86 a.u.` and `w=0.014` inferred from
  `G_parallel approximately 2-2.002319+w` [Eq. (25)]. This reproduces its
  HfF+ scale, but the paper identifies the admixed Pi component as chiefly
  `|5s 5d pi>`; it gives no ThF+ `6d` radial matrix element and no ThF+
  single-channel validation of Eq. (25). The project interval
  `<r^-3>_6d/<r^-3>_5d=0.6...1.0` is therefore a scenario range, not a
  source-calibrated bound.

Ng et al. measure only `|g_F=3/2|=0.0149(3)` [Table I]. Their conversion gives
`G_parallel=-0.042(2)` if `g_F>0` and `+0.048(2)` otherwise; they say the
spectroscopy is not sensitive to the sign. Using the positive branch gives
`w=0.0499` and a raw Petrov-style scale of roughly `+216...+360 MHz` for the
assumed radial interval. Substituting the negative branch into Eq. (25) would
instead give `w` near `-0.040`, which is not a physical result because Petrov
defines `w` as the weight of the admixed configuration. It shows that the
single-channel ansatz cannot be inverted on that branch; it does **not** imply
a negative admixture or a sign-flipped quadrupole constant. The positive raw
scale is therefore conditional on the positive-`G_parallel` branch, while the
other branch leaves this route without a quantitative prediction. Finally, the
positive raw value is not yet the B&C Eq. (9.52) parameter because the factor
and sign conversion remains open.

[claim 4] The quoted `+/-1000 MHz` and `+/-100 MHz` are calibrated uncertainty
estimates.

- **verdict: STOP**
- **classification: unsupported**
- **source:** Petrov 2018, printed pp. 3-4, Sec. IV and Eqs. (24)-(25); Ng et
  al. 2022, printed p. 5.
- **reason:** Petrov describes the CCSD(T) calculation behind HfF+
  `eQq0=-2100 MHz` and notes an approximately 88 MHz triples contribution, but
  does not attach an uncertainty to that constant. Its `eQq2` is a one-channel
  model, also without an uncertainty. Ng's experimental error bars apply to
  `|g_F|` and to each conditional `G_parallel` branch, not to the ThF+ electric
  quadrupole constants. No inspected source assigns probability coverage to the
  assumed `R_el` or radial intervals. Their half-widths therefore cannot be
  propagated as ordinary standard uncertainties.

## Permitted use and exact missing bridge

The retained values may be used to test numerical stability and to display how
spectra respond to GHz-scale `q=0` and few-hundred-MHz raw `q=2` inputs. Any
result using them must identify them as sensitivity placeholders, sweep the
unresolved sign and normalization branches when those affect the conclusion,
and avoid presenting the output as a prediction for 229ThF+.

A quantitative replacement requires molecule-specific ThF+ electronic EFG
matrix elements `T^2_0(grad E)` and `T^2_{+/-2}(grad E)` in the same molecular
axis, phase, and Racah-tensor convention as the B&C Eqs. (9.52)-(9.53) kernel.
For reuse of Petrov's definitions, an author or erratum must instead resolve
both whether Eq. (19) intended the scalar-product `q,-q` pairing and whether
`sqrt(2 pi/5)` in Eqs. (19), (22), and (23) is intentional. A ThF+-specific
calculation must also replace the magnetic-HFS proxy `R_el`, the assumed 6d/5d
radial interval, and the single-channel mapping from `G_parallel` to `w` before
an uncertainty bar can be calibrated.

summary: 4 total claims; 1 OK, 0 WRONG, 0 NEEDS-CONTEXT, 3 STOP
