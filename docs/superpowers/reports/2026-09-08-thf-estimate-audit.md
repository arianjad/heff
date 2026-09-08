# ThF+ parameter estimate audit

Status: in progress, 2026-09-08. Start 5db8df7.

Requested work: repair missing 232Th TOML provenance and assess validity of
current estimates/placeholders against primary sources. No new physics kernels.

Metadata repair: regression first failed for missing native Param metadata;
TOML now copies thf_v1 value/unit/status/uncertainty/source/isotope/convention/note.
Numeric values and conventions are unchanged.

Audit scope: F spin-rotation scale/uncertainty; inherited isotope constants;
229Th hyperfine and moment inputs; transferred quadrupole estimates/normalization;
227Th Schmidt placeholders and zero-held Th spin-rotation. Primary sources are
being reread; earlier reports are navigation, not acceptance evidence.

## Shared inputs: verified findings

The canonical Zotero B&C PDF (storage CKZKCGXY, 1045 pages, 2003 title/copyright)
was checked after reading its contents. Relevant pages were rendered and read.
Table8.12, printed481/PDF513, gives CsF F-spin rotation15.1kHz and B=.183782cm-1.
Scaling with the current ThF B gives19.936385kHz. This verifies the arithmetic,
not transferability. Printed421/PDF453 discusses opposite-sign first/second-order
contributions; the earlier p453 reference had omitted that this was PDF numbering.
No quantitative factor-three bound or ThF sign follows. The coefficient remains
20kHz for exploratory sensitivity; its note no longer claims that error bound.

B&C Eq7.199, printed345/PDF377, gives Y_kl ~ mu^(-(k+2l)/2) U_kl.
Using mass numbers only to size the effect (not an isotope fit), B proportional
mu^-1 gives +7.213712MHz for229 and +12.128782MHz for227 relative to232.
D proportional mu^-2 gives +.007733/.013006kHz. Vibrational averaging and
Born-Oppenheimer breakdown limit this approximation. We did not replace B0 with
these diagnostic values. Unscaled shared B0,D0,omega_ef,A_par,d_mf,G_par in odd
isotopes are now explicitly estimates, with source error bars kept in notes and
no claimed target-isotope uncertainty. Numerical parameters remain unchanged.

Source discovery included read-only Zotero title/attachment queries and current
web searches for ThF nuclear spin rotation and isotope-dependent rotational
constants. Ng2022 (https://arxiv.org/abs/2202.01346) and Petrov/Skripnikov2025
(https://arxiv.org/abs/2503.02840) are the relevant primary ThF records found;
this search is not proof that no later calculation exists. Th-specific reviews
are linked below when complete.

Metadata verification:32 model/isotope checks passed after232repair;52 focused
model/isotope/parameter checks passed after transfer-label clarification.
