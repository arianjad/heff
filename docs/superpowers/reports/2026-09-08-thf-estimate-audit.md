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
