# Statement of Significance

**Problem:** Integrating single-cohort-per-condition multi-omics datasets across conditions is routine in biomedical informatics, but when the biological contrast of interest is perfectly aliased with cohort and batch the integration is non-identifiable, and a standard harmonization can manufacture large, plausible, but false differential lists.

**What is already known:** Batch-free anchor and negative-control checks (removing unwanted variation, surrogate variable analysis, guided principal-component analysis, and the routine mRNA-protein concordance check) exist, but they are routinely omitted in the single-cohort-per-condition setting, and how they degrade under perfect aliasing has not been calibrated.

**What this paper adds:** A formal non-identifiability framing for this integration design, a ground-truth calibration of how a standard anchor check degrades with the degree of aliasing, and a reproducible two-disease benchmark, demonstrated where the merge is most tempting (liver and kidney cancer proteogenomics).
