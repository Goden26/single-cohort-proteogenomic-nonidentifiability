Anh Tuan Quan, MD, MSc
International Graduate Program in Medicine, College of Medicine, Taipei Medical University, Taipei, Taiwan
ORCID: https://orcid.org/0009-0008-2281-970X

Dear Editor-in-Chief,

I am pleased to submit the manuscript "Merging single-cohort proteogenomic datasets is non-identifiable for cross-tumor contrasts: a calibrated benchmark in liver and kidney cancer" for consideration as a Research Paper in the Journal of Biomedical Informatics. This is an original methodological contribution on the identifiability of multi-omics data integration, and it offers a reusable evaluation procedure for the common single-cohort-per-condition setting.

The work fits the scope of JBI, which centers on methodological rigor, reproducibility, and the critical evaluation of informatics approaches. Its primary contribution is a formal non-identifiability result for a widespread integration design: when each condition is represented by a single cohort, the biological contrast is perfectly aliased with cohort and batch, so no procedure can separate biology from technical variation. We make this concrete with a harmonization trap. Standard quantile harmonization of two merged proteomes flagged 63 percent of liver proteins and 76 percent of kidney proteins as significant, none of which validated against a batch-free RNA anchor. We further show that a dedicated batch corrector cannot rescue this: ComBat without a protected covariate erases the contrast, while protecting the condition leaves a rank-deficient design that standard implementations reject, an impossibility that follows directly from the aliasing.

We are explicit about what is and is not new. The batch-free anchor check is established practice, drawing on negative-control and orthogonal-validation methods and the routine mRNA-protein concordance check. Our contribution is the framing, a ground-truth calibration of how that check degrades across aliasing regimes via a controlled titration, and a reproducible two-disease benchmark (a powered liver case and a corroborating kidney case) that gives practitioners a concrete way to judge when a cross-cohort contrast is trustworthy. The worked example also yields surviving positive deliverables computed only from batch-controlled layers: a cross-validated lineage axis, a lesson on choosing the correct comparator for dependency analysis, and a clinically anchored amino-acid-transport (LAT1) dependency that aligns with a positive randomized phase-2 trial.

The study is fully reproducible: all data and code are public, and every reported number is regenerable from the released scripts. The manuscript is original, has not been published, and is not under consideration elsewhere. All author declarations, the data and code availability statement, and the disclosure of AI-assisted tool use are provided within the manuscript.

Thank you for your consideration.

Sincerely,

Dr. Anh Tuan Quan
International Graduate Program in Medicine, College of Medicine, Taipei Medical University
No. 250 Wuxing Street, Xinyi District, Taipei 11031, Taiwan
Tel.: +84 987 750 775
Email: tuan.quan.md@outlook.com
