# Supplementary material

**Figure S1. Ground-truth calibration of the anchor check.** Using real TCGA-LIHC expression as substrate, a known 200-gene signature was spiked into a synthetic biological contrast, a synthetic batch effect was added to a synthetic cohort, and the aliasing of the contrast with the cohort was swept across 11 levels from 0 to 1 (20 repeats per level). Concordance is measured against a separate noisy orthogonal anchor (concordance to the truth about 0.20), not the ground truth. As aliasing increases, true-signal recall stays high while the false-discovery proportion rises to 0.93 and noisy-anchor concordance shifts from 0.105 to 0.042. Pooling the 220 runs, noisy-anchor concordance predicts false-discovery proportion (Spearman -0.68, level-block bootstrap 95% CI -0.72 to -0.43; -0.68 under a multiplicative batch model), while the within-level association is near zero.

**Table S1. Genes passing the direct iCCA versus HCC dependency contrast (FDR below 0.1).** Negative dChronos means more essential in iCCA lines. TFRC and TRPM7 remain pan-essential, so they are preferential rather than iCCA-exclusive. FDR reflects rank separation rather than effect magnitude; NHEJ1 has no transport rationale and is likely an incidental survivor.

| Gene | dChronos (iCCA minus HCC) | FDR | Chronos iCCA | Chronos HCC | RNA log2FC (iCCA minus HCC) |
|---|---:|---:|---:|---:|---:|
| TFRC | -0.548 | 0.063 | -0.989 | -0.441 | 0.642 |
| TRPM7 | -0.521 | 0.075 | -1.052 | -0.531 | 0.633 |
| SLC3A2 | -0.369 | 0.063 | -0.725 | -0.355 | 0.759 |
| NHEJ1 | -0.159 | 0.045 | -0.171 | -0.012 | 0.460 |
