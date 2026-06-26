# c1 — Real-data analysis results
**Pan-liver lineage axis + batch-controlled druggable-kinase map across HCC / iCCA / cHCC-CCA**
Analysis date: 2026-06-22. All data are public/open. Dry-lab only (Python 3.13 venv; no GPU).

> ⚠️ **SUPERSEDED NOTICE — read before quoting.** The "Headline findings" sections 1–3 below were the FIRST-DRAFT analysis. They did **not** survive peer review and are **retracted/corrected** in the "Stage 4 revision" and "Deep-dive" sections lower in this file (and in `manuscript_c1_revised.md`, the authoritative artifact). Specifically: the in-sample "AUC = 1.000" → cross-validated **0.964**; "genuine intermediate state" (cHCC-CCA) → "consistent with, n=8, preliminary"; the "iCCA-selective druggable kinase axis (EGFR/GRB2/CRKL/CCND1/LAT1)" was an **all-other-lineages comparator artifact** and is withdrawn — only an amino-acid-transport (LAT1) signal survives, mainly on batch-free expression. Do not cite sections 1–3 as findings.

## Datasets actually used (all verified, openly downloadable, no login)
| Cohort | Source / accession | N (used) | Data type |
|---|---|---|---|
| Gao 2019 HBV-HCC | Cell 2019, PMID 31585088; suppl. Table S1 (Elsevier CDN; raw at NODE OEP000321) | 159 tumors | proteome (6,478), phosphosite (26,418) |
| Dong/Lu 2022 iCCA | Cancer Cell 2022, PMID 34971568; suppl. Table S1 (raw at NODE OEP001105) | 214 (prot/phos), 255 (mRNA) | proteome (8,296), phosphosite (18,347), mRNA |
| TCGA-LIHC | UCSC Xena GDC hub (STAR TPM) | 371 primary tumors | RNA-seq (batch-free vs CHOL: same pipeline) |
| TCGA-CHOL | UCSC Xena GDC hub (STAR TPM) | 35 primary tumors | RNA-seq |
| cHCC-CCA (i-CHC) | GEO **GSE241466** (mixed series: 7 normal + 5 HCC + 5 iCCA + **8 i-CHC**) | **8 i-CHC used** | RNA-seq |
| DepMap 24Q4 | figshare (CRISPRGeneEffect, Model) | 24 HCC + 31 iCCA(IHCH) + 1 combined (KMCH1) + 1,122 other | CRISPR Chronos gene effect |
| Kinase–substrate network | OmniPath (Enzsub) | 1,657 kinases / 40,109 phospho-relations | KSEA reference |

## Headline findings (real numbers)

### 1. A batch-free HCC↔iCCA lineage axis is robust; cHCC-CCA is positioned by an ordered trend
- TCGA-LIHC vs TCGA-CHOL (same GDC STAR pipeline → **no cohort/batch confound**): the lineage axis separates HCC from iCCA **perfectly (AUC = 1.000)**; PC1 = 26% variance.
- cHCC-CCA (i-CHC, GSE241466; **8 true i-CHC samples**, after excluding the cohort's 7 normal + 5 HCC + 5 iCCA) on the cholangiocyte−hepatocyte single-sample marker score (platform-robust): median **HCC −0.165 < cHCC-CCA −0.082 < iCCA +0.126**; Kruskal p = 8.7e-21; **ordered trend p = 4.0e-23** (Spearman ρ = 0.46). i-CHC is **significantly intermediate — distinct from BOTH HCC (p = 0.045) AND iCCA (p = 2.3e-5)** — and **75%** of i-CHC samples fall in the HCC–iCCA intermediate band. → i-CHC is transcriptomically a genuine intermediate state. *(NOTE: an earlier run used all 25 GSE241466 samples and was contaminated by the cohort's normal/HCC/iCCA samples, which falsely flattened cHCC-CCA toward HCC; corrected at the Stage 2.5 integrity gate.)*

### 2. CENTRAL NEGATIVE RESULT — open single-cohort proteome/phospho are too confounded for a cross-tumor manifold
The Gao (HCC) and Dong (iCCA) proteome/phospho datasets are each from one cohort, so **cohort is perfectly confounded with tumor type**. Tested honestly:
- **Proteome:** within-sample percentile direction (iCCA−HCC) vs the batch-free TCGA-RNA direction → Spearman **ρ = −0.15** (i.e., NOT concordant); lineage-marker score iCCA vs HCC p = 0.12 (NS); sign agreement on strong genes 42% (< chance). Per-cohort z-scoring instead **erased** the contrast (LDA AUC 0.506).
- **Phospho kinase activity (KSEA, 237 kinases scored in both):** cross-tumor activity difference vs batch-free RNA → Spearman **ρ = −0.045, p = 0.49, sign agreement 45%** (≈ chance). The 134 "FDR-significant" cross-tumor kinases are dominated by batch.
→ **The fashionable "cross-cohort proteogenomic kinase manifold spanning HCC↔iCCA" is NOT achievable in a demonstrably biological (batch-robust) way from these open cohorts.** This reframes c1 from a positive manifold claim into a methods caution + a batch-controlled alternative.

### 3. A batch-controlled druggable-kinase map DOES yield coherent, actionable lineage biology
Using only single-platform / batch-controlled layers (TCGA RNA + DepMap):
- **TCGA RNA (batch-free):** 571 kinases lineage-differential (FDR < 0.05 & |log2FC| > 1) between HCC and iCCA.
- **DepMap (single platform):** **11 iCCA-selective kinase dependencies** (vs all other lineages; Δ Chronos < −0.2, p < 0.1), of which **9/11 are RNA-concordant** (iCCA-up in expression AND iCCA-selective dependency): **EGFR, GRB2, CRKL, CCND1, WNK1, TRPM7, SLC7A5, SLC3A2, ILK** (+ CDK1, PPP1R12A). This is a coherent **RTK→GRB2/CRKL/PTPN11→MAPK + LAT1 (SLC7A5/SLC3A2) amino-acid-transport** program — the canonical cholangiocarcinoma actionable axis (EGFR, SHP2, FGFR all clinically pursued in iCCA).
- **Direct HCC-line vs iCCA-line contrast:** HCC lines are preferentially dependent on **FGFR1, PTK2 (FAK), MTOR, PLK1, ITGB5**; note FGFR1 RNA is iCCA-high (log2FC +3.16) yet dependency is HCC-higher — a real expression/dependency discordance worth flagging.
- Concordance contrast that validates the approach: RNA+DepMap concordance for iCCA-selective kinases = **9/11 (82%)** vs proteome/phospho concordance ≈ chance (42–45%).

## Figures
- `figures/fig1_rna_lineage_manifold.png` — TCGA RNA manifold (AUC 1.0) + cHCC-CCA intermediate-trend violin.
- `figures/fig2_proteome_confound_check.png` — proteome cohort confound (PCA) + failed RNA concordance (ρ=−0.15).
- `figures/fig3_depmap_rna_kinases.png` — iCCA-selective druggable kinases (RNA log2FC vs DepMap selectivity) + per-line Chronos bars.

## Output tables
`results/kinase_lineage_integrated.csv` (1,606 kinases, all metrics), `results/icca_selective_kinases.csv`,
`results/hcc_selective_kinases.csv`, `results/direct_icca_gt_hcc_kinases.csv`, `results/direct_hcc_gt_icca_kinases.csv`,
`results/robust_phospho_kinases.csv`.

## Honest limitations
1. Cross-cohort proteome/phospho confound is the dominant limitation (see finding 2) — we do NOT use those layers for cross-tumor claims.
2. TCGA-CHOL (n=35) and cHCC-CCA (n=8 true i-CHC) are small; iCCA includes extrahepatic in some resources (we restricted to IHCH for DepMap).
3. cHCC-CCA evidence is RNA-only (no open phospho/proteome) → kinase ACTIVITY in cHCC-CCA cannot be measured directly; positioning is by RNA proxy.
4. DepMap lineage-selectivity vs "all other lineages" conflates liver-lineage essentiality; the direct HCC-vs-iCCA line contrast (n=24 vs 31) is the cleaner but lower-powered test.
5. KSEA used per-sample z of substrate sites (≥5 substrates); OmniPath KSN is literature-biased toward well-studied kinases.

## One-line thesis (reframed by the data)
*Open, batch-controlled transcriptomic + functional-genomic data define a robust HCC↔iCCA lineage axis along which cHCC-CCA (i-CHC) sits as a significant intermediate state (distinct from both HCC and iCCA), and nominate a concordant set of iCCA-selective druggable kinase dependencies (EGFR/GRB2/CRKL/CCND1/LAT1) — while single-cohort proteome/phosphoproteome data are shown to be too cohort-confounded to support the cross-tumor kinase-activity manifold that has been proposed.*

---

## Stage 4 revision (post peer-review) — honest re-analysis findings

Peer review (4 reviewers + EIC, decision = **major revision**) caught a factual error and construct-validity problems. Re-analysis (`scripts/70_revision_analyses.py`, `results/revision_stats.json`):

- **Gene-list error fixed:** the true RNA-concordant set (per the stated rule) is {GRB2, TRPM7, CCND1, CRKL, WNK1, PPP1R12A, SLC3A2, CDK1, SLC7A5}; EGFR (log2FC 0.38) and ILK (0.44) FAIL the 0.5 threshold and were wrongly named in the first draft.
- **A — "iCCA-selective" largely collapses under a fair test.** Redefined as the DIRECT iCCA-line-vs-HCC-line contrast with BH-FDR: of the original 11 "all-other-lineage" hits, only **2 survive FDR<0.1**, **7 are pan-liver-essential** (Chronos_other < −0.5), and **3 (GRB2, CCND1, CRKL) are flagged selective for BOTH lineages** → liver-essential, not iCCA-specific. The genuinely iCCA-vs-HCC-preferential set is small: **TFRC, TRPM7, SLC3A2, NHEJ1** (FDR 0.04–0.08) — a **nutrient/iron (TFRC), amino-acid-transport (SLC3A2/LAT1), and ion-channel-kinase (TRPM7)** theme, NOT the textbook EGFR axis. EGFR/GRB2/CRKL/CCND1 are NOT iCCA-vs-HCC selective.
- **B — lineage axis is robust, NOT marker-circular.** Cross-validated AUC = **0.964** (all genes) and **0.966 on non-marker genes** (permutation p≈0). The in-sample AUC 1.000 was optimistic; the honest CV value stands and survives marker hold-out.
- **C — cHCC-CCA (i-CHC, n=8): "consistent with intermediate."** Strongly distinct from iCCA (Cliff's δ=−0.86, p=2.3e-5) but only moderately/marginally from HCC (δ=0.41, p=0.045); median −0.082, bootstrap CI [−0.21, 0.00]. The trend p=4e-23 mainly reflects the pure-lineage gradient. Preliminary (n=8, single cohort, RNA-only).
- **D — negative result reframed as design identifiability.** Proteome ρ=−0.15 (CI [−0.18,−0.13]; sign-agreement 42%, binom p=2.3e-12 = significantly DISCORDANT); phospho-kinase ρ=−0.045 (sign 45%, binom p=0.17, NS). A non-aliased positive control (DepMap-line vs TCGA iCCA−HCC direction) gave ρ=0.015 (~0), showing per-gene cross-dataset concordance is intrinsically weak — so the negative result is best framed as **perfect-aliasing non-identifiability** (no harmonization can separate tumor type from cohort), with the ρ values as illustrations, not as a clean batch-specific empirical claim. (HBV etiology is a second variable aliased with cohort, since Gao is HBV-HCC only.)
- **E — gene-class corrected:** the OmniPath "enzyme" universe includes non-kinases — SLC7A5/SLC3A2 (LAT1 transporter), GRB2/CRKL (adaptors), CCND1 (cyclin), PPP1R12A (phosphatase subunit), ILK (pseudokinase). The "kinase" headline is inaccurate; correct framing = "signaling/druggable dependencies."

**Net effect on the thesis:** the robust, defensible contributions are (1) the methods/identifiability caution and (2) the CV-validated lineage axis + preliminary cHCC-CCA positioning. The originally headline "iCCA-selective druggable kinase axis (EGFR/GRB2/CRKL/CCND1/LAT1)" does NOT survive a fair direct-contrast test and must be retracted/deflated; only a small, exploratory nutrient/ion-transport iCCA-preferential signal (TFRC/SLC3A2/TRPM7) remains.

---

## Deep-dive (post-revision) — the surviving positive: amino-acid transport (LAT1) in iCCA

To strengthen (or refute) the exploratory TFRC/SLC3A2 signal, we tested iron-handling vs amino-acid-transport as PROGRAMS (`scripts/80_transport_deepdive.py`, `results/transport_deepdive_stats.json`):

- **Iron-handling program: does NOT generalize** — competitive permutation p = 0.14, batch-free TCGA expression p = 0.55 (47% genes iCCA-up). Only TFRC holds individually; iron-addiction is NOT supported as a program.
- **Amino-acid-transport program (LAT1 + glutamine + xCT): ROBUST across 3 independent axes:**
  - DepMap dependency preferentially essential in iCCA lines: competitive permutation **p = 0.0012**, line-level MWU p = 0.045 (iCCA mean Chronos −0.23 vs HCC −0.17); top genes TFRC, SLC3A2, SLC7A5, GLS, SLC2A1.
  - Batch-free TCGA expression coordinately iCCA-up: median log2FC **+0.79**, set-score MWU **p = 1.2×10⁻¹⁰**, 86% of genes iCCA-up.
  - Co-dependency module: TFRC–SLC3A2 Spearman ρ = 0.56; TFRC/SLC3A2/SLC7A5 cluster.
- **Clinical anchor (verified, real):** LAT1 (SLC7A5/SLC3A2) inhibitor **nanvuranlat (JPH203)** improved PFS in a randomized placebo-controlled **phase-2 BTC trial** (Furuse et al., Clin Cancer Res 2024, PMID 39058429; HR≈0.56); LAT1 overexpression is prognostic in BTC (Kaira 2013 PMID 24131658; Yanagisawa 2014 PMID 24890221). Caveats: phase-2 PFS (not OS/phase-3), BTC broadly (not iCCA-exclusive); SLC7A11/GLS ferroptosis strands are preclinical in CCA.

**Net:** the deep-dive converts the weak 4-gene exploratory signal into ONE genuine, triangulated, clinically-anchored positive — **iCCA preferentially depends on amino-acid transport (LAT1)** — while honestly rejecting the iron-program over-generalization. This is the "constructive" half of the comparator-rigor message: the same direct-contrast discipline that kills the EGFR artifact recovers LAT1, the one axis with a positive randomized trial. New: Figure 4; refs 16–18.
