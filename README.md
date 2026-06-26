# Liver cancer genes study: START HERE (project README)

Read this one file to understand the whole project: what it is, what it found, where everything is, how to reproduce it, and what is left to do. Last updated 2026-06-25.

---

## 1. TL;DR

A dry-lab, open-data-only computational study (project codename "c1") that began as a search for a max-novelty liver-cancer gene study and, after real analysis and four rounds of simulated peer review, became an honest methods/benchmark paper. Final one-line thesis:

> Merging single-cohort proteogenomic datasets is non-identifiable for cross-tumor contrasts. We demonstrate this in liver (HCC vs iCCA) and kidney (ccRCC vs pRCC), calibrate on ground truth how a standard batch-free anchor check degrades with aliasing, and provide a cross-validated HCC-iCCA lineage axis plus a small clinically-anchored LAT1 dependency as the surviving positive results.

Status: manuscript finalized and submission-ready. Honest tier estimate: applied / methods-caution journal (Q1-applied or Q2), roughly 55 to 70 percent acceptance odds. It is NOT a top-methods novelty paper (see Section 8). It IS a degree-qualifying first-author SCI paper.

Author: Anh Tuan Quan, MD, MSc. International Ph.D. Program in Medicine, Taipei Medical University. Email tuan.quan.md@outlook.com. ORCID 0009-0008-2281-970X.

---

## 2. Folder structure

```
Liver cancer genes study/
├── README.md                 <- this file
├── Data/        (3.8 GB)      raw input datasets (26 files), all public
├── Study/       (881 MB)      the analysis pipeline:
│     ├── scripts/             Python modules 20..93 + shell scripts (the pipeline)
│     ├── results/             parsed matrices (parquet), stat JSONs, result CSVs
│     ├── figures/             fig1..fig6 (PNG, the manuscript figures)
│     ├── logs/                download / install logs
│     └── venv/                Python 3.13 env (STALE PATHS after the folder move;
│                              recreate with scripts/setup_venv.sh if rerunning)
├── Submission package/ (1.7M) THE FINAL manuscript, ready to submit:
│     ├── manuscript.md / .docx
│     ├── title_page.md / .docx
│     ├── cover_letter.md / .docx
│     ├── declarations.md / .docx
│     └── submission_package.pdf   (combined: title+cover+manuscript+declarations)
└── Materials/   (40 KB)       supporting records:
      ├── RESULTS.md           lab-notebook results writeup (NOTE: its top
      │                        sections 1-3 are SUPERSEDED, see its own banner)
      ├── response_to_reviewers.md   point-by-point responses (review rounds 1-2)
      └── _tnr_ref.docx        Times New Roman pandoc reference template
```

All intermediate and superseded manuscript versions were deleted; only the final is kept (in Submission package/).

---

## 3. What the paper claims (genuine contributions, honest)

1. Formal non-identifiability: when two tumor types each come from one proteomic cohort from a different lab, tumor type is perfectly aliased with batch, so the between-type contrast cannot be separated from technical variation by any normalization.
2. A concrete harmonization-trap demonstration: a standard quantile harmonization of the two largest open liver proteogenomic cohorts flags 3,796 of 6,067 proteins (63 percent, including 407 named signaling enzymes) as "significant" iCCA-vs-HCC, none of which validate against a batch-free RNA axis. The same failure recurs in kidney (76 percent significant, anchor concordance rho = 0.004).
3. A ground-truth calibration: using real TCGA-LIHC expression with a spiked signature and a noisy orthogonal anchor, anchor concordance predicts the false-discovery load across 11 aliasing levels (per-replicate Spearman -0.68, 95 percent CI -0.74 to -0.61, robust to additive and multiplicative batch). We calibrate a known check (RUV, gPCA, sva, the routine mRNA-protein concordance check); we do NOT claim a new method.
4. A reusable, cross-validated HCC-iCCA lineage axis (nested CV-AUC 0.966; 0.966 without marker genes), along which intermediate-cell combined HCC-CCA (i-CHC) sits as an intermediate phenotype, replicated in an independent cohort (Cliff's delta 0.58, P = 1.7e-5).
5. A functional-genomics caution plus the one surviving positive: an apparent iCCA-selective kinase program (EGFR/GRB2/CRKL) is a comparator artifact (7 of 11 hits pan-essential); only LAT1 amino-acid transport survives a fair direct iCCA-vs-HCC test, is over-expressed in iCCA tumors, and aligns with the positive randomized phase-2 trial of the LAT1 inhibitor nanvuranlat.

---

## 4. Datasets used (Data/, all open)

- Gao 2019 HBV-HCC proteogenomics: NODE OEP000321 (159 tumors; proteome, phosphosite). PMID 31585088.
- Dong 2022 iCCA proteogenomics: NODE OEP001105 (214 proteome/phospho, 255 mRNA). PMID 34971568.
- TCGA-LIHC (371) and TCGA-CHOL (35) RNA-seq: UCSC Xena GDC hub (batch-free anchor; same pipeline).
- cHCC-CCA: GEO GSE241466 (8 true i-CHC of a 25-sample mixed series); independent replication GSE231854 (25 cHCC-CCA, 56 HCC).
- DepMap Public 24Q4: CRISPR Chronos + expression (24 HCC, 31 iCCA lines).
- Kidney generality: ccRCC iProX IPX0001962000 (Qu 2022, PMID 35440542), pRCC PRIDE PXD013523 (Al Ahmad 2019, PMID 31484429), with TCGA-KIRC/KIRP anchor.
- OmniPath kinase-substrate network (KSEA).
- NOTE: there is no "CPTAC-HCC" proteome; the real open liver proteomes are the NODE deposits above.

---

## 5. How to reproduce

```
cd "Liver cancer genes study/Study"
bash scripts/setup_venv.sh                  # recreate the Python 3.13 environment
./venv/bin/python3 scripts/20_parse.py      # parse proteome/phospho/clinical/DepMap
bash scripts/run_all.sh                     # regenerate stats, figures, and submission files
```
Script map (Study/scripts/): 20 parse; 30 RNA lineage axis; 40 proteome confound check; 50 KSEA; 60 DepMap+RNA kinases; 70 revision stats; 71 fig3; 80 transport deep-dive (LAT1); 81 leave-out; 90 independent replications (GSE231854, DepMap line expr); 91 harmonization trap (liver); 92 ground-truth titration v2; 93 kidney generality. Every reported number traces to a JSON/CSV in Study/results/.

The submission docs are rebuilt with pandoc + tectonic using Materials/_tnr_ref.docx as the Times New Roman reference. Helper to extract a single file from a huge remote zip via HTTP range reads is inside scripts/93 (used for the 7.5 GB pRCC zip).

---

## 6. Status and how to submit (what is left)

The manuscript is final. To actually submit, the author (not the analysis) must:
1. Pick a journal. Suggested order: GigaScience or PLOS Computational Biology (reproducible benchmark / negative-with-content friendly) -> Briefings in Bioinformatics (applied track) -> BMC Bioinformatics or NAR Genomics and Bioinformatics -> backup Scientific Reports. All are SCI-indexed (degree-qualifying).
2. Confirm the final journal submission form selected the subscription / non-open-access route for Journal of Biomedical Informatics. The code repository and Zenodo archive are already public: https://github.com/Goden26/single-cohort-proteogenomic-nonidentifiability and DOI 10.5281/zenodo.20888892.
3. Trim the abstract to the chosen journal limit (currently 456 words; some journals want 250 to 350).
4. Frame the cover letter around the general methods caution (already done) to avoid a desk-reject on "narrow/negative" grounds.

---

## 7. Project history (six simulated review rounds; the rigor is the point)

- Round 1: major revision. Caught a factual gene-list error and collapsed the original positive headline; the "iCCA-selective EGFR kinase axis" was a comparator artifact.
- Round 2: minor revision. Caught LAT1 dependency circularity (disclosed via leave-out test).
- Round 3: minor revision. All numbers re-verified; added the harmonization-trap demonstration.
- Round 4: major revision. Caught my novelty over-claim: the "diagnostic" is established practice; the titration was partly circular and the -1.0 statistic was near-mechanical. Fixed by walking back the framing, citing prior literature, and rebuilding the titration with a noisy separate-modality anchor.
- Round 6 (full re-review of the round-5 revision: EIC + 3 reviewers + Devil's Advocate): ACCEPT WITH MINOR REVISION, 84/100 (up from 68). Integrity pass clean (all 38 load-bearing numbers reproduce). Implemented the two new must-fixes: (a) disclosed the forced-pinv ComBat protect-tumour-type anchor rho +0.18 in-text and explained it as the protected covariate being re-imprinted (circular, not recovery); (b) made ALL scripts path-portable (the round-5 "all portable" claim was false; scripts 10/20/30/50/60/71/90/93 still had dead paths and 93 crashed) and verified script 93 regenerates the kidney numbers (76%, rho 0.004); added scripts/run_all.sh as the top-level reproduction runner. Minors: abstract carries the "across aliasing regimes / not per-replicate" qualifier; EGFR BH-FDR 0.625 now emitted to manuscript_inline_numbers.json; report 0.975 (not 0.98) for z-score full-LDA; softened "implementations abort" to the sva/pyComBat confounded-design guards.
- Round 5 (full panel: EIC + 3 reviewers + Devil's Advocate, major revision): implemented all eight major points. Code: nested the lineage-axis CV (selection+scaling re-fit per fold; AUC 0.964->0.966, robust); added ComBat baselines (no-covariate erases the contrast to AUC 0.50; protect-tumour-type is rank-deficient and standard ComBat aborts); added a four-way normalization sweep showing the negative proteome rho is harmonization-dependent (-0.21 raw to -0.004 per-cohort z); made the titration cluster-aware (pooled -0.68, level-block CI -0.72 to -0.43, within-level ~ -0.02, shared-truth caveat disclosed); unified the LAT1 leave-out seed so the 0.0024->0.20 collapse comes from one run; fixed the Figure 2 title and deleted the stale v1 titration CSV; made all scripts path-portable. Text: de-strawmaned the framing (prospective caution), added the etiology/grade confound caveat, bounded the kidney "across two cancers" claim. Regenerated figures + docx/pdf via scripts/build_submission.py.

Integrity catches along the way (all fixed): a contaminated cHCC-CCA cohort (used all 25 GSE241466 samples instead of the 8 true i-CHC); two fabricated reference PMIDs (corrected after PubMed verification); the novelty over-claim. This is why "run the real analysis and let it be adversarially reviewed" matters.

---

## 8. Honest caveats and gotchas (read before extending)

- The core idea (validate a cross-cohort contrast against a batch-free anchor) is ESTABLISHED practice, not novel. Do not re-inflate the framing to "a new diagnostic"; four reviewers rejected that unanimously. The novelty ceiling is real and incremental.
- "iCCA-selective" must be defined by the DIRECT iCCA-vs-HCC contrast, never "vs all other lineages" (which manufactures tissue-essentiality false positives). This is a load-bearing lesson of the paper.
- The kidney generality is asymmetric (ccRCC TMT n=232 vs pRCC label-free n=16); it shows the failure mode recurs, not a powered second calibration. Stated honestly in the text.
- LAT1 over-expression is a tumor-level phenomenon (it is NOT elevated in DepMap cell lines), so the CRISPR dependency leg is weak; the claim rests on tumor expression plus the nanvuranlat trial.
- Study/venv has stale absolute paths after the folder move; recreate it before rerunning.
- Formatting rule for any document for this author: Times New Roman, black, NO strange symbols (write rho, delta, e-notation, "approximately"; no Greek/superscripts/multiplication signs), and NO dashes (no em or en dash); keep ordinary hyphens only in identifiers like TCGA-LIHC and cHCC-CCA. Run the humanizer pass.

---

## 9. If you want a genuinely higher-novelty paper

Do NOT try to push this one further; the honest ceiling is reached. A higher-novelty, positive-result option from the same research line is project "c18" (germline liver-tissue cis-pQTL -> Mendelian randomization for HCC/iCCA risk), feasibility about 7/10, all open data. See the ranked topic menu at ~/Downloads/liver-cancer-gene-novelty-menu.md for c18 and other candidates.
