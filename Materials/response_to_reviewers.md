# Response to Reviewers — c1 manuscript (major revision)

We thank the four reviewers and the Editor. The panel identified a factual error and several construct-validity problems; addressing them substantially changed our conclusions, and we have rewritten the manuscript accordingly (now a methods/cautionary-benchmark paper). All new analyses are in `scripts/70_revision_analyses.py` → `results/revision_stats.json`. Point-by-point below; reviewer-facing changes are reflected in `manuscript_c1_revised.md`.

## Must-fix 1 — Headline gene-list factual error (R1/R2/R3, critical)
**Agreed; this was a real error.** The first draft named EGFR and ILK as RNA-concordant; both fail our own prespecified rule (log2FC > 0.5): EGFR log2FC = 0.379, ILK = 0.440, and both are flagged `rna_concordant_iCCA = False` in `results/icca_selective_kinases.csv`. **More importantly, the entire "iCCA-selective" framing was a comparator artifact (Must-fix 2), so we no longer present any "concordant 9/11" headline.** The revised paper reports the corrected analysis (Must-fix 2) instead. The CSV-vs-text discrepancy is gone because the claim it supported has been withdrawn.

## Must-fix 2 — "iCCA-selective" must be HCC-vs-iCCA, with FDR (R4 critical; R1/R2/R3 major)
**Agreed and implemented.** We recomputed selectivity as the **direct iCCA-line-vs-HCC-line contrast** with BH-FDR (`scripts/70`, section A). Of the original 11 "vs-all-other-lineages" hits: **only 2 survive FDR < 0.1**, **7 are pan-lineage-essential** (Chronos in other lineages < −0.5), and **3 (GRB2, CCND1, CRKL) are selective for *both* liver lineages**. We now state plainly that the "vs-all-other-lineages" comparator measures liver-vs-non-liver essentiality, and we present the result as a **cautionary artifact**, not a target map. The genuinely iCCA-vs-HCC-preferential set (direct FDR < 0.1) is small and exploratory: TFRC, SLC3A2, TRPM7, NHEJ1 (new Figure 3, `results/rev_*.csv`).

## Must-fix 3 — Absolute Chronos / essentiality floor; soften "druggable dependency" (R1/R3/R4)
**Implemented.** We annotate pan-essential and absolute essentiality per gene, report absolute Chronos for every hit, and removed "druggable dependency" language. CDK1 (Chronos ≈ −2.5 in both lineages) is now explicitly named as pan-essential; EGFR's marginal absolute dependency (−0.56) is stated. The exploratory survivors are framed as hypotheses, with TFRC/TRPM7 noted as preferential-but-pan-essential.

## Must-fix 4 — In-sample AUC = 1.000 is circular (R3 major; R1 marker circularity)
**Agreed and implemented.** We report the in-sample 0.99985 only as an acknowledged-optimistic sanity check, and now lead with a **stratified 5-fold cross-validated AUC = 0.964** (200-fold permutation null, p < 0.005) and a **non-marker hold-out AUC = 0.966** (markers removed). This directly answers the marker-circularity concern: the axis separates HCC from iCCA at CV-AUC 0.966 even without the lineage-marker genes.

## Must-fix 5 — n=8 i-CHC over-interpreted (R2/R3/R4)
**Agreed and implemented.** We de-emphasize the ordered-trend p = 4×10⁻²³ (now explicitly attributed to the pure-lineage arms, not the 8 middle samples), report **Cliff's δ** (−0.86 vs iCCA, p = 2.3×10⁻⁵; 0.41 vs HCC, p = 0.045) and a **bootstrap CI** for the i-CHC median (−0.21 to 0.00), and **disclose in the main text** the 25→8 sample-selection sensitivity that reversed the earlier (contaminated) result. We downgrade "genuine intermediate state" to "consistent with an intermediate transcriptomic phenotype, pending replication," and use "i-CHC" rather than "cHCC-CCA" for the headline.

## Must-fix 6 — Gene-class mislabeling; EGFR clinical overreach (R2 major; R1/R4)
**Agreed and implemented.** We state that SLC7A5/SLC3A2 are the LAT1 transporter, GRB2/CRKL adaptors, CCND1 a cyclin, PPP1R12A a phosphatase subunit, and ILK a pseudokinase, and we note the OmniPath "enzyme" universe includes non-kinases. The set is no longer called "kinases." EGFR is no longer presented as iCCA-selective or as "the canonical actionable axis"; we note the negative EGFR-inhibitor BTC trial history and reserve "actionable" for FGFR2/IDH1. The FGFR1 example is restricted to a generic expression≠dependency lesson with an explicit FGFR1≠FGFR2 statement (the actionable fusion program is not assayed in DepMap).

## Must-fix 7 — Reframe the negative result; positive control; CIs (R1/R3/R4)
**Agreed and implemented.** We reframe the proteome/phospho result as **non-identifiability under perfect aliasing** (no harmonization can separate aliased contrasts) and demote the tautological z-scoring AUC 0.506. We add a **non-aliased positive control**: DepMap-line vs TCGA iCCA−HCC direction (ρ = 0.015, n = 19,193), showing per-gene cross-dataset concordance is intrinsically weak — so a near-zero aliased ρ is uninformative by construction. We add a **bootstrap CI** for the proteome ρ (−0.18 to −0.13) and a **binomial sign test** (42%, p = 2×10⁻¹²; phospho 45%, p = 0.17, NS). We softened "anti-correlated" and, lacking a specific published manifold to refute, replaced "debunking a fashionable program" with "pre-empting a tempting but flawed analysis." HBV etiology is now named as a second variable aliased with cohort.

## Must-fix 8 — References, cohort N, reproducibility (all reviewers, minor)
**Implemented.** All 15 references were independently verified (PMIDs added); the placeholder line and the unsupported in-text "Comms Biol 2023" cite were removed. The Abou-Alfa/FIGHT-202 citation is now used only for the FGFR2 program. The Dong cohort N is reconciled in Methods (262 patients; 214 proteome/phospho used; 255 mRNA). We note the TCGA-CHOL histology caveat and state multiple-testing handling for every selection. The released `scripts/` reproduce all reported numbers (`70_revision_analyses.py` regenerates `revision_stats.json`).

## Editorial (titles/tone)
Title shortened and de-editorialized; "the honest route" and competitor-dishonesty framing removed; AUC reported as 0.96 (CV) with the in-sample 0.99985 labeled optimistic.

---
**Summary of changed conclusions:** the original positive headline (an iCCA-selective druggable kinase axis) did not survive a fair test and has been withdrawn. The paper now stands on (1) a cross-validated, marker-independent lineage axis with preliminary i-CHC positioning and (2) a cautionary identifiability benchmark for single-cohort proteogenomic integration, with a small exploratory transporter signal offered as a hypothesis.

---

# Response to Re-review (Stage 3', minor revision — scores 79/86/82)

We thank the panel; the re-review confirmed the round-1 fixes landed and that the LAT1 finding survives on its non-circular legs. All residual items are addressed (no new data needed beyond the leave-out test, `results/leaveout_stats.json`).

- **Circularity of the dependency leg (BIOSTAT, major).** Agreed and now disclosed in text: the dependency competitive-permutation signal (p=0.0012) is carried by the LAT1 discovery genes and collapses to **p=0.18** when SLC3A2+SLC7A5 are removed (ΔChronos −0.059→−0.011). We no longer call it "three independent axes." We rest the claim on the **selection-independent expression axis**, which survives leave-discovery-genes-out (p=1.2×10⁻¹⁰→1.6×10⁻¹⁰, 83% of remaining genes iCCA-up), plus the clinical anchor. Methods, Results, abstract, Figure 4 legend and Conclusions all updated.
- **Make expression load-bearing (moderate).** Done — stated explicitly as the selection-independent evidence.
- **Post-hoc programs / multiplicity (moderate).** Disclosed: programs curated *a priori* from canonical membership (not data-derived); no multiplicity correction; analysis labeled exploratory.
- **Dependency is the weakest/novel leg (HBP-ONC, moderate).** Added: expression and clinical legs corroborate established LAT1 biology, while the lineage-preferential dependency is the study's own modest observation (line MWU p=0.045; Chronos −0.23 vs −0.17, above the −0.5 floor).
- **Figure 4 MWU 0.045 vs 0.096 (minor, multiple reviewers).** Fixed — each MWU labeled by gene set (AA-alone 0.045; combined 0.096); combined permutation 0.0024 defined.
- **2-of-11 vs 4-survivors count (VERIFY-ED, major).** Reconciled: 2 of the 11 all-other hits survive the direct test (SLC3A2, TRPM7); the direct contrast independently flags 4 (adding TFRC, NHEJ1).
- **Promotional headers (HBP-ONC, minor).** Softened ("survives triangulation"→"recovered by the comparator-correct pipeline"; "the one signal that strengthens"→"corroborated rather than dissolved").
- **9/11 artifact phrasing (VERIFY-ED, minor).** Marked explicitly as the refuted first-draft headline.
- **RESULTS.md superseded banner (HBP-ONC, minor).** Added at the top so the released artifact cannot be quoted out of context.
