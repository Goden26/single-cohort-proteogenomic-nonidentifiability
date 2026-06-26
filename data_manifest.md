# Public Data Manifest

This repository does not version the raw public input files because they total several gigabytes and are already available from public archives. The analysis scripts expect the raw files in `Data/` at the project root. The parsed matrices and result tables needed to verify the reported numbers are included under `Study/results/`.

## Primary liver proteogenomics

- Gao 2019 HBV-related HCC proteogenomics: `gao_S1.xlsx`, `gao_S2_cnv.xlsx`; source Cell supplementary files for Gao et al. 2019.
- Dong 2022 iCCA proteogenomics: `dong_S1.xlsx`, `dong_S2_cnv.xlsx`, `dong_S5_sub.xlsx`; source Cancer Cell supplementary files for Dong et al. 2022.

## TCGA and independent RNA cohorts

- TCGA-LIHC and TCGA-CHOL TPM, clinical, and survival tables from the UCSC Xena GDC hub.
- GSE241466 raw expression table for intermediate-cell cHCC-CCA positioning.
- GSE231854 FPKM table for independent cHCC-CCA replication.

## DepMap and pathway/network resources

- DepMap 24Q4 CRISPRGeneEffect, ExpressionTPM, and Model annotation files from figshare.
- OmniPath kinase-substrate network snapshot stored as `omnipath_ksn.parquet`.

## Kidney generality check

- ccRCC proteome matrix from iProX IPX0001962000.
- pRCC MaxQuant proteinGroups table from PRIDE PXD013523.
- TCGA-KIRC and TCGA-KIRP TPM tables from UCSC Xena GDC hub.

## Reproduction

From the project root:

```bash
cd Study
bash scripts/setup_venv.sh
bash scripts/download_data.sh
bash scripts/dl_kidney.sh
./venv/bin/python3 scripts/20_parse.py
bash scripts/run_all.sh
```

`Study/results/parsed/manifest.json` records the dimensions of the parsed matrices included in this release.
