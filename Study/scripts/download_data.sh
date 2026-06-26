#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY="$(dirname "$HERE")"
BASE="$(dirname "$STUDY")"
DATA="$BASE/Data"
mkdir -p "$DATA" "$STUDY/logs"
cd "$DATA"
LOG="$STUDY/logs/download.log"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)"
echo "download run start" > "$LOG"
dl() { # name url
  local name="$1" url="$2"
  if curl -sSL -A "$UA" --retry 3 --retry-delay 2 -o "$name" "$url" 2>>"$LOG"; then
    local sz=$(stat -f%z "$name" 2>/dev/null); local ft=$(file -b "$name" | cut -c1-40)
    echo "OK   $name  ${sz}B  [$ft]" | tee -a "$LOG"
  else
    echo "FAIL $name  ($url)" | tee -a "$LOG"
  fi
}
# Gao 2019 HBV-HCC
dl gao_S1.xlsx       "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419310037-mmc1.xlsx"
dl gao_S2_cnv.xlsx   "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419310037-mmc2.xlsx"
# Dong 2022 iCCA
dl dong_S1.xlsx      "https://ars.els-cdn.com/content/image/1-s2.0-S1535610821006590-mmc2.xlsx"
dl dong_S2_cnv.xlsx  "https://ars.els-cdn.com/content/image/1-s2.0-S1535610821006590-mmc3.xlsx"
dl dong_S5_sub.xlsx  "https://ars.els-cdn.com/content/image/1-s2.0-S1535610821006590-mmc6.xlsx"
# DepMap 24Q4
dl depmap_CRISPRGeneEffect.csv "https://ndownloader.figshare.com/files/51064667"
dl depmap_Model.csv            "https://ndownloader.figshare.com/files/51065297"
dl depmap_ExpressionTPM.csv    "https://ndownloader.figshare.com/files/51065489"
# TCGA Xena
dl TCGA-LIHC.tpm.tsv.gz       "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_tpm.tsv.gz"
dl TCGA-LIHC.clinical.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.clinical.tsv.gz"
dl TCGA-LIHC.survival.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.survival.tsv.gz"
dl TCGA-CHOL.tpm.tsv.gz       "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CHOL.star_tpm.tsv.gz"
dl TCGA-CHOL.clinical.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CHOL.clinical.tsv.gz"
dl TCGA-CHOL.survival.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CHOL.survival.tsv.gz"
# cHCC-CCA
dl GSE241466_iCHC.txt.gz      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241466/suppl/GSE241466_Raw_data.txt.gz"
dl GSE231854_fpkm.txt.gz      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231854/suppl/GSE231854_gene_fpkm.txt.gz"
echo "DOWNLOAD_DONE total:" | tee -a "$LOG"; du -sh . | tee -a "$LOG"
