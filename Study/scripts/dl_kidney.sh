#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY="$(dirname "$HERE")"
BASE="$(dirname "$STUDY")"
DATA="$BASE/Data"
mkdir -p "$DATA" "$STUDY/logs"
cd "$DATA"
UA="Mozilla/5.0"; LOG="$STUDY/logs/kidney_dl.log"; echo "kidney dl start" > "$LOG"
dl(){ curl -sSL -A "$UA" --retry 3 -o "$1" "$2" 2>>"$LOG" && echo "OK $1 $(stat -f%z "$1" 2>/dev/null)B [$(file -b "$1"|cut -c1-30)]" | tee -a "$LOG" || echo "FAIL $1" | tee -a "$LOG"; }
dl ccRCC_proteome.xlsx       "https://ndownloader.figshare.com/files/31799282"
dl pRCC_maxquant.zip         "https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD013523/txt.zip"
dl TCGA-KIRC.tpm.tsv.gz      "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.star_tpm.tsv.gz"
dl TCGA-KIRP.tpm.tsv.gz      "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRP.star_tpm.tsv.gz"
echo "KIDNEY_DL_DONE" | tee -a "$LOG"
