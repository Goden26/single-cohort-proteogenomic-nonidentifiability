#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY="$(dirname "$HERE")"
BASE="$(dirname "$STUDY")"
DATA="$BASE/Data"
mkdir -p "$DATA" "$STUDY/logs"
cd "$DATA"
UA="Mozilla/5.0"; LOG="$STUDY/logs/kidney2_dl.log"; echo "start" > "$LOG"
dl(){ curl -sSL -A "$UA" --retry 3 -o "$1" "$2" 2>>"$LOG" && echo "OK $1 $(stat -f%z "$1")B" | tee -a "$LOG" || echo "FAIL $1" | tee -a "$LOG"; }
dl ccRCC_proteome.xlsx  "https://ndownloader.figshare.com/files/31799282"
dl TCGA-KIRC.tpm.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.star_tpm.tsv.gz"
dl TCGA-KIRP.tpm.tsv.gz  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRP.star_tpm.tsv.gz"
echo "K2_DONE" | tee -a "$LOG"
