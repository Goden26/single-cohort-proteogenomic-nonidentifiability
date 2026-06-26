#!/usr/bin/env bash
# Top-level reproduction script: regenerates every released stat file, Table 2, the figures,
# and the inline-number JSON, then rebuilds the submission docx/pdf. Portable (no hardcoded paths):
# resolves the Study dir from this script's own location. Assumes results/parsed/ exists
# (created by 20_parse.py) and the raw inputs are in ../Data.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"          # .../Study/scripts
STUDY="$(dirname "$HERE")"                                     # .../Study
PY="$STUDY/venv/bin/python3"
cd "$STUDY"
for s in 30_rna_manifold 40_proteome_manifold 50_kinase_activity 60_depmap_rna_kinases \
         70_revision_analyses 71_fig3_revised 80_transport_deepdive 81_leaveout \
         90_improvements 91_harmonization_trap 92_identifiability_diagnostic \
         93_kidney_generality 94_inline_numbers; do
  echo ">>> running $s"
  "$PY" "scripts/$s.py"
done
echo ">>> rebuilding submission docx/pdf"
"$PY" "scripts/build_submission.py"
echo "ALL_DONE — stat files, Table 2, figures, and submission package regenerated."
