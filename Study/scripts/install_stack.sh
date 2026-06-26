#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY="$(dirname "$HERE")"
mkdir -p "$STUDY/logs"
LOG="$STUDY/logs/pip.log"
echo "pip install run starting" > "$LOG"
python3 -m pip install --upgrade pip >>"$LOG" 2>&1
for p in scipy scikit-learn seaborn pycombat decoupler mofapy2 anndata gseapy upsetplot adjustText openpyxl; do
  echo "=== installing $p ===" >>"$LOG"
  if python3 -m pip install "$p" >>"$LOG" 2>&1; then
    echo "RESULT $p OK" | tee -a "$LOG"
  else
    echo "RESULT $p FAIL" | tee -a "$LOG"
  fi
done
echo "=== final versions ===" >>"$LOG"
python3 -c "import importlib.util as u;[print(p, 'OK' if u.find_spec(p) else 'MISSING') for p in ['scipy','sklearn','seaborn','pycombat','decoupler','mofapy2','anndata','gseapy','openpyxl']]" 2>>"$LOG" | tee -a "$LOG"
echo "INSTALL_DONE" | tee -a "$LOG"
