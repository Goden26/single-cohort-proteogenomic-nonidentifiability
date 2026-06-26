#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY="$(dirname "$HERE")"
cd "$STUDY"
mkdir -p logs
LOG=logs/venv.log
echo "creating venv with python3.13" > "$LOG"
/usr/local/bin/python3.13 -m venv venv >>"$LOG" 2>&1
./venv/bin/python -m pip install --upgrade pip >>"$LOG" 2>&1
echo "=== core install ===" >>"$LOG"
./venv/bin/pip install numpy pandas scipy scikit-learn matplotlib seaborn statsmodels openpyxl adjustText >>"$LOG" 2>&1 && echo "CORE_OK" | tee -a "$LOG" || echo "CORE_FAIL" | tee -a "$LOG"
echo "=== bio install (decoupler/omnipath) ===" >>"$LOG"
./venv/bin/pip install decoupler omnipath >>"$LOG" 2>&1 && echo "DECOUPLER_OK" | tee -a "$LOG" || echo "DECOUPLER_FAIL" | tee -a "$LOG"
echo "=== mofapy2 (optional) ===" >>"$LOG"
./venv/bin/pip install mofapy2 >>"$LOG" 2>&1 && echo "MOFA_OK" | tee -a "$LOG" || echo "MOFA_FAIL" | tee -a "$LOG"
echo "=== versions ===" | tee -a "$LOG"
./venv/bin/python -c "import importlib.metadata as m
for p in ['numpy','pandas','scipy','scikit-learn','matplotlib','seaborn','statsmodels','openpyxl','decoupler','omnipath','mofapy2']:
    try: print(p, m.version(p))
    except Exception as e: print(p,'MISSING')" 2>&1 | tee -a "$LOG"
echo "VENV_DONE" | tee -a "$LOG"
