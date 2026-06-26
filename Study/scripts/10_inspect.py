#!/usr/bin/env python
"""Inspect downloaded data: Excel sheet names/shapes + first rows; CSV/TSV heads.
Run with the venv python once downloads complete."""
import os, gzip, sys
import pandas as pd
DATA=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"Data")

def inspect_xlsx(fn, max_sheets=20):
    path = os.path.join(DATA, fn)
    if not os.path.exists(path):
        print(f"  !! MISSING {fn}"); return
    try:
        xl = pd.ExcelFile(path, engine="openpyxl")
    except Exception as e:
        print(f"  !! ERROR opening {fn}: {e}"); return
    print(f"\n### {fn}  ({os.path.getsize(path)/1e6:.1f} MB) — {len(xl.sheet_names)} sheets")
    for s in xl.sheet_names[:max_sheets]:
        try:
            df = xl.parse(s, nrows=4)
            cols = list(df.columns)[:8]
            print(f"  - [{s}]  cols~{df.shape[1]}  first_cols={cols}")
        except Exception as e:
            print(f"  - [{s}]  (parse error: {e})")

def head_text(fn, n=3, gz=False, sep=None):
    path = os.path.join(DATA, fn)
    if not os.path.exists(path):
        print(f"  !! MISSING {fn}"); return
    op = gzip.open if gz else open
    print(f"\n### {fn}  ({os.path.getsize(path)/1e6:.1f} MB)")
    with op(path, "rt", errors="replace") as f:
        for i, line in enumerate(f):
            if i >= n: break
            print("   ", line.rstrip()[:240])

if __name__ == "__main__":
    print("===== EXCEL FILES =====")
    for fn in ["gao_S1.xlsx","gao_S2_cnv.xlsx","dong_S1.xlsx","dong_S2_cnv.xlsx","dong_S5_sub.xlsx"]:
        inspect_xlsx(fn)
    print("\n===== CSV/TSV FILES =====")
    head_text("depmap_Model.csv", 2)
    head_text("depmap_CRISPRGeneEffect.csv", 1)
    head_text("TCGA-LIHC.tpm.tsv.gz", 2, gz=True)
    head_text("TCGA-CHOL.tpm.tsv.gz", 2, gz=True)
    head_text("TCGA-LIHC.survival.tsv.gz", 3, gz=True)
    head_text("GSE241466_iCHC.txt.gz", 2, gz=True)
    head_text("GSE231854_fpkm.txt.gz", 2, gz=True)
