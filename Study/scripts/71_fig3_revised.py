#!/usr/bin/env python
"""Revised Figure 3: the comparator artifact. Left = the 11 'all-other-lineage' iCCA-selective hits
are mostly essential in ALL groups (pan-essential); right = only 4 survive the direct iCCA-vs-HCC test."""
import pandas as pd, numpy as np
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"
a=pd.read_csv(f"{OUT}/rev_allother11_audit.csv").sort_values("icca_selective")
d=pd.read_csv(f"{OUT}/rev_direct_icca_selective.csv").sort_values("dep_diff_iCCA_minus_HCC")

fig,ax=plt.subplots(1,2,figsize=(14,5.4))
# LEFT: the 11 all-other hits, Chronos in iCCA / HCC / other lineages
x=np.arange(len(a)); w=0.27
ax[0].bar(x-w,a["chronos_iCCA"],w,label="iCCA lines",color="#d95f0e")
ax[0].bar(x,  a["chronos_HCC"], w,label="HCC lines", color="#2c7fb8")
ax[0].bar(x+w,a["chronos_other"],w,label="other lineages",color="#999999")
ax[0].axhline(-0.5,ls=":",c="k",lw=.8,label="essentiality")
for i,(_,r) in enumerate(a.iterrows()):
    tag=""
    if r["common_essential"]: tag+="P"          # pan-essential
    if r["selective_both"]: tag+="B"             # selective for both liver lineages
    if tag: ax[0].text(i,0.12,tag,ha="center",fontsize=8,color="firebrick",fontweight="bold")
ax[0].set_xticks(x); ax[0].set_xticklabels(a["kinase"],rotation=60,ha="right",fontsize=8)
ax[0].set_ylabel("Mean Chronos gene effect"); ax[0].legend(frameon=False,fontsize=8,loc="lower left")
ax[0].set_title("'iCCA-selective vs all-other-lineages' (n=11)\nP=pan-essential, B=selective for BOTH liver lineages")

# RIGHT: the 4 that survive the DIRECT iCCA-vs-HCC test
x2=np.arange(len(d)); w2=0.38
ax[1].bar(x2-w2/2,d["chronos_iCCA"],w2,label="iCCA lines",color="#d95f0e")
ax[1].bar(x2+w2/2,d["chronos_HCC"], w2,label="HCC lines", color="#2c7fb8")
ax[1].axhline(-0.5,ls=":",c="k",lw=.8)
for i,(_,r) in enumerate(d.iterrows()):
    ax[1].text(i,min(r["chronos_iCCA"],r["chronos_HCC"])-0.08,f"FDR {r['dep_fdr_icca_vs_hcc']:.02f}",ha="center",fontsize=8)
ax[1].set_xticks(x2); ax[1].set_xticklabels(d["kinase"],fontsize=10,fontweight="bold")
ax[1].set_ylabel("Mean Chronos gene effect"); ax[1].legend(frameon=False,fontsize=9)
ax[1].set_title("Survive the DIRECT iCCA-vs-HCC test (FDR<0.1)\nnutrient/ion transport: TFRC, SLC3A2, TRPM7, NHEJ1")
plt.tight_layout(); plt.savefig(f"{FIG}/fig3_depmap_rna_kinases.png",dpi=160)
print("saved revised fig3")
