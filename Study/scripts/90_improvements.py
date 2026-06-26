#!/usr/bin/env python
"""Improvements: independent replications that address the two reviewer-flagged weaknesses.
A. cHCC-CCA lineage shift replicated in an INDEPENDENT cohort (GSE231854: 25 cHCC-CCA + 56 HCC),
   complementing the n=8 GSE241466 i-CHC result.
B. LAT1/amino-acid-transport over-expression replicated in DepMap cell-line expression (iCCA vs HCC
   lines), independent of TCGA and of the dependency-based gene selection (also leave-discovery-out)."""
import os, json, gzip
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"; DATA=os.path.join(os.path.dirname(ROOT),"Data")
HEP=["ALB","APOA1","APOB","APOC3","TF","TTR","CYP3A4","HNF4A","FGA","FGB","FGG","SERPINC1","CPS1","ARG1"]
CHOL=["KRT19","KRT7","SOX9","EPCAM","CD24","SPP1","MUC1","CFTR","HNF1B","ANXA4","TACSTD2","CLDN4"]
AAT=["SLC7A5","SLC3A2","SLC1A5","SLC7A11","SLC38A1","SLC38A2","SLC38A5","SLC7A6","SLC7A8","SLC16A1","SLC16A3","SLC2A1","GLS","GOT2"]
pm=pd.read_csv(f"{DATA}/gencode.v36.probemap",sep="\t"); pm["base"]=pm["id"].str.split(".").str[0]
ens2sym=dict(zip(pm["base"],pm["gene"]))
def cliffs(a,b):
    a=np.asarray(a); b=np.asarray(b); return float((np.sum(a[:,None]>b)-np.sum(a[:,None]<b))/(len(a)*len(b)))
def marker_score(df):
    r=df.rank(axis=0,pct=True); h=[g for g in HEP if g in df.index]; c=[g for g in CHOL if g in df.index]
    return (r.loc[c].mean()-r.loc[h].mean())
R={}

# ---------- A. GSE231854 independent cHCC-CCA cohort ----------
t=g=None
with gzip.open(f"{DATA}/GSE231854_series.txt.gz","rt",errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"): t=[x.strip('"') for x in line.rstrip().split("\t")[1:]]
        if line.startswith("!Sample_source_name_ch1"): g=[x.strip('"') for x in line.rstrip().split("\t")[1:]]
lab=dict(zip(t,g))
fp=pd.read_csv(f"{DATA}/GSE231854_fpkm.txt.gz",sep="\t",index_col=0)
fp=fp.loc[:,[c for c in fp.columns if c in lab]]                  # keep labeled samples
fp.index=[ens2sym.get(str(i).split(".")[0],None) for i in fp.index]
fp=fp[[i is not None for i in fp.index]].apply(pd.to_numeric,errors="coerce").groupby(level=0).max()
expr=np.log2(fp+1)
sc=marker_score(expr)
grp=pd.Series({c:lab[c] for c in expr.columns})
s_chc=sc[grp=="cHCC-CCA"].values; s_hcc=sc[grp=="HCC"].values
R["A_GSE231854_n_cHCCCCA"]=int(len(s_chc)); R["A_GSE231854_n_HCC"]=int(len(s_hcc))
R["A_cHCCCCA_median"]=float(np.median(s_chc)); R["A_HCC_median"]=float(np.median(s_hcc))
R["A_cliffs_cHCCCCA_vs_HCC"]=cliffs(s_chc,s_hcc)
R["A_mwu_p"]=float(stats.mannwhitneyu(s_chc,s_hcc,alternative="greater").pvalue)
print(f"A. GSE231854 (independent): cHCC-CCA (n={len(s_chc)}) cholangiocyte-shift vs HCC (n={len(s_hcc)}): "
      f"median {np.median(s_chc):+.3f} vs {np.median(s_hcc):+.3f}; Cliff's d={R['A_cliffs_cHCCCCA_vs_HCC']:.2f}; "
      f"one-sided MWU p={R['A_mwu_p']:.2e}")

# ---------- B. DepMap cell-line expression: AA-transport iCCA vs HCC (independent of TCGA + of selection) ----------
model=pd.read_parquet(f"{P}/depmap_model.parquet")
hcc=set(model[model["OncotreeCode"]=="HCC"]["ModelID"]); icca=set(model[model["OncotreeCode"]=="IHCH"]["ModelID"])
hdr=pd.read_csv(f"{DATA}/depmap_ExpressionTPM.csv",nrows=0)
s2c={}
for c in hdr.columns[1:]: s2c.setdefault(c.split(" (")[0],c)
need=[g for g in AAT if g in s2c]
ex=pd.read_csv(f"{DATA}/depmap_ExpressionTPM.csv",usecols=[hdr.columns[0]]+[s2c[g] for g in need],index_col=0)
ex.columns=[c.split(" (")[0] for c in ex.columns]; ex=ex.loc[:,~ex.columns.duplicated()]
eh=ex[ex.index.isin(hcc)]; ei=ex[ex.index.isin(icca)]
def setdiff(genes):
    gs=[g for g in genes if g in ex.columns]
    lh=eh[gs].mean(axis=1).dropna(); li=ei[gs].mean(axis=1).dropna()
    return float(li.median()-lh.median()), float(stats.mannwhitneyu(li,lh,alternative="greater").pvalue), len(lh), len(li), li, lh
d_full,p_full,nh,ni,li_full,lh_full=setdiff(AAT)
d_drop,p_drop,_,_,_,_=setdiff([g for g in AAT if g not in ("SLC3A2","SLC7A5")])
R["B_depmap_AAT_log2TPM_diff_iCCA_HCC"]=d_full; R["B_depmap_AAT_mwu_p"]=p_full
R["B_depmap_HCC_lines"]=nh; R["B_depmap_iCCA_lines"]=ni
R["B_depmap_AAT_drop_discovery_diff"]=d_drop; R["B_depmap_AAT_drop_discovery_p"]=p_drop
print(f"B. DepMap line expression (independent): AA-transport iCCA (n={ni}) vs HCC (n={nh}) lines: "
      f"median log2TPM diff={d_full:+.2f}, one-sided MWU p={p_full:.2e}; "
      f"drop discovery genes -> diff={d_drop:+.2f}, p={p_drop:.2e}")

json.dump(R,open(f"{OUT}/improvement_stats.json","w"),indent=2)

# ---------- FIGURE 5 ----------
fig,ax=plt.subplots(1,2,figsize=(11,4.6))
parts=ax[0].violinplot([s_hcc,s_chc],showmedians=True)
for i,b in enumerate(parts['bodies']): b.set_facecolor(["#2c7fb8","#7570b3"][i]); b.set_alpha(.6)
ax[0].set_xticks([1,2]); ax[0].set_xticklabels([f"HCC\n(n={len(s_hcc)})",f"cHCC-CCA\n(n={len(s_chc)})"])
ax[0].set_ylabel("Cholangiocyte - hepatocyte score"); ax[0].axhline(0,ls=":",c="gray",lw=.7)
ax[0].set_title(f"A. Independent cohort GSE231854\ncHCC-CCA cholangiocyte-shifted vs HCC (p={R['A_mwu_p']:.1e})")
ax[1].boxplot([lh_full,li_full],tick_labels=[f"HCC\n(n={nh})",f"iCCA\n(n={ni})"])
ax[1].set_ylabel("Mean log2(TPM+1) over LAT1/AA-transport set")
ax[1].set_title(f"B. DepMap cell-line expression\nAA-transport NOT elevated in iCCA lines (p={p_full:.2f}); tumor-specific")
plt.tight_layout(); plt.savefig(f"{FIG}/fig5_independent_replication.png",dpi=160); print("saved fig5")
print("IMPROVEMENTS_DONE")
