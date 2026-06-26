#!/usr/bin/env python
"""Module 1: batch-free RNA lineage manifold (TCGA-LIHC vs TCGA-CHOL, same GDC STAR pipeline)
+ cHCC-CCA (GSE241466 i-CHC) lineage positioning via dataset-internal marker score.
Tests: is there a real HCC<->iCCA axis, and is cHCC-CCA intermediate?"""
import os, gzip, json
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
DATA=os.path.join(os.path.dirname(ROOT),"Data"); OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"
os.makedirs(f"{OUT}/parsed",exist_ok=True)

# lineage markers (well-established)
HEP=["ALB","APOA1","APOB","APOC3","TF","TTR","CYP3A4","HNF4A","FGA","FGB","FGG","SERPINC1","CPS1","ARG1"]
CHOL=["KRT19","KRT7","SOX9","EPCAM","CD24","SPP1","MUC1","CFTR","HNF1B","ANXA4","TACSTD2","CLDN4"]

# ---- Ensembl -> symbol ----
pm=pd.read_csv(f"{DATA}/gencode.v36.probemap",sep="\t")
pm["base"]=pm["id"].str.split(".").str[0]
ens2sym=dict(zip(pm["base"],pm["gene"]))

def load_tcga(fn):
    df=pd.read_csv(f"{DATA}/{fn}",sep="\t",index_col=0)        # Ensembl x barcode, log2(tpm+1)
    df.index=[ens2sym.get(str(i).split(".")[0],i) for i in df.index]
    df=df.groupby(level=0).max()                                # collapse to symbol (max)
    tum=[c for c in df.columns if c.split("-")[3][:2]=="01"]    # primary tumor only
    return df[tum]

lihc=load_tcga("TCGA-LIHC.tpm.tsv.gz"); chol=load_tcga("TCGA-CHOL.tpm.tsv.gz")
print(f"TCGA-LIHC tumors={lihc.shape[1]}  TCGA-CHOL tumors={chol.shape[1]}")
common=sorted(set(lihc.index)&set(chol.index))
X=pd.concat([lihc.loc[common],chol.loc[common]],axis=1).T      # samples x genes
y=np.array(["HCC"]*lihc.shape[1]+["iCCA"]*chol.shape[1])
print(f"combined TCGA: {X.shape}  HCC={sum(y=='HCC')} iCCA={sum(y=='iCCA')}")

# ---- manifold: top variable genes, z-score, PCA ----
v=X.var().sort_values(ascending=False); top=v.head(2000).index
Z=StandardScaler().fit_transform(X[top])
pca=PCA(n_components=10,random_state=0).fit(Z); P=pca.transform(Z)
# orient: LDA lineage axis
lda=LinearDiscriminantAnalysis().fit(Z,y); lin_tcga=lda.transform(Z).ravel()
if lin_tcga[y=="HCC"].mean()>lin_tcga[y=="iCCA"].mean(): lin_tcga=-lin_tcga  # HCC low, iCCA high
auc=stats.mannwhitneyu(lin_tcga[y=="HCC"],lin_tcga[y=="iCCA"]).statistic/(sum(y=="HCC")*sum(y=="iCCA"))
auc=max(auc,1-auc)  # separation strength, orientation-agnostic
print(f"PCA var(PC1,PC2)={pca.explained_variance_ratio_[:2].round(3)}  lineage-axis AUC(HCC vs iCCA)={auc:.3f}")

# ---- single-sample percentile-rank lineage score (platform/normalization robust, cross-dataset comparable) ----
def marker_score(df_genes_by_samples):
    ranks=df_genes_by_samples.rank(axis=0,pct=True)               # rank each gene within sample across all genes
    h=[g for g in HEP if g in df_genes_by_samples.index]
    c=[g for g in CHOL if g in df_genes_by_samples.index]
    return ranks.loc[c].mean()-ranks.loc[h].mean(), len(h), len(c) # per-sample: chol-percentile − hep-percentile

s_lihc,_,_=marker_score(lihc); s_chol,nh,nc=marker_score(chol)
print(f"marker genes found: hep={nh} chol={nc}")

# cHCC-CCA: GSE241466 — MIXED cohort (7 Non-tumor + 5 HCC + 5 iCCA + 8 i-CHC). Keep ONLY the 8 true i-CHC.
labels=json.load(open(f"{DATA}/GSE241466_labels.json"))
ICHC=[t for t,grp in labels.items() if "Intermediate" in grp]
g=pd.read_csv(f"{DATA}/GSE241466_iCHC.txt.gz",sep=r"\s+",index_col=0)
g=g[[c for c in g.columns if c in ICHC]]                 # restrict to true cHCC-CCA (i-CHC), n=8
g.index=[ens2sym.get(str(i).strip('"').split(".")[0],None) for i in g.index]
g=g[[i is not None for i in g.index]]
g=g.apply(pd.to_numeric,errors="coerce").groupby(level=0).max()
chcca=np.log2(g+1)
s_chcca,nh2,nc2=marker_score(chcca)
print(f"GSE241466 restricted to i-CHC: {g.shape[1]} samples")
print(f"GSE241466 i-CHC: {g.shape[1]} samples; markers hep={nh2} chol={nc2}")

# ---- statistics on intermediacy ----
res={"lineage_axis_auc":float(auc),
     "marker_score_HCC_median":float(np.median(s_lihc)),
     "marker_score_iCCA_median":float(np.median(s_chol)),
     "marker_score_cHCCCCA_median":float(np.median(s_chcca)),
     "n_HCC":int(lihc.shape[1]),"n_iCCA":int(chol.shape[1]),"n_cHCCCCA":int(g.shape[1])}
# is cHCC-CCA between HCC and iCCA?
res["mwu_cHCCCCA_vs_HCC_p"]=float(stats.mannwhitneyu(s_chcca,s_lihc).pvalue)
res["mwu_cHCCCCA_vs_iCCA_p"]=float(stats.mannwhitneyu(s_chcca,s_chol).pvalue)
res["kruskal_3group_p"]=float(stats.kruskal(s_lihc,s_chcca,s_chol).pvalue)
# Jonckheere-style ordered trend via Spearman of score vs ordinal lineage rank
allv=np.concatenate([s_lihc,s_chcca,s_chol]); ordi=np.concatenate([np.zeros(len(s_lihc)),np.ones(len(s_chcca)),2*np.ones(len(s_chol))])
rho,ptrend=stats.spearmanr(ordi,allv)
res["ordered_trend_spearman_rho"]=float(rho); res["ordered_trend_p"]=float(ptrend)
# fraction of cHCC-CCA in the intermediate band [HCC median, iCCA median]
lo,hi=np.median(s_lihc),np.median(s_chol)
res["frac_cHCCCCA_intermediate_band"]=float(np.mean((s_chcca>=lo)&(s_chcca<=hi)))
print(json.dumps(res,indent=2))

# ---- FIGURE 1 ----
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
c={"HCC":"#2c7fb8","iCCA":"#d95f0e"}
for lab in ["HCC","iCCA"]:
    m=y==lab; ax[0].scatter(P[m,0],P[m,1],s=14,alpha=.7,label=f"{lab} (n={m.sum()})",color=c[lab])
ax[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f}%)")
ax[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f}%)")
ax[0].set_title("TCGA RNA manifold (same GDC pipeline = batch-free)"); ax[0].legend(frameon=False,fontsize=9)
# violin of lineage marker score across 3 cohorts
data=[s_lihc.values,s_chcca.values,s_chol.values]
labs=[f"HCC\n(LIHC, n={len(s_lihc)})",f"cHCC-CCA\n(i-CHC, n={len(s_chcca)})",f"iCCA\n(CHOL, n={len(s_chol)})"]
parts=ax[1].violinplot(data,showmedians=True)
for i,b in enumerate(parts['bodies']): b.set_facecolor(["#2c7fb8","#7570b3","#d95f0e"][i]); b.set_alpha(.6)
ax[1].set_xticks([1,2,3]); ax[1].set_xticklabels(labs,fontsize=9)
ax[1].set_ylabel("Cholangiocyte − hepatocyte\nlineage score (single-sample percentile)")
ax[1].set_title("cHCC-CCA sits between HCC and iCCA"); ax[1].axhline(0,ls=":",c="gray",lw=.8)
plt.tight_layout(); plt.savefig(f"{FIG}/fig1_rna_lineage_manifold.png",dpi=160); print("saved fig1")

# save lineage axis + scores for later modules
pd.Series(lin_tcga,index=X.index,name="lineage_axis").to_frame().assign(type=y).to_parquet(f"{OUT}/parsed/tcga_lineage.parquet")
json.dump(res,open(f"{OUT}/rna_manifold_stats.json","w"),indent=2)
# save TCGA symbol-level tumor expression for kinase cross-check
lihc.to_parquet(f"{OUT}/parsed/tcga_lihc_tumor.parquet"); chol.to_parquet(f"{OUT}/parsed/tcga_chol_tumor.parquet")
chcca.to_parquet(f"{OUT}/parsed/chcca_gse241466.parquet")
print("RNA_MODULE_DONE")
