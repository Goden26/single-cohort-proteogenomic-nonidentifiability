#!/usr/bin/env python
"""Module 2 (reframed, confound-aware): the Gao(HCC) and Dong(iCCA) proteomes are perfectly
cohort-confounded with tumor type, so unsupervised separation is NOT used as evidence.
Instead we test whether the proteome carries genuine LINEAGE biology that is robust to the
global reference offset, using within-sample percentile normalization, and validate the
per-protein iCCA-vs-HCC direction against the BATCH-FREE TCGA RNA axis."""
import os, json
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"; DATA=os.path.join(os.path.dirname(ROOT),"Data")
HEP=["ALB","APOA1","APOB","APOC3","TF","TTR","CYP3A4","HNF4A","FGA","FGB","FGG","SERPINC1","CPS1","ARG1"]
CHOL=["KRT19","KRT7","SOX9","EPCAM","CD24","SPP1","MUC1","CFTR","HNF1B","ANXA4","TACSTD2","CLDN4"]

gao=pd.read_parquet(f"{P}/gao_prot_T.parquet")     # proteins x 159 HCC
dong=pd.read_parquet(f"{P}/dong_prot.parquet")     # proteins x 214 iCCA
common=sorted(set(gao.index)&set(dong.index))
print(f"common proteins={len(common)}  HCC={gao.shape[1]}  iCCA={dong.shape[1]}")

# within-sample percentile rank (cancels global TMT reference offset; batch-robust)
gaoR=gao.rank(axis=0,pct=True); dongR=dong.rank(axis=0,pct=True)

def pct_marker(R):
    h=[g for g in HEP if g in R.index]; c=[g for g in CHOL if g in R.index]
    return (R.loc[c].mean()-R.loc[h].mean()), len(h), len(c)
sH,nh,nc=pct_marker(gaoR); sI,_,_=pct_marker(dongR)
mw=stats.mannwhitneyu(sI,sH,alternative="greater")
print(f"proteome lineage marker score: HCC median={sH.median():.3f}  iCCA median={sI.median():.3f}  (iCCA>HCC p={mw.pvalue:.1e})")

# per-protein iCCA-vs-HCC direction (batch-robust) vs batch-free TCGA RNA direction
prot_dir=(dongR.loc[common].mean(axis=1)-gaoR.loc[common].mean(axis=1))
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); cholrna=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
g2=sorted(set(common)&set(lihc.index)&set(cholrna.index))
rna_dir=(cholrna.loc[g2].mean(axis=1)-lihc.loc[g2].mean(axis=1))
rho,pr=stats.spearmanr(prot_dir.loc[g2],rna_dir)
# agreement on strongly-directional genes
strong=rna_dir[abs(rna_dir)>1.0].index
agree=np.mean(np.sign(prot_dir.loc[strong])==np.sign(rna_dir.loc[strong]))
print(f"per-protein iCCA-vs-HCC direction vs batch-free TCGA RNA: Spearman rho={rho:.3f} p={pr:.1e} (n={len(g2)})")
print(f"sign agreement on |RNA log2FC|>1 genes (n={len(strong)}): {agree*100:.0f}%")

# --- M3: normalization sweep — is the iCCA-vs-HCC direction vs RNA sign harmonization-DEPENDENT? ---
# Compute the proteome direction under four standard normalizations and correlate each with the
# batch-free RNA axis. If the sign/magnitude moves with the normalization choice, the negative rho is
# a property of the harmonization, not a fixed biological anti-correlation.
def qn_merge(H,I):
    m=pd.concat([H,I],axis=1); a=m.values.copy()
    order=np.argsort(a,axis=0); ranks=np.argsort(order,axis=0)
    qn=pd.DataFrame(np.sort(a,axis=0).mean(axis=1)[ranks],index=m.index,columns=m.columns)
    return qn.iloc[:,:H.shape[1]],qn.iloc[:,H.shape[1]:]
def rho_vs_rna(direction):
    gg=sorted(set(direction.dropna().index)&set(rna_dir.index))
    return float(stats.spearmanr(direction.loc[gg],rna_dir.loc[gg]).correlation),len(gg)
Hc,Ic=gao.loc[common],dong.loc[common]
sweep={}
# (a) raw mean-abundance difference (no normalization) — dominated by the global TMT reference offset
sweep["raw_abundance"]=dict(zip(["rho","n"],rho_vs_rna(Ic.mean(axis=1)-Hc.mean(axis=1))))
# (b) within-sample percentile rank (the paper's main choice)
sweep["within_sample_percentile"]=dict(zip(["rho","n"],rho_vs_rna(dongR.loc[common].mean(axis=1)-gaoR.loc[common].mean(axis=1))))
# (c) per-cohort z-score (removes the between-cohort contrast by construction)
Hz=Hc.sub(Hc.mean(axis=1),axis=0).div(Hc.std(axis=1)+1e-9,axis=0)
Iz=Ic.sub(Ic.mean(axis=1),axis=0).div(Ic.std(axis=1)+1e-9,axis=0)
sweep["per_cohort_zscore"]=dict(zip(["rho","n"],rho_vs_rna(Iz.mean(axis=1)-Hz.mean(axis=1))))
# (d) quantile-normalized merge (the naive-analyst harmonization)
Hq,Iq=qn_merge(Hc,Ic)
sweep["quantile_merge"]=dict(zip(["rho","n"],rho_vs_rna(Iq.mean(axis=1)-Hq.mean(axis=1))))
print("M3 normalization sweep (proteome iCCA-HCC direction vs batch-free RNA rho):")
for k,v in sweep.items(): print(f"   {k:26s} rho={v['rho']:+.3f} (n={v['n']})")

# PCA on rank-normalized common proteins — shown ONLY to illustrate the confound (cohort==type)
X=pd.concat([gaoR.loc[common],dongR.loc[common]],axis=1).T.fillna(0.5)
y=np.array(["HCC"]*gao.shape[1]+["iCCA"]*dong.shape[1])
pca=PCA(n_components=5,random_state=0).fit(X.values); PC=pca.transform(X.values)

res={"common_proteins":len(common),"n_HCC":int(gao.shape[1]),"n_iCCA":int(dong.shape[1]),
 "marker_score_HCC_median":float(sH.median()),"marker_score_iCCA_median":float(sI.median()),
 "marker_iCCA_gt_HCC_p":float(mw.pvalue),
 "proteome_vs_RNA_direction_spearman":float(rho),"proteome_vs_RNA_direction_p":float(pr),
 "sign_agreement_strong_genes":float(agree),"n_strong_genes":int(len(strong)),
 "M3_normalization_sweep_rho_vs_RNA":sweep}
json.dump(res,open(f"{ROOT}/results/proteome_manifold_stats.json","w"),indent=2)
print(json.dumps(res,indent=2))

fig,ax=plt.subplots(1,2,figsize=(12,4.8))
c={"HCC":"#2c7fb8","iCCA":"#d95f0e"}
for lab in ["HCC","iCCA"]:
    m=y==lab; ax[0].scatter(PC[m,0],PC[m,1],s=14,alpha=.6,label=f"{lab} (n={m.sum()})",color=c[lab])
ax[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f}%)"); ax[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f}%)")
ax[0].set_title("Proteome PCA — cohort==tumor-type\n(confounded; not used as separation evidence)"); ax[0].legend(frameon=False,fontsize=9)
ax[1].scatter(rna_dir,prot_dir.loc[g2],s=6,alpha=.25,color="gray")
for g in [x for x in HEP+CHOL if x in g2]:
    ax[1].annotate(g,(rna_dir[g],prot_dir.loc[g]),fontsize=7,color=("#d95f0e" if g in CHOL else "#2c7fb8"))
ax[1].axhline(0,ls=":",c="k",lw=.6); ax[1].axvline(0,ls=":",c="k",lw=.6)
ax[1].set_xlabel("TCGA RNA iCCA−HCC (batch-free, log2)"); ax[1].set_ylabel("Proteome iCCA−HCC (within-sample percentile)")
ax[1].set_title(f"Proteome direction is discordant with the batch-free axis\nSpearman ρ={rho:.2f}, p={pr:.0e}; only {agree*100:.0f}% sign-agree (harmonization-dependent)")
plt.tight_layout(); plt.savefig(f"{FIG}/fig2_proteome_confound_check.png",dpi=160); print("saved fig2")
print("PROTEOME_MODULE_DONE")
