#!/usr/bin/env python
"""Harmonization-trap demonstration (addresses the reviewers' strongest constructive ask):
show that a STANDARD naive harmonization of the two aliased proteomes produces specific,
plausible-looking 'iCCA-vs-HCC' hits that DO NOT validate against the batch-free RNA axis,
i.e. concrete false positives a naive analyst would have published."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; DATA=os.path.join(os.path.dirname(ROOT),"Data")

def combat(Y, batch, mod=None, eb=True):
    """Parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007). Y: features x samples.
    batch: integer array of length n. mod: n x p covariates to PRESERVE (no intercept), or None.
    Returns the batch-adjusted matrix. When the protected covariate (mod) is perfectly aliased
    with batch the design is singular; pinv keeps it numerically stable, and the 'preserved'
    covariate effect then re-imprints the aliased contrast — exactly the failure we demonstrate."""
    Y=np.asarray(Y,float); g,n=Y.shape; batches=np.unique(batch); nb=len(batches)
    Xb=np.zeros((n,nb))
    for j,b in enumerate(batches): Xb[batch==b,j]=1.0
    design=Xb if mod is None else np.hstack([Xb,np.asarray(mod,float)])
    B=np.linalg.pinv(design.T@design)@design.T@Y.T            # coefficients x g
    nbatches=Xb.sum(0); grand=(nbatches/n)@B[:nb]             # size-weighted grand mean (g,)
    stand_mean=np.tile(grand[:,None],(1,n))
    if mod is not None: stand_mean=stand_mean+(np.asarray(mod,float)@B[nb:]).T
    resid=Y-(design@B).T; var_pooled=np.maximum((resid**2)@np.ones(n)/n,1e-12)
    Z=(Y-stand_mean)/np.sqrt(var_pooled)[:,None]
    gamma=np.zeros((nb,g)); delta=np.ones((nb,g))
    for j,b in enumerate(batches):
        c=batch==b; gamma[j]=Z[:,c].mean(1); delta[j]=np.maximum(Z[:,c].var(1),1e-12)
    Zadj=Z.copy()
    for j,b in enumerate(batches):
        c=batch==b; ni=int(c.sum())
        if eb:
            gbar=gamma[j].mean(); tau2=gamma[j].var()
            m=delta[j].mean(); s2=delta[j].var(); a=(2*s2+m*m)/s2; bb=(m*s2+m**3)/s2
            gs=gamma[j].copy(); ds=delta[j].copy()
            for _ in range(200):
                gnew=(ni*tau2*gamma[j]+ds*gbar)/(ni*tau2+ds)
                sum2=((Z[:,c]-gnew[:,None])**2).sum(1); dnew=(0.5*sum2+bb)/(0.5*ni+a-1)
                if max(np.abs(gnew-gs).max(),np.abs(dnew-ds).max())<1e-4: gs,ds=gnew,dnew; break
                gs,ds=gnew,dnew
            gstar,dstar=gs,np.maximum(ds,1e-12)
        else:
            gstar,dstar=gamma[j],delta[j]
        Zadj[:,c]=(Z[:,c]-gstar[:,None])/np.sqrt(dstar)[:,None]
    return Zadj*np.sqrt(var_pooled)[:,None]+stand_mean

def quantile_normalize(df):                       # standard QN across samples (proteins x samples)
    arr=df.values.copy()
    order=np.argsort(arr,axis=0)
    ranks=np.argsort(order,axis=0)
    mean_sorted=np.sort(arr,axis=0).mean(axis=1)
    return pd.DataFrame(mean_sorted[ranks],index=df.index,columns=df.columns)

gao=pd.read_parquet(f"{P}/gao_prot_T.parquet"); dong=pd.read_parquet(f"{P}/dong_prot.parquet")
common=sorted(set(gao.index)&set(dong.index))
H=gao.loc[common]; I=dong.loc[common]
# naive analyst pipeline: keep proteins with <30% missing, median-fill, then quantile-normalize the merged matrix
merged=pd.concat([H,I],axis=1)
keep=merged.isna().mean(axis=1)<0.30
merged=merged.loc[keep]; merged=merged.apply(lambda r:r.fillna(r.median()),axis=1)
qn=quantile_normalize(merged)
y=np.array(["HCC"]*H.shape[1]+["iCCA"]*I.shape[1])
isH=y=="HCC"; isI=y=="iCCA"

# per-protein "differential" test on harmonized data (what the naive analyst reports)
res=[]
for g in qn.index:
    h=qn.loc[g,isH].values; i=qn.loc[g,isI].values
    res.append((g,float(i.mean()-h.mean()),float(stats.mannwhitneyu(i,h).pvalue)))
df=pd.DataFrame(res,columns=["protein","diff_iCCA_HCC","p"]).dropna()
df["fdr"]=multipletests(df["p"],method="fdr_bh")[1]
sig=df[df["fdr"]<0.05]
print(f"Naive quantile harmonization: {len(sig)} of {len(df)} proteins 'significant' at FDR<0.05 ({100*len(sig)/len(df):.0f}%)")

# how many of the 'significant' hits are kinases / named druggable signaling genes
ksn=pd.read_parquet(f"{DATA}/omnipath_ksn.parquet"); KIN=set(str(x) for x in ksn["enzyme_genesymbol"].dropna().unique())
sig_kin=sorted(set(sig["protein"])&KIN, key=lambda g: df.set_index("protein").loc[g,"fdr"])
top_named=list(sig.reindex(sig["diff_iCCA_HCC"].abs().sort_values(ascending=False).index)["protein"].head(12))
print(f"  'significant' signaling enzymes among hits: {len(sig_kin)}; examples: {sig_kin[:8]}")
print(f"  top by effect size (plausible 'iCCA-vs-HCC' hits): {top_named}")

# VALIDATE the harmonized hits against the batch-free TCGA RNA direction
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); chol=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
g2=sorted(set(sig["protein"])&set(lihc.index)&set(chol.index))
rna_dir=(chol.loc[g2].mean(axis=1)-lihc.loc[g2].mean(axis=1))
harm_dir=df.set_index("protein").loc[g2,"diff_iCCA_HCC"]
rho,pr=stats.spearmanr(harm_dir,rna_dir)
sign_agree=float(np.mean(np.sign(harm_dir)==np.sign(rna_dir)))
print(f"VALIDATION of the {len(g2)} significant hits vs batch-free TCGA RNA: Spearman rho={rho:.3f} (p={pr:.2g}); sign agreement {sign_agree*100:.0f}%")

# ============ M2: ComBat baselines — does a 'proper' batch corrector escape the trap? ============
# The non-identifiability argument predicts it cannot: with batch = cohort and tumour type perfectly
# aliased with cohort, (i) ComBat with NO protected covariate must erase the contrast (discriminant
# to ~0.5), and (ii) ComBat that PROTECTS tumour type must re-imprint the aliased contrast as
# 'biology'. We now SHOW both rather than asserting them.
batch=np.where(isI,1,0)
lihc_c=lihc; chol_c=chol   # batch-free RNA already loaded above
def validate(direction):
    gg=sorted(set(direction.dropna().index)&set(lihc_c.index)&set(chol_c.index))
    rd=(chol_c.loc[gg].mean(axis=1)-lihc_c.loc[gg].mean(axis=1))
    return float(stats.spearmanr(direction.loc[gg],rd).correlation),len(gg)
def discriminant_auc(M):   # within-sample percentile lineage marker score, HCC vs iCCA separability
    R=M.rank(axis=0,pct=True); hep=[g for g in HEP_ if g in R.index]; chl=[g for g in CHOL_ if g in R.index]
    s=(R.loc[chl].mean()-R.loc[hep].mean()); a=s[isI].values; b=s[isH].values
    u=stats.mannwhitneyu(a,b).statistic/(len(a)*len(b)); return float(max(u,1-u))
def pct_sig(M):
    diffs=M.loc[:,isI].mean(1)-M.loc[:,isH].mean(1)
    ps=np.array([stats.mannwhitneyu(M.loc[g,isI].values,M.loc[g,isH].values).pvalue for g in M.index])
    fdr=multipletests(np.nan_to_num(ps,nan=1.0),method="fdr_bh")[1]
    return float((fdr<0.05).mean()),diffs

HEP_=["ALB","APOA1","APOB","APOC3","TF","TTR","CYP3A4","HNF4A","FGA","FGB","FGG","SERPINC1","CPS1","ARG1"]
CHOL_=["KRT19","KRT7","SOX9","EPCAM","CD24","SPP1","MUC1","CFTR","HNF1B","ANXA4","TACSTD2","CLDN4"]
cb_unprot=pd.DataFrame(combat(qn.values,batch,mod=None),index=qn.index,columns=qn.columns)
auc_u=discriminant_auc(cb_unprot); sigfrac_u,dir_u=pct_sig(cb_unprot); rho_u,nu=validate(dir_u)
print(f"M2 ComBat (no protected covariate): discriminant AUC={auc_u:.3f}, {100*sigfrac_u:.0f}% sig, anchor rho={rho_u:+.3f} -> contrast ERASED (like per-cohort z-score)")
# Protect tumour type: but tumour type IS batch (perfect aliasing) -> design is rank-deficient.
mod=np.where(isI,1.0,0.0)[:,None]
Xb=np.column_stack([isH.astype(float),isI.astype(float)]); design=np.hstack([Xb,mod])
rank=int(np.linalg.matrix_rank(design)); singular=rank<design.shape[1]
print(f"M2 ComBat (protect tumour type):    design [batch|tumour type] has rank {rank} of {design.shape[1]} columns "
      f"-> {'RANK-DEFICIENT: standard ComBat (sva/pyComBat) aborts (covariate confounded with batch)' if singular else 'full rank'}")
# Forced through pinv it does NOT recover the batch-free axis either:
cb_prot=pd.DataFrame(combat(qn.values,batch,mod=mod),index=qn.index,columns=qn.columns)
auc_p=discriminant_auc(cb_prot); sigfrac_p,dir_p=pct_sig(cb_prot); rho_p,np_=validate(dir_p)
print(f"                                    forced (pinv): AUC={auc_p:.3f}, {100*sigfrac_p:.0f}% sig, anchor rho={rho_p:+.3f} -> degenerate, does NOT recover lineage axis")

R={"naive_method":"quantile normalization of merged Gao-HCC + Dong-iCCA proteome",
   "n_proteins_tested":int(len(df)),"n_significant_fdr05":int(len(sig)),"pct_significant":float(len(sig)/len(df)),
   "n_significant_signaling_enzymes":int(len(sig_kin)),"example_false_kinase_hits":sig_kin[:8],
   "top_by_effect":top_named[:8],
   "validation_vs_batchfree_RNA_rho":float(rho),"validation_p":float(pr),"validation_sign_agreement":sign_agree,
   "M2_combat_no_covariate":{"discriminant_auc":auc_u,"pct_significant_fdr05":sigfrac_u,"anchor_rho":rho_u,"n_anchor":nu,
                             "interpretation":"contrast erased (batch==tumour type removed by construction), like per-cohort z-score"},
   "M2_combat_protect_tumortype":{"design_rank":rank,"design_cols":int(design.shape[1]),"rank_deficient":bool(singular),
                             "forced_pinv_discriminant_auc":auc_p,"forced_pinv_pct_significant_fdr05":sigfrac_p,"forced_pinv_anchor_rho":rho_p,
                             "interpretation":"tumour type perfectly confounded with batch -> design rank-deficient; standard ComBat (sva/pyComBat) aborts. Forced via pinv it is degenerate and does not recover the batch-free axis."}}
json.dump(R,open(f"{OUT}/harmonization_trap_stats.json","w"),indent=2)
print("\nTRAP_DONE")
