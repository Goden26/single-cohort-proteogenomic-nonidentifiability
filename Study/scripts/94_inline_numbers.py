#!/usr/bin/env python
"""Emit the small inline numbers quoted in the manuscript text that were previously not captured in
any released stat file (reviewer minor: 'emit currently-unstored numbers in the regeneration script').
Writes results/manuscript_inline_numbers.json. All values are recomputed from the released artifacts."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import StratifiedKFold, cross_val_predict
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; DATA=os.path.join(os.path.dirname(ROOT),"Data")

class TopKVariance(BaseEstimator,TransformerMixin):
    def __init__(self,k=2000): self.k=k
    def fit(self,X,y=None): v=np.var(np.asarray(X),axis=0); self.idx_=np.argsort(v)[::-1][:min(self.k,X.shape[1])]; return self
    def transform(self,X): return np.asarray(X)[:,self.idx_]

R={}

# 1) per-cohort z-scoring and the lineage-MARKER discriminant: removes the mean/marker contrast (~0.5).
#    We also report a full-LDA CV-AUC, which stays high because per-cohort z-scoring removes the mean
#    shift but not residual multivariate (covariance) structure; that residual separability is the
#    batch artifact the paper explicitly declines to read as biology.
HEP=["ALB","APOA1","APOB","APOC3","TF","TTR","CYP3A4","HNF4A","FGA","FGB","FGG","SERPINC1","CPS1","ARG1"]
CHOL=["KRT19","KRT7","SOX9","EPCAM","CD24","SPP1","MUC1","CFTR","HNF1B","ANXA4","TACSTD2","CLDN4"]
gao=pd.read_parquet(f"{P}/gao_prot_T.parquet"); dong=pd.read_parquet(f"{P}/dong_prot.parquet")
common=sorted(set(gao.index)&set(dong.index))
H=gao.loc[common]; I=dong.loc[common]
Hz=H.sub(H.mean(axis=1),axis=0).div(H.std(axis=1)+1e-9,axis=0)        # per-cohort z-score
Iz=I.sub(I.mean(axis=1),axis=0).div(I.std(axis=1)+1e-9,axis=0)
hep=[g for g in HEP if g in common]; chl=[g for g in CHOL if g in common]
sH=(Hz.loc[chl].mean()-Hz.loc[hep].mean()); sI=(Iz.loc[chl].mean()-Iz.loc[hep].mean())  # per-sample marker score
u=stats.mannwhitneyu(sI,sH).statistic/(len(sI)*len(sH))
R["per_cohort_zscore_marker_discriminant_auc"]=float(max(u,1-u))
X=pd.concat([Hz,Iz],axis=1).T.fillna(0.0).values
y=np.array([0]*H.shape[1]+[1]*I.shape[1])
pipe=Pipeline([("var",TopKVariance(2000)),("sc",StandardScaler()),("lda",LinearDiscriminantAnalysis())])
sc=cross_val_predict(pipe,X,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),method="decision_function")
auc=stats.mannwhitneyu(sc[y==1],sc[y==0]).statistic/((y==1).sum()*(y==0).sum())
R["per_cohort_zscore_fullLDA_cv_auc"]=float(max(auc,1-auc))

# 2,3) EGFR raw direct-contrast P and FGFR1 log2FC, from the released integrated table
ki=pd.read_csv(f"{OUT}/kinase_lineage_integrated.csv").set_index("kinase")
R["EGFR_raw_dep_p_icca_vs_hcc"]=float(ki.loc["EGFR","dep_p_icca_vs_hcc"])
audit=pd.read_csv(f"{OUT}/rev_allother11_audit.csv").set_index("kinase")   # BH-FDR over the direct contrast
R["EGFR_dep_fdr_icca_vs_hcc"]=float(audit.loc["EGFR","dep_fdr_icca_vs_hcc"])
R["EGFR_log2FC_iCCA_HCC"]=float(ki.loc["EGFR","log2FC_iCCA_HCC"])
R["FGFR1_log2FC_iCCA_HCC"]=float(ki.loc["FGFR1","log2FC_iCCA_HCC"])
R["FGFR1_dep_diff_iCCA_minus_HCC"]=float(ki.loc["FGFR1","dep_diff_iCCA_minus_HCC"])  # >0: more essential in HCC

# 4) OmniPath kinase-substrate network size
ksn=pd.read_parquet(f"{DATA}/omnipath_ksn.parquet")
R["omnipath_phosphorylation_relations"]=int(len(ksn))
R["omnipath_enzymes"]=int(ksn["enzyme_genesymbol"].nunique())

json.dump(R,open(f"{OUT}/manuscript_inline_numbers.json","w"),indent=2)
print(json.dumps(R,indent=2)); print("INLINE_NUMBERS_DONE")
