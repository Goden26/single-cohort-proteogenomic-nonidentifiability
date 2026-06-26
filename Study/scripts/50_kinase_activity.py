#!/usr/bin/env python
"""Module 3: per-sample KSEA kinase activity from phosphoproteomes (Gao HCC, Dong iCCA),
using the OmniPath kinase-substrate network. Cross-tumor HCC-vs-iCCA activity is reported but
FLAGGED as cohort-confounded; we quantify its concordance with batch-free TCGA RNA and keep
only concordant kinases as robust phospho-supported lineage kinases."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"

ksn=pd.read_parquet(f"{ROOT}/data/omnipath_ksn.parquet")
ksn=ksn.dropna(subset=["enzyme_genesymbol","substrate_genesymbol","residue_type","residue_offset"])
for c in ["enzyme_genesymbol","substrate_genesymbol","residue_type"]: ksn[c]=ksn[c].astype(str)
ksn["target"]=ksn["substrate_genesymbol"]+"_"+ksn["residue_type"]+ksn["residue_offset"].astype(int).astype(str)
k2t={k:set(g["target"]) for k,g in ksn.groupby("enzyme_genesymbol")}
print(f"KSN kinases={len(k2t)}  median substrates/kinase={int(np.median([len(v) for v in k2t.values()]))}")

def ksea(phos, min_sub=5):
    """phos: sites x samples (site id 'GENE:S123'); returns kinase x sample activity (per-sample z of substrate sites)."""
    ph=phos.copy(); ph.index=[str(s).replace(":","_") for s in ph.index]
    ph=ph[~ph.index.duplicated()]
    Z=ph.sub(ph.mean(axis=0),axis=1).div(ph.std(axis=0)+1e-9,axis=1)   # z within each sample across sites
    acts={}
    for k,targets in k2t.items():
        sites=[t for t in targets if t in Z.index]
        if len(sites)>=min_sub:
            acts[k]=Z.loc[sites].mean(axis=0)
    return pd.DataFrame(acts).T   # kinase x sample

gao=pd.read_parquet(f"{P}/gao_phos_T.parquet"); dong=pd.read_parquet(f"{P}/dong_phos.parquet")
aH=ksea(gao); aI=ksea(dong)
print(f"HCC kinase-activity matrix={aH.shape}  iCCA={aI.shape}")
ks=sorted(set(aH.index)&set(aI.index))
print(f"kinases scored in both cohorts={len(ks)}")

# cross-tumor activity difference (CONFOUNDED) + stats
rows=[]
for k in ks:
    h=aH.loc[k].dropna(); i=aI.loc[k].dropna()
    if len(h)<10 or len(i)<10: continue
    u=stats.mannwhitneyu(i,h)
    rows.append({"kinase":k,"act_diff_iCCA_minus_HCC":float(i.median()-h.median()),
                 "p":float(u.pvalue),"n_sub_HCC":int((gao.index.str.split(":").str[0].isin([k])).sum())})
ph_df=pd.DataFrame(rows)
ph_df["fdr"]=multipletests(ph_df["p"],method="fdr_bh")[1]
print(f"cross-tumor kinase tests={len(ph_df)}  FDR<0.05={(ph_df['fdr']<0.05).sum()}")

# concordance with batch-free TCGA RNA kinase expression (LIHC vs CHOL)
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); chol=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
rna_diff={}
for k in ph_df["kinase"]:
    if k in lihc.index and k in chol.index:
        rna_diff[k]=float(chol.loc[k].mean()-lihc.loc[k].mean())
ph_df["rna_diff_iCCA_minus_HCC"]=ph_df["kinase"].map(rna_diff)
m=ph_df.dropna(subset=["rna_diff_iCCA_minus_HCC"])
rho,pr=stats.spearmanr(m["act_diff_iCCA_minus_HCC"],m["rna_diff_iCCA_minus_HCC"])
sign_agree=float(np.mean(np.sign(m["act_diff_iCCA_minus_HCC"])==np.sign(m["rna_diff_iCCA_minus_HCC"])))
print(f"phospho-activity vs batch-free RNA (kinase level): Spearman rho={rho:.3f} p={pr:.1e}  sign-agree={sign_agree*100:.0f}%  (n={len(m)})")

# robust phospho-supported kinases: significant cross-tumor AND sign-concordant with batch-free RNA
m=m.assign(concordant=(np.sign(m["act_diff_iCCA_minus_HCC"])==np.sign(m["rna_diff_iCCA_minus_HCC"])))
robust=m[(m["fdr"]<0.05)&(m["concordant"])].sort_values("act_diff_iCCA_minus_HCC")
print(f"robust phospho-supported lineage kinases (FDR<0.05 & RNA-concordant)={len(robust)}")
print("  iCCA-high:",list(robust[robust['act_diff_iCCA_minus_HCC']>0]['kinase'].tail(12)))
print("  HCC-high :",list(robust[robust['act_diff_iCCA_minus_HCC']<0]['kinase'].head(12)))

aH.to_parquet(f"{P}/kinase_act_HCC.parquet"); aI.to_parquet(f"{P}/kinase_act_iCCA.parquet")
ph_df.to_parquet(f"{P}/kinase_crosstumor.parquet")
json.dump({"kinases_both":len(ks),"crosstumor_FDR05":int((ph_df['fdr']<0.05).sum()),
 "phospho_vs_RNA_spearman":float(rho),"phospho_vs_RNA_p":float(pr),"sign_agree":sign_agree,
 "n_robust_concordant":int(len(robust))},open(f"{OUT}/kinase_activity_stats.json","w"),indent=2)
robust.to_csv(f"{OUT}/robust_phospho_kinases.csv",index=False)
print("KINASE_MODULE_DONE")
