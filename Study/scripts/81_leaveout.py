#!/usr/bin/env python
"""Anti-circularity robustness: does the amino-acid-transport program stay iCCA-preferential
AFTER removing the discovery genes (SLC3A2, SLC7A5, TFRC)? If yes, the program signal is not
merely the cherry-picked genes. Tests dependency (competitive permutation) + batch-free expression."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; DATA=os.path.join(os.path.dirname(ROOT),"Data")
# M7: identical universe + seeds to 80_transport_deepdive.py so AAT_full perm p matches byte-for-byte.
BG_SEED=2024
IRON=["TFRC","TFR2","SLC11A2","STEAP1","STEAP2","STEAP3","ACO1","IREB2","FTH1","FTL","SLC40A1",
      "HMOX1","CYBRD1","SLC25A37","SLC25A28","HAMP","SLC46A1","NCOA4","FXN"]
AAT=["SLC7A5","SLC3A2","SLC1A5","SLC7A11","SLC38A1","SLC38A2","SLC38A5","SLC7A6","SLC7A8","SLC16A1","SLC16A3","SLC2A1","GLS","GOT2"]
DISCOVERY=["TFRC","SLC3A2","SLC7A5"]   # genes that came out of the direct-contrast selection (+ LAT1 partner)
PERM_SEEDS={"AAT_full":202,"AAT_minus_SLC3A2":203,"AAT_minus_LAT1(SLC3A2+SLC7A5)":204}  # AAT_full == aa_transport in 80

model=pd.read_parquet(f"{P}/depmap_model.parquet")
hcc=set(model[model["OncotreeCode"]=="HCC"]["ModelID"]); icca=set(model[model["OncotreeCode"]=="IHCH"]["ModelID"])
hdr=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",nrows=0)
sym2col={}
for c in hdr.columns[1:]: sym2col.setdefault(c.split(" (")[0],c)
allsyms=sorted(sym2col.keys()); bg=list(np.random.default_rng(BG_SEED).choice(allsyms,3000,replace=False))
UNIVERSE=sorted(set(IRON+AAT+DISCOVERY+bg)&set(sym2col))   # identical to 80_transport_deepdive.py
ce=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",usecols=[hdr.columns[0]]+[sym2col[g] for g in UNIVERSE],index_col=0)
ce.columns=[c.split(" (")[0] for c in ce.columns]; ce=ce.loc[:,~ce.columns.duplicated()]
inH=ce.index.isin(hcc); inI=ce.index.isin(icca)
diff=pd.Series({g:(ce[g][inI].dropna().mean()-ce[g][inH].dropna().mean()) for g in ce.columns
                if ce[g][inI].notna().sum()>=3 and ce[g][inH].notna().sum()>=3}).dropna()

lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); chol=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
def expr_p(genes):
    gs=[g for g in genes if g in lihc.index and g in chol.index]
    allg=pd.concat([lihc.loc[gs],chol.loc[gs]],axis=1)
    z=allg.sub(allg.mean(axis=1),axis=0).div(allg.std(axis=1)+1e-9,axis=0)
    sH=z.iloc[:,:lihc.shape[1]].mean(); sI=z.iloc[:,lihc.shape[1]:].mean()
    dl=chol.loc[gs].mean(axis=1)-lihc.loc[gs].mean(axis=1)
    return float(stats.mannwhitneyu(sI,sH).pvalue), len(gs), float((dl>0).mean())
def dep_perm(genes,seed,B=5000):                     # deterministic competitive permutation (M7)
    gs=[g for g in genes if g in diff.index]; obs=float(diff[gs].mean()); rg=np.random.default_rng(seed)
    null=np.array([diff.sample(len(gs),random_state=int(rg.integers(1_000_000_000))).mean() for _ in range(B)])
    return float((null<=obs).mean()), obs, len(gs)

variants={"AAT_full":AAT,
          "AAT_minus_SLC3A2":[g for g in AAT if g!="SLC3A2"],
          "AAT_minus_LAT1(SLC3A2+SLC7A5)":[g for g in AAT if g not in ("SLC3A2","SLC7A5")]}
R={}
for name,genes in variants.items():
    pp,obs,nd=dep_perm(genes,PERM_SEEDS[name]); pe,ne,frac_up=expr_p(genes)
    R[name]={"n_genes_dep":nd,"dep_competitive_perm_p":pp,"dep_mean_diff":obs,"tcga_expr_p":pe,
             "tcga_expr_n_genes":ne,"tcga_frac_genes_iCCA_up":frac_up}
    print(f"{name:32s} (n={nd}): dependency perm p={pp:.4f} (Δ={obs:+.3f}) | TCGA expr p={pe:.2e} ({frac_up*100:.0f}% iCCA-up)")
json.dump(R,open(f"{OUT}/leaveout_stats.json","w"),indent=2)
print("\nLEAVEOUT_DONE")
