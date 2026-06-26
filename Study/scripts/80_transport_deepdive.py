#!/usr/bin/env python
"""Deep-dive on the iron / amino-acid-transport iCCA-preferential signal (TFRC, SLC3A2, ...).
Triangulate across independent batch-controlled layers:
A. DepMap gene-SET dependency: iron & AA-transport programs, iCCA vs HCC lines + competitive permutation null.
B. Batch-free TCGA expression corroboration (CHOL vs LIHC) for the same programs.
C. Co-dependency module structure across liver+biliary lines.
All honest: report effect sizes, p, and whether the program-level claim holds."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"; DATA=os.path.join(os.path.dirname(ROOT),"Data")
rng=np.random.default_rng(7)
# M7: shared deterministic background + per-set permutation seeds so the amino-acid-transport (AAT)
# competitive-permutation p is BYTE-IDENTICAL here and in 81_leaveout.py (same universe, same seed).
BG_SEED=2024; DISCOVERY=["TFRC","SLC3A2","SLC7A5"]
PERM_SEEDS={"iron_handling":101,"aa_transport":202,"transport_combined":303}

IRON=["TFRC","TFR2","SLC11A2","STEAP1","STEAP2","STEAP3","ACO1","IREB2","FTH1","FTL","SLC40A1",
      "HMOX1","CYBRD1","SLC25A37","SLC25A28","HAMP","SLC46A1","NCOA4","FXN"]
AAT=["SLC7A5","SLC3A2","SLC1A5","SLC7A11","SLC38A1","SLC38A2","SLC38A5","SLC7A6","SLC7A8",
     "SLC16A1","SLC16A3","SLC2A1","GLS","GOT2"]
SETS={"iron_handling":IRON,"aa_transport":AAT,"transport_combined":sorted(set(IRON+AAT))}
def comp_perm(diff, genes, seed, B=5000):   # competitive permutation null with a deterministic seed
    gs=[g for g in genes if g in diff.index]; obs=float(diff[gs].mean()); rg=np.random.default_rng(seed)
    null=np.array([diff.sample(len(gs),random_state=int(rg.integers(1_000_000_000))).mean() for _ in range(B)])
    return float((null<=obs).mean()), obs, len(gs)

# ---- DepMap line groups ----
model=pd.read_parquet(f"{P}/depmap_model.parquet")
hcc=set(model[model["OncotreeCode"]=="HCC"]["ModelID"]); icca=set(model[model["OncotreeCode"]=="IHCH"]["ModelID"])

# ---- read CRISPRGeneEffect: set genes + random background for competitive null ----
hdr=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",nrows=0)
sym2col={}
for c in hdr.columns[1:]:
    sym2col.setdefault(c.split(" (")[0],c)
allsyms=sorted(sym2col.keys())                       # deterministic order
bg=list(np.random.default_rng(BG_SEED).choice(allsyms,size=min(3000,len(allsyms)),replace=False))
UNIVERSE=sorted(set(IRON+AAT+DISCOVERY+bg)&set(sym2col))   # identical universe to 81_leaveout.py
cols=[hdr.columns[0]]+[sym2col[g] for g in UNIVERSE]
ce=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",usecols=cols,index_col=0)
ce.columns=[c.split(" (")[0] for c in ce.columns]; ce=ce.loc[:,~ce.columns.duplicated()]
inH=ce.index.isin(hcc); inI=ce.index.isin(icca)
print(f"DepMap: HCC lines={inH.sum()}  iCCA lines={inI.sum()}  genes loaded={ce.shape[1]}")

# per-gene iCCA-minus-HCC dependency difference (negative = more essential in iCCA)
def gene_diff(g):
    if g not in ce.columns: return np.nan
    i=ce[g][inI].dropna(); h=ce[g][inH].dropna()
    return i.mean()-h.mean() if len(i)>=3 and len(h)>=3 else np.nan
diff_all=pd.Series({g:gene_diff(g) for g in ce.columns}).dropna()
bg_diffs=diff_all  # background distribution of iCCA-HCC dep differences

A={}
for name,genes in SETS.items():
    # competitive permutation with a deterministic per-set seed (M7: matches 81_leaveout.py exactly)
    pcomp,obs,_=comp_perm(diff_all,genes,PERM_SEEDS[name])
    gs=[g for g in genes if g in diff_all.index]
    # line-level MWU: per-line mean Chronos over set
    lineH=ce.loc[inH,gs].mean(axis=1).dropna(); lineI=ce.loc[inI,gs].mean(axis=1).dropna()
    mwu=stats.mannwhitneyu(lineI,lineH)
    A[name]={"n_genes":len(gs),"mean_iCCA_minus_HCC_dep":float(obs),"competitive_perm_p":pcomp,
             "line_MWU_p":float(mwu.pvalue),"line_iCCA_meanChronos":float(lineI.mean()),"line_HCC_meanChronos":float(lineH.mean())}
    print(f"A. [{name}] n={len(gs)} mean ΔdepiCCA-HCC={obs:+.3f} competitive p={pcomp:.4f} | line MWU p={mwu.pvalue:.3f} (iCCA {lineI.mean():.2f} vs HCC {lineH.mean():.2f})")
# per-gene table for the combined set
pg=pd.DataFrame({"gene":SETS["transport_combined"]})
pg["dep_iCCA_minus_HCC"]=pg["gene"].map(lambda g:diff_all.get(g,np.nan))
pg=pg.dropna().sort_values("dep_iCCA_minus_HCC")
pg.to_csv(f"{OUT}/transport_pergene_dep.csv",index=False)
print("   most iCCA-preferential transport genes:",list(pg.head(8)["gene"]))

# ---- B. batch-free TCGA expression corroboration ----
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); chol=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
B={}
for name,genes in SETS.items():
    gs=[g for g in genes if g in lihc.index and g in chol.index]
    dl=chol.loc[gs].mean(axis=1)-lihc.loc[gs].mean(axis=1)  # per-gene log2FC iCCA-HCC
    # gene-set expression: per-sample mean z, compare
    allg=pd.concat([lihc.loc[gs],chol.loc[gs]],axis=1)
    z=allg.sub(allg.mean(axis=1),axis=0).div(allg.std(axis=1)+1e-9,axis=0)
    sH=z.iloc[:,:lihc.shape[1]].mean(); sI=z.iloc[:,lihc.shape[1]:].mean()
    mwu=stats.mannwhitneyu(sI,sH)
    B[name]={"n_genes":len(gs),"median_log2FC_iCCA_HCC":float(dl.median()),"set_expr_MWU_p":float(mwu.pvalue),
             "frac_genes_iCCA_up":float((dl>0).mean())}
    print(f"B. [{name}] TCGA expr: median log2FC(iCCA-HCC)={dl.median():+.2f}  set-score MWU p={mwu.pvalue:.1e}  {(dl>0).mean()*100:.0f}% genes iCCA-up")

# ---- C. co-dependency module across liver+biliary lines ----
liv=ce[inH|inI]
core=[g for g in ["TFRC","SLC3A2","SLC7A5","SLC1A5","SLC11A2","STEAP3","SLC7A11","SLC2A1","ACO1","FTH1"] if g in liv.columns]
cc=liv[core].corr(method="spearman")
# mean off-diagonal corr among the iron/AA core
od=cc.values[np.triu_indices_from(cc.values,1)]
C={"core_genes":core,"mean_codependency_spearman":float(np.nanmean(od)),
   "TFRC_SLC3A2_corr":float(cc.loc["TFRC","SLC3A2"]) if "TFRC" in cc.index and "SLC3A2" in cc.columns else None}
print(f"C. co-dependency module (n={len(core)} genes, {inH.sum()+inI.sum()} liver+biliary lines): mean Spearman={C['mean_codependency_spearman']:.2f}; TFRC-SLC3A2 r={C['TFRC_SLC3A2_corr']}")

json.dump({"A_depmap_geneset":A,"B_tcga_expr":B,"C_codependency":C},open(f"{OUT}/transport_deepdive_stats.json","w"),indent=2)

# ---- FIGURE 4 ----
fig,ax=plt.subplots(1,3,figsize=(16,4.8))
# panel 1: per-gene iCCA-HCC dependency diff (combined set)
pgs=pg.copy(); colors=["#d95f0e" if v<0 else "#2c7fb8" for v in pgs["dep_iCCA_minus_HCC"]]
ax[0].barh(range(len(pgs)),pgs["dep_iCCA_minus_HCC"],color=colors)
ax[0].set_yticks(range(len(pgs))); ax[0].set_yticklabels(pgs["gene"],fontsize=7)
ax[0].axvline(0,c="k",lw=.6); ax[0].set_xlabel("DepMap ΔChronos (iCCA − HCC)\n← more essential in iCCA")
ax[0].set_title("Iron/AA-transport gene dependency\n(iCCA vs HCC lines)")
# panel 2: line-level set Chronos iCCA vs HCC for combined set
gs=[g for g in SETS["transport_combined"] if g in ce.columns]
lineH=ce.loc[inH,gs].mean(axis=1).dropna(); lineI=ce.loc[inI,gs].mean(axis=1).dropna()
ax[1].boxplot([lineH,lineI],tick_labels=[f"HCC\n(n={len(lineH)})",f"iCCA\n(n={len(lineI)})"])
ax[1].set_ylabel("Mean Chronos over transport program")
ax[1].set_title(f"Program-level dependency\nMWU p={A['transport_combined']['line_MWU_p']:.3f}, perm p={A['transport_combined']['competitive_perm_p']:.3f}")
# panel 3: co-dependency heatmap
im=ax[2].imshow(cc.values,cmap="RdBu_r",vmin=-1,vmax=1)
ax[2].set_xticks(range(len(core))); ax[2].set_xticklabels(core,rotation=90,fontsize=7)
ax[2].set_yticks(range(len(core))); ax[2].set_yticklabels(core,fontsize=7)
ax[2].set_title("Co-dependency (liver+biliary lines)")
plt.colorbar(im,ax=ax[2],fraction=0.046)
plt.tight_layout(); plt.savefig(f"{FIG}/fig4_transport_deepdive.png",dpi=160); print("saved fig4")
print("TRANSPORT_DEEPDIVE_DONE")
