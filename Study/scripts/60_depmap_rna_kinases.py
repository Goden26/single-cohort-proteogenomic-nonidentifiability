#!/usr/bin/env python
"""Module 4 (robust core): lineage-differential + selectively-essential kinases across HCC/iCCA
from BATCH-CONTROLLED layers only — TCGA RNA (LIHC vs CHOL, same GDC pipeline) and DepMap CRISPR
(liver vs biliary cell lines, single platform) — then position each kinase in cHCC-CCA (i-CHC RNA)."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"; DATA=os.path.join(os.path.dirname(ROOT),"Data")

ksn=pd.read_parquet(f"{DATA}/omnipath_ksn.parquet")
KIN=sorted({str(x) for x in ksn["enzyme_genesymbol"].dropna().unique()})
print(f"kinases (KSN phospho-enzymes)={len(KIN)}")

# ---------- TCGA RNA differential (batch-free) ----------
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet"); chol=pd.read_parquet(f"{P}/tcga_chol_tumor.parquet")
rows=[]
for k in KIN:
    if k in lihc.index and k in chol.index:
        h=lihc.loc[k].values; i=chol.loc[k].values
        rows.append({"kinase":k,"log2FC_iCCA_HCC":float(i.mean()-h.mean()),
                     "rna_p":float(stats.mannwhitneyu(i,h).pvalue),
                     "hcc_tpm":float(h.mean()),"icca_tpm":float(i.mean())})
rna=pd.DataFrame(rows); rna["rna_fdr"]=np.nan
vp=rna["rna_p"].notna()
rna.loc[vp,"rna_fdr"]=multipletests(rna.loc[vp,"rna_p"],method="fdr_bh")[1]
print(f"TCGA kinase DE: tested={len(rna)}  FDR<0.05 & |log2FC|>1 = {((rna.rna_fdr<0.05)&(rna.log2FC_iCCA_HCC.abs()>1)).sum()}")

# ---------- DepMap dependency (single platform) ----------
model=pd.read_parquet(f"{P}/depmap_model.parquet")
def codes(df,col,vals): return set(df[df[col].isin(vals)]["ModelID"])
hcc_ids=codes(model,"OncotreeCode",["HCC"])
icca_ids=codes(model,"OncotreeCode",["IHCH"])               # intrahepatic cholangiocarcinoma
comb_ids=codes(model,"OncotreeCode",["HCCIHCH"])            # combined HCC-CC (KMCH1)
print(f"DepMap lines: HCC={len(hcc_ids)} iCCA(IHCH)={len(icca_ids)} combined={len(comb_ids)}")

hdr=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",nrows=0)
col2sym={c:c.split(" (")[0] for c in hdr.columns[1:]}
kin_cols=[c for c,s in col2sym.items() if s in set(KIN)]
ce=pd.read_csv(f"{DATA}/depmap_CRISPRGeneEffect.csv",usecols=[hdr.columns[0]]+kin_cols,index_col=0)
ce.columns=[col2sym[c] for c in ce.columns]; ce=ce.loc[:,~ce.columns.duplicated()]
print(f"DepMap CRISPR loaded: {ce.shape[0]} lines x {ce.shape[1]} kinases")
inHCC=ce.index.isin(hcc_ids); inICCA=ce.index.isin(icca_ids); inOther=~(inHCC|inICCA|ce.index.isin(comb_ids))

dep=[]
for k in ce.columns:
    v=ce[k]; h=v[inHCC].dropna(); i=v[inICCA].dropna(); o=v[inOther].dropna()
    if len(h)<3 or len(i)<3: continue
    dep.append({"kinase":k,"chronos_HCC":float(h.mean()),"chronos_iCCA":float(i.mean()),
                "chronos_other":float(o.mean()),
                "icca_selective":float(i.mean()-o.mean()),"hcc_selective":float(h.mean()-o.mean()),
                "dep_p_icca_vs_other":float(stats.mannwhitneyu(i,o).pvalue),
                "dep_p_hcc_vs_other":float(stats.mannwhitneyu(h,o).pvalue),
                "dep_diff_iCCA_minus_HCC":float(i.mean()-h.mean()),               # direct liver-lineage contrast
                "dep_p_icca_vs_hcc":float(stats.mannwhitneyu(i,h).pvalue)})
dep=pd.DataFrame(dep)
print(f"DepMap kinase dependency tested={len(dep)}")

# ---------- integrate ----------
M=rna.merge(dep,on="kinase",how="inner")
# robust selectively-essential lineage kinases
M["iCCA_essential"]=(M["chronos_iCCA"]<-0.4)
M["HCC_essential"]=(M["chronos_HCC"]<-0.4)
M["iCCA_selective_dep"]=(M["icca_selective"]<-0.2)&(M["dep_p_icca_vs_other"]<0.1)
M["HCC_selective_dep"]=(M["hcc_selective"]<-0.2)&(M["dep_p_hcc_vs_other"]<0.1)
M["rna_iCCA_up"]=(M["rna_fdr"]<0.05)&(M["log2FC_iCCA_HCC"]>1)
M["rna_HCC_up"]=(M["rna_fdr"]<0.05)&(M["log2FC_iCCA_HCC"]<-1)

# cHCC-CCA positioning of each kinase (i-CHC RNA percentile vs HCC/iCCA medians)
chcca=pd.read_parquet(f"{P}/chcca_gse241466.parquet")
def chcca_pos(k):
    if k not in chcca.index or k not in lihc.index or k not in chol.index: return np.nan
    h=np.median(lihc.loc[k]); i=np.median(chol.loc[k])
    if abs(i-h)<1.0: return np.nan                       # ill-defined when HCC~iCCA expression
    med=np.median(chcca.loc[k]); lo,hi=sorted([h,i])
    return float(np.clip((med-lo)/(hi-lo),-0.5,1.5))     # 0=HCC-like,1=iCCA-like,~0.5 intermediate
M["chcca_position"]=M["kinase"].map(chcca_pos)
# direct liver-lineage contrast hits: more essential in iCCA than HCC lines AND iCCA-up in RNA
M["icca_over_hcc_dep"]=(M["dep_diff_iCCA_minus_HCC"]<-0.2)&(M["dep_p_icca_vs_hcc"]<0.1)
M["hcc_over_icca_dep"]=(M["dep_diff_iCCA_minus_HCC"]>0.2)&(M["dep_p_icca_vs_hcc"]<0.1)
direct_icca=M[M["icca_over_hcc_dep"]].sort_values("dep_diff_iCCA_minus_HCC")
direct_hcc =M[M["hcc_over_icca_dep"]].sort_values("dep_diff_iCCA_minus_HCC",ascending=False)
print(f"\nDIRECT contrast — kinases more essential in iCCA than HCC lines (Δ<-0.2,p<0.1)={len(direct_icca)}:")
print(direct_icca[["kinase","dep_diff_iCCA_minus_HCC","chronos_iCCA","chronos_HCC","dep_p_icca_vs_hcc","log2FC_iCCA_HCC","rna_fdr"]].head(12).to_string(index=False))
print(f"\nDIRECT contrast — kinases more essential in HCC than iCCA lines={len(direct_hcc)}:")
print(direct_hcc[["kinase","dep_diff_iCCA_minus_HCC","chronos_HCC","chronos_iCCA","dep_p_icca_vs_hcc","log2FC_iCCA_HCC","rna_fdr"]].head(12).to_string(index=False))
direct_icca.to_csv(f"{OUT}/direct_icca_gt_hcc_kinases.csv",index=False); direct_hcc.to_csv(f"{OUT}/direct_hcc_gt_icca_kinases.csv",index=False)

# PRIMARY robust hits = DepMap lineage-selective dependency (single platform, tested vs other lineages)
M["rna_concordant_iCCA"]=(M["rna_fdr"]<0.05)&(M["log2FC_iCCA_HCC"]>0.5)
M["rna_concordant_HCC"] =(M["rna_fdr"]<0.05)&(M["log2FC_iCCA_HCC"]<-0.5)
icca_hits=M[M["iCCA_selective_dep"]].sort_values("icca_selective")
hcc_hits =M[M["HCC_selective_dep"]].sort_values("hcc_selective")
print(f"\niCCA-selective kinase dependencies (DepMap, Δ<-0.2 & p<0.1)={len(icca_hits)};  RNA-concordant={int(icca_hits['rna_concordant_iCCA'].sum())}")
print(icca_hits[["kinase","icca_selective","chronos_iCCA","chronos_other","dep_p_icca_vs_other","log2FC_iCCA_HCC","rna_fdr","chcca_position"]].head(15).to_string(index=False))
print(f"\nHCC-selective kinase dependencies (DepMap)={len(hcc_hits)};  RNA-concordant={int(hcc_hits['rna_concordant_HCC'].sum())}")
print(hcc_hits[["kinase","hcc_selective","chronos_HCC","chronos_other","dep_p_hcc_vs_other","log2FC_iCCA_HCC","rna_fdr","chcca_position"]].head(12).to_string(index=False))

M.to_csv(f"{OUT}/kinase_lineage_integrated.csv",index=False)
icca_hits.to_csv(f"{OUT}/icca_selective_kinases.csv",index=False); hcc_hits.to_csv(f"{OUT}/hcc_selective_kinases.csv",index=False)
json.dump({"n_kinases_integrated":int(len(M)),
 "rna_DE_sig_fdr05_lfc1":int(((rna.rna_fdr<0.05)&(rna.log2FC_iCCA_HCC.abs()>1)).sum()),
 "icca_selective_dep":int(len(icca_hits)),"icca_selective_dep_rna_concordant":int(icca_hits['rna_concordant_iCCA'].sum()),
 "hcc_selective_dep":int(len(hcc_hits)),"hcc_selective_dep_rna_concordant":int(hcc_hits['rna_concordant_HCC'].sum()),
 "depmap_HCC_lines":len(hcc_ids),"depmap_iCCA_lines":len(icca_ids),
 "top_icca_kinases":list(icca_hits["kinase"].head(12)),"top_hcc_kinases":list(hcc_hits["kinase"].head(8))},
 open(f"{OUT}/depmap_rna_kinase_stats.json","w"),indent=2)

# ---------- FIGURE 3: RNA log2FC vs iCCA-selective dependency ----------
from adjustText import adjust_text
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
ax[0].scatter(M["log2FC_iCCA_HCC"],M["icca_selective"],s=12,alpha=.35,color="lightgray")
hi=icca_hits.head(11)
ax[0].scatter(hi["log2FC_iCCA_HCC"],hi["icca_selective"],s=45,color="#d95f0e",edgecolor="k",lw=.4,zorder=3)
txt=[ax[0].text(r["log2FC_iCCA_HCC"],r["icca_selective"],r["kinase"],fontsize=9,fontweight="bold") for _,r in hi.iterrows()]
adjust_text(txt,ax=ax[0],arrowprops=dict(arrowstyle="-",color="gray",lw=.5))
ax[0].axhline(-0.2,ls=":",c="r",lw=.8); ax[0].axvline(0.5,ls=":",c="gray",lw=.8)
ax[0].set_xlabel("TCGA RNA log2FC (iCCA − HCC), batch-free"); ax[0].set_ylabel("DepMap iCCA-selective dependency\n(ΔChronos: iCCA − other lineages)")
ax[0].set_title("iCCA-selective druggable kinases\n(lower-right = iCCA-up RNA + iCCA-essential)")
# top kinases chronos by group
top=pd.concat([icca_hits.head(6),hcc_hits.head(6)])
xb=np.arange(len(top)); w=0.4
ax[1].bar(xb-w/2,top["chronos_HCC"],w,label="HCC lines",color="#2c7fb8")
ax[1].bar(xb+w/2,top["chronos_iCCA"],w,label="iCCA lines",color="#d95f0e")
ax[1].axhline(-0.4,ls=":",c="k",lw=.8,label="essentiality")
ax[1].set_xticks(xb); ax[1].set_xticklabels(top["kinase"],rotation=60,ha="right",fontsize=8)
ax[1].set_ylabel("Mean Chronos gene effect"); ax[1].legend(frameon=False,fontsize=8)
ax[1].set_title("Lineage-selective kinase dependencies (DepMap)")
plt.tight_layout(); plt.savefig(f"{FIG}/fig3_depmap_rna_kinases.png",dpi=160); print("saved fig3")
print("DEPMAP_RNA_MODULE_DONE")
