#!/usr/bin/env python
"""Cross-organ generality: the harmonization trap in a SECOND cancer (kidney).
Two single-cohort, different-lab proteomes (ccRCC, China, TMT/FOT; pRCC, Germany, label-free) are
naively harmonized and tested ccRCC-vs-pRCC; the 'significant' hits are validated against the
batch-free TCGA-KIRC vs KIRP RNA anchor. If the trap generalizes, many hits will be significant yet
fail to validate (low/negative anchor concordance), as in liver."""
import os, json, gzip
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); DATA=os.path.join(os.path.dirname(ROOT),"Data"); OUT=f"{ROOT}/results"

# ---- ccRCC tumor proteome (Matrix MBR: Gene Symbol + _TA tumor columns) ----
cc=pd.read_excel(f"{DATA}/ccRCC_proteome.xlsx",sheet_name="Matrix MBR",engine="openpyxl")
cc=cc.rename(columns={cc.columns[1]:"gene"})
tcols=[c for c in cc.columns if str(c).endswith("_TA")]
ccm=cc[["gene"]+tcols].dropna(subset=["gene"]).set_index("gene")
ccm=ccm.apply(pd.to_numeric,errors="coerce"); ccm=ccm[~ccm.index.duplicated()]
ccm=np.log2(ccm.clip(lower=0)+1)
print(f"ccRCC tumor proteome: {ccm.shape[0]} genes x {ccm.shape[1]} tumors")

# ---- pRCC tumor proteome (proteinGroups.txt: Gene names + Intensity TT* columns) ----
pr=pd.read_csv(f"{DATA}/pRCC_proteinGroups.txt",sep="\t",low_memory=False)
ttcols=[c for c in pr.columns if c.startswith("Intensity TT")]
pr["gene"]=pr["Gene names"].astype(str).str.split(";").str[0]
prm=pr[["gene"]+ttcols].replace(0,np.nan).dropna(subset=["gene"])
prm=prm[prm["gene"]!="nan"].set_index("gene"); prm=prm[~prm.index.duplicated()]
prm=np.log2(prm.apply(pd.to_numeric,errors="coerce")+1)
print(f"pRCC tumor proteome: {prm.shape[0]} genes x {prm.shape[1]} tumors")

common=sorted(set(ccm.index)&set(prm.index))
print(f"common proteins ccRCC vs pRCC = {len(common)}")
H=ccm.loc[common]; I=prm.loc[common]
merged=pd.concat([H,I],axis=1); merged=merged.apply(lambda r:r.fillna(r.median()),axis=1)

def quantile_normalize(df):
    a=df.values.copy(); order=np.argsort(a,axis=0); ranks=np.argsort(order,axis=0)
    ms=np.sort(a,axis=0).mean(axis=1); return pd.DataFrame(ms[ranks],index=df.index,columns=df.columns)
qn=quantile_normalize(merged)
isC=np.array([True]*H.shape[1]+[False]*I.shape[1])     # True = ccRCC
res=[]
for g in qn.index:
    a=qn.loc[g,isC].values; b=qn.loc[g,~isC].values
    res.append((g,float(a.mean()-b.mean()),float(stats.mannwhitneyu(a,b).pvalue)))
df=pd.DataFrame(res,columns=["gene","diff_cc_pr","p"]).dropna()
df["fdr"]=multipletests(df["p"],method="fdr_bh")[1]
sig=df[df["fdr"]<0.05]
print(f"\nNaive harmonization (kidney): {len(sig)} of {len(df)} proteins 'significant' ccRCC-vs-pRCC at FDR<0.05 ({100*len(sig)/len(df):.0f}%)")

# ---- batch-free anchor: TCGA-KIRC vs KIRP RNA ----
pm=pd.read_csv(f"{DATA}/gencode.v36.probemap",sep="\t"); pm["base"]=pm["id"].str.split(".").str[0]
e2s=dict(zip(pm["base"],pm["gene"]))
def load_tcga(fn):
    d=pd.read_csv(f"{DATA}/{fn}",sep="\t",index_col=0)
    d.index=[e2s.get(str(i).split(".")[0],i) for i in d.index]; d=d.groupby(level=0).max()
    return d[[c for c in d.columns if c.split("-")[3][:2]=="01"]]
kirc=load_tcga("TCGA-KIRC.tpm.tsv.gz"); kirp=load_tcga("TCGA-KIRP.tpm.tsv.gz")
print(f"TCGA-KIRC tumors={kirc.shape[1]} KIRP tumors={kirp.shape[1]}")
g2=sorted(set(sig["gene"])&set(kirc.index)&set(kirp.index))
rna_dir=(kirc.loc[g2].mean(axis=1)-kirp.loc[g2].mean(axis=1))     # ccRCC(KIRC)-pRCC(KIRP), batch-free
prot_dir=df.set_index("gene").loc[g2,"diff_cc_pr"]
rho,pr_p=stats.spearmanr(prot_dir,rna_dir)
sign=float(np.mean(np.sign(prot_dir)==np.sign(rna_dir)))
top=list(sig.reindex(sig["diff_cc_pr"].abs().sort_values(ascending=False).index)["gene"].head(10))
print(f"VALIDATION of {len(g2)} 'significant' kidney hits vs batch-free TCGA-KIRC/KIRP RNA: Spearman rho={rho:.3f} (p={pr_p:.1g}); sign agreement {sign*100:.0f}%")
print(f"  example plausible-but-false hits: {top}")

R={"organ":"kidney","types":"ccRCC vs pRCC","ccRCC_tumors":int(H.shape[1]),"pRCC_tumors":int(I.shape[1]),
   "common_proteins":len(common),"n_tested":int(len(df)),"n_significant_fdr05":int(len(sig)),
   "pct_significant":float(len(sig)/len(df)),"anchor_concordance_rho":float(rho),"anchor_sign_agreement":sign,
   "example_false_hits":top,"tcga_KIRC_n":int(kirc.shape[1]),"tcga_KIRP_n":int(kirp.shape[1])}
json.dump(R,open(f"{OUT}/kidney_generality_stats.json","w"),indent=2)
print("\nKIDNEY_GENERALITY_DONE")
