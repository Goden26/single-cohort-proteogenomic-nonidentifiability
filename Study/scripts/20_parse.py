#!/usr/bin/env python
"""Parse Gao HCC + Dong iCCA proteome/phospho/clinical and DepMap into clean matrices.
Output: results/parsed/*.parquet  (features x samples), plus a manifest."""
import os, re, json
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
DATA=os.path.join(ROOT,"data"); OUT=os.path.join(ROOT,"results","parsed")
os.makedirs(OUT,exist_ok=True)
manifest={}

def dedup_index(df):
    return df[~df.index.duplicated(keep='first')]

# ---------- GAO HBV-HCC ----------
gp=pd.read_excel(f"{DATA}/gao_S1.xlsx",sheet_name="3. 6,478 proteins at gene level",engine="openpyxl")
gp=gp.rename(columns={gp.columns[0]:"gene"}); gp["gene"]=gp["gene"].astype(str).str.strip()
gp=gp[gp["gene"]!="nan"].set_index("gene"); gp=dedup_index(gp)
gp=gp.apply(pd.to_numeric,errors="coerce")
gao_prot_T=gp[[c for c in gp.columns if str(c).startswith("T")]]
print(f"Gao proteome: {gp.shape}  tumor(T) cols={gao_prot_T.shape[1]}")

gph=pd.read_excel(f"{DATA}/gao_S1.xlsx",sheet_name="4. 26,418 phosphosites",engine="openpyxl")
gph=gph.rename(columns={gph.columns[0]:"site"}); gph["site"]=gph["site"].astype(str).str.strip()
gph=gph[gph["site"]!="nan"].set_index("site"); gph=dedup_index(gph)
gph=gph.apply(pd.to_numeric,errors="coerce")
gao_phos_T=gph[[c for c in gph.columns if str(c).startswith("T")]]
print(f"Gao phospho: {gph.shape}  tumor(T) cols={gao_phos_T.shape[1]}")

gcl=pd.read_excel(f"{DATA}/gao_S1.xlsx",sheet_name="1. Clinical and molecular data",engine="openpyxl")
gcl.columns=[str(c).strip() for c in gcl.columns]
gcl=gcl.rename(columns={"Tumor (T) sample ID":"sample","Proteomic subtype":"prot_subtype","mRNA subtype":"mrna_subtype"})
print(f"Gao clinical: {gcl.shape}  prot_subtype values={gcl['prot_subtype'].value_counts().to_dict()}")

# ---------- DONG iCCA ----------
def parse_dong_matrix(sheet, id_row, feat_col, val_start_col):
    raw=pd.read_excel(f"{DATA}/dong_S1.xlsx",sheet_name=sheet,header=None,engine="openpyxl")
    sids=raw.iloc[id_row, val_start_col:].tolist()
    sids=[(str(int(s)) if (isinstance(s,float) and not np.isnan(s)) else (None if (isinstance(s,float) and np.isnan(s)) else str(s))) for s in sids]
    feats=raw.iloc[id_row+1:, feat_col].astype(str).str.strip().tolist()
    vals=raw.iloc[id_row+1:, val_start_col:].copy()
    vals.index=feats; vals.columns=sids
    vals=vals.loc[:, [c for c in vals.columns if c is not None]]
    vals=vals[vals.index!="nan"]
    vals=vals.apply(pd.to_numeric,errors="coerce")
    return dedup_index(vals)

dong_prot=parse_dong_matrix("S1D. proteins expression", id_row=2, feat_col=0, val_start_col=2)
print(f"Dong proteome: {dong_prot.shape}")
dong_phos=parse_dong_matrix("S1E. 18,347 phosphosites", id_row=2, feat_col=0, val_start_col=1)
print(f"Dong phospho: {dong_phos.shape}")
dong_mrna=parse_dong_matrix("S1C. mRNA expression", id_row=1, feat_col=0, val_start_col=1)
print(f"Dong mRNA: {dong_mrna.shape}")

dcl=pd.read_excel(f"{DATA}/dong_S1.xlsx",sheet_name="S1A. Clinical info",header=1,engine="openpyxl")
dcl=dcl.rename(columns={"Patient_ID":"patient"})
print(f"Dong clinical: {dcl.shape}")

# Dong subtype S5E
sraw=pd.read_excel(f"{DATA}/dong_S5_sub.xlsx",sheet_name="S5E.ssGSEA results",header=None,engine="openpyxl")
dong_sub=pd.DataFrame({"patient":sraw.iloc[3:,0].astype(str).str.replace(r"\.0$","",regex=True),
                       "prot_subgroup":sraw.iloc[3:,1].astype(str)}).dropna()
dong_sub=dong_sub[dong_sub["patient"]!="nan"]
print(f"Dong subtype: {dong_sub.shape}  values={dong_sub['prot_subgroup'].value_counts().to_dict()}")

# ---------- DepMap ----------
model=pd.read_csv(f"{DATA}/depmap_Model.csv")
liver=model[model["OncotreeLineage"]=="Liver"]
biliary=model[model["OncotreeLineage"]=="Biliary Tract"]
print(f"DepMap liver lines={len(liver)}  biliary lines={len(biliary)}")
print("  Liver OncotreeCode:",liver["OncotreeCode"].value_counts().to_dict())
print("  Biliary OncotreeCode:",biliary["OncotreeCode"].value_counts().to_dict())

# ---------- SAVE ----------
gao_prot_T.to_parquet(f"{OUT}/gao_prot_T.parquet")
gao_phos_T.to_parquet(f"{OUT}/gao_phos_T.parquet")
gcl.to_parquet(f"{OUT}/gao_clinical.parquet")
dong_prot.to_parquet(f"{OUT}/dong_prot.parquet")
dong_phos.to_parquet(f"{OUT}/dong_phos.parquet")
dong_mrna.to_parquet(f"{OUT}/dong_mrna.parquet")
dcl.to_parquet(f"{OUT}/dong_clinical.parquet")
dong_sub.to_parquet(f"{OUT}/dong_subtype.parquet")
model.to_parquet(f"{OUT}/depmap_model.parquet")

manifest={"gao_prot_T":list(gao_prot_T.shape),"gao_phos_T":list(gao_phos_T.shape),
 "dong_prot":list(dong_prot.shape),"dong_phos":list(dong_phos.shape),"dong_mrna":list(dong_mrna.shape),
 "depmap_liver":int(len(liver)),"depmap_biliary":int(len(biliary)),
 "common_proteins_HCC_iCCA":int(len(set(gao_prot_T.index)&set(dong_prot.index))),
 "common_phosphosites_HCC_iCCA":int(len(set(gao_phos_T.index)&set(dong_phos.index)))}
json.dump(manifest,open(f"{OUT}/manifest.json","w"),indent=2)
print("\n=== MANIFEST ==="); print(json.dumps(manifest,indent=2))
print("PARSE_DONE")
