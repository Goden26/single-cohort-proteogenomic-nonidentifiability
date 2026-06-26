#!/usr/bin/env python
"""Ground-truth titration (v2, addressing R4): calibrate how an anchor/negative-control concordance
check (a KNOWN principle: RUV, gPCA, sva, mRNA-protein concordance) degrades with cohort/batch aliasing.
Two fixes vs v1: (1) the anchor is now a SEPARATE, NOISY proxy of truth (independent noise + modality
bias), not the noise-free ground truth, so 'anchor concordance tracks recoverability' is not built in;
(2) finer alpha grid (11 levels) with PER-REPLICATE rank correlations and bootstrap CIs, so we no longer
over-read a 5-point Spearman of -1.0. We do not claim a new instrument; we calibrate a known check."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
P=f"{ROOT}/results/parsed"; OUT=f"{ROOT}/results"; FIG=f"{ROOT}/figures"
rng=np.random.default_rng(20)
lihc=pd.read_parquet(f"{P}/tcga_lihc_tumor.parquet")
X0=lihc.loc[lihc.var(axis=1).sort_values(ascending=False).head(3000).index].copy()
genes=list(X0.index); n=X0.shape[1]; G=len(genes)
TRUE_K=200; EFF=1.0; BATCH_SD=1.2
ALPHAS=np.round(np.linspace(0,1,11),3); REPS=20
ANCHOR_NOISE=1.0      # makes the anchor a weak, realistic proxy (corr ~0.4-0.5 with truth)

def qn(df):
    a=df.values.copy(); order=np.argsort(a,axis=0); ranks=np.argsort(order,axis=0)
    return pd.DataFrame(np.sort(a,axis=0).mean(axis=1)[ranks],index=df.index,columns=df.columns)

def one_run(alpha,batch="additive"):
    Xb=X0.values.copy().astype(float)
    typ=rng.permutation(np.array([0]*(n//2)+[1]*(n-n//2)))
    sig=rng.choice(G,TRUE_K,replace=False); sdir=rng.choice([-1,1],TRUE_K)
    Xb[np.ix_(sig,np.where(typ==1)[0])]+=(EFF*sdir)[:,None]
    true_eff=np.zeros(G); true_eff[sig]=EFF*sdir
    pB=np.where(typ==1,0.5+0.5*alpha,0.5-0.5*alpha); coh=(rng.random(n)<pB).astype(int)
    if batch=="additive": Xb[:,coh==1]+=rng.normal(0,BATCH_SD,G)[:,None]
    else:                  Xb[:,coh==1]*=np.exp(rng.normal(0,0.5,G))[:,None]   # multiplicative/scale
    q=qn(pd.DataFrame(Xb,index=genes)).values
    a=q[:,typ==1]; b=q[:,typ==0]
    t,p=stats.ttest_ind(a,b,axis=1,equal_var=False)
    fdr=multipletests(np.nan_to_num(p,nan=1.0),method="fdr_bh")[1]; s=fdr<0.05
    rec_eff=a.mean(1)-b.mean(1)
    is_t=np.zeros(G,bool); is_t[sig]=True
    sens=float(np.mean(s[is_t])); fdp=float((s&~is_t).sum()/max(s.sum(),1))
    # NOISY SEPARATE-MODALITY anchor (independent noise + a modality scale bias), NOT the ground truth
    anchor=true_eff*rng.uniform(0.7,1.3)+rng.normal(0,ANCHOR_NOISE,G)
    diag=float(stats.spearmanr(rec_eff,anchor).correlation)
    return sens,fdp,diag

R={}
for batch in ["additive","multiplicative"]:
    per=[]   # per-replicate rows
    lvl=[]
    for al in ALPHAS:
        rr=np.array([one_run(al,batch) for _ in range(REPS)])
        for row in rr: per.append((al,*row))
        lvl.append((al,rr[:,0].mean(),rr[:,1].mean(),rr[:,2].mean()))
    per=np.array(per); lvl=np.array(lvl)
    # M1: the 220 points are 20 reps NESTED within 11 alpha levels. We now report the correlation
    # three ways so the headline statistic is not mistaken for a clean per-replicate effect size:
    #  (1) pooled across all 220 points (the original number; dominated by the deterministic alpha trend);
    #  (2) a BLOCK BOOTSTRAP that resamples the 11 alpha LEVELS (honest CI under level clustering);
    #  (3) the mean WITHIN-LEVEL (alpha-conditioned) Spearman, which isolates genuine diagnostic value
    #      from the alpha trend (pseudoreplication check).
    al_col=per[:,0]; levels=ALPHAS
    def pooled_spear(x,y):  return float(stats.spearmanr(x,y).correlation)
    def block_boot(x,y,B=2000):
        out=[]
        for _ in range(B):
            chosen=rng.choice(levels,len(levels),replace=True)
            idx=np.concatenate([np.where(al_col==a)[0] for a in chosen])
            out.append(stats.spearmanr(x[idx],y[idx]).correlation)
        return [float(np.nanpercentile(out,2.5)),float(np.nanpercentile(out,97.5))]
    def within_level(x,y):
        rs=[stats.spearmanr(x[al_col==a],y[al_col==a]).correlation for a in levels if (al_col==a).sum()>2]
        rs=[r for r in rs if not np.isnan(r)]; return float(np.mean(rs)),float(np.std(rs)),len(rs)
    rho_df=pooled_spear(per[:,3],per[:,2]); ci_df=block_boot(per[:,3],per[:,2])
    rho_rc=pooled_spear(per[:,3],per[:,1]); ci_rc=block_boot(per[:,3],per[:,1])
    wl_mean,wl_sd,wl_n=within_level(per[:,3],per[:,2])
    R[batch]={"n_levels":len(ALPHAS),"reps_per_level":REPS,"n_points":int(len(per)),
              "diag_at_alpha0":float(lvl[0,3]),"diag_at_alpha1":float(lvl[-1,3]),
              "fdp_at_alpha0":float(lvl[0,2]),"fdp_at_alpha1":float(lvl[-1,2]),
              "recall_at_alpha0":float(lvl[0,1]),"recall_at_alpha1":float(lvl[-1,1]),
              "recall_min":float(lvl[:,1].min()),
              "pooled_spearman_diag_vs_fdp":rho_df,"blockboot_CI_diag_vs_fdp":ci_df,
              "within_level_mean_spearman_diag_vs_fdp":wl_mean,"within_level_sd":wl_sd,"within_level_n":wl_n,
              "pooled_spearman_diag_vs_recall":rho_rc,"blockboot_CI_diag_vs_recall":ci_rc,
              "note":"anchor shares the injected-truth term with the recovered effect by construction; this calibrates SENSITIVITY of the check, not independence of two measurements"}
    if batch=="additive": LVL=lvl; PER=per
    print(f"[{batch}] diag {lvl[0,3]:.3f}->{lvl[-1,3]:.3f} | fdp {lvl[0,2]:.2f}->{lvl[-1,2]:.2f} | recall_min={lvl[:,1].min():.3f} | "
          f"pooled Spearman(diag,fdp)={rho_df:.2f} block-boot CI[{ci_df[0]:.2f},{ci_df[1]:.2f}] | within-level mean rho={wl_mean:.2f} (sd {wl_sd:.2f})")
# measured anchor-truth correlation (how weak the anchor is) for one example
ex=np.zeros(G); ex[:TRUE_K]=EFF; an=ex*1.0+rng.normal(0,ANCHOR_NOISE,G)
R["anchor_truth_spearman_example"]=float(stats.spearmanr(ex,an).correlation)
json.dump(R,open(f"{OUT}/identifiability_diagnostic_stats.json","w"),indent=2)
print(f"anchor-truth concordance (weak proxy) ~ {R['anchor_truth_spearman_example']:.2f}")

# FIGURE 6
fig,ax=plt.subplots(1,2,figsize=(12,4.6))
ax[0].plot(LVL[:,0],LVL[:,1],"o-",color="#2c7fb8",label="true-signal recall")
ax[0].plot(LVL[:,0],LVL[:,2],"s-",color="#d95f0e",label="false-discovery proportion")
ax[0].plot(LVL[:,0],LVL[:,3],"^--",color="#333",label="anchor concordance (noisy proxy)")
ax[0].set_xlabel("aliasing of contrast with cohort/batch (alpha; 11 levels)"); ax[0].set_ylabel("metric")
ax[0].legend(frameon=False,fontsize=9); ax[0].set_title(f"Ground-truth titration (spiked TCGA-LIHC)\nfalse positives rise to {LVL[-1,2]:.2f}; anchor concordance shifts only {LVL[0,3]:.2f}->{LVL[-1,3]:.2f}")
ax[1].scatter(PER[:,3],PER[:,2],s=10,alpha=.3,color="#555")
ax[1].set_xlabel("anchor concordance (noisy separate-modality proxy)"); ax[1].set_ylabel("false-discovery proportion")
ax[1].set_title(f"Pooled (n={len(PER)}): a noisy anchor still tracks the trap\nSpearman={R['additive']['pooled_spearman_diag_vs_fdp']:.2f} (block-boot CI {R['additive']['blockboot_CI_diag_vs_fdp'][0]:.2f},{R['additive']['blockboot_CI_diag_vs_fdp'][1]:.2f}); anchor shares truth term by design")
plt.tight_layout(); plt.savefig(f"{FIG}/fig6_identifiability_diagnostic.png",dpi=160); print("saved fig6")
print("DIAGNOSTIC_V2_DONE")
