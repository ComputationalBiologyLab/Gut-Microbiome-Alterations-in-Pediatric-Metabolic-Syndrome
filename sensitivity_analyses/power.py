"""Power and precision assessment (primary DerSimonian-Laird analysis of 61 genera)."""
import numpy as np, pandas as pd
from scipy import stats
import robustness
Z=stats.norm.ppf
def run(seed=2026,nsim=4000):
    C,g61=robustness.load()
    S2=pd.read_csv(robustness.B+'meta_analysis/three_study_random_effects.tsv',sep='\t').set_index('taxon').loc[g61]
    se=S2.pooled_SE.values; eff=S2.pooled_effect.values; lo=S2.lower_CI.values; hi=S2.upper_CI.values
    out={}
    out['max_abs_eff']=float(np.abs(eff).max()); out['n']=len(g61); out['se_gt_effect']=int((se>np.abs(eff)).sum())
    out['median_se']=float(np.median(se)); out['median_abs_eff']=float(np.median(np.abs(eff)))
    # (a) minimum detectable effect at 80% power
    a05=Z(.975); aB=Z(1-.05/61/2); z80=Z(.8)
    mde05=(a05+z80)*se; mdeB=(aB+z80)*se
    out['mde']=dict(nominal_median=float(np.median(mde05)),nominal_q=np.percentile(mde05,[25,75]).tolist(),bonf_median=float(np.median(mdeB)),bonf_q=np.percentile(mdeB,[25,75]).tolist(),
                    fold_nominal=float(np.exp(np.median(mde05))),fold_bonf=float(np.exp(np.median(mdeB))))
    # (b) power at assumed true effects
    rows=[]
    for d in (0.25,np.log(1.5),np.log(2),1.0):
        for lab,alpha in (('P < 0.05 (nominal)',.05),('P < 0.05/61 (Bonferroni)',.05/61)):
            za=Z(1-alpha/2); pw=stats.norm.cdf(d/se-za)+stats.norm.cdf(-d/se-za)
            rows.append({'True effect (ln fold-change)':round(float(d),3),'Fold change':round(float(np.exp(d)),2),'Test threshold':lab,'Mean power across 61 genera':float(pw.mean()),
                         'Genera with power >= 80% (n)':int((pw>=.8).sum()),'Genera with power >= 80% (%)':float(100*(pw>=.8).mean())})
    out['power']=pd.DataFrame(rows)
    # (c) interval-based evidence for or against effects of a given size (DL z-interval, and Hartung-Knapp-like wider t(2) interval)
    rows=[]
    t2=stats.t.ppf(.95,2); z90=Z(.95)
    for d in (np.log(1.5),np.log(2),1.0):
        for lab,m in (('DerSimonian-Laird (normal)',z90),('Wider interval (t, 2 df)',t2)):
            ci_lo=eff-m*se; ci_hi=eff+m*se; within=((ci_lo>-d)&(ci_hi<d))
            rows.append({'Effect size margin (ln fold-change)':round(float(d),3),'Fold change':round(float(np.exp(d)),2),'Interval':lab,'Genera whose 90% interval lies within the margin (n)':int(within.sum()),'Genera whose 90% interval lies within the margin (%)':float(100*within.mean())})
    out['equiv']=pd.DataFrame(rows)
    out['ci_halfwidth_median']=float(np.median((hi-lo)/2)); out['max_abs_bound_median']=float(np.median(np.maximum(np.abs(lo),np.abs(hi))))
    out['ci_contains_zero']=int(((lo<0)&(hi>0)).sum())
    # (d) concordance: interval for the observed share and implied upper bound on the share of reproducible genera
    k,n=16,61; cp=stats.beta.ppf([.025,.975],[k,k+1],[n-k+1,n-k])
    out['conc']=dict(k=k,n=n,pct=100*k/n,ci=(100*cp[0],100*cp[1]),pi_max=float((cp[1]-.25)/.75))
    kp,npair=93,183; cp2=stats.beta.ppf([.025,.975],[kp,kp+1],[npair-kp+1,npair-kp]); out['pair']=dict(k=kp,n=npair,pct=100*kp/npair,ci=(100*cp2[0],100*cp2[1]))
    # (e) simulation: power of the cross-cohort direction-concordance test
    rng=np.random.default_rng(seed)
    SE3=np.column_stack([C[c].loc[g61,'standard_error'].values for c in ('US','China','Mexico')])
    PV=np.array([stats.binomtest(k_,61,.25).pvalue for k_ in range(62)])
    pvals=lambda kk: PV[int(kk)]
    rows=[]
    for d in (0.25,0.5,0.75,1.0):
        pr_same=np.prod(stats.norm.cdf(d/SE3),axis=1)+np.prod(1-stats.norm.cdf(d/SE3),axis=1)   # P(all three signs agree | true shared effect d)
        for pi in (0.05,0.10,0.20,0.30,0.50):
            m=int(round(pi*61)); hits=0
            for _ in range(nsim):
                idx=rng.permutation(61)[:m]; p=np.full(61,.25); p[idx]=pr_same[idx]
                kk=(rng.random(61)<p).sum()
                if kk>.25*61 and pvals(kk)<.05: hits+=1
            rows.append({'True shared effect (ln fold-change)':d,'Fold change':round(float(np.exp(d)),2),'Share of genera with a shared effect':pi,'Genera (n)':m,'Mean probability that the three cohorts agree in direction':float(pr_same.mean()),'Power of the concordance test':hits/nsim})
    out['conc_power']=pd.DataFrame(rows)
    return out
if __name__=='__main__':
    o=run(nsim=1500)
    for k,v in o.items():
        print('\n==',k); print(v.round(3).to_string(index=False) if isinstance(v,pd.DataFrame) else v)
