"""Assumptions behind the 25% chance benchmark for cross-cohort direction concordance."""
import numpy as np, pandas as pd
from scipy import stats
import robustness
B=robustness.B
COH={'US':'PRJNA1016449','China':'PRJNA838778','Mexico':'PRJNA819279'}
def counts(c,g61):
    d=pd.read_csv(f'{B}{COH[c]}_analysis/results/genus_abundance.tsv',sep='\t',skiprows=1,index_col=0)
    d.index=[i.split('g__')[-1] for i in d.index]
    d=d.groupby(level=0).sum()
    md=pd.read_csv(f'{B}{COH[c]}_analysis/00_metadata/sample_metadata.tsv',sep='\t')
    idc='sample_id'; md=md.drop_duplicates(idc).set_index(idc)
    lab=md.reindex(d.columns)['mets_status']
    keep=lab.isin(['MetS']) | lab.astype(str).str.lower().str.contains('control|non')
    d=d.loc[:,keep.values]; y=(lab[keep]=='MetS').values.astype(int)
    return d.reindex(g61).fillna(0).values, y
def eff(X,y):
    # CLR on the 61 shared genera (pseudocount 0.5), mean difference MetS - comparator
    L=np.log(X+0.5); L=L-L.mean(axis=0)
    a=L[:,y==1]; b=L[:,y==0]
    d=a.mean(1)-b.mean(1); se=np.sqrt(a.var(1,ddof=1)/a.shape[1]+b.var(1,ddof=1)/b.shape[1])
    return d,se
def run(nperm=5000,seed=7):
    C,g61=robustness.load(); rng=np.random.default_rng(seed); out={}
    E=np.column_stack([C[c].loc[g61,'effect_size'].values for c in C]); SE=np.column_stack([C[c].loc[g61,'standard_error'].values for c in C])
    n=len(g61); conc=lambda M: ((M>0).all(1)|(M<0).all(1))
    k_obs=int(conc(E).sum()); out['k_obs']=k_obs; out['n']=n
    # (a) null probability using each cohort's observed share of positive effects
    pp=(E>0).mean(0); p_marg=float(np.prod(pp)+np.prod(1-pp)); out['pos_share']=dict(zip(C,pp.tolist())); out['p_marg']=p_marg
    out['binom_marg']=stats.binomtest(k_obs,n,p_marg).pvalue
    # (b) precision strata: genera whose SE is below / above the cross-cohort median in every cohort (by rank sum)
    z=np.abs(E)/SE; rk=(SE/np.nanmedian(SE,axis=0)).mean(1); hi=rk<=np.median(rk)
    rows=[]
    for lab,m in (('More precise half of genera',hi),('Less precise half of genera',~hi)):
        k=int(conc(E[m]).sum()); nn=int(m.sum()); rows.append({'Stratum':lab,'Genera (n)':nn,'Concordant (n)':k,'Concordant (%)':100*k/nn,'Binomial P (vs 25%)':stats.binomtest(k,nn,.25).pvalue,'Median SE':float(np.median(SE[m]))})
    prev=np.column_stack([C[c].loc[g61,'prevalence'].values for c in C]).mean(1); hp=prev>=np.median(prev)
    for lab,m in (('More prevalent half of genera',hp),('Less prevalent half of genera',~hp)):
        k=int(conc(E[m]).sum()); nn=int(m.sum()); rows.append({'Stratum':lab,'Genera (n)':nn,'Concordant (n)':k,'Concordant (%)':100*k/nn,'Binomial P (vs 25%)':stats.binomtest(k,nn,.25).pvalue,'Median SE':float(np.median(SE[m]))})
    out['strata']=pd.DataFrame(rows)
    # correlated cohorts: equicorrelated Gaussian latent signs, P(all same sign)=1/4+3/(4pi)*arcsin(r)
    rr=[0,.05,.1,.2,.3,.5]; rows=[]
    for r in rr:
        p0=.25+3/(4*np.pi)*np.arcsin(r); rows.append({'Latent between-cohort correlation under the null':r,'Expected share concordant':p0,'Expected genera (n)':p0*n,'Binomial P for 16 of 61 (upper tail)':stats.binomtest(k_obs,n,p0,alternative='greater').pvalue})
    out['corr']=pd.DataFrame(rows)
    # (c) sign-flip null on the reported effects: random sign per cohort and genus, magnitudes retained -> 25% by construction; check spread of k
    # (d) sample-label permutation within each cohort (keeps precision, zero-inflation, genus-genus correlation)
    data={c:counts(c,g61) for c in C}
    obs=[eff(*data[c]) for c in C]; E2=np.column_stack([o[0] for o in obs]); out['k_obs_clr']=int(conc(E2).sum())
    out['clr_agree_with_ancom']=[float((np.sign(E2[:,i])==np.sign(E[:,i])).mean()) for i in range(3)]
    ks=np.empty(nperm,int); ks_pos=[]
    for b in range(nperm):
        cols=[]
        for c in C:
            X,y=data[c]; cols.append(eff(X,rng.permutation(y))[0])
        ks[b]=conc(np.column_stack(cols)).sum()
    out['perm_mean_k']=float(ks.mean()); out['perm_mean_pct']=float(100*ks.mean()/n); out['perm_sd']=float(ks.std())
    out['perm_q']=np.percentile(ks,[2.5,50,97.5]).tolist()
    out['perm_p_ge']=float(((ks>=out['k_obs_clr']).sum()+1)/(nperm+1)); out['perm_p_ge_ancom']=float(((ks>=k_obs).sum()+1)/(nperm+1))
    out['perm_p_le']=float(((ks<=out['k_obs_clr']).sum()+1)/(nperm+1))
    out['perm_dispersion']=float(ks.var()/(n*out['perm_mean_k']/n*(1-out['perm_mean_k']/n)))  # variance relative to binomial
    out['ks']=ks
    # (e) cohort-level dependence: correlation between cohorts' permuted-null effects is zero by design; report observed cross-cohort correlation of effects and of SEs
    out['cor_eff']={f'{a}-{b}':float(stats.spearmanr(E[:,i],E[:,j])[0]) for (i,a),(j,b) in [((0,'US'),(1,'China')),((0,'US'),(2,'Mexico')),((1,'China'),(2,'Mexico'))]}
    out['cor_se']={f'{a}-{b}':float(stats.spearmanr(SE[:,i],SE[:,j])[0]) for (i,a),(j,b) in [((0,'US'),(1,'China')),((0,'US'),(2,'Mexico')),((1,'China'),(2,'Mexico'))]}
    # (f) dependence-robust null: exact P if the cohorts' signs were correlated with pairwise sign-agreement q (exchangeable-ish): bound
    return out
if __name__=='__main__':
    o=run(nperm=2000); ks=o.pop('ks')
    for k,v in o.items(): print(k,':',v.round(3).to_string(index=False) if isinstance(v,pd.DataFrame) else v)
