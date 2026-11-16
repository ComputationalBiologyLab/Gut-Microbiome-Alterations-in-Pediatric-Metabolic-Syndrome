"""Transparent calculation of the attenuation of effects outside the discovery cohort."""
import numpy as np, pandas as pd
from scipy import stats
import robustness
def run(seed=11,nboot=10000):
    C,_=robustness.load(); names=list(C)
    sig=[]
    for c in names:
        for t,r in C[c].iterrows():
            if r['FDR']<0.1: sig.append((t,c))
    out={'n_sig':len(sig)}
    rows=[]
    for t,c in sig:
        others=[o for o in names if o!=c]
        if all(t in C[o].index for o in others):
            d=C[c].loc[t]; r={'Genus':t,'Discovery cohort':c,'Discovery effect':d.effect_size,'Discovery SE':d.standard_error,'Discovery FDR':d.FDR}
            for o in others: r[f'{o} effect']=C[o].loc[t,'effect_size']; r[f'{o} SE']=C[o].loc[t,'standard_error']
            r['_o']=others; rows.append(r)
    out['n_meas']=len(rows)
    T=[]
    for r in rows:
        o=r['_o']; e_o=np.array([r[f'{x} effect'] for x in o]); s_o=np.array([r[f'{x} SE'] for x in o]); sd=np.sign(r['Discovery effect'])
        T.append({'Genus':r['Genus'],'Discovery cohort':r['Discovery cohort'],'Discovery effect':r['Discovery effect'],'Discovery SE':r['Discovery SE'],'Discovery FDR':r['Discovery FDR'],
          'Other cohort 1':o[0],'Effect, other cohort 1':e_o[0],'SE, other cohort 1':s_o[0],'Other cohort 2':o[1],'Effect, other cohort 2':e_o[1],'SE, other cohort 2':s_o[1],
          'Mean |effect| in other two cohorts':np.abs(e_o).mean(),'Mean effect in other two, oriented to discovery direction':(sd*e_o).mean(),
          'Same direction in both others (n of 2)':int((np.sign(e_o)==sd).sum()),'Expected mean |effect| if true effect were zero':0.798*s_o.mean()})
    T=pd.DataFrame(T); out['table']=T
    T['|Mean effect| in other two cohorts']=np.abs((T['Effect, other cohort 1']+T['Effect, other cohort 2'])/2)
    dm=np.abs(T['Discovery effect']); om=T['|Mean effect| in other two cohorts']; om_sep=T['Mean |effect| in other two cohorts']
    out['disc_mean']=dm.mean(); out['oth_mean']=om.mean(); out['red']=100*(1-om.mean()/dm.mean())
    out['sep_mean']=om_sep.mean(); out['sep_red']=100*(1-om_sep.mean()/dm.mean())
    out['red_median_pergene']=float(100*(1-(om/dm)).median())
    rng=np.random.default_rng(seed); n=len(T); bs=[]
    for _ in range(nboot):
        i=rng.integers(0,n,n); bs.append(100*(1-om.values[i].mean()/dm.values[i].mean()))
    out['red_ci']=np.percentile(bs,[2.5,97.5]).tolist()
    bs2=[];bs3=[];ov=T['Mean effect in other two, oriented to discovery direction'].values
    for _ in range(nboot):
        i=rng.integers(0,n,n); bs2.append(100*(1-om_sep.values[i].mean()/dm.values[i].mean())); bs3.append(100*(1-ov[i].mean()/dm.values[i].mean()))
    out['sep_ci']=np.percentile(bs2,[2.5,97.5]).tolist(); out['or_ci']=np.percentile(bs3,[2.5,97.5]).tolist(); out['or_red']=100*(1-ov.mean()/dm.mean())
    out['oriented']=float(T['Mean effect in other two, oriented to discovery direction'].mean())
    out['oriented_ci']=stats.t.interval(.95,n-1,loc=T['Mean effect in other two, oriented to discovery direction'].mean(),scale=stats.sem(T['Mean effect in other two, oriented to discovery direction'])) 
    out['null_exp']=T['Expected mean |effect| if true effect were zero'].mean()
    # direction in the other cohorts: 20 cohort-genus pairs
    same=int(T['Same direction in both others (n of 2)'].sum()); out['dir_pairs']=(same,2*n); out['dir_p']=stats.binomtest(same,2*n,.5).pvalue
    # first-other cohort direction count as reported in the manuscript (average of the two)
    sd=np.sign(T['Discovery effect']); avg=T['Mean effect in other two, oriented to discovery direction']
    out['dir_avg']=int((avg>0).sum())
    # empirical-Bayes (normal prior, mean zero, variance from all genera of the cohort) shrinkage of discovery estimates -> winner's-curse-corrected
    eb=[]
    for _,r in T.iterrows():
        c=r['Discovery cohort']; e=C[c].effect_size.values; s=C[c].standard_error.values
        t2=max(0,np.var(e,ddof=1)-np.mean(s**2)); f=t2/(t2+r['Discovery SE']**2); eb.append((f,f*r['Discovery effect']))
    T['EB shrinkage factor']=[x[0] for x in eb]; T['Empirical-Bayes discovery effect']=[x[1] for x in eb]
    out['eb_mean_abs']=float(np.abs(T['Empirical-Bayes discovery effect']).mean()); out['eb_f_mean']=float(np.mean([x[0] for x in eb]))
    out['eb_red']=100*(1-om.mean()/np.abs(T['Empirical-Bayes discovery effect']).mean()); out['eb_red_sep']=100*(1-om_sep.mean()/np.abs(T['Empirical-Bayes discovery effect']).mean())
    # what selection alone would produce: genera with true effect zero, selected at |z| >= z_thr in discovery cohort (z at FDR<0.1 boundary), expected |estimate| = SE*E|Z| given |Z|>z_thr
    zthr=np.array([stats.norm.isf(.5*0.1*0+1e-300)]) if False else None
    sel=[]
    for _,r in T.iterrows():
        zc=abs(r['Discovery effect'])/r['Discovery SE']
        # expected |estimate| given selection at the observed z (lower bound on the winner's-curse inflation under no true effect)
        sel.append(r['Discovery SE']*zc)
    out['zmean']=float(np.mean([abs(r['Discovery effect'])/r['Discovery SE'] for _,r in T.iterrows()]))
    # share of the discovery cohorts' own genera with |z| as large: fraction of tested genera selected
    out['n_tested']={c:int(len(C[c])) for c in names}
    # simulation: all true effects zero, same SEs; how many of tested genera would be nominally 'significant' and how large would selected effects be
    out['sim']=None
    out['table']=T
    return out
if __name__=='__main__':
    o=run(); T=o.pop('table')
    pd.set_option('display.width',250); print(T.round(3).to_string())
    for k,v in o.items(): print(k,v)
