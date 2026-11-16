"""Alternative random-effects estimators, interval-based reporting and heterogeneity uncertainty (k = 3)."""
import numpy as np, pandas as pd
from scipy import stats, optimize
import robustness
def _w(v,t2): return 1/(v+t2)
def tau_dl(y,v):
    w=1/v; mu=(w*y).sum()/w.sum(); Q=(w*(y-mu)**2).sum(); k=len(y)
    return max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
def tau_pm(y,v):
    k=len(y)
    def f(t2):
        w=_w(v,t2); mu=(w*y).sum()/w.sum(); return (w*(y-mu)**2).sum()-(k-1)
    if f(0)<=0: return 0.0
    hi=1.0
    while f(hi)>0 and hi<1e4: hi*=2
    return optimize.brentq(f,0,hi)
def tau_reml(y,v):
    k=len(y)
    def nll(t2):
        w=_w(v,t2); mu=(w*y).sum()/w.sum(); return 0.5*(np.log(v+t2).sum()+np.log(w.sum())+(w*(y-mu)**2).sum())
    r=optimize.minimize_scalar(nll,bounds=(0,50),method='bounded',options={'xatol':1e-8}); 
    return 0.0 if nll(0)<=r.fun+1e-12 else r.x
def tau_sj(y,v):
    k=len(y); t0=max(((y-y.mean())**2).sum()/k,1e-9)
    w=1/(v/t0+1); mu=(w*y).sum()/w.sum(); return (w*(y-mu)**2).sum()/(k-1)
EST={'DerSimonian-Laird':tau_dl,'REML':tau_reml,'Paule-Mandel':tau_pm,'Sidik-Jonkman':tau_sj}
def pool(y,v,t2,hk):
    k=len(y); w=_w(v,t2); mu=(w*y).sum()/w.sum()
    if hk:
        q=(w*(y-mu)**2).sum()/(k-1)   # Hartung-Knapp variance, as in Section 3.7 (no truncation at 1)
        se=np.sqrt(q/w.sum()); crit=stats.t.ppf(.975,k-1); p=2*stats.t.sf(abs(mu/se),k-1)
    else:
        se=np.sqrt(1/w.sum()); crit=stats.norm.ppf(.975); p=2*stats.norm.sf(abs(mu/se))
    return mu,se,mu-crit*se,mu+crit*se,p
def qprofile(y,v,level=.95):
    """Q-profile CI for tau^2 (Viechtbauer 2007) and implied I^2 CI using the typical within-study variance."""
    k=len(y); w0=1/v; sv=(k-1)*w0.sum()/(w0.sum()**2-(w0**2).sum())
    def Q(t2):
        w=_w(v,t2); mu=(w*y).sum()/w.sum(); return (w*(y-mu)**2).sum()
    a=(1-level)/2; qu=stats.chi2.ppf(1-a,k-1); ql=stats.chi2.ppf(a,k-1)
    def solve(target):
        if Q(0)<target: return 0.0
        hi=1.0
        while Q(hi)>target and hi<1e6: hi*=2
        return optimize.brentq(lambda t:Q(t)-target,0,hi)
    lo=solve(qu); up=solve(ql)
    i2=lambda t: 100*t/(t+sv)
    return lo,up,i2(lo),i2(up),sv
def bayes(y,v,scale=0.5,n=400):
    """Posterior for the pooled effect with flat prior on mu and half-normal(scale) prior on tau (grid over tau, mu integrated analytically)."""
    taus=np.linspace(0,3,n); lp=[]; mus=[]; vs=[]
    for t in taus:
        w=1/(v+t*t); W=w.sum(); mu=(w*y).sum()/W
        lp.append(-0.5*np.log(v+t*t).sum()-0.5*np.log(W)-0.5*(w*(y-mu)**2).sum()-0.5*(t/scale)**2); mus.append(mu); vs.append(1/W)
    lp=np.array(lp); pw=np.exp(lp-lp.max()); pw/=pw.sum(); mus=np.array(mus); vs=np.array(vs)
    xs=np.linspace(y.min()-6*np.sqrt(v.max()+9),y.max()+6*np.sqrt(v.max()+9),4001)
    dens=(pw[:,None]*stats.norm.pdf(xs[None,:],mus[:,None],np.sqrt(vs)[:,None])).sum(0); dens/=dens.sum()
    cdf=np.cumsum(dens); q=lambda a: xs[np.searchsorted(cdf,a)]
    return float((xs*dens).sum()),q(.025),q(.975),float(dens[xs>0].sum()),float((pw*taus**2).sum())
def run():
    C,g61=robustness.load(); rows=[]; cis=[]
    for t in g61:
        y=np.array([C[c].loc[t,'effect_size'] for c in C]); v=np.array([C[c].loc[t,'standard_error']**2 for c in C])
        r={'taxon':t}
        for nm,f in EST.items():
            t2=f(y,v); r[f'tau2 {nm}']=t2
            for hk in (False,True):
                mu,se,lo,hi,p=pool(y,v,t2,hk); k_=f'{nm}{" + HK" if hk else ""}'
                r[f'{k_} effect']=mu; r[f'{k_} lower']=lo; r[f'{k_} upper']=hi; r[f'{k_} P']=p
        lo,up,il,iu,sv=qprofile(y,v); r.update({'tau2 CI lower':lo,'tau2 CI upper':up,'I2 CI lower (%)':il,'I2 CI upper (%)':iu,'typical within-study variance':sv})
        w=1/v; mu=(w*y).sum()/w.sum(); Q=(w*(y-mu)**2).sum(); r['Q']=Q; r['Q P']=stats.chi2.sf(Q,2); r['I2 (Q-based) (%)']=max(0,(Q-2)/Q)*100
        m,l,u,pp,t2b=bayes(y,v); r.update({'Bayes mean':m,'Bayes lower':l,'Bayes upper':u,'Bayes P(effect>0)':pp,'Bayes mean tau2':t2b})
        rows.append(r)
    D=pd.DataFrame(rows).set_index('taxon')
    S=[]
    for nm in EST:
        for hk in (False,True):
            k_=f'{nm}{" + HK" if hk else ""}'; p=D[f'{k_} P'].values; fdr=robustness.bh(p)
            w=D[f'{k_} upper']-D[f'{k_} lower']
            S.append({'Estimator':k_,'Median tau2':D[f'tau2 {nm}'].median(),'Genera with tau2 > 0':int((D[f'tau2 {nm}']>0).sum()),'Median 95% interval width (ln)':w.median(),'Genera with P < 0.05':int((p<.05).sum()),'Expected at P < 0.05':3.05,
                      'Genera whose 95% interval excludes 0':int(((D[f'{k_} lower']>0)|(D[f'{k_} upper']<0)).sum()),'Smallest FDR':fdr.min(),'Genera with FDR < 0.2':int((fdr<.2).sum())})
    k_='Bayesian (half-normal 0.5)'; S.append({'Estimator':k_,'Median tau2':D['Bayes mean tau2'].median(),'Genera with tau2 > 0':np.nan,'Median 95% interval width (ln)':(D['Bayes upper']-D['Bayes lower']).median(),'Genera with P < 0.05':np.nan,'Expected at P < 0.05':np.nan,
        'Genera whose 95% interval excludes 0':int(((D['Bayes lower']>0)|(D['Bayes upper']<0)).sum()),'Smallest FDR':np.nan,'Genera with FDR < 0.2':np.nan})
    return D,pd.DataFrame(S)
if __name__=='__main__':
    D,S=run(); pd.set_option('display.width',250); print(S.round(3).to_string(index=False))
    print(D[['tau2 CI lower','tau2 CI upper','I2 CI lower (%)','I2 CI upper (%)','I2 (Q-based) (%)']].describe().round(2))
    print((D['I2 CI upper (%)']>75).sum(),(D['I2 CI lower (%)']>0).sum())
    print(D.loc[['Terrisporobacter','Alistipes','[Eubacterium]_ruminantium_group'],[c for c in D.columns if 'effect' in c or 'lower' in c or 'upper' in c]].round(2).T)
