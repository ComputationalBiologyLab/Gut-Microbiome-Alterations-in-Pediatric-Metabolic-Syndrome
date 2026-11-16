"""Robustness analyses (diagnostic definitions / comparator groups).
Recomputed from the per-cohort ANCOM-BC2 outputs so the supplementary table is reproducible."""
import itertools, numpy as np, pandas as pd
from scipy import stats
import os
B=os.environ.get('MS_ROOT','.').rstrip('/')+'/'
def load():
    us=pd.read_csv(B+'PRJNA1016449_analysis/08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv',sep='\t').set_index('taxon')
    cn=pd.read_csv(B+'PRJNA838778_analysis/08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv',sep='\t').set_index('taxon')
    mx=pd.read_csv(B+'PRJNA819279_analysis/08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv',sep='\t').set_index('taxon')
    g61=list(pd.read_csv(B+'meta_analysis/three_study_random_effects.tsv',sep='\t').taxon)
    ok=lambda d: d[np.isfinite(d.effect_size)&np.isfinite(d.standard_error)&(d.standard_error>0)]
    return {'US':ok(us),'China':ok(cn),'Mexico':ok(mx)}, g61
def bh(p):
    p=np.array(p);n=len(p);o=np.argsort(p);r=np.empty(n);prev=1
    for i,idx in enumerate(o[::-1]):
        k=n-i;prev=min(prev,p[idx]*n/k);r[idx]=prev
    return r
def ivw(y,v):
    y=np.array(y);v=np.array(v);w=1/v;mu=(w*y).sum()/w.sum();se=np.sqrt(1/w.sum());return mu,se,2*(1-stats.norm.cdf(abs(mu/se)))
def dl(y,v,hk=False):
    y=np.array(y);v=np.array(v);w=1/v;mu=(w*y).sum()/w.sum();Q=(w*(y-mu)**2).sum();k=len(y)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()));ws=1/(v+tau2);mr=(ws*y).sum()/ws.sum()
    if hk:
        q=(ws*(y-mr)**2).sum()/(k-1);se=np.sqrt(q/ws.sum());p=2*(1-stats.t.cdf(abs(mr/se),k-1))
    else:
        se=np.sqrt(1/ws.sum());p=2*(1-stats.norm.cdf(abs(mr/se)))
    return mr,se,p
def pool(C,g,cohorts,fn):
    r=[]
    for t in g:
        y=[C[c].loc[t,'effect_size'] for c in cohorts]; v=[C[c].loc[t,'standard_error']**2 for c in cohorts]
        mu,se,p=fn(y,v); r.append((t,mu,se,p))
    d=pd.DataFrame(r,columns=['taxon','pooled','se','p']); d['FDR']=bh(d.p); return d
def run():
    C,g61=load(); out={}
    # A. pairwise agreement
    rows=[]
    for a,b in itertools.combinations(C,2):
        for label,genes in (('61 genera shared by all three cohorts',g61),('all genera shared by the pair',sorted(set(C[a].index)&set(C[b].index)))):
            x=C[a].loc[genes,'effect_size'];y=C[b].loc[genes,'effect_size'];n=len(genes);k=int((np.sign(x)==np.sign(y)).sum())
            rho,pr=stats.spearmanr(x,y)
            rows.append({'Cohort pair':f'{a} vs {b}','Genus set':label,'Genera (n)':n,'Same direction (n)':k,'Same direction (%)':round(100*k/n,1),
                         'Binomial P (vs 50%)':stats.binomtest(k,n,.5).pvalue,'Spearman rho':rho,'Spearman P':pr})
    out['pairwise']=pd.DataFrame(rows)
    # B. estimators and leave-one-out
    rows=[]
    def add(name,cohorts,fn):
        d=pool(C,g61,cohorts,fn)
        conc=int(sum(abs(np.sign([C[c].loc[t,'effect_size'] for c in cohorts]).sum())==len(cohorts) for t in g61))
        rows.append({'Analysis':name,'Cohorts':'+'.join(cohorts),'Smallest FDR':d.FDR.min(),'Genera with FDR < 0.2':int((d.FDR<.2).sum()),
                     'Genera with P < 0.05':int((d.p<.05).sum()),'Expected at P < 0.05':round(.05*len(d),2),'Same direction in all pooled cohorts':f'{conc}/{len(d)}'})
        return d
    all3=list(C)
    add('DerSimonian-Laird random effects (primary analysis)',all3,dl)
    add('Hartung-Knapp random effects',all3,lambda y,v:dl(y,v,True))
    fe=add('Fixed effect (inverse variance)',all3,ivw)
    for left in all3:
        co=[c for c in all3 if c!=left]; add(f'Leave out {left} (fixed effect, k = 2)',co,ivw)
    out['estimators']=pd.DataFrame(rows)
    # C. BMI-balanced pair vs Mexico
    pc=pool(C,g61,['US','China'],ivw).set_index('taxon'); x=pc.pooled; y=C['Mexico'].loc[g61,'effect_size']
    k=int((np.sign(x)==np.sign(y)).sum()); rho,pr=stats.spearmanr(x,y)
    out['balanced']=pd.DataFrame([{'Comparison':'US + China pooled (fixed effect) vs Mexico','Genera (n)':61,'Same direction (n)':k,'Same direction (%)':round(100*k/61,1),'Binomial P (vs 50%)':stats.binomtest(k,61,.5).pvalue,'Spearman rho':rho,'Spearman P':pr,
                                   'Smallest FDR of US + China pooled':pc.FDR.min()}])
    # D. fixed-effect hits and heterogeneity
    rows=[]
    for t in g61:
        y=np.array([C[c].loc[t,'effect_size'] for c in all3]);v=np.array([C[c].loc[t,'standard_error']**2 for c in all3]);w=1/v
        mu=(w*y).sum()/w.sum();Q=(w*(y-mu)**2).sum()
        rows.append({'taxon':t,'US':y[0],'China':y[1],'Mexico':y[2],'Fixed-effect pooled':mu,'Weight, Mexico':w[2]/w.sum(),'Cochran Q P':1-stats.chi2.cdf(Q,2),'I2 (%)':max(0,(Q-2)/Q)*100 if Q>0 else 0})
    h=pd.DataFrame(rows); h['Fixed-effect FDR']=fe.set_index('taxon').loc[h.taxon,'FDR'].values
    h['Same direction in all three']=[bool(abs(np.sign([r.US,r.China,r.Mexico]).sum())==3) for r in h.itertuples()]
    out['fe_hits']=h[h['Fixed-effect FDR']<.2].sort_values('Fixed-effect FDR').reset_index(drop=True)
    out['het']=dict(q_lt05=int((h['Cochran Q P']<.05).sum()),n=len(h),median_I2=float(h['I2 (%)'].median()),
                    fe_hits=len(out['fe_hits']),fe_hits_q=int((out['fe_hits']['Cochran Q P']<.05).sum()),fe_hits_same=int(out['fe_hits']['Same direction in all three'].sum()))
    return out
if __name__=='__main__':
    o=run()
    for k,v in o.items(): print('\n==',k); print(v if not isinstance(v,pd.DataFrame) else v.round(4).to_string(index=False))
