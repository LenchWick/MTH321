# %% Setup and reproducibility
from pathlib import Path
import sys, platform, time, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit, ndtr, ndtri
OUT = Path(__file__).resolve().parents[1] / 'outputs'; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10, 'axes.grid':True, 'grid.alpha':.25,
                     'figure.dpi':160, 'savefig.dpi':250, 'savefig.bbox':'tight'})
SEED=20260916
M_GBM=65536; M_SV=32768; BATCH=1024
RESULTS={'seed':SEED, 'python':sys.version.split()[0], 'numpy':np.__version__,
         'platform':platform.platform(), 'M_GBM':M_GBM,'M_SV':M_SV}
def summary(x):
    x=np.asarray(x); m=float(x.mean()); se=float(x.std(ddof=1)/np.sqrt(x.size))
    return {'mean':m,'se':se,'lo':m-1.96*se,'hi':m+1.96*se}
def slope_ci(batch_errors, hs, seed=17):
    # Resample independent path batches; retain dependence across levels.
    a=np.asarray(batch_errors); hs=np.asarray(hs)
    fit=lambda v: np.polyfit(np.log(hs),np.log(v),1)[0]
    rng=np.random.default_rng(seed)
    slopes=[fit(a[rng.integers(len(a),size=len(a))].mean(axis=0)) for _ in range(500)]
    return [float(fit(a.mean(axis=0))), *np.quantile(slopes,[.025,.975]).tolist()]
def savefig(name):
    plt.tight_layout(); plt.savefig(OUT/(name+'.png')); plt.savefig(OUT/(name+'.pdf'))
    show_figure(plt.gcf()) if 'show_figure' in globals() else None
    plt.close()
def gbm(dW,h,method,mu=.2,sigma=.3,s0=100.):
    # Rows are time, columns independent paths; no Python loop over paths.
    if method=='exact': return s0*np.exp((mu-.5*sigma*sigma)*dW.shape[0]*h+sigma*dW.sum(axis=0))
    f=1+mu*h+sigma*dW
    if method=='Milstein': f=f+.5*sigma*sigma*(dW*dW-h)
    return s0*np.prod(f,axis=0)
def sv(d1,d2,h,mu=.05,kappa=2.,theta=-.2,xi=.6, balanced=False):
    x=np.full(d1.shape[1],np.log(100.)); y=np.zeros_like(x)
    decay=np.exp(-kappa*h)
    noise=xi*np.sqrt(-np.expm1(-2*kappa*h)/(2*kappa*h))
    for a,b in zip(d1,d2):
        g=.1+.4*expit(y)
        x += (mu-.5*g*g)*h+g*a
        if balanced:
            y=theta+decay*(y-theta)+noise*np.sqrt(1+y*y)*b
        else:
            y += kappa*(theta-y)*h+xi*np.sqrt(1+y*y)*b
        assert np.all(np.isfinite(y)) and np.all(np.abs(x)<700)
    s=np.exp(x)
    assert np.all(np.isfinite(s)) and np.all(s>0) and np.all(np.isfinite(y))
    return s,y
print('Reproducibility:', RESULTS)
# %% GBM: nested paths, strong and weak errors
ng=np.array([4,8,16,32,64,128,256,512,1024]); hg=1/ng
rng=np.random.default_rng(np.random.SeedSequence([SEED,1]))
gerr={k:[] for k in ['EM','Milstein']}; gdiff={k:[] for k in gerr}
gterm={k:[] for k in ['exact','EM','Milstein']}; batches={k:[] for k in gerr}
start=time.perf_counter()
for b in range(M_GBM//BATCH):
    fine=rng.normal(size=(ng[-1],BATCH))*np.sqrt(hg[-1]); exact=gbm(fine,hg[-1],'exact')
    errors={k:[] for k in gerr}; diffs={k:[] for k in gerr}
    for n,h in zip(ng,hg):
        dw=fine.reshape(n,ng[-1]//n,BATCH).sum(axis=1)
        for method in gerr:
            s=gbm(dw,h,method); errors[method].append(abs(s-exact)); diffs[method].append(s-exact)
            if n==32: gterm[method].append(s)
    gterm['exact'].append(exact)
    for method in gerr:
        gerr[method].append(np.array(errors[method])); gdiff[method].append(np.array(diffs[method]))
        batches[method].append(np.array(errors[method]).mean(axis=1))
G={}; analytic_mean=100*np.exp(.2); analytic_var=10000*np.exp(.4)*np.expm1(.09)
for method in gerr:
    gerr[method]=np.concatenate(gerr[method],axis=1); gdiff[method]=np.concatenate(gdiff[method],axis=1)
    G[method]={'strong':[summary(v) for v in gerr[method]],'weak_signed':[summary(v) for v in gdiff[method]],
               'slope':slope_ci(np.array(batches[method])[:,2:],hg[2:])}
for k in gterm: gterm[k]=np.concatenate(gterm[k])
bias_exact=100*np.exp(ng*np.log1p(.2/ng))-analytic_mean
G.update({'n':ng.tolist(),'analytic_mean':analytic_mean,'analytic_variance':analytic_var,
          'exact_mean_check':summary(gterm['exact']), 'exact_sample_variance':float(gterm['exact'].var(ddof=1)),
          'analytic_bias':bias_exact.tolist(),'elapsed_seconds':time.perf_counter()-start})
G['variance_check']=summary((gterm['exact']-analytic_mean)**2)
for method in ['EM','Milstein']:
    k=7 if method=='EM' else 9
    wb=-gdiff[method].reshape(9,M_GBM//BATCH,BATCH).mean(axis=2).T[:,:k]
    G[method]['weak_slope']=slope_ci(wb,hg[:k])
G['analytic_weak_slope']=float(np.polyfit(np.log(hg),np.log(abs(bias_exact)),1)[0])
for method in ['EM','Milstein']:
    print(method,'strong slope, bootstrap 95% CI:',G[method]['slope'])
    print('Weak slope and bootstrap CI:',G[method]['weak_slope'])
    print(' n     strong L1 (SE)        signed weak bias (SE)      analytic bias')
    for j,n in enumerate(ng):
        a=G[method]['strong'][j]; d=G[method]['weak_signed'][j]
        print(f'{n:4}  {a["mean"]:.6f} ({a["se"]:.6f})  {d["mean"]:+.6f} ({d["se"]:.6f})  {bias_exact[j]:+.6f}')
print('Exact-law moments:', {k:G[k] for k in ['analytic_mean','analytic_variance','exact_mean_check','exact_sample_variance']})
fig,ax=plt.subplots(1,2,figsize=(10,4))
for method in ['EM','Milstein']:
    ss=G[method]['strong']; ww=G[method]['weak_signed']
    ax[0].errorbar(hg,[v['mean'] for v in ss],yerr=[1.96*v['se'] for v in ss],marker='o',label=method)
    ax[1].errorbar(hg,-np.array([v['mean'] for v in ww]),yerr=[1.96*v['se'] for v in ww],marker='o',label=method)
ax[1].plot(hg,-bias_exact,'k--',label='analytic bias magnitude')
for a in ax: a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('step h'); a.legend()
ax[1].set_yscale('symlog',linthresh=.003); ax[1].axhline(0,c='gray',ls=':')
ax[0].set_ylabel('mean absolute terminal error'); ax[1].set_ylabel('negative signed mean bias'); savefig('gbm_convergence')
q=np.array([.01,.05,.5,.95,.99]); qe=100*np.exp(.155+.3*ndtri(q))
bins=np.linspace(0,400,81)
G['histogram']={'centers':(.5*(bins[1:]+bins[:-1])).tolist(),**{k:(np.histogram(v,bins)[0]/(len(v)*5)).tolist() for k,v in gterm.items()}}
G['quantiles']={'p':q.tolist(),'analytic':qe.tolist(),**{k:np.quantile(v,q).tolist() for k,v in gterm.items()}}
fig,ax=plt.subplots(1,2,figsize=(10,4)); xx=np.linspace(25,300,500)
density=np.exp(-.5*((np.log(xx/100)-.155)/.3)**2)/(xx*.3*np.sqrt(2*np.pi))
for k in ['exact','EM','Milstein']: ax[0].hist(gterm[k],bins=80,density=True,histtype='step',label=k)
ax[0].plot(xx,density,'k--',label='lognormal density'); ax[0].set_xlim(25,300); ax[0].legend(); ax[0].set_xlabel('terminal price'); ax[0].set_ylabel('density')
for k in ['exact','EM','Milstein']: ax[1].plot(qe,G['quantiles'][k],'o-',label=k)
ax[1].plot(qe,qe,'k--'); ax[1].set_xlabel('analytic quantile'); ax[1].set_ylabel('sample quantile'); ax[1].legend(); savefig('gbm_distribution')
RESULTS['gbm']=G
# %% Non-affine model: all reference grids use the same Brownian path
ns=np.array([8,16,32,64,128,256,1024,2048,4096]); hs=1/ns
rng=np.random.default_rng(np.random.SeedSequence([SEED,2])); st=[]; yt=[]; bt=[]; moment=np.zeros(5); ninc=0
start=time.perf_counter()
for b in range(M_SV//BATCH):
    z1=rng.normal(size=(4096,BATCH)); z2=rng.normal(size=z1.shape)
    d1=z1/64; d2=(-.7*z1+np.sqrt(.51)*z2)/64
    moment += [d1.sum(),d2.sum(),(d1*d1).sum(),(d2*d2).sum(),(d1*d2).sum()]; ninc+=d1.size
    ss=[]; yy=[]; bb=[]
    for n,h in zip(ns,hs):
        a=d1.reshape(n,4096//n,BATCH).sum(axis=1); c=d2.reshape(n,4096//n,BATCH).sum(axis=1)
        s,y=sv(a,c,h); ss.append(s); yy.append(y)
        if n<=256: bb.append(sv(a,c,h,balanced=True)[0])
    st.append(ss); yt.append(yy); bt.append(bb)
st=np.concatenate(st,axis=1); yt=np.concatenate(yt,axis=1); bt=np.concatenate(bt,axis=1)
V={'n':ns.tolist(),'elapsed_seconds':time.perf_counter()-start,'reference_comparisons':{}}
means=moment[:2]/ninc; variances=moment[2:4]/ninc-means**2
V['brownian']={'count':ninc,'means':means.tolist(),'variances':variances.tolist(),
               'correlation':float((moment[4]/ninc-means.prod())/np.sqrt(variances.prod()))}
for ref in [1024,2048,4096]:
    j=int(np.where(ns==ref)[0][0]); err=abs(st[:6]-st[j]); payoff=np.maximum(st[:6]-100,0)-np.maximum(st[j]-100,0)
    batch=err.reshape(6,M_SV//BATCH,BATCH).mean(axis=2).T
    V['reference_comparisons'][str(ref)]={'strong':[summary(a) for a in err], 'payoff_bias':[summary(a) for a in payoff], 'slope':slope_ci(batch[:,:5],hs[:5])}
V['refinement']={}
for a,b in [(1024,2048),(2048,4096),(1024,4096)]:
    i=int(np.where(ns==a)[0][0]); j=int(np.where(ns==b)[0][0])
    V['refinement'][f'{a}-{b}']={'strong':summary(abs(st[i]-st[j])),
      'payoff':summary(np.maximum(st[i]-100,0)-np.maximum(st[j]-100,0))}
V['terminal']={'mean':summary(st[-1]),'payoff':summary(np.maximum(st[-1]-100,0)),
 'min_price':float(st.min()),'max_price':float(st.max()),'min_y':float(yt.min()),'max_y':float(yt.max()),
 'analytic_mean':100*np.exp(.05)}
# An exact discrete martingale identity validates E[S_n] for log EM.
print('Brownian validation:',V['brownian']); print('Reference terminal diagnostics:',V['terminal'])
for ref,r in V['reference_comparisons'].items():
 print('Reference',ref,'empirical slope and CI',r['slope'])
 for n,a,p in zip(ns[:6],r['strong'],r['payoff_bias']): print(f'n={n:3} L1={a["mean"]:.6f} +/- {1.96*a["se"]:.6f}; payoff bias={p["mean"]:+.6f} +/- {1.96*p["se"]:.6f}')
print('Reference refinement:',V['refinement'])
fig,ax=plt.subplots(1,2,figsize=(10,4))
for ref,r in V['reference_comparisons'].items():
 ax[0].errorbar(hs[:6],[v['mean'] for v in r['strong']],yerr=[1.96*v['se'] for v in r['strong']],marker='o',label='ref 1/'+ref)
 ax[1].errorbar(hs[:6],[v['mean'] for v in r['payoff_bias']],yerr=[1.96*v['se'] for v in r['payoff_bias']],marker='o',label='ref 1/'+ref)
ax[0].set_yscale('log'); ax[0].set_ylabel('mean absolute terminal difference'); ax[1].set_ylabel('signed call-payoff difference'); ax[1].axhline(0,c='k',ls='--')
for a in ax: a.set_xscale('log'); a.set_xlabel('step h'); a.legend()
savefig('nonaffine_convergence'); RESULTS['nonaffine']=V
# %% Cost at a common GBM pathwise accuracy target
# Time only integration on pre-generated increments, 7 repeats after warm-up.
target=.5; cost={}; rng=np.random.default_rng(np.random.SeedSequence([SEED,3]))
for method in ['EM','Milstein']:
    acceptable=[j for j,a in enumerate(G[method]['strong']) if a['hi']<=target]
    assert acceptable, 'Extend grid for matched target'
    j=acceptable[0]; n=int(ng[j]); dw=rng.normal(size=(n,8192))/np.sqrt(n)
    gbm(dw,1/n,method); ts=[]
    for _ in range(7):
        tic=time.perf_counter(); gbm(dw,1/n,method); ts.append(time.perf_counter()-tic)
    cost[method]={'n':n,'error':G[method]['strong'][j], 'seconds_median':float(np.median(ts)), 'seconds_range':[min(ts),max(ts)], 'path_updates':n*8192}
dw=rng.normal(size=(1,8192)); gbm(dw,1,'exact'); ts=[]
for _ in range(7):
 tic=time.perf_counter(); gbm(dw,1,'exact'); ts.append(time.perf_counter()-tic)
cost['exact']={'n':1,'seconds_median':float(np.median(ts)),'path_updates':8192}
RESULTS['cost']={'target':target,'paths':8192,'methods':cost}; print('Matched-accuracy cost:',RESULTS['cost'])
fig,ax=plt.subplots(figsize=(6,4)); labels=list(cost); ax.bar(labels,[cost[k]['seconds_median']*1000 for k in labels]); ax.set_ylabel('median integration time (ms), 8192 paths'); ax.set_title('GBM terminal L1 target <= 0.5'); savefig('cost')
# %% Positivity: analytic one-step probabilities and independent stress samples
rng=np.random.default_rng(np.random.SeedSequence([SEED,4])); positivity=[]
for label,mu,sigma in [('baseline',.2,.3),('stress',.05,1.5)]:
 for h in [1.,.5,.25,.125]:
    z=rng.normal(size=200000); em=1+mu*h+sigma*np.sqrt(h)*z
    mi=em+.5*sigma*sigma*h*(z*z-1)
    minimum=.5+(mu-.5*sigma*sigma)*h
    p_em=float(ndtr(-(1+mu*h)/(sigma*np.sqrt(h))))
    if minimum<0:
      root=np.sqrt(-2*minimum/(sigma*sigma*h)); center=-1/(sigma*np.sqrt(h))
      p_mi=float(ndtr(center+root)-ndtr(center-root))
    else:p_mi=0.
    row={'case':label,'h':h,'EM_probability':p_em,'EM_count':int((em<0).sum()),
         'Milstein_probability':p_mi,'Milstein_count':int((mi<0).sum()),'samples':len(z)}
    positivity.append(row); print(row)
RESULTS['positivity']=positivity
# %% Derived extension: moment-balanced factor update and payoff control variate
E={'strong':[summary(abs(v-st[-1])) for v in bt],
   'payoff_bias':[summary(np.maximum(v-100,0)-np.maximum(st[-1]-100,0)) for v in bt]}
batch=abs(bt-st[-1]).reshape(6,M_SV//BATCH,BATCH).mean(axis=2).T
E['slope']=slope_ci(batch[:,:5],hs[:5]); print('Moment-balanced method:',E)
# Reference payoff CV: beta uses a disjoint pilot, evaluation uses remaining paths.
pilot=4096; stock=st[-1]; call=np.maximum(stock-100,0)
beta=float(np.cov(stock[:pilot],call[:pilot],ddof=1)[0,1]/stock[:pilot].var(ddof=1))
raw=call[pilot:]; controlled=raw-beta*(stock[pilot:]-100*np.exp(.05))
E['control_variate']={'pilot':pilot,'evaluation':len(raw),'beta':beta,'raw':summary(raw),
 'controlled':summary(controlled),'variance_reduction':float(raw.var(ddof=1)/controlled.var(ddof=1))}
# The model and log-EM both have exactly known first price moment.
print('Independent-pilot payoff control variate:',E['control_variate'])
# Stability / stationary second moment stress experiment (T=40).
rng=np.random.default_rng(np.random.SeedSequence([SEED,5])); stationary=[]
vstar=.36*(1+.04)/(4-.36)
for h in [1.,.5,.25,.125]:
 n=round(40/h); dw=rng.normal(size=(n,32768))*np.sqrt(h)
 for balanced in [False,True]:
  y=np.full(32768,-.2); a=np.exp(-2*h); c=.6*np.sqrt(-np.expm1(-4*h)/4)
  for w in dw:
   y=(-.2+a*(y+.2)+c*np.sqrt(1+y*y)*w/np.sqrt(h)) if balanced else y+2*(-.2-y)*h+.6*np.sqrt(1+y*y)*w
  second=summary((y+.2)**2)
  A=a*a+.36*(-np.expm1(-4*h))/4 if balanced else (1-2*h)**2+.36*h
  theoretical=vstar if balanced else (.36*(1+.04)/(4-.36-4*h) if A<1 else None)
  stationary.append({'h':h,'method':'balanced' if balanced else 'EM','second_about_theta':second,
   'A':float(A),'stationary_variance':theoretical,'finite_T_expected_second':float(.36*(1+.04)*(1-a*a)/4*(1-A**n)/(1-A)) if balanced else float(.36*(1+.04)*h*(1-A**n)/(1-A))})
E['stationary']=stationary; print('Factor stability stress:',stationary)
# Match a common empirical reference-error ceiling for the two SV schemes.
cs={}; rng=np.random.default_rng(np.random.SeedSequence([SEED,6]))
for name,errs,balanced in [('EM',V['reference_comparisons']['4096']['strong'],False),('balanced',E['strong'],True)]:
 j=next(j for j,e in enumerate(errs) if e['hi']<.5); n=int(ns[j]); h=1/n
 w1=rng.normal(size=(n,8192))*np.sqrt(h); w2=-.7*w1+np.sqrt(.51*h)*rng.normal(size=w1.shape)
 sv(w1,w2,h,balanced=balanced); ts=[]
 for _ in range(7):
  tic=time.perf_counter(); sv(w1,w2,h,balanced=balanced); ts.append(time.perf_counter()-tic)
 cs[name]={'n':n,'error':errs[j],'seconds_median':float(np.median(ts)),'path_updates':n*8192}
E['cost']={'target':.5,'paths':8192,'methods':cs}; print('SV matched reference-accuracy cost:',E['cost'])
fig,ax=plt.subplots(1,2,figsize=(10,4))
for name,vals in [('EM',V['reference_comparisons']['4096']['strong']),('balanced',E['strong'])]:
 ax[0].errorbar(hs[:6],[v['mean'] for v in vals],yerr=[1.96*v['se'] for v in vals],marker='o',label=name)
ax[0].set_xscale('log');ax[0].set_yscale('log');ax[0].set_xlabel('step h');ax[0].set_ylabel('mean absolute terminal difference');ax[0].legend()
for name in ['EM','balanced']:
 rows=[r for r in stationary if r['method']==name and r['h']<1]
 ax[1].plot([r['h'] for r in rows],[r['second_about_theta']['mean'] for r in rows],'o-',label=name)
ax[1].axhline(vstar,c='k',ls='--',label='exact stationary variance');ax[1].set_xlabel('step h');ax[1].set_ylabel('E[(Y(40)-theta)^2]');ax[1].legend();savefig('extension')
RESULTS['extension']=E
# %% Transient animation: prices, volatility and the factor on shared paths
from matplotlib.animation import FuncAnimation, PillowWriter
rng=np.random.default_rng(np.random.SeedSequence([SEED,7])); n=256; h=1/n
w1=rng.normal(size=(n,6))*np.sqrt(h);w2=-.7*w1+np.sqrt(.51*h)*rng.normal(size=w1.shape)
xx=np.full((n+1,6),np.log(100.)); yy=np.zeros((n+1,6))
for i in range(n):
 g=.1+.4*expit(yy[i]);xx[i+1]=xx[i]+(.05-.5*g*g)*h+g*w1[i]
 yy[i+1]=yy[i]+2*(-.2-yy[i])*h+.6*np.sqrt(1+yy[i]**2)*w2[i]
prices=np.exp(xx);vol=.1+.4*expit(yy); tt=np.arange(n+1)*h
fig,axes=plt.subplots(2,1,figsize=(7,5),sharex=True);lines=[]
for ax,data,label in zip(axes,[prices,vol],['price S','volatility g(Y)']):
 for j in range(6): lines.append((ax.plot([],[],lw=1)[0],data[:,j]))
 ax.set_xlim(0,1);ax.set_ylim(data.min()*.95,data.max()*1.05);ax.set_ylabel(label)
axes[-1].set_xlabel('time (years)')
def animate(k):
 for line,data in lines: line.set_data(tt[:k+1],data[:k+1])
 return [item[0] for item in lines]
ani=FuncAnimation(fig,animate,frames=range(0,n+1,4),interval=70,blit=True)
ani.save(OUT/'paths.gif',writer=PillowWriter(fps=14));plt.close(fig)
if 'show_animation' in globals(): show_animation(OUT/'paths.gif')
print('Animation generated: outputs/paths.gif (six simulated paths, fixed seed).')
# %% Independent checks and export
checks={}
rng=np.random.default_rng(44); d=rng.normal(size=(16,100))/4
checks['exact_transition_composition']=bool(np.allclose(100*np.prod(np.exp(.155/16+.3*d),axis=0),gbm(d,1/16,'exact'),rtol=1e-13))
checks['zero_noise_EM_equals_Milstein']=bool(np.array_equal(gbm(d,1/16,'EM',sigma=0),gbm(d,1/16,'Milstein',sigma=0)))
x=100.; w=.12; h=.1
checks['scalar_Milstein_formula']=bool(np.isclose(gbm(np.array([[w]]),h,'Milstein')[0],x*(1+.2*h+.3*w+.5*.09*(w*w-h))))
# xi=0, theta=Y0=0 freezes Y; log EM must equal exact GBM with sigma=g(0)=.3.
s,y=sv(d,d,1/16,mu=.2,theta=0,xi=0)
checks['nonaffine_frozen_factor_reduces_to_GBM']=bool(np.allclose(s,gbm(d,1/16,'exact'),rtol=1e-13))
checks['Brownian_variance']=bool(np.all(abs(np.array(V['brownian']['variances'])*4096-1)<.002))
checks['Brownian_correlation']=bool(abs(V['brownian']['correlation']+.7)<.002)
checks['GBM_moment_within_4SE']=bool(abs(G['exact_mean_check']['mean']-analytic_mean)<4*G['exact_mean_check']['se'])
checks['SV_martingale_within_4SE']=bool(abs(V['terminal']['mean']['mean']-100*np.exp(.05))<4*V['terminal']['mean']['se'])
checks['EM_order']=.45<G['EM']['slope'][0]<.6
checks['Milstein_order']=.9<G['Milstein']['slope'][0]<1.1
checks['balanced_frozen_factor']=bool(np.allclose(sv(d,d,1/16,mu=.2,theta=0,xi=0,balanced=True)[0],gbm(d,1/16,'exact'),rtol=1e-13))
checks['variance_identity']=bool(np.isclose((np.exp(-.4)**2+.36*(-np.expm1(-.8))/4)*vstar+.36*(-np.expm1(-.8))/4*(1+.04),vstar))
assert all(checks.values()), checks
RESULTS['checks']=checks
(OUT/'results.json').write_text(json.dumps(RESULTS,indent=2))
print('Validation checks:',json.dumps(checks,indent=2))
print('All outputs written to',OUT.resolve())

# Submission visualization checks: bootstrap whole independent batches.
from submission_figures import supplementary_figures
supplementary_figures(OUT, G, V, batches, gerr, st, hg, hs, M_SV, BATCH)

# Independent whole-experiment seeds required by the visualization guide.
from repeat_seeds import independent_seed_check
independent_seed_check(OUT)
