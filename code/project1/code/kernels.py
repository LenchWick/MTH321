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
