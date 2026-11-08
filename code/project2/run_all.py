"""Variable diffusion project: run all verification, comparisons and plots.

Run: python run_all.py. All files are written beside this script in outputs/.
No network, input data, FEM package, or hidden helper script is needed.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import csv
import json
import platform
import time
import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import splu, cg, eigsh, LinearOperator
import sympy as sy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits, threadpool_info

OUT = Path(os.environ.get('MTH321_OUTPUT_DIR', Path(__file__).resolve().parent / 'outputs')).resolve()
OUT.mkdir(exist_ok=True)
LOG = []
TIMINGS = []


def say(s):
    print(s, flush=True)
    LOG.append(str(s))


def save(name, rows):
    with (OUT / (name + '.csv')).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def exact(x, y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)


def forcing(x, y):
    return 2*np.pi**2*(1+x)*exact(x, y)-np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)


def layered(x, y, contrast=10):
    flux = 2*contrast/(contrast+1)
    return np.where(x <= .5, flux*x, flux*.5+flux*(x-.5)/contrast)+0*y


def assemble(N, coefficient, source, boundary, average='arithmetic'):
    """Shared edge coefficients; interior index k=(i-1)*(N-1)+(j-1)."""
    if not isinstance(N, (int, np.integer)) or N < 2:
        raise ValueError('N must be an integer at least 2')
    h = 1/N
    t = np.linspace(0, 1, N+1)
    X, Y = np.meshgrid(t, t, indexing='ij')
    a = np.broadcast_to(np.asarray(coefficient(X, Y), dtype=float), X.shape)
    if not np.all(np.isfinite(a)) or np.min(a) <= 0:
        raise ValueError('Conductivity must be finite and positive')
    if average == 'arithmetic':
        mean = lambda l, r: (l+r)/2
    elif average == 'harmonic':
        mean = lambda l, r: 2*l*r/(l+r)
    else:
        raise ValueError(average)
    ax, ay = mean(a[:-1], a[1:]), mean(a[:, :-1], a[:, 1:])
    E, W = ax[1:, 1:-1], ax[:-1, 1:-1]
    T, B = ay[1:-1, 1:], ay[1:-1, :-1]
    m = N-1
    ids = np.arange(m*m).reshape(m, m)
    rows = [ids.ravel()]
    cols = [ids.ravel()]
    vals = [((E+W+T+B)/h**2).ravel()]
    # Each undirected edge is inserted in both rows with the same weight.
    for p, q, weight in [(ids[:-1], ids[1:], ax[1:-1, 1:-1]),
                          (ids[:, :-1], ids[:, 1:], ay[1:-1, 1:-1])]:
        rows.extend([p.ravel(), q.ravel()])
        cols.extend([q.ravel(), p.ravel()])
        vals.extend([-weight.ravel()/h**2, -weight.ravel()/h**2])
    A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows),
                     np.concatenate(cols))), shape=(m*m, m*m)).tocsc()
    b = np.broadcast_to(np.asarray(source(X[1:-1, 1:-1], Y[1:-1, 1:-1]), dtype=float), (m, m)).copy()
    g = np.broadcast_to(np.asarray(boundary(X, Y), dtype=float), X.shape).copy()
    if not np.all(np.isfinite(b)) or not np.all(np.isfinite(g)):
        raise ValueError('Source and boundary must be finite')
    b[0] += W[0]*g[0, 1:-1]/h**2
    b[-1] += E[-1]*g[-1, 1:-1]/h**2
    b[:, 0] += B[:, 0]*g[1:-1, 0]/h**2
    b[:, -1] += T[:, -1]*g[1:-1, -1]/h**2
    return A, b.ravel(), (X, Y, g, ax, ay)


def metrics(A, b, u):
    r = b-A@u
    rel = np.linalg.norm(r)/max(np.linalg.norm(b), 1.)
    anorm = np.max(np.asarray(abs(A).sum(axis=1)))
    denominator = anorm*np.max(abs(u))+np.max(abs(b))
    back = np.max(abs(r))/denominator if denominator else 0.
    return float(rel), float(back)


def direct(A, b, case='unspecified'):
    # One untimed warm-up; five independent factor-and-solve repetitions.
    splu(A, permc_spec='COLAMD').solve(b)
    factors, solves, totals = [], [], []
    for repeat in range(5):
        t = time.perf_counter()
        lu = splu(A, permc_spec='COLAMD')
        factors.append(time.perf_counter()-t)
        t = time.perf_counter()
        u = lu.solve(b)
        solves.append(time.perf_counter()-t)
        totals.append(factors[-1]+solves[-1])
        TIMINGS.append(dict(case=case,method='LU',dof=len(b),tolerance=0,repeat=repeat+1,
            setup_s=0.,factor_s=factors[-1],solve_s=solves[-1],total_s=totals[-1]))
    direct.stats = dict(median=float(np.median(totals)),minimum=min(totals),maximum=max(totals))
    return u, float(np.median(factors)), float(np.median(solves)), lu.L.nnz+lu.U.nnz


def pcg(A, b, tol, case='unspecified'):
    totals = []
    for repeat in range(6):  # repeat zero is an unrecorded warm-up
        t = time.perf_counter()
        diagonal = A.diagonal()
        M = LinearOperator(A.shape, matvec=lambda v: v/diagonal, dtype=float)
        setup = time.perf_counter()-t
        count = [0]
        def callback(u):
            count[0] += 1
        t = time.perf_counter()
        u, info = cg(A, b, M=M, rtol=tol, atol=0., maxiter=20000, callback=callback)
        elapsed = time.perf_counter()-t
        if info != 0:
            raise RuntimeError(f'CG did not converge: {info}')
        if repeat:
            totals.append(setup+elapsed)
            TIMINGS.append(dict(case=case,method='PCG',dof=len(b),tolerance=tol,repeat=repeat,
                setup_s=setup,factor_s=0.,solve_s=elapsed,total_s=setup+elapsed))
    pcg.stats = dict(minimum=min(totals),maximum=max(totals))
    return u, count[0], float(np.median(totals))


def error(u, ref, h):
    return float(h*np.linalg.norm(u-ref)), float(np.max(abs(u-ref)))


def graph_neumann(N):
    """Cell-centred conservative operator with zero boundary flux."""
    h = 1/N
    t = (np.arange(N)+.5)*h
    X, Y = np.meshgrid(t, t, indexing='ij')
    ids = np.arange(N*N).reshape(N, N)
    p = np.concatenate([ids[:-1].ravel(), ids[:, :-1].ravel()])
    q = np.concatenate([ids[1:].ravel(), ids[:, 1:].ravel()])
    ax = np.broadcast_to(1+np.arange(1, N)[:, None]*h, (N-1, N))
    ay = np.broadcast_to(1+X[:, :-1], (N, N-1))
    w = np.concatenate([ax.ravel(), ay.ravel()])/h**2
    A = sp.coo_matrix((np.concatenate([w, w, -w, -w]),
                     (np.concatenate([p, q, p, q]), np.concatenate([p, q, q, p]))),
                     shape=(N*N, N*N)).tocsc()
    ref = np.cos(np.pi*X)*np.cos(np.pi*Y)
    b = (2*np.pi**2*(1+X)*ref+np.pi*np.sin(np.pi*X)*np.cos(np.pi*Y)).ravel()
    return A, b, ref.ravel()


def neumann_solve(A, b):
    # Check compatibility BEFORE any projection; reject substantial mismatch.
    defect = abs(np.mean(b))
    if defect > 1e-12*max(1., np.max(abs(b))):
        raise ValueError('Incompatible Neumann forcing; no zero-flux solution.')
    bp = b-np.mean(b)  # remove only quadrature/round-off imbalance
    z = sp.csc_matrix(np.ones((len(b), 1))/np.sqrt(len(b)))
    K = sp.bmat([[A, z], [z.T, None]], format='csc')
    sol = splu(K).solve(np.r_[bp, 0.])
    return sol[:-1], float(sol[-1]), float(defect)


def figure(name):
    plt.tight_layout()
    plt.savefig(OUT / (name+'.png'), dpi=180)
    plt.savefig(OUT / (name+'.pdf'))
    plt.close()


def main():
    say('Predictions: smooth order 2; harmonic interface nodal exactness; '
        'arithmetic interface resistance error; Neumann constant nullspace.')
    x, y = sy.symbols('x y')
    us = sy.sin(sy.pi*x)*sy.sin(sy.pi*y)
    fs = -sy.diff((1+x)*sy.diff(us, x), x)-sy.diff((1+x)*sy.diff(us, y), y)
    assert sy.simplify(fs-(2*sy.pi**2*(1+x)*us-sy.pi*sy.cos(sy.pi*x)*sy.sin(sy.pi*y))) == 0
    un = sy.cos(sy.pi*x)*sy.cos(sy.pi*y)
    fn = -sy.diff((1+x)*sy.diff(un, x), x)-sy.diff((1+x)*sy.diff(un, y), y)
    assert sy.simplify(fn-(2*sy.pi**2*(1+x)*un+sy.pi*sy.sin(sy.pi*x)*sy.cos(sy.pi*y))) == 0
    say('Independent symbolic forcing checks: PASS (Dirichlet and Neumann).')
    smooth, comparisons = [], []
    for N in [16, 32, 64, 128]:
        t0 = time.perf_counter()
        A, b, data = assemble(N, lambda x,y: 1+x, forcing, exact)
        assembly = time.perf_counter()-t0
        u, factor, solve, fill = direct(A, b, f'smooth_N{N}')
        direct_stats = direct.stats.copy()
        X, Y, g, ax, ay = data
        e2, ei = error(u, exact(X[1:-1, 1:-1], Y[1:-1, 1:-1]).ravel(), 1/N)
        rel, back = metrics(A, b, u)
        asym = sp.linalg.norm(A-A.T)/sp.linalg.norm(A)
        row = dict(N=N, h=1/N, dof=A.shape[0], nnz=A.nnz, L2=e2, Linf=ei,
                   residual=rel, backward=back, asymmetry=asym, assembly_s=assembly,
                   factor_s=factor, solve_s=solve, LU_nnz=fill)
        smooth.append(row)
        v, it, elapsed = pcg(A, b, 1e-12, f'smooth_N{N}')
        gap = error(v, u, 1/N)[0]
        assert gap < e2*1e-5
        comparisons.append(dict(N=N, tolerance=1e-12, iterations=it, pcg_s=elapsed,
                                direct_s=direct_stats['median'], gap_L2=gap,
                                direct_min_s=direct_stats['minimum'],direct_max_s=direct_stats['maximum'],
                                pcg_min_s=pcg.stats['minimum'],pcg_max_s=pcg.stats['maximum'],
                                residual=metrics(A,b,v)[0]))
        if N == 64:
            U = g.copy(); U[1:-1, 1:-1] = u.reshape(N-1,N-1)
            np.savez(OUT/'smooth_field.npz',X=X,Y=Y,U=U,error=U-exact(X,Y))
            plt.figure(figsize=(5,5)); plt.spy(A,markersize=.2)
            plt.title('Conservative five-point matrix, N=64'); figure('matrix_spy')
            plt.figure(figsize=(9,4))
            plt.subplot(121); plt.pcolormesh(X,Y,U,shading='auto'); plt.colorbar(); plt.title('Smooth solution')
            plt.subplot(122); plt.pcolormesh(X,Y,U-exact(X,Y),shading='auto'); plt.colorbar(); plt.title('Nodal error')
            figure('smooth_field')
    save('smooth', smooth); save('solver_comparison', comparisons)
    for key in ['L2', 'Linf']:
        p = np.polyfit(np.log([r['h'] for r in smooth[-3:]]), np.log([r[key] for r in smooth[-3:]]),1)[0]
        say(f'Smooth fitted {key} order (last three meshes): {p:.8f}')
        assert 1.95 < p < 2.05
    plt.figure()
    for key in ['L2','Linf']:
        plt.loglog([r['h'] for r in smooth], [r[key] for r in smooth], 'o-', label=key)
    plt.xlabel('h'); plt.ylabel('error'); plt.legend(); plt.grid(True); figure('smooth_convergence')
    # Structural and nonzero-boundary audit independent of the sine case.
    A, b, _ = assemble(8, lambda x,y: 1+x, lambda x,y: -1+0*x*y, lambda x,y: x+2*y)
    X,Y = np.meshgrid(np.arange(1,8)/8,np.arange(1,8)/8,indexing='ij')
    affine_error = np.max(abs(splu(A).solve(b)-(X+2*Y).ravel()))
    assert A[0,0] == 288 and A[0,1] == -72 and A[0,7] == -76 and b[0] == 25
    off = A-sp.diags(A.diagonal())
    assert affine_error < 1e-12 and np.max(off.data) <= 0
    eigen_min = float(eigsh(A,k=1,which='SA',return_eigenvectors=False)[0])
    say(f'Affine boundary audit error {affine_error:.3e}; smallest eigenvalue {eigen_min:.6g}')
    # Tolerance sweep separates linear-solve and discretisation errors.
    A,b,data = assemble(64,lambda x,y: 1+x,forcing,exact)
    ud = splu(A).solve(b)
    ref = exact(data[0][1:-1,1:-1],data[1][1:-1,1:-1]).ravel()
    sweep=[]
    for tol in [1e-4,1e-6,1e-8,1e-10,1e-12]:
        u,it,t = pcg(A,b,tol, f'tolerance_N64_{tol}')
        sweep.append(dict(tolerance=tol, iterations=it,residual=metrics(A,b,u)[0],
                          gap_L2=error(u,ud,1/64)[0],L2=error(u,ref,1/64)[0],
                          QoI=u.reshape(63,63)[31,31]))
    save('tolerance_sweep',sweep)
    interfaces=[]
    for N in [31,63,127,255]:
        for avg in ['arithmetic','harmonic']:
            A,b,data=assemble(N,lambda x,y: np.where(x<.5,1.,10.),lambda x,y:0*x*y,layered,avg)
            u,fac,sol,fill=direct(A,b, f'interface_N{N}_{avg}')
            X,Y,U,ax,ay=data; U[1:-1,1:-1]=u.reshape(N-1,N-1)
            e2,ei=error(u,layered(X[1:-1,1:-1],Y[1:-1,1:-1]).ravel(),1/N)
            mid=N//2
            q=-ax[:,mid]*(U[1:,mid]-U[:-1,mid])*N
            # Adjacent faces on the central horizontal cut, not a false 1D assumption.
            jump=abs(q[mid-1]-q[mid+1])
            # Independent constitutive reconstruction on the cut edge.
            qphysical=-(20/11)*(U[mid+1,mid]-U[mid,mid])*N
            physical_defect=max(abs(qphysical-q[mid-1]),abs(qphysical-q[mid+1]))
            interfaces.append(dict(N=N,average=avg,L2=e2,Linf=ei,
                flux_max_error=np.max(abs(q+20/11)),flux_jump=jump,
                resistance_flux=float(qphysical),physical_flux_defect=float(physical_defect),
                minimum=np.min(U),maximum=np.max(U),residual=metrics(A,b,u)[0],
                total_s=fac+sol))
            assert np.min(U)>-1e-10 and np.max(U)<1+1e-10
            if avg=='harmonic': assert ei<1e-9
            if N==63:
                save('profile_'+avg,[dict(x=float(X[i,mid]),u=float(U[i,mid])) for i in range(N+1)])
                save('flux_'+avg,[dict(x=(i+.5)/N,q=float(q[i])) for i in range(N)])
                if avg=='arithmetic':
                    plt.figure(50,figsize=(9,4))
                else:
                    plt.figure(50)
                plt.subplot(121)
                plt.plot(X[:,mid],U[:,mid],label=avg)
                plt.subplot(122); plt.plot((np.arange(N)+.5)/N,q,label=avg)
    plt.figure(50);plt.subplot(121);plt.xlabel('x');plt.ylabel('u(x,y_mid)');plt.legend()
    plt.subplot(122);plt.axhline(-20/11,color='k',ls=':',label='exact');plt.xlabel('x face');plt.ylabel('q_x');plt.legend()
    figure('interface_profiles'); save('interface',interfaces)
    contrasts=[]
    for N in [31,63,127]:
        for c in [1,10,100,1000]:
            A,b,_=assemble(N,lambda x,y:np.where(x<.5,1.,float(c)),lambda x,y:1+0*x*y,lambda x,y:0*x*y,'harmonic')
            u,it,t=pcg(A,b,1e-10, f'contrast_N{N}_c{c}')
            contrasts.append(dict(N=N,contrast=c,iterations=it,pcg_s=t,residual=metrics(A,b,u)[0]))
    save('contrast',contrasts)
    plt.figure()
    for c in [1,10,100,1000]:
        rr=[r for r in contrasts if r['contrast']==c]
        plt.plot([r['N'] for r in rr],[r['iterations'] for r in rr],'o-',label=f'contrast {c}')
    plt.xlabel('N');plt.ylabel('Jacobi-PCG iterations');plt.legend();figure('contrast_iterations')
    neumann=[]
    for N in [16,32,64,128]:
        t=time.perf_counter(); A,b,ref=graph_neumann(N); u,lag,defect=neumann_solve(A,b)
        elapsed=time.perf_counter()-t
        e2,ei=error(u,ref,1/N)
        null=float(np.max(abs(A@np.ones(N*N))))
        neumann.append(dict(N=N,h=1/N,L2=e2,Linf=ei,residual=metrics(A,b,u)[0],
            mean_u=np.mean(u),compatibility_defect=defect,null_defect=null,multiplier=lag,total_s=elapsed))
        assert null<1e-8 and abs(np.mean(u))<1e-10
    save('neumann',neumann)
    A,b,ref=graph_neumann(16)
    try:
        neumann_solve(A,b+1)
    except ValueError:
        say('Incompatible Neumann source b+1: correctly REJECTED.')
    else:
        raise AssertionError('Incompatible source was silently accepted')
    p=float(np.polyfit(np.log([r['h'] for r in neumann[-3:]]),np.log([r['L2'] for r in neumann[-3:]]),1)[0])
    assert 1.9<p<2.1
    say(f'Neumann fitted L2 order: {p:.8f}')
    plt.figure();plt.loglog([r['h'] for r in neumann],[r['L2'] for r in neumann],'o-',label='Neumann L2')
    plt.xlabel('h');plt.ylabel('error');plt.legend();figure('neumann_convergence')
    plt.figure()
    plt.loglog([r['factor_s']+r['solve_s'] for r in smooth],[r['L2'] for r in smooth],'o-',label='sparse LU')
    plt.loglog([r['pcg_s'] for r in comparisons],[r['L2'] for r in smooth],'s-',label='Jacobi-PCG')
    plt.xlabel('solver time (seconds)');plt.ylabel('L2 error');plt.legend();figure('cost_accuracy')
    meta=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
              matplotlib=matplotlib.__version__,sympy=sy.__version__,platform=platform.platform(),
              machine=platform.machine(),repeats_direct=5,repeats_pcg=5,warmups=1,
              threads_requested=1,threadpools=threadpool_info(),
              affine_error=float(affine_error),eigen_min=eigen_min,neumann_order=p)
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
    save('timing_samples', TIMINGS)
    # Spatial scan of the actual steady field, not a transient simulation.
    from matplotlib.animation import FuncAnimation, PillowWriter
    field=np.load(OUT/'smooth_field.npz')
    fig,axes=plt.subplots(2,1,figsize=(6,5),sharex=True)
    numerical,=axes[0].plot([],[],label='computed N=64')
    analytic,=axes[0].plot([],[],'--',label='analytic')
    errline,=axes[1].plot([],[])
    axes[0].set_ylim(-.05,1.05);axes[0].set_ylabel('u');axes[0].legend()
    axes[1].set_xlim(0,1);axes[1].set_ylim(-1.1*np.max(abs(field['error'])),1.1*np.max(abs(field['error'])))
    axes[1].set_xlabel('x');axes[1].set_ylabel('nodal error')
    def scan(j):
        numerical.set_data(field['X'][:,j],field['U'][:,j])
        analytic.set_data(field['X'][:,j],exact(field['X'][:,j],field['Y'][:,j]))
        errline.set_data(field['X'][:,j],field['error'][:,j])
        axes[0].set_title(f'Spatial slice, y={field["Y"][0,j]:.3f} (steady field)')
        return numerical,analytic,errline
    fig.tight_layout()
    animation=FuncAnimation(fig,scan,frames=range(0,65,2),interval=100)
    animation.save(OUT/'spatial_scan.gif',writer=PillowWriter(fps=10));plt.close(fig)
    say('ALL AUTOMATED AUDITS PASSED. Human contribution and viva checks remain team-owned.')
    (OUT/'console.txt').write_text('\n'.join(LOG)+'\n')


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
