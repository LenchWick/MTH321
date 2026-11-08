"""Independent whole-experiment seeds for the visualization guide.

These smaller supplementary runs do not replace the main convergence table.
Student-t bands quantify the mean fitted log-error across independent seeds.
They do not cover numerical-reference bias or model-fit error.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
from kernels import gbm, sv


def independent_seed_check(out):
    grids = np.array([8, 16, 32, 64, 128, 256])
    finest, paths, batch = 4096, 4096, 512
    errors = {'EM': [], 'Milstein': [], 'nonaffine': []}
    seeds = [20261001, 20261002, 20261003, 20261004, 20261005]
    for seed in seeds:
        sums = {name: np.zeros(len(grids)) for name in errors}
        rng = np.random.default_rng(seed)
        for _ in range(paths // batch):
            d1 = rng.normal(size=(finest, batch)) / np.sqrt(finest)
            d2 = -.7*d1 + np.sqrt(1-.7**2)*rng.normal(size=d1.shape)/np.sqrt(finest)
            exact = gbm(d1, 1/finest, 'exact')
            ref = sv(d1, d2, 1/finest)[0]
            for j, n in enumerate(grids):
                w1 = d1.reshape(n, finest//n, batch).sum(axis=1)
                w2 = d2.reshape(n, finest//n, batch).sum(axis=1)
                for method in ['EM', 'Milstein']:
                    sums[method][j] += abs(gbm(w1, 1/n, method)-exact).sum()
                sums['nonaffine'][j] += abs(sv(w1, w2, 1/n)[0]-ref).sum()
        for name in errors:
            errors[name].append((sums[name]/paths).tolist())
    h = 1/grids
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    fits = {}
    for name, ax, color, order, window in [
        ('EM', axes[0], 'tab:blue', .5, slice(1, None)),
        ('Milstein', axes[0], 'tab:orange', 1., slice(1, None)),
        ('nonaffine', axes[1], 'tab:blue', .5, slice(None, 5))
    ]:
        e = np.asarray(errors[name])
        coefficients = np.array([np.polyfit(np.log(h[window]), np.log(row[window]), 1)
                                 for row in e])
        lx = np.linspace(np.log(h[window]).min(), np.log(h[window]).max(), 100)
        lines = coefficients[:, :1]*lx + coefficients[:, 1:]
        center = lines.mean(axis=0)
        half = t.ppf(.975, len(seeds)-1)*lines.std(axis=0, ddof=1)/np.sqrt(len(seeds))
        for row in e:
            ax.loglog(h, row, color=color, alpha=.22, lw=.9)
        ax.loglog(h, e.mean(axis=0), 'o-', color=color, label=name+' mean error')
        ax.loglog(np.exp(lx), np.exp(center), '--', color=color, label=name+' fit')
        ax.fill_between(np.exp(lx), np.exp(center-half), np.exp(center+half),
                        color=color, alpha=.18, label='95% seed-mean fitted band')
        ax.loglog(np.exp(lx), e.mean(axis=0)[-1]*(np.exp(lx)/h[-1])**order,
                  ':', color=color, label=f'reference order {order:g}')
        p = coefficients[:, 0]
        hw = t.ppf(.975, len(seeds)-1)*p.std(ddof=1)/np.sqrt(len(seeds))
        fits[name] = {'seed_slopes': p.tolist(), 'mean_slope': float(p.mean()),
                      't95_interval_mean_slope': [float(p.mean()-hw), float(p.mean()+hw)]}
    for ax, title in zip(axes, ['GBM: independent seeds agree with theoretical orders',
                               'Non-affine EM: finite-reference half-order regime']):
        ax.set_title(title)
        ax.set_xlabel('step h [year]')
        ax.set_ylabel('terminal L1 price error [price units]')
        ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(out/('independent_seed_bands.'+ext), dpi=250, bbox_inches='tight')
    plt.close(fig)
    record = {'seeds': seeds, 'paths_per_seed': paths, 'finest_grid': finest,
              'grids': grids.tolist(), 'errors': errors, 'fits': fits,
              'interpretation': 'approximate Student-t band for mean fitted log-error; no reference-bias coverage'}
    (out/'independent_seeds.json').write_text(json.dumps(record, indent=2))
    print('Independent seed check:', json.dumps(fits, indent=2))


if __name__ == '__main__':
    independent_seed_check(Path(__file__).resolve().parents[1]/'outputs')
