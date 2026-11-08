"""Additional evidence requested by the course visualization guide."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def fit_band(ax, hs, batch_errors, theoretical_order, label, color, seed):
    """Resample independent batches, keeping the whole grid vector together."""
    x = np.log(np.asarray(hs))
    batch_errors = np.asarray(batch_errors)
    rng = np.random.default_rng(seed)
    grid = np.linspace(x.min(), x.max(), 100)
    coefs = np.array([
        np.polyfit(x, np.log(batch_errors[
            rng.integers(len(batch_errors), size=len(batch_errors))
        ].mean(axis=0)), 1) for _ in range(500)
    ])
    predicted = coefs[:, :1] * grid + coefs[:, 1:]
    lo, hi = np.quantile(predicted, [.025, .975], axis=0)
    central = np.polyfit(x, np.log(batch_errors.mean(axis=0)), 1)
    ax.fill_between(np.exp(grid), np.exp(lo), np.exp(hi), color=color,
                    alpha=.22, label=label + ' batch-bootstrap 95% band')
    ax.plot(np.exp(grid), np.exp(central[0]*grid+central[1]), color=color,
            ls='--', label=f'{label} fit p={central[0]:.4f}')
    theory = batch_errors.mean(axis=0)[-1] * (np.exp(grid)/hs[-1])**theoretical_order
    ax.plot(np.exp(grid), theory, color=color, ls=':',
            label=f'{label} reference order {theoretical_order:g}')
    return dict(slope=float(central[0]),
                slope_ci=np.quantile(coefs[:, 0], [.025, .975]).tolist(),
                h=np.exp(grid).tolist(), fit_lower=np.exp(lo).tolist(),
                fit_upper=np.exp(hi).tolist())


def supplementary_figures(out, gbm, sv, batches, gerr, st, hg, hs, m_sv, batch):
    fits = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for method, color in [('EM', 'tab:blue'), ('Milstein', 'tab:orange')]:
        rows = gbm[method]['strong']
        axes[0].errorbar(hg, [v['mean'] for v in rows],
                         yerr=[1.96*v['se'] for v in rows], marker='o',
                         color=color, label=method + ' data')
        fits[method] = fit_band(axes[0], hg[2:], np.asarray(batches[method])[:, 2:],
                                .5 if method == 'EM' else 1., method, color, 17)
    ref = sv['reference_comparisons']['4096']
    axes[1].errorbar(hs[:6], [v['mean'] for v in ref['strong']],
                     yerr=[1.96*v['se'] for v in ref['strong']], marker='s',
                     color='tab:blue', label='EM vs 4096-step reference')
    e = abs(st[:6]-st[-1]).reshape(6, m_sv//batch, batch).mean(axis=2).T
    fits['nonaffine'] = fit_band(axes[1], hs[:5], e[:, :5], .5, 'EM', 'tab:blue', 17)
    for ax in axes:
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('step h [year]')
        ax.set_ylabel('terminal L1 price error [price units]')
        ax.legend(fontsize=7)
    axes[0].set_title('GBM: measured rates follow orders 1/2 and 1')
    axes[1].set_title('Non-affine EM: empirical half order, finite reference')
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(out / ('convergence_bands.'+ext), dpi=250, bbox_inches='tight')
    plt.close(fig)
    (out/'fit_bands.json').write_text(json.dumps(fits, indent=2))

    # Use the archived measurement quoted by the report, not a newly timed run.
    recorded = json.loads((Path(out).parent/'archived_outputs/results.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, cost, title in [
        (axes[0], recorded['cost'], 'GBM: Milstein needs less integration work'),
        (axes[1], recorded['extension']['cost'], 'Non-affine: no general speed advantage')
    ]:
        for method, color, marker in [('EM', 'tab:blue', 'o'),
                                     ('Milstein', 'tab:orange', 's'),
                                     ('balanced', 'tab:orange', 's')]:
            if method not in cost['methods']:
                continue
            v = cost['methods'][method]
            x = v['seconds_median']*1000
            ax.errorbar(x, v['error']['mean'], yerr=1.96*v['error']['se'],
                        color=color, marker=marker, ms=8, label=f'{method}, N={v["n"]}')
        ax.axhline(cost['target'], color='black', ls='--', label='common error ceiling 0.5')
        ax.set_xlabel('recorded median integration wall time [ms]')
        ax.set_ylabel('terminal L1 price discrepancy [price units]')
        ax.set_title(title); ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(out / ('error_vs_cost.'+ext), dpi=250, bbox_inches='tight')
    plt.close(fig)
