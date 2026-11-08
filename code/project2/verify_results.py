"""Check numerical evidence and figure inventory without requiring LaTeX or slides tooling."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def verify(out, figures):
    def rows(name):
        with (out/(name+'.csv')).open() as stream:
            return [{k: (v if k in {'average', 'case', 'method'} else float(v))
                     for k, v in row.items()} for row in csv.DictReader(stream)]

    smooth, interface, solvers, samples, contrast, neumann, sweep = map(
        rows, ['smooth', 'interface', 'solver_comparison', 'timing_samples', 'contrast', 'neumann', 'tolerance_sweep'])
    for table in [smooth, interface, solvers, samples, contrast, neumann, sweep]:
        assert table, 'Empty output table'
        assert all(np.isfinite(v) for row in table for v in row.values() if isinstance(v, (float, int)))
    assert [r['N'] for r in smooth] == [16, 32, 64, 128]

    def order(table, column):
        last = table[-3:]
        return float(np.polyfit(np.log([r.get('h', 1/r['N']) for r in last]),
                                np.log([r[column] for r in last]), 1)[0])

    p2, pi = order(smooth, 'L2'), order(smooth, 'Linf')
    assert 1.95 < p2 < 2.05 and 1.95 < pi < 2.05
    assert all(r['asymmetry'] == 0 and r['backward'] < 1e-12 for r in smooth)
    field = np.load(out/'smooth_field.npz')
    ref = np.sin(np.pi*field['X'])*np.sin(np.pi*field['Y'])
    np.testing.assert_allclose(field['error'], field['U']-ref, atol=1e-14)
    e = field['error'][1:-1, 1:-1]
    row64 = next(r for r in smooth if r['N'] == 64)
    np.testing.assert_allclose([np.linalg.norm(e)/64, np.max(abs(e))],
                               [row64['L2'], row64['Linf']], rtol=1e-9, atol=1e-14)

    harmonic = [r for r in interface if r['average'] == 'harmonic']
    arithmetic = [r for r in interface if r['average'] == 'arithmetic']
    assert len(harmonic) == len(arithmetic) == 4
    assert all(r['Linf'] < 1e-9 and r['flux_max_error'] < 1e-8 for r in harmonic)
    assert all(r['minimum'] >= -1e-10 and r['maximum'] <= 1+1e-10 for r in interface)
    pa = order(arithmetic, 'L2')
    assert .9 < pa < 1.2
    for average in ['arithmetic', 'harmonic']:
        profile, flux = rows('profile_'+average), rows('flux_'+average)
        assert len(profile) == 64 and len(flux) == 63
        u = np.array([r['u'] for r in profile])
        q = np.array([r['q'] for r in flux])
        row = next(r for r in interface if r['N'] == 63 and r['average'] == average)
        j = 31
        reconstructed = -(20/11)*(u[j+1]-u[j])*63
        jump = abs(q[j-1]-q[j+1])
        defect = max(abs(reconstructed-q[j-1]), abs(reconstructed-q[j+1]))
        np.testing.assert_allclose([reconstructed, jump, defect],
                                   [row['resistance_flux'], row['flux_jump'], row['physical_flux_defect']], atol=1e-10)
        np.testing.assert_allclose(np.max(abs(q+20/11)), row['flux_max_error'], atol=1e-11)

    for row in samples:
        np.testing.assert_allclose(row['total_s'], row['setup_s']+row['factor_s']+row['solve_s'], rtol=1e-10, atol=1e-14)
    for row in solvers:
        spatial = next(r['L2'] for r in smooth if r['N'] == row['N'])
        assert row['gap_L2'] < spatial*1e-5
        for method, keys in [('LU', ['direct_s', 'direct_min_s', 'direct_max_s']),
                             ('PCG', ['pcg_s', 'pcg_min_s', 'pcg_max_s'])]:
            values = [r['total_s'] for r in samples if r['case'] == f"smooth_N{int(row['N'])}" and r['method'] == method]
            assert len(values) == 5
            np.testing.assert_allclose([np.median(values), min(values), max(values)], [row[k] for k in keys], rtol=1e-12)
    assert len(contrast) == 12 and all(r['residual'] < 2e-9 for r in contrast)
    assert all(r['null_defect'] < 1e-8 and abs(r['mean_u']) < 1e-10 for r in neumann)
    pn = order(neumann, 'L2')
    assert 1.9 < pn < 2.1
    assert sweep[-1]['gap_L2'] < sweep[0]['gap_L2'] and max(r['L2'] for r in sweep)/min(r['L2'] for r in sweep) < 1.001
    metadata = json.loads((out/'metadata.json').read_text())
    assert metadata['repeats_direct'] == metadata['repeats_pcg'] == 5
    assert metadata['warmups'] == 1 and metadata['threads_requested'] == 1
    assert metadata['threadpools'] and all(p['num_threads'] == 1 for p in metadata['threadpools'])
    expected = ['smooth_convergence', 'smooth_field', 'matrix_spy', 'interface_profiles',
                'cost_accuracy', 'contrast_iterations', 'neumann_convergence']
    for stem in expected:
        for extension in ['.png', '.pdf']:
            file = figures/(stem+extension)
            assert file.is_file() and file.stat().st_size > 0, file
    assert (figures/'spatial_scan.gif').is_file()
    result = dict(status='PASS', smooth_L2_order=p2, smooth_Linf_order=pi,
                  arithmetic_L2_order=pa, neumann_L2_order=pn, timed_samples_consistent=True,
                  field_norms_consistent=True, layered_flux_diagnostics_consistent=True,
                  figure_files_checked=15)
    (out/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'code/project2/outputs')
    parser.add_argument('--figure-dir', type=Path, default=ROOT/'figures/project2')
    args = parser.parse_args()
    print(json.dumps(verify(args.output_dir, args.figure_dir), indent=2))
