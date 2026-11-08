# MTH321 Project 1 — Financial SDEs

## Install and run

Use Python 3.10 or newer. From this extracted folder:

```sh
python -m pip install -r requirements.txt
python code/validation_tests.py
python code/run_all.py
```

The full run uses 65,536 GBM paths and 32,768 non-affine paths, batched in
groups of 1,024. Allow several minutes and approximately 1–2 GB of available
memory. No data downloads or network access are used by the experiment.
Outputs are written beside this README in `outputs/`, regardless of the
current working directory. An assertion failure stops the run.

## Files and ownership

- `code/run_all.py`: directly implemented price-space EM, scalar Milstein,
  exact GBM, joint log-state EM, moment-balanced factor update, all experiments,
  12 integrated checks, animation, and JSON export.
- `code/kernels.py`: importable copies of the core definitions for tests.
- `code/validation_tests.py`: eight independent analytic/special-case tests.
- `code/submission_figures.py`: reference-rate overlays, independent-batch
  bootstrap fitted-line bands, and a measured error-versus-cost comparison.
- `code/repeat_seeds.py`: five independent whole-experiment seeds, with
  Student-t bands for the mean fitted log-error.
- `requirements.txt`: portable dependency bounds.
- `archived_outputs/`: the original recorded simulation used in report tables.
- `outputs/`: freshly regenerated results, including `fit_bands.json`.
- `verification_run.log`: actual complete run output for this preparation.
- `validation_tests.log`: actual independent test results.

Numerical methods and analytic checks are implemented here, rather than
delegated to an SDE solver. NumPy supplies array arithmetic and random draws;
SciPy `expit`, `ndtr`, and `ndtri` supply sigmoid/normal scalar functions;
Matplotlib and Pillow produce static figures and GIFs. AI assistance is
disclosed in the report. Personal authorship must be confirmed by the team.
The GBM-specific product and moment-balanced factor update deliberately use
the model's structure; this package is not advertised as a general ODE solver.

## Reproduce figures

One full run regenerates:

| File | Evidence |
|---|---|
| `gbm_convergence` | Strong and signed mean bias with pointwise intervals |
| `gbm_distribution` | Histograms and exact-law quantiles |
| `nonaffine_convergence` | Three nested reference grids; unresolved payoff bias |
| `extension` | Balanced update: price-error comparison and factor variance |
| `cost` | Recorded-run integration timings at the matched target |
| `convergence_bands` | Theoretical orders and 500 whole-batch bootstrap bands |
| `independent_seed_bands` | Five separate seeds and mean fitted log-error bands |
| `error_vs_cost` | Archived wall time versus achieved error at ceiling 0.5 |
| `paths.gif` | Six evolving price/volatility paths, not a statistical test |

All static figures have PNG (250 DPI) and vector PDF versions. Core report
figures are also encoded as inline numerical data in the LaTeX report.
Supplementary figure captions state their interpretation and source.
Signed weak biases use a linear/symmetric-log vertical axis where necessary
because negative estimates and intervals crossing zero cannot be logged.

## Statistical and timing conventions

Fixed base seed: 20260916, with separate experiment streams. Strong comparisons
share Brownian trajectories between methods and grids. Coarse increments sum
fine increments. Distinct path batches are independent; resampling whole
batches preserves cross-grid dependence. The displayed fitted-line band is
a pointwise bootstrap sampling band, not a bound on model/reference error.
The main table uses independent batches as allowed by the problem pack.
In addition, five smaller whole-experiment repetitions use independent seeds,
4,096 paths per seed and a 4,096-step reference. Their approximate Student-t
band concerns the mean fitted log-error across seeds; it does not cover
reference or model-fit bias. `independent_seeds.log` records those checks.

Integration timings are seven-repeat medians after warm-up on 8,192 paths.
They exclude RNG, reference solves, imports, input construction, plotting,
animation, and bootstrap calculations. The common terminal L1 error ceiling
is 0.5, using the upper sampling endpoint for grid selection.
`error_vs_cost` intentionally reads `archived_outputs/results.json` to match
the timing measurement quoted in the report; fresh timings are separately
stored in `outputs/results.json`. Exact GBM has zero discretisation error
and is discussed separately from nonzero-error comparisons.

Platform/library versions and wall times may differ. Do not replace archived
report tables with new timings without updating all associated statements.
No historical Git activity or personal participation is asserted by this ZIP.
