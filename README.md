# MTH321 — Group D

## Project 1 — Financial SDEs

[Report PDF](report/project1/project1_report.pdf) · [Slides PDF](report/project1/project1_slides.pdf) · [LaTeX sources](report/project1/) · [Code and recorded outputs](code/project1/) · [Figures](figures/project1/) · [Speaker script and provenance](notes/project1/)

Exact GBM benchmarking, Euler–Maruyama/Milstein convergence and non-affine stochastic volatility. The latest report has 39 pages, with the existing signed AI log and ICS appendix; slides have 17 pages.

From the repository root:

```sh
python -m pip install -r code/project1/requirements.txt
python code/project1/code/validation_tests.py
python code/project1/code/run_all.py
```

The full Monte Carlo run takes several minutes and writes to `code/project1/outputs/`. See [Project 1 instructions](notes/project1/README.md) for archive preservation, computational cost and evidence limits. Project 1 has its own dependency file and does not change the Project 2 root entry point.

---

## Project 2 — Variable diffusion


Project 2: conservative finite differences for variable diffusion,
`-div(a grad u) = f`, on the unit square.

**Start here:** [numerical implementation](code/project2/run_all.py),
[independent tests](code/project2/test_solver.py),
[figures and raw-data index](figures/project2/README.md),
[model and algorithm](notes/project2/METHOD.md),
[verification and measured results](notes/project2/VERIFICATION.md).

## Reproduce the project

Python 3.10 or newer is required. Create and activate a virtual environment,
then run these commands from the repository root:

```sh
python -m pip install -r requirements.txt
python run_all.py
```

The command runs the ten regression tests, all numerical experiments, exports
the figures to `figures/project2/`, and checks the resulting evidence. No
LaTeX, browser, API key or private data is required.

To validate the checked-in results without running new timing experiments:

```sh
python run_all.py --check-only
```

To reproduce into a separate directory and preserve the published reference:

```sh
python run_all.py --output-dir build/project2-rerun
```

The separate run contains its own CSV tables, figures, environment metadata,
test log and `verification.json`. Failures raise an exception and return a
nonzero exit code.

## Repository layout

```text
code/project2/       conservative assembly, solvers, tests and raw outputs
figures/project2/    PNG/PDF figures and the stationary spatial-slice GIF
notes/project2/      model, verification, AI audit and data provenance
report/project2/     presentation files and final-report status
run_all.py          one-command reproduction entry point
requirements.txt    numerical and plotting dependencies
```

The original `code / figures / notes / report` framework is retained. Project 2
files use their own subdirectories so other course projects can coexist.

## Experiments

| Requirement | Implemented evidence |
| --- | --- |
| Conservative FDM, coefficients and Dirichlet elimination | Shared positive face coefficients; explicit interior assembly; affine boundary test |
| Symmetry and positive definiteness | Energy identity; zero Frobenius asymmetry; small-grid eigenvalue check |
| Four-mesh smooth verification | N=16,32,64,128; discrete E2/E-infinity; fitted orders and residuals |
| Interface comparison | N=31,63,127,255; arithmetic versus harmonic/resistance coefficients; temperature and flux diagnostics |
| Maximum principle | Full layered fields stay in the boundary range [0,1] |
| Solver difficulty | Jacobi-PCG iterations versus mesh and conductivity contrast 1,10,100,1000 |
| Cost versus accuracy | Matched LU/PCG solutions; one warm-up and five measured repetitions; raw samples retained |
| Challenge extension | Pure Neumann compatibility rejection, constant nullspace and a mean-zero gauge |

## Reference outputs and interpretation

The checked-in numerical data were measured on **7 October 2026**, and support
the current presentation. Their environment is recorded in
[metadata.json](code/project2/outputs/metadata.json). The numerical evidence
was audited again on **8 November 2026**. Times depend on the machine and
software environment; new runs need not reproduce the historical times or the
last floating-point digits. A rerun does not automatically rewrite an existing
presentation or report.

Core numerical methods, edge assembly, boundary handling and diagnostics are
explicit in our code. SciPy supplies sparse LU, CG and the eigenvalue solver;
SymPy independently checks manufactured forcing. See
[method details and library boundaries](notes/project2/METHOD.md).

AI assistance and the limits of the evidence are described in
[AI_AUDIT.md](notes/project2/AI_AUDIT.md). This upload records a present-day
import and repository integration. It does not reconstruct earlier development
history or claim personal contributions. Each member must supply their own
actual ICS and be able to explain and rerun the implementation.

## Presentation and report

The [slides PDF](report/project2/Project2_Presentation.pdf) and
[editable presentation](report/project2/Project2_Presentation.pptx) accompany
the code. The report directory explains the status of the separate final report.
The course also requires a report PDF and runnable code ZIP; a GitHub link does
not replace those submissions.

