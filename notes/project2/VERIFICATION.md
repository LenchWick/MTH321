# Verification and reference results

The implementation uses discrete nodal norms

$$E_2=h\sqrt{\sum_{i,j}|U_{ij}-u_*(x_i,y_j)|^2},\qquad
E_\infty=\max_{i,j}|U_{ij}-u_*(x_i,y_j)|.$$

For `r=b-AU`, the reported residual is
`rho=||r||2/max(||b||2,1)`. The backward-error measure is
`eta=||r||inf/(||A||inf*||U||inf+||b||inf)`, with zero defined for the exact
zero solution. Residual and solution error are different quantities.

## Smooth refinement

| N | E2 | E-infinity | Relative residual |
| --- | --- | --- | --- |
| 16 | 1.5975864283e-3 | 3.2050306080e-3 | 7.9733254340e-15 |
| 32 | 3.9879758022e-4 | 8.0094066771e-4 | 3.3104872374e-14 |
| 64 | 9.9662020102e-5 | 2.0029767177e-4 | 1.3730634780e-13 |
| 128 | 2.4913170093e-5 | 5.0069739232e-5 | 5.5229879464e-13 |

Source: [smooth.csv](../../code/project2/outputs/smooth.csv). Fitting log(error)
against log(h) on the last three meshes gives E2 order `2.00033807044` and
maximum-norm order `1.99984226369`. Pairwise orders and fit assumptions are
discussed in the report; round-off-level harmonic errors are not fitted.

## Independent checks

The ten tests in [test_solver.py](../../code/project2/test_solver.py) check:

- Constant-coefficient assembly against a separately constructed Kronecker matrix.
- The zero solution, including integer-valued zero input data.
- Nonzero affine boundary data and the hand-derived first corner row.
- Callable forcing against an independently obtained symbolic derivative.
- Exact harmonic layered nodes and a falsified arithmetic closure.
- Neumann constant nullspace and mean-zero gauge.
- Rejection of an incompatible Neumann source.
- The smallest allowed grid and invalid grids, coefficients and data.

At N=8, the affine audit gives `A00=288`, `A01=-72`, `A07=-76`, `b0=25` for
`u=x+2y`. The recorded smallest eigenvalue is `28.33254682968626`. The
Frobenius asymmetry is zero on all smooth meshes. The energy identity is the
SPD proof; the numerical eigenvalue is corroboration.

## Interface physics

The source-free reference has `q_exact=-20/11`. Harmonic averaging reproduces
the analytic nodal field to floating-point solve accuracy. Arithmetic E2 errors
fit order `1.01058164051` on the last three odd meshes. Both closures remain
in `[0,1]`, so the maximum principle alone does not validate interface accuracy.

Let the interface cut-edge index be `m=floor(N/2)`. The one-sided jump uses
the entirely-left and entirely-right adjacent faces:

`J = |q[m-1] - q[m+1]|`.

On this particular symmetric test, J is tiny for both closures. We therefore
independently reconstruct the cut-edge flux using the physical series resistance:

`q_res = -(20/11) * (U[m+1,j]-U[m,j]) / h`.

The constitutive mismatch is
`D_phys=max(|q_res-q[m-1]|,|q_res-q[m+1]|)`.
It compares reconstruction with adjacent numerical fluxes, rather than with
the exact flux. At N=63, arithmetic gives `q_res=-0.608635182421163` and
`D_phys=1.2324797290530936`. Small algebraic residuals, symmetry and a tiny J
can therefore conceal a wrong interface resistance.

## Solver accuracy and cost

The N=64 tolerance sweep tightens PCG from `1e-4` to `1e-12`; the gap to LU
falls while the spatial E2 stays near `9.9662e-5`. Both LU and PCG operate on
the same matrix, so this cross-check cannot independently establish assembly.

For the main solver comparison, LU includes factorisation plus triangular
solution; PCG includes Jacobi setup and iterations. Each has one excluded
warm-up and five measured repetitions. The total-time median and sample
min/max come from paired raw totals, not from adding separate medians.
Assembly and plotting are excluded. Loaded BLAS pools are limited to one thread.

At N=128, the reference medians are LU `36.496375 ms` and PCG `19.847083 ms`.
These are historical machine measurements, not universal complexity or speed
claims. Sample ranges are not confidence intervals. An independently
recomputed finite-precision CG residual can slightly exceed the requested
recursive stopping tolerance; the tight gap to LU remains negligible compared
with spatial error.

The contrast experiment uses the separate load `f=1`, `g=0`, harmonic faces,
and PCG `rtol=1e-10`. Its iteration counts must not be directly compared with
the sine benchmark as if the right-hand sides were identical.

## Neumann extension and limits

The compatible cosine reference yields observed E2 order `2.00035224442`.
The graph has a constant null vector; the computed temperature has mean zero.
Adding one to the source must raise an exception before projection. The
observed rate concerns the selected smooth data and is not a general theorem
for arbitrary Neumann problems.

Other limits are the square geometry, controlled planar interface, selected
loads, one horizontal flux cut and one historical timing environment. The
project does not establish curved-interface accuracy or a universal solver ranking.
