# Model and numerical method

## Physical model

We model steady heat conduction in a heterogeneous medium on `Omega=(0,1)^2`:

$$-\nabla\cdot(a\nabla u)=f,\qquad u=g\text{ on }\partial\Omega,\qquad a>0.$$

Temperature is `u`, conductivity is `a`, and heat flux is `q=-a grad u`.
The divergence form preserves conservation. Ideal material contact without a
singular interface source imposes continuous temperature and continuous normal
flux; the temperature gradient may jump.

## Exact benchmarks

**Smooth manufactured solution:** `a=1+x`, `u*=sin(pi x) sin(pi y)`, zero
Dirichlet temperature, and

$$f=2\pi^2(1+x)\sin(\pi x)\sin(\pi y)-\pi\cos(\pi x)\sin(\pi y).$$

The implementation checks this expression independently using symbolic
differentiation before solving. Heat flux is continuous but varies in space.

**Layered reference:** `a=1` for `x<1/2`, `a=10` for `x>1/2`, `f=0`, and full
analytic boundary data from

$$u_*(x)=\begin{cases}20x/11,&x\leq1/2,\\(2x+9)/11,&x\geq1/2.\end{cases}$$

Both sides give `u*(1/2)=10/11` and `q_x=-20/11`. The top and bottom boundaries
also carry this analytic profile. The arithmetic numerical solution is
generally two-dimensional and is not assumed to be a perturbed 1D solution.
Odd interval counts place the interface halfway between nodal columns.

## Assembly and boundary treatment

With `h=1/N`, interior index `k=(i-1)(N-1)+(j-1)`, the conservative stencil is

$$ (A_hU)_{ij}=h^{-2}\left[a_E(U_{ij}-U_{i+1,j})+a_W(U_{ij}-U_{i-1,j})
+a_N(U_{ij}-U_{i,j+1})+a_S(U_{ij}-U_{i,j-1})\right]. $$

The same positive coefficient is shared by adjacent rows. A known boundary
neighbour adds `a_face*g/h^2` to the right-hand side; its weight stays in the
diagonal. Only interior values enter the unknown vector.

For smooth `a`, the arithmetic face average is a second-order midpoint
approximation. At an interface, segment resistances give

$$a_{\rm face}=\frac{h}{d_L/a_L+d_R/a_R}.$$

For a bisected edge this becomes `2*aL*aR/(aL+aR)`. Our layered test compares
the harmonic value `20/11` with the deliberately inappropriate arithmetic
value `11/2`.

For a nonzero interior vector extended by zero to the boundary,

$$h^2v^TA_hv=\sum_{\{p,r\}\text{ interior}}a_{pr}(v_p-v_r)^2+
\sum_{p,b\text{ boundary}}a_{pb}v_p^2>0.$$

Here each undirected edge is counted once. Positive weights, shared edges
and connection to Dirichlet values prove symmetry and positive definiteness.
The matrix sign pattern also gives an M-matrix and the discrete maximum
principle for the source-free layered problem. Boundedness alone does not
establish accuracy. A sparse eigenvalue estimate corroborates the energy
proof; it is not a replacement for it.

## Solvers and software boundaries

The project explicitly implements sparse COO edge assembly, conductivity
averaging, boundary elimination, norms and Neumann compatibility treatment.
SciPy provides SuperLU (`splu`), conjugate gradients (`cg`) and ARPACK
(`eigsh`). Jacobi-PCG uses the cached inverse diagonal, zero initial guess,
`atol=0` and a stated relative tolerance. Nonzero solver status is rejected.
The independently recomputed residual is also recorded.

The pure-Neumann extension uses a conservative cell-centred graph with
`A*1=0`. A substantial nonzero mean of the right-hand side is rejected before
removing floating-point imbalance. A normalized constant-vector constraint
fixes `mean(U)=0`. The augmented system is symmetric indefinite and is solved
by sparse LU, not ordinary SPD CG.

## Algorithm

1. Predict signs, expected orders, interface flux and nullspace behaviour.
2. Independently verify manufactured forcing.
3. Evaluate conductivity, source and full boundary data on each mesh.
4. Assemble shared interior edges and add known boundary contributions.
5. Solve with sparse LU or justified Jacobi-PCG and recompute the residual.
6. Compare against analytic references in named norms.
7. Repeat meshes, closure choices, tolerances and contrasts under matched settings.
8. Record raw tables, timings, environment and figures; run the evidence checks.

The root runner adds file organization and verification. The imported numerical
core differs from the audited source only by an optional output-directory
environment variable; mathematical assembly and solvers are unchanged.
