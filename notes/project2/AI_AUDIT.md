# AI assistance and claim audit

Codex assisted with mathematical exposition, implementation, experimental
design, numerical execution, verification, presentation layout and repository
integration. The course problem pack supplied the core benchmarks. The team
must review, explain and reproduce the work and document their actual changes.
This file does not assert human signatures, personal ICS claims or historical
development activity.

| Claim and assumptions | Risk / possible falsifier | Independent evidence | Decision and limit |
| --- | --- | --- | --- |
| Smooth manufactured forcing with a=1+x | Omitting the derivative of a solves a different PDE | Symbolic differentiation, exact-reference refinement and the independent affine boundary test | Formula and second-order results pass for the selected smooth case |
| Shared positive face assembly is SPD after Dirichlet elimination | Transposed indices or boundary row replacement can hide behind a plausible sine field | Kronecker operator, nonsymmetric affine field, hand corner row, zero asymmetry and the energy identity | Algebraic structure passes; a positive random Rayleigh quotient is not a proof |
| Harmonic resistance reproduces the layered interface | A symmetric matrix, bounded temperature and small scheme-flux jump can hide wrong resistance | Analytic piecewise-linear reference and independent cut-edge physical reconstruction | Harmonic passes; arithmetic is a deliberate comparison. One cut is not a global flux norm |
| A pure-Neumann solution needs compatible source and a gauge | Projecting an incompatible source would silently change the problem | Constant null vector, mean-zero solve, independent residual and explicit rejection of b+1 | Reject incompatibility first; fix uniqueness separately |

## Five audit gates

1. **Transcription:** state source, coefficient, full boundary data, flux sign and independent forcing checks.
2. **Structure and indexing:** use shared edges, the explicit index map, the corner row and independent operator tests.
3. **Linearisation:** the selected diffusion problems are linear. The indefinite Neumann saddle system is distinguished from the SPD Dirichlet system.
4. **Independent reproduction:** the one-command runner executes tests and experiments, records environment/raw timings and verifies output consistency. Team members still need to run it themselves.
5. **Interrogation:** explain how a tiny residual can coexist with wrong interface physics, why arithmetic can be bounded yet inaccurate, and why compatibility differs from gauge fixing.

## Repository provenance

The reference measurements date from 7 October 2026. A later numerical/slides
audit on 8 November 2026 re-solved 29 linear systems and checked chart data and
raw timing medians. It clarified the profile cut and the definition of D_phys.
Repository integration then adds output-directory support, a root runner,
evidence verification and automation. These integration changes are tested
as actual present work. No earlier Git history is manufactured.
