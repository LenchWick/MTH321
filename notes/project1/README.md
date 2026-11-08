# Project 1 provenance and verification

This is the Group D financial SDE project: exact GBM, price-space Euler-Maruyama and scalar Milstein, followed by the specified non-affine stochastic-volatility model and a moment-balanced factor update.

The [signed report](../../report/project1/project1_report.pdf) contains the mathematical derivations, uncertainty analysis, limitations, AI transparency log and ICS records. The [speaker script](speaker_script.md) accompanies the approved presentation.

## Reproduce from the repository root

```sh
python -m pip install -r code/project1/requirements.txt
python code/project1/code/validation_tests.py
python code/project1/code/run_all.py
```

The full run uses 65,536 GBM and 32,768 non-affine paths. It can take several minutes and around 1-2 GB of memory. It writes to `code/project1/outputs/`; copy that package elsewhere before running if the checked-in outputs should remain unchanged. The `archived_outputs/` measurements remain separate and are used for the report's historical timing comparison. The independent analytic tests are fast and do not require a full Monte Carlo rerun.

## Evidence boundaries

The checked-in tables, figures and logs are October 2026 records, not a new November full experiment. GBM has an exact pathwise benchmark; non-affine errors use finite numerical references. Reported intervals describe sampling uncertainty, not certified continuum errors. In the non-affine payoff experiment a reliable empirical weak order was not established. Preserving selected factor moments does not prove universally better price accuracy or lower cost.

AI assistance is disclosed in Appendix A of the report. This GitHub commit imports the existing project on 8 November 2026; it does not reconstruct earlier Git development history or fabricate individual contributions. The existing signed personal statements are retained.
