# ReSuS — MATLAB

## Overview
Relaxed subset simulation estimates failure probabilities with hybrid LSF / Gaussian-tail intermediate domains. This revision restores the **original sequential probability-product estimator**, Eqs. (21)–(24) of the cited paper, in place of the previous historical-sample MIS estimator. The default sampler remains `ES`; the four unified kernels remain unchanged.

## Structure
- `drive_resus.m` — example selection, settings and execution.
- `SeMCS.m` — iteration and stored levels; `Support-Fun/ReSuSProduct.m` — native probability-product postprocessing.
- `ReSuS.m` — hybrid distributions, thresholds, seed selection and convergence.
- `Limit-state-functions/` — `LSF_2D_PiecewiseLinear.m`, `LSF_2D_ChangingTopological.m`, `LSF_2D_LinearLog.m`, examples A–C.
- `MCMC-samplers/` — `MCMCBase.m`, `CW_aMH.m`, `aCS.m`, `HMC.m`, `ES.m`, the locally included sampling interface and kernels.
- `Operations on Logarithmic values/logGauss.m` — standard-normal log density used by the framework and model.

## Quick start
Run from this repository's root:

```matlab
out = drive_resus('A',1000,411,'ES','paper');
out_c = drive_resus('C',1000,411,'aCS','paper');
```
The returned result includes `pf`, `log_pf`, `conditional_probabilities`, actual `n_calls`, `converged`, and sampled levels. `converged = false` is an unfinished run, not a validated estimate. `estimator` identifies this revision as `product`.

## Examples and samplers
| Example | Model / dimension | Driver selector |
|---|---|---|
| A | Piecewise-linear LSF; 2 independent standard normals | `A` |
| B | Nonlinear changing-topology LSF; 2 independent standard normals | `B` |
| C | Linear / logistic LSF, paper Eq. (26); 2 independent standard normals | `C` |

- Available samplers: `CW_aMH`, `aCS`, `HMC`, `ES`; default `ES`. The quick start shows explicit selection.
- Defaults are unchanged: example A, N=5000, seed 1, p0=0.1 and `legacy` profile. Explicit `paper` runs above use N=1000; the level budget defaults to 100.
- `legacy` retains midpoint thresholds and strict seed eligibility. `paper` uses order-statistic thresholds, inclusive eligibility and Algorithm 1 chain-count rounding. **Both profiles now use the product estimator**; selecting `paper` does not replace the unified ES with a historical ESS kernel.
- C uses the paper's unscaled `5-x1` branch, not the `0.03` scale present in an earlier supplied file. No D–H models are bundled.
- The `SeMCS` constructor accepts a model and its dimension. Models return one finite LSF value per column of a dimension-by-sample input. Custom physical/nonnormal models must be mapped to Gaussian/Nataf reference coordinates; common kernels are not restricted to 2D.

## Dependencies
- MATLAB with `classdef` and implicit expansion (R2016b or later); no Statistics Toolbox is needed for the bundled examples.
- Each repository runs independently, with no sibling-toolbox import or external example data.

## Method notes
- **Implementation:** each transition contributes its empirical hybrid-event probability on the **current sampled level**, divided by the actual number of samples (which can differ from requested N). The last factor is the true-failure proportion `g <= 0`; no extra unsampled hybrid factor is included. Factors are accumulated in log space. A zero estimate is retained as `pf=0`, `log_pf=-Inf`.
- **Sampler policy:** `CW_aMH` uses seed-spread-scaled **uniform random-walk increments in whitened Gaussian-reference coordinates**. `ES` is the coordinatewise initial-angle `[0, pi]` version, not a claim of identical historical scalar-angle ESS. `aCS` starts with `sigma0=1`, `lambda0=0.6` and carries the final lambda between levels (five adaptation groups, target acceptance 0.44). `HMC` uses analytic Gaussian-reference rotation; its adaptive trajectory time is carried between levels. It is not generic physical-coordinate gradient HMC.
- **Shared interface:** `SamGen(seed_u, seed_y, Nsteps)` returns `(u, y, Ncal, BayStc)` (MATLAB uses multiple outputs); arrays are dimension/function × chain × step. `G` carries thresholds and adaptation. Local copies of common kernels keep each repository independently runnable; no sibling repository is required.
- **Outputs:** `conditional_probabilities` contains the factors, `log_normalizers` their cumulative level normalizers, and `sample_counts` actual level sizes. `n_calls` counts all LSF evaluations, including rejected proposals. Standard deviation across independent runs is not the standard error of their mean. Random seeds do not identify equal streams across languages.

## Reference & license
Binbin Li, Weili Xia and Zihan Liao. **Relaxed subset simulation for reliability estimation.** *Reliability Engineering & System Safety* 264 (2025), 111302. [DOI](https://doi.org/10.1016/j.ress.2025.111302).

Licensed under the [MIT License](LICENSE). Original copyright and source notices are retained. [Other language](https://github.com/ZihanLiao-ZJU/Relaxed-Subset-Simulation-Python).
