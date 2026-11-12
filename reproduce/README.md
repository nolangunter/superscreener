# Reproducing the simulation study

These scripts regenerate the simulation results of the chapter: Table 2 (fixed budget),
appendix Tables A1 to A3, and Figures 1 to 3. Run them from the repository root after
installing the package (or with `pkgload` installed, without installing).

| Step | Command | Output |
|---|---|---|
| 1. Smoke test (small designs, light learner; minutes) | `Rscript reproduce/run_simulation.R --quick` | `reproduce/output/sim_quick/` |
| 2. Replicates | `Rscript reproduce/run_simulation.R --reps 1000 --cores 50` | `reproduce/output/sim/reps/` (one `.rds` per replicate) |
| 3. Summaries and LaTeX tables | `Rscript reproduce/summarize_simulation.R` | `reproduce/output/sim/summary/` |
| 4. Figures | `Rscript reproduce/make_figures.R` | `reproduce/output/sim/summary/figures/` |

Step 2 can be resumed: replicates whose file exists are skipped, and a larger `--reps` adds
replicates. `--scenarios n500,small_n` runs a subset. Steps 3 and 4 need `dplyr`, `tidyr`,
`readr` and `ggplot2`; to summarize the smoke test, add `--outdir reproduce/output/sim_quick`.

## Scenarios

| Scenario | n | p | Proxy correlation | Surface |
|---|---|---|---|---|
| `n500` | 500 | 300 | U(0.5, 0.85) | base (nonlinear) |
| `prox_high` | 1000 | 300 | U(0.90, 0.97) | base (nonlinear) |
| `small_n` | 100 | 2000 | U(0.5, 0.85) | base (nonlinear) |
| `strong_nl` | 200 | 500 | U(0.5, 0.85) | strongly nonlinear |

Every scenario has 10 true covariates, 2 proxies per true covariate, p - 30 noise covariates,
prevalence 0.30 and a true-risk AUC of 0.80. See `?make_design`.

## Cost

Minutes per replicate on one core of the server used for the chapter: `n500` about 10,
`prox_high` about 40, `small_n` about 5, `strong_nl` about 4. The chapter's 1000 replicates of
all four scenarios took about 1000 core-hours. `--reps 50` gives stable means for most
quantities in about 50 core-hours.

## Exactness

Every replicate has a fixed seed, so a rerun on the same machine with the same package
versions gives identical results. Before release, these scripts were checked in two ways:

- `summarize_simulation.R`, run on the chapter's 4000 stored replicate files, reproduces the
  chapter's summary CSVs and appendix tables exactly.
- `run_simulation.R` reproduces the original research code exactly (every prediction, risk,
  choice and stability index) when both run on the same machine, after removing from the
  original the methods the chapter does not report (see below).

Results on another machine can differ in the random forest ranking, because `ranger`'s
forests depend on its version. The LRT, lasso and Spearman rankings and the self-sized
screeners do not. The checks above used R 4.4.2 with SuperLearner 2.0.29, glmnet 4.1.10,
ranger 0.17.0 and xgboost 1.7.8.1. Means over many replicates do not depend on these details.

## Differences from the research code

The research code also computed a few methods that the chapter does not report (a risk-weighted
consensus, fold-aggregated rankings and a partially linear ranking). They are left out here.
The Spearman ranking is called `spearman` here and `corrank` in the research code.

## Case study

The Cameroon Trauma Registry data are not public, so the case study cannot be rerun from this
repository. [`examples/apply_to_your_data.R`](../examples/apply_to_your_data.R) runs the same
analysis (price path, nested cvAUC with intervals, selection frequencies) on any dataset with a
binary outcome. The registry analysis differed from this template in two ways stated in the
chapter: it selected on the log-likelihood (`loss = "neg_loglik"`) and used a smaller inner
Super Learner.
