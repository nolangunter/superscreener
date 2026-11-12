# superscreener

SuperScreener chooses both the **screening rule** and the **number of variables** for a
binary-outcome prediction model. It adds a price for each variable to the Super Learner's
cross-validation selector, so a larger set is chosen only when it lowers the cross-validated
risk by enough to pay for the extra variables. Changing the price shows how much accuracy is
given up as fewer variables are kept, which matters when every variable has to be collected.

This repository has the R package and the code for the simulation study in

> Gunter N, Hubbard A. *Choosing how many variables to collect: cross-validated screening
> with an explicit price per variable.* Chapter 5 of N. Gunter, *Algorithmic Choices in
> Targeted Learning*, PhD dissertation, University of California, Berkeley, 2026.

## Installation

```r
# install.packages("remotes")
remotes::install_github("nolangunter/superscreener")
```

Requires R >= 4.1 and the packages `SuperLearner`, `glmnet`, `ranger` and `xgboost`.

## Quick start

```r
library(superscreener)

# simulated data: 10 true predictors, 2 correlated proxies of each, 80 noise columns
design <- make_design(p = 110)
train  <- simulate_data(n = 400, design, seed = 1)

# a light learner inside the cross-validation keeps this to a few minutes;
# the final model uses the default Super Learner (random forest, lasso, boosting, mean)
fit <- superscreener(train$X, train$Y, price = 0.01,
                     selection_library = c("SL.glm", "SL.glmnet", "SL.mean"), seed = 1)
fit                  # chosen screener, number of variables, selected variables
price_path(fit)      # what each price would choose (no refitting)
predict(fit, newdata = train$X[1:5, ])

# honest accuracy of the whole procedure, including the choice of screener and size
ncv <- nested_cv_superscreener(train$X, train$Y, prices = c(0, 0.01, 0.04),
                               selection_library = c("SL.glm", "SL.glmnet", "SL.mean"),
                               sl_library = c("SL.glm", "SL.glmnet", "SL.mean"), cores = 4, seed = 1)
ncv
```

With the default library everywhere, each fit refits 34 candidates in 5 folds with a random
forest and boosting, which takes much longer (see *Speed* below).

[`examples/quickstart.R`](examples/quickstart.R) runs these steps with comments (about 10
minutes of CPU time). [`examples/apply_to_your_data.R`](examples/apply_to_your_data.R) is a
command-line template for your own CSV file.

## How it works

A **candidate** is a screener applied to the data, followed by the same Super Learner. The
default library has 34 candidates:

| Screener | Type | Description |
|---|---|---|
| `lrt` | ranking | univariate logistic likelihood-ratio statistic |
| `lasso` | ranking | order of entry along the lasso path |
| `rf` | ranking | random forest permutation importance (`ranger`) |
| `consensus` | ranking | mean of the `lrt`, `lasso` and `rf` ranks |
| `spearman` | ranking | absolute Spearman correlation (`screen.corRank`) |
| `sl_glmnet` | self-sized | `SuperLearner::screen.glmnet` |
| `sl_corP` | self-sized | `SuperLearner::screen.corP` |
| `lasso_1se` | self-sized | nonzero lasso coefficients at `lambda.1se` |
| `stabsel` | self-sized | stability selection, 25 complementary pairs |

Each ranking is used at budgets `K = 5, 10, ..., 30` columns. Self-sized screeners choose
their own number of columns.

`superscreener()` then runs six steps:

1. Split the data into `V` folds, stratified by outcome.
2. In each training sample, run every screener and fit the Super Learner to each candidate's set.
3. Predict the held-out fold. The cross-validated risk R(j) of candidate j is the mean loss
   (Brier score by default) of its held-out predictions. R(0) is the risk of the
   intercept-only model.
4. Run every screener on the whole sample. s(j) is the number of variables candidate j selects.
5. Set the price per variable to `lambda = c * (R(0) - min_j R(j))` and choose the candidate
   with the smallest `R(j) + lambda * s(j)`. Ties go to the smaller set.
6. Fit the Super Learner to the chosen set on the whole sample.

The relative price `c` (`price`) is the share of the best candidate's risk reduction that one
more variable has to buy. With `c = 0.01`, a variable is kept only if it lowers the
cross-validated risk by 1% of what the best candidate achieves over the intercept-only model.
`c = 0` is the ordinary discrete Super Learner over the candidates.

## Main functions

| Function | Purpose |
|---|---|
| `superscreener()` | fit SuperScreener (steps 1 to 6) |
| `price_path()` | choice at each price, from the stored risks |
| `set_price()` | refit the final model at another price |
| `predict()` | predictions from the final model |
| `nested_cv_superscreener()` | nested cvAUC with 95% influence-curve intervals, difference from the full model, selection frequencies and stability over outer folds |
| `rank_columns()`, `self_sized_sets()` | the screeners on their own |
| `choose_candidate()` | step 5 for any table of candidates |
| `expand_factors()` | categorical variables to indicator columns, plus the column-to-variable map |
| `cv_auc()`, `brier()`, `auc()` | metrics |
| `selection_frequency()`, `nogueira_index()` | stability of selected sets |
| `make_design()`, `simulate_data()`, `chapter_scenarios()` | the simulation design of the chapter |

## Practical notes

**Variables and columns.** The price is charged per *variable*. A categorical variable coded
as several indicator columns costs one unit:

```r
ex  <- expand_factors(my_data_frame)          # numeric columns + map from column to variable
fit <- superscreener(ex$X, y, groups = ex$groups)
```

Variables can also have different costs (`cost = c(age = 1, lab_test = 5, ...)`). The
chapter uses equal costs.

**Report the nested cvAUC.** The cross-validated AUC that `superscreener()` stores for the
chosen candidate is optimistic, because the same folds chose it. `nested_cv_superscreener()`
reruns the whole procedure inside each outer training sample.

**Weak signal.** When no candidate beats the intercept-only model by much (in the chapter's
simulation, at n = 100 to 200 with a weak or strongly nonlinear signal), the priced rule can
choose very small or empty sets. `min_size = 5` drops candidates with fewer than five
variables. In the chapter this floor was added after seeing the results, so it is reported
there as a post hoc rule.

**Speed.** Every candidate is refit in every fold, so cost grows with the number of
candidates and the size of the Super Learner. A smaller `selection_library` (for example
`c("SL.glm", "SL.glmnet", "SL.mean")`) speeds up the selection while the final model still uses
`sl_library`. `nested_cv_superscreener()` repeats the whole procedure in every outer fold,
so it costs about `V_outer` times one fit. `cores =` runs outer folds in parallel
(any operating system; the package must be installed).

**New data with categorical variables.** `predict()` needs the same columns the model was
fit on. If you use `expand_factors()`, encode the training and new data together (stack
them, expand, then split) so both get the same indicator columns.

**Parallel runs and your own learners.** With `cores > 1`, learner wrappers you defined in
your session (e.g. `SL.mine`) are copied to the workers. Anything they depend on must be
available there too.

**Outcome.** Only binary outcomes are supported (0/1, or TRUE/FALSE, or a two-level factor
with levels "0" and "1"). Missing values must be imputed first. Each class needs at least as
many observations as there are folds.

## Reproducing the simulation study

See [`reproduce/README.md`](reproduce/README.md). The registry data of the case study are
not public and are not in this repository.

## Tests

```r
devtools::test()
```

## License

MIT. See [LICENSE.md](LICENSE.md).
