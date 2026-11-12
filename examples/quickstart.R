# quickstart.R: SuperScreener on simulated data, end to end (a few minutes on a laptop).
# Run from the repository root:  Rscript examples/quickstart.R

library(superscreener)

## 1. Data: binary outcome, 10 true predictors, 2 correlated proxies of each, 80 noise columns
design <- make_design(p = 110, surface = "nonlinear")
train <- simulate_data(n = 400, design, seed = 1)
test <- simulate_data(n = 5000, design, seed = 2)

## 2. SuperScreener at the chapter's pre-specified price c = 1%.
##    `selection_library` is the learner fit to every candidate inside the cross-validation;
##    a small one keeps this example fast. The final model uses `sl_library`.
fit <- superscreener(train$X, train$Y,
                     K = c(5, 10, 15, 20, 25, 30),
                     price = 0.01,
                     selection_library = c("SL.glm", "SL.glmnet", "SL.mean"),
                     sl_library = sl_library_default(),
                     seed = 1)
fit
head(fit$candidates[order(fit$candidates$cv_risk), ], 10)   # best candidates by CV risk

## 3. What a higher or lower price would choose (no refitting)
price_path(fit, prices = c(0, 0.0025, 0.005, 0.01, 0.02, 0.04))

## 4. Predict new data, and check against the truth (simulation only)
p_hat <- predict(fit, test$X)
cat(sprintf("test AUC %.3f (true risk: %.3f)\n", auc(test$Y, p_hat), auc(test$Y, test$p)))
table(design$roles[fit$selected])                       # true / proxy / noise among the selected

## 5. Refit at another price without rerunning the cross-validation
fit4 <- set_price(fit, 0.04)
fit4$selected

## 6. Honest accuracy of the whole procedure (screener and K chosen inside each outer fold)
ncv <- nested_cv_superscreener(train$X, train$Y, prices = c(0, 0.01, 0.04),
                               K = c(5, 10, 20, 30), selection_library = c("SL.glm", "SL.glmnet", "SL.mean"),
                               sl_library = c("SL.glmnet", "SL.ranger_300", "SL.mean"),
                               V_outer = 5, seed = 1)
ncv
## variables selected in at least 80% of outer training folds at c = 1%
names(which(ncv$frequency[, "c=0.01"] >= 0.8))
