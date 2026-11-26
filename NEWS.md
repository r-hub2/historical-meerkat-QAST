# QAST 0.1.0

First release.

## Models

* `semi_aft_fg()` fits the bivariate semiparametric AFT model of Fine and
  Gelber (2001) for survival and quality-adjusted survival, with binary
  covariates and truncation at `L`. As in the article, the estimating
  function of `theta` sums over all pairs of subjects (the estimates do not
  depend on the order of the rows) and weights each subject by the
  Kaplan-Meier censoring survival of its censoring group. The new argument
  `cens_group` sets the censoring partition; the default does not
  stratify, as in the application of the article, and a warning flags
  groups with no censoring survival left at `L`. This differs from the
  code used in the thesis, which paired subjects by their position within
  each covariate group and estimated the censoring separately for each
  covariate.
* `semi_aft_prop()` fits the bivariate semiparametric AFT model proposed in
  the doctoral thesis of da Silva (2026, Institute of Mathematics and
  Statistics, University of São Paulo), with IPCW weights from a Cox model
  for the quality-adjusted component. `citation("QAST")` lists the thesis
  and Fine and Gelber (2001).
* Fits store the intercepts, which the rank-based estimating equations do
  not identify: the Kaplan-Meier mean of the survival time residuals
  (Ding and Nan, 2015) and the IPCW mean of the quality-adjusted residuals
  (Zhao and Tsiatis, 2000), in `fit$intercept`, shown by `print()` and
  `summary()`.
* `predict()` method for fits: linear predictors (`type = "lp"`),
  survival probabilities given the covariates (`type = "survival"`) and
  mean times restricted to `tau` (`type = "rmean"`, computed exactly on
  the step curves), from the Kaplan-Meier curve of the survival time
  residuals and the Zhao and Tsiatis (1997) estimator of the
  quality-adjusted residuals; accepts `newdata` for formula and matrix
  fits. With `se.fit = TRUE` (fits with `var_comp = TRUE`), standard
  errors (and, for the survival probabilities, pointwise confidence
  limits) from the jackknife or bootstrap replicates of the fit, with the
  intercept and the residual distribution re-estimated in each replicate.
* `qast_pp()`: P-P plots of the residuals by covariate group, the
  lack-of-fit diagnostic of Fine and Gelber (2001, Section 4), with the
  trimmed conditional quality-adjusted distribution for `semi_aft_fg()`
  fits; `print()` and `plot()` methods.
* `qast_metrics()`: Brier score over time and integrated Brier score
  (with those of the model without covariates as a reference) and IPCW
  mean absolute error of the restricted mean, for a fit or for a matrix of
  survival probabilities from any model. For the quality-adjusted time,
  whose censoring is induced informative, the weights use the censoring
  survival on the survival time scale, as the estimator of Zhao and
  Tsiatis (1997), instead of the censoring curve of the quality-adjusted
  time used by the usual tools. The total IPCW mass is reported as a
  positivity diagnostic.
* `surv_qt_fun()` implements the estimator of Zhao and Tsiatis (1997) for
  the survival curve of the quality-adjusted time, evaluated as
  P(Q > q) at each subject (the convention of `predict()` and of the
  Kaplan-Meier curve).

## Interface

* Both models accept a formula with the new `Qsurv(time, qtime, status)`
  response, e.g. `semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
  L = 5)`; factors enter as indicators, without an intercept. Vectors and a
  covariate matrix are still accepted.
* Fits return an object of class `qast_fit` with `print()`, `summary()`,
  `coef()`, `vcov()` and `confint()` methods. `vcov()` returns the full
  jackknife or bootstrap covariance matrix of each component.
* Quality per unit of survival: `coef(fit, "diff")` is `theta - beta`,
  whose exponential (with log links) is the change in quality of life per
  unit of survival and the ratio of mean quality-adjusted survival per
  unit of mean survival (Fine and Gelber, 2001). Its variance uses the
  covariance of beta and theta over the same replicates; `summary()` shows
  the table and the normal interval of `exp(theta - beta)`, and
  `confint(fit, component = "diff")` gives the normal or, for bootstrap
  fits, the percentile interval.
* Only the data (and `L` for `semi_aft_fg()`) are required; the defaults are
  `h1 = h2 = "log"`, `method = "Nelder-Mead"`, `norm_min = "L2"` and
  `var_comp = FALSE`. `M` is required only with `var_method = "bootstrap"`.
* Inputs are validated: `NA`, non-positive times, `qt > t`, `status` outside
  {0, 1}, mismatched lengths and invalid `L` or `M` stop with a message that
  points to the offending observations.
* Constant covariate columns (an intercept) are rejected, since the
  estimating equations are invariant to shifts; linearly dependent columns
  trigger a warning.
* With a single covariate, optimizers that do not support one coefficient
  or perform poorly with it are replaced by `"subplex"`, with a message.

## Computation

* The C++ kernels are exported with `rng = false`: they do not use the R
  random number generator, and the default `RNGScope` of the wrappers
  reloaded `.Random.seed` at every call, which restarted the random walk
  of `method = "SANN"` at each evaluation of the objective.
* `method = "bobyqa"` uses the BOBYQA of NLopt (package nloptr) instead of
  the one of minqa, whose RESCUE step evaluates the function at an
  uninitialized point (minqa 1.2.8, `src/rescue.f`).
* Estimating equations evaluated by C++ kernels (Rcpp/RcppArmadillo).
* In `semi_aft_fg()` with `phi = "piecewise"`, the objective is evaluated
  without storing the pairs: O(n) memory and O(n log^2 n) per evaluation,
  with the same estimates.
* In `semi_aft_fg()`, `theta` is obtained in closed form for
  `phi = "identity"` and `norm_min = "L2"` (`fast_path = TRUE`), and the
  pair marginals of the U-statistic are computed in O(N log^2 N).
* Jackknife and bootstrap variances run in parallel (fork) on the physical
  cores; set `options(QAST.cores = k)` to change the number of processes.
* Each jackknife or bootstrap replicate refits the whole estimator on the
  same subsample or resample: beta, the censoring weights and theta, with
  the beta of the replicate in the theta equation of `semi_aft_fg()` (which
  depends on beta through the trimming at `L`). Keeping the full-sample
  beta, as the code used in the thesis did, underestimates the standard
  errors of theta. The bootstrap draws one set of resamples for both
  components, as `boot::boot()`, which is no longer a dependency.
* Replicates whose covariates lose rank in the resample (e.g. a rare binary
  covariate that disappears) or whose optimizer does not converge are
  discarded, with a warning that counts them; warnings of the replicates
  are reported once per distinct message.
