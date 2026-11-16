# QAST 0.1.0

First release.

## Models

* `semi_aft_fg()` fits the bivariate semiparametric AFT model of Fine and
  Gelber (2001) for survival and quality-adjusted survival, with binary
  covariates and truncation at `L`.
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
* `surv_qt_fun()` implements the estimator of Zhao and Tsiatis (1997) for
  the survival curve of the quality-adjusted time.

## Interface

* Both models accept a formula with the new `Qsurv(time, qtime, status)`
  response, e.g. `semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
  L = 5)`; factors enter as indicators, without an intercept. Vectors and a
  covariate matrix are still accepted.
* Fits return an object of class `qast_fit` with `print()`, `summary()`,
  `coef()` and `vcov()` methods. `vcov()` returns the full jackknife or
  bootstrap covariance matrix of each component.
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
  are replaced by `"subplex"`, with a message.

## Computation

* Estimating equations evaluated by C++ kernels (Rcpp/RcppArmadillo).
* In `semi_aft_fg()`, `theta` is obtained in closed form for
  `phi = "identity"` and `norm_min = "L2"` (`fast_path = TRUE`), and the
  pair marginals of the U-statistic are computed in O(N log^2 N).
* Jackknife and bootstrap variances run in parallel (fork) on the physical
  cores; set `options(QAST.cores = k)` to change the number of processes.
