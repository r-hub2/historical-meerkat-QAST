#' Bivariate semiparametric AFT regression model (proposed model)
#'
#' Fits the package's proposed bivariate semiparametric accelerated failure
#' time model for the joint analysis of survival time and quality-adjusted
#' survival time, with inverse probability of censoring weights (IPCW) from
#' a Cox model for the quality-adjusted component.
#'
#' @section Reference: this is the model proposed in the doctoral thesis of
#'   da Silva (2026), Institute of Mathematics and Statistics, University of
#'   Sao Paulo (see References). It is distinct from the Fine and Gelber
#'   (2001) model, implemented in \code{\link{semi_aft_fg}}.
#'
#' @section Numerical kernel: the estimating equations are evaluated by the
#'   C++ kernels \code{U1_obj_cpp}/\code{U2_obj_cpp} in src/aft_core.cpp.
#'   The beta equation is the same in both models, so \code{U1_obj_cpp} also
#'   serves the beta stage of \code{semi_aft_fg}.
#'
#' @param t Observed survival times, or a formula
#'   \code{Qsurv(t, qt, status) ~ covariates} (see \code{\link{Qsurv}}); with
#'   a formula, \code{qt}, \code{status} and \code{z} come from it and from
#'   \code{data}
#' @param qt Quality-adjusted survival times (\code{0 <= qt <= t})
#' @param status Event indicator (1 = event, 0 = censored)
#' @param z Covariate matrix (n x p), WITHOUT an intercept. The estimating
#'   equations are invariant to shifts, so the model does not estimate an
#'   intercept; constant columns are rejected with an error
#' @param h1 Link function for \code{t} (default \code{"log"})
#' @param h2 Link function for \code{qt} (default \code{"log"})
#' @param method Optimizer: "Nelder-Mead" (default), "uobyqa", "mads",
#'   "subplex", "BB", "SANN", "bobyqa", "newuoa", "hjkb", "BFGS". With a
#'   single covariate, methods that do not support one coefficient
#'   (Nelder-Mead, uobyqa, mads, bobyqa, hjkb) are replaced by "subplex",
#'   with a message
#' @param norm_min Minimization norm, \code{"L2"} (default) or \code{"L1"}
#' @param cens_weight If TRUE (default), IPCW weights for informative
#'   censoring
#' @param M Number of bootstrap replicates; required only with
#'   \code{var_method = "bootstrap"}
#' @param var_comp If TRUE, estimates the variances (default FALSE)
#' @param var_method \code{"jackknife"} (default, delete-1) or
#'   \code{"bootstrap"} (resampling of subjects). Each replicate
#'   re-estimates beta, the Cox censoring weights and theta on the same subsample or
#'   resample, which also gives the joint covariance of the replicates.
#'   Replicates whose covariates lose rank (e.g. a binary covariate
#'   that became constant) or whose optimizer does not converge are
#'   discarded, with a warning.
#' @param data Data frame with the formula variables, when \code{t} is a
#'   formula
#' @return Object of class \code{"qast_fit"} (see \code{\link{qast_fit}}),
#'   with \code{beta}, \code{theta}, \code{intercept} and, if \code{var_comp = TRUE},
#'   \code{var_beta}, \code{var_theta}, \code{vcov_beta} and
#'   \code{vcov_theta}
#' @references da Silva, A. (2026). \emph{Regression analysis of
#'   quality-adjusted survival time in cancer patients admitted to intensive
#'   care units: a new estimation method}. Doctoral thesis, Institute of
#'   Mathematics and Statistics, University of Sao Paulo. Advisor: Gisela
#'   Tunes da Silva.
#' @seealso \code{\link{semi_aft_fg}} for the Fine and Gelber (2001) model;
#'   \code{\link{qast_fit}} for \code{print}, \code{summary}, \code{coef}
#'   and \code{vcov}
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(age = rnorm(n), trt = rbinom(n, 1, 0.5))
#' d$t <- rexp(n, rate = exp(-0.3 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#'
#' fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d)
#' fit
#' @export
#' @importFrom stats make.link lm runif optim coef vcov
#' @importFrom survival coxph Surv
#' @importFrom minqa uobyqa
#' @importFrom subplex subplex
#' @importFrom dfoptim mads nmk
#' @importFrom optimx optimr
#' @importFrom BB BBoptim
semi_aft_prop <- function(t, qt, status, z, h1 = "log", h2 = "log",
                          method = "Nelder-Mead", norm_min = "L2",
                          cens_weight = TRUE, M = NULL, var_comp = FALSE,
                          var_method = "jackknife", data = NULL) {

  cl <- match.call()
  formula_info <- NULL

  # --- Formula --------------------------------------------------------------
  if (inherits(t, "formula")) {
    names(cl)[names(cl) == "t"] <- ""   # Call: shows the formula without "t ="
    md <- .qast_model_data(t, data, parent.frame())
    t <- md$t; qt <- md$qt; status <- md$status; z <- md$z
    formula_info <- md$info
  }

  # --- Validation -----------------------------------------------------------
  method <- match.arg(method,
                      c("uobyqa", "mads", "Nelder-Mead", "subplex", "BB", "SANN",
                        "bobyqa", "newuoa", "hjkb", "BFGS"))
  norm_min <- match.arg(norm_min, c("L1", "L2"))
  var_method <- match.arg(var_method, c("jackknife", "bootstrap"))
  .validate_var_args(var_comp, var_method, M)
  if (var_comp && !cens_weight)
    stop("var_comp = TRUE with cens_weight = FALSE is not supported: ",
         "the variance of theta requires the IPCW weights.", call. = FALSE)

  h1_name <- h1; h2_name <- h2
  h1 <- make.link(h1)
  h2 <- make.link(h2)
  z <- as.matrix(z)
  .validate_qast_data(t, qt, status, z, h1, h2)
  method <- .qast_method_for_p(method, ncol(z))

  # --- Precomputation (constants in the optimization) ----------------------
  lt <- h1$linkfun(t)
  lq <- h2$linkfun(qt)
  n <- length(t)
  if (is.null(colnames(z))) colnames(z) <- paste0("x", seq_len(ncol(z)))
  # The model does not estimate an intercept: constant columns are not
  # identifiable (equations invariant to shifts)
  validate_no_constant_z(z)
  warn_collinear_z(z)          # warning, not error
  # As double once: with integers RcppArmadillo would convert z and status
  # at every evaluation of the kernels.
  storage.mode(z) <- "double"
  status <- as.numeric(status)

  # --- Helpers ---------------------------------------------------------------
  run_optimizer <- function(par, fn) {
    switch(method,
           "uobyqa"       = uobyqa(par = par, fn = fn),
           "mads"         = mads(par = par, fn = fn, control = list(trace = FALSE)),
           "Nelder-Mead"  = nmk(par = par, fn = fn),
           "subplex"      = subplex(par = par, fn = fn),
           "BB"           = BBoptim(par = par, fn = fn,
                                    control = list(trace = FALSE), quiet = TRUE),
           "SANN"         = optim(par = par, fn = fn, method = "SANN"),
           "bobyqa" = ,
           "newuoa" = ,
           "hjkb"   = ,
           "BFGS"   = optimr(par = par, fn = fn, method = method,
                             control = list(trace = FALSE)),
           stop("unknown method: ", method)
    )
  }

  # Convergence code: 0 = converged
  extract_conv <- function(res) {
    if (method == "uobyqa") return(res$ierr)
    cc <- res$convergence
    if (method == "mads") cc <- if (cc < 1e-06) 0 else -1
    cc
  }

  # Starting coefficients by lm; NA (collinearity) -> runif
  init_coef <- function(y, zmat) {
    b <- tryCatch(lm(y ~ zmat)$coef[-1], error = function(e) rep(0, ncol(zmat)))
    if (anyNA(b)) b[is.na(b)] <- runif(sum(is.na(b)))
    b
  }

  # ---------------------------------------------------------------------------
  # Stage 1: estimation of beta
  # ---------------------------------------------------------------------------
  U1min <- run_optimizer(init_coef(lt, z),
                         function(beta) U1_obj_cpp(beta, lt, z, status, norm_min))
  stopifnot("The beta estimation algorithm did not converge" =
              extract_conv(U1min) == 0)
  betahat <- U1min$par

  # ---------------------------------------------------------------------------
  # Stage 2: estimation of theta
  # ---------------------------------------------------------------------------
  if (!cens_weight) {
    # Without IPCW weights, theta uses the UNweighted equation (U1_obj_cpp) on lq
    U2min <- run_optimizer(init_coef(lq, z),
                           function(theta) U1_obj_cpp(theta, lq, z, status,
                                                      norm_min))
  } else {
    # IPCW weights from a Cox model (survival::coxph)
    HG_surv <- HG(X = z, t = t, status = status)
    U2min <- run_optimizer(init_coef(lq, z),
                           function(theta) U2_obj_cpp(theta, lq, z, status,
                                                      HG_surv, norm_min))
  }
  stopifnot("The theta estimation algorithm did not converge" =
              extract_conv(U2min) == 0)
  thetahat <- U2min$par

  # ---------------------------------------------------------------------------
  # Variance (var_comp = TRUE, which requires cens_weight = TRUE)
  #
  # Each replicate refits beta and theta on the same subsample (jackknife)
  # or resample (bootstrap), with the Cox censoring weights recomputed on
  # it. A replicate whose covariates lose rank, or whose optimizer does not
  # converge, gives NA and is discarded (with a warning).
  # ---------------------------------------------------------------------------
  vcov_beta <- vcov_theta <- NULL
  replicates <- list(beta = NULL, theta = NULL)
  if (var_comp) {
    p <- ncol(z)
    rank_z <- qr(z)$rank
    na_p <- rep(NA_real_, p)

    refit <- function(rows) {
      z_s <- z[rows, , drop = FALSE]
      if (.qast_rank_deficient(z_s, rank_z)) return(c(na_p, na_p))
      t_s <- t[rows]; status_s <- status[rows]
      lt_s <- lt[rows]; lq_s <- lq[rows]
      rb <- run_optimizer(betahat, function(b)
        U1_obj_cpp(b, lt_s, z_s, status_s, norm_min))
      b <- if (extract_conv(rb) == 0) rb$par else na_p
      HG_s <- HG(X = z_s, t = t_s, status = status_s)
      rt <- run_optimizer(thetahat, function(th)
        U2_obj_cpp(th, lq_s, z_s, status_s, HG_s, norm_min))
      th <- if (extract_conv(rt) == 0) rt$par else na_p
      c(b, th)
    }
    blocks <- list(beta = seq_len(p), theta = p + seq_len(p))

    if (var_method == "jackknife") {
      vr <- .jackknife_se(n, refit = function(i) refit(-i), blocks = blocks,
                          cov = TRUE)
      .qast_warn_discarded(vr, n, "jackknife")
      replicates <- lapply(vr, .qast_jack_reps)
    } else {
      vr <- .qast_bootstrap(n, M, refit, blocks)
      .qast_warn_discarded(vr, M, "bootstrap")
      replicates <- lapply(vr, .qast_boot_reps)
    }
    vcov_beta  <- vr$beta$cov
    vcov_theta <- vr$theta$cov
  }

  .new_qast_fit(
    model = "prop", call = cl, z = z, status = status,
    betahat = betahat, thetahat = thetahat,
    vcov_beta = vcov_beta, vcov_theta = vcov_theta,
    var_comp = var_comp, var_method = var_method, M = M,
    lt = lt, lq = lq, t = t, qt = qt, formula_info = formula_info,
    replicates = replicates,
    settings = list(h1 = h1_name, h2 = h2_name, method = method,
                    norm_min = norm_min, cens_weight = cens_weight)
  )
}
