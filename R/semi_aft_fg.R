#' Bivariate semiparametric AFT regression model (Fine and Gelber, 2001)
#'
#' Fits the bivariate semiparametric accelerated failure time model of Fine
#' and Gelber (2001) for the joint analysis of survival time and
#' quality-adjusted survival time. The estimating function of the
#' quality-adjusted coefficients is a U-statistic over all pairs of
#' subjects, truncated at \code{L} and weighted by the inverse of the
#' Kaplan-Meier censoring survival in the censoring group of each subject.
#'
#' @param t Observed survival times, or a formula
#'   \code{Qsurv(t, qt, status) ~ covariates} (see \code{\link{Qsurv}}); with
#'   a formula, \code{qt}, \code{status} and \code{z} come from it and from
#'   \code{data}, and factors enter as indicators (without an intercept)
#' @param qt Quality-adjusted survival times (\code{0 <= qt <= t})
#' @param status Event indicator (1 = event, 0 = censored)
#' @param z Binary covariate matrix (n x p), WITHOUT an intercept. Binary
#'   covariates are central to this model; multi-level factors must be
#'   coded as indicators (the formula interface does it). Constant columns
#'   are rejected, since the model does not estimate an intercept
#' @param h1 Link function for \code{t} (default \code{"log"})
#' @param h2 Link function for \code{qt} (default \code{"log"})
#' @param method Optimizer: "Nelder-Mead" (default), "uobyqa", "mads",
#'   "subplex", "BB", "SANN", "bobyqa", "newuoa", "hjkb", "BFGS". With a
#'   single covariate, methods that do not support one coefficient
#'   (Nelder-Mead, uobyqa, mads, bobyqa, hjkb) are replaced by "subplex",
#'   with a message
#' @param norm_min Minimization norm, \code{"L2"} (default) or \code{"L1"}
#' @param L Truncation time for the quality-adjusted parameters (required,
#'   positive)
#' @param phi Odd function: \code{"identity"} (default) or
#'   \code{"piecewise"}
#' @param M Number of bootstrap replicates; required only with
#'   \code{var_method = "bootstrap"}
#' @param var_comp If TRUE, estimates the variances (default FALSE)
#' @param var_method \code{"jackknife"} (default, delete-1) or
#'   \code{"bootstrap"} (resampling of subjects). Each replicate
#'   re-estimates beta, the Kaplan-Meier censoring weights and theta (with the beta of
#'   the replicate, since the theta equation depends on beta) on the same subsample or
#'   resample, which also gives the joint covariance of the replicates.
#'   Replicates whose covariates lose rank (e.g. a binary covariate
#'   that became constant) or whose optimizer does not converge are
#'   discarded, with a warning.
#' @param fast_path If TRUE (default) and \code{phi = "identity"},
#'   \code{norm_min = "L2"}, theta is estimated in closed form: least
#'   squares on the affine representation, which is the EXACT minimum of the
#'   L2 objective (the optimizer only approximates it). In that case
#'   \code{method} is used for beta only; theta does not go through the
#'   optimizer, also in the variance replicates. If the system is
#'   rank-deficient or the residual gradient is large: a warning and a
#'   fallback to the optimizer. With other phi/norm_min combinations, or
#'   with \code{fast_path = FALSE}, theta is always estimated by
#'   \code{method}
#' @param data Data frame with the formula variables, when \code{t} is a
#'   formula
#' @param cens_group Censoring partition. The model assumes that the
#'   censoring depends on the covariates only through a discrete grouping
#'   G chosen by the analyst, and weights each subject by the Kaplan-Meier
#'   censoring survival of its group, \eqn{\hat H_G(X)}{H_G(X)}.
#'   \code{NULL} (default) does not stratify, as in the application of
#'   Fine and Gelber (2001). Otherwise, a vector with one value per
#'   observation defining the groups, e.g. a covariate associated with the
#'   censoring, or a discretized one (Fine and Gelber stratified by age
#'   dichotomized at the median as a sensitivity analysis). To choose it,
#'   one can check the effect of the covariates on the censoring, as in the
#'   article, e.g. with \code{survival::coxph(Surv(t, 1 - status) ~ z)}.
#'   \code{L} must satisfy \eqn{\hat H_G(L) > 0}{H_G(L) > 0} in every
#'   group (a warning flags the groups where it does not)
#'
#' @return Object of class \code{"qast_fit"} (see \code{\link{qast_fit}}),
#'   with \code{beta}, \code{theta}, \code{intercept} and, if
#'   \code{var_comp = TRUE},
#'   \code{var_beta}, \code{var_theta}, \code{vcov_beta} and
#'   \code{vcov_theta}
#' @references Fine, J. P. and Gelber, R. D. (2001). Joint regression
#'   analysis of survival and quality-adjusted survival. \emph{Biometrics},
#'   57(2), 376-382.
#' @seealso \code{\link{semi_aft_prop}}; \code{\link{qast_fit}} for
#'   \code{print}, \code{summary}, \code{coef} and \code{vcov}
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
#' d$t <- rexp(n, rate = exp(-0.4 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#'
#' fit <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                    L = quantile(d$t, 0.9))
#' fit
#'
#' \donttest{
#' # With jackknife standard errors
#' summary(semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                     L = quantile(d$t, 0.9), var_comp = TRUE))
#' }
#' @export
#' @importFrom stats make.link lm runif optim coef vcov confint
#' @importFrom minqa uobyqa
#' @importFrom subplex subplex
#' @importFrom dfoptim mads nmk
#' @importFrom optimx optimr
#' @importFrom BB BBoptim
semi_aft_fg <- function(t, qt, status, z, h1 = "log", h2 = "log",
                        method = "Nelder-Mead", norm_min = "L2",
                        L, phi = "identity", M = NULL, var_comp = FALSE,
                        var_method = "jackknife", fast_path = TRUE,
                        data = NULL, cens_group = NULL) {
  
  cl <- match.call()
  formula_info <- NULL
  
  # --- Formula --------------------------------------------------------------
  if (inherits(t, "formula")) {
    names(cl)[names(cl) == "t"] <- ""   # Call: shows the formula without "t ="
    md <- .qast_model_data(t, data, parent.frame())
    t <- md$t; qt <- md$qt; status <- md$status; z <- md$z
    formula_info <- md$info
  }
  if (missing(L)) stop("L (truncation time) is required", call. = FALSE)
  L <- unname(L)
  
  # --- Validation -----------------------------------------------------------
  method <- match.arg(method,
                      c("uobyqa", "mads", "Nelder-Mead", "subplex", "BB", "SANN",
                        "bobyqa", "newuoa", "hjkb", "BFGS"))
  phi <- match.arg(phi, c("identity", "piecewise"))
  norm_min <- match.arg(norm_min, c("L1", "L2"))
  var_method <- match.arg(var_method, c("jackknife", "bootstrap"))
  .validate_var_args(var_comp, var_method, M)
  
  h1_name <- h1; h2_name <- h2
  h1 <- make.link(h1)
  h2 <- make.link(h2)
  .validate_qast_data(t, qt, status, z, h1, h2, L = L)
  method <- .qast_method_for_p(method, NCOL(z))
  
  # Discrete covariates: central assumption of the Fine and Gelber model
  z <- validate_discrete_z(z)
  validate_no_constant_z(z)   # no intercept / empty group
  z_full_rank <- warn_collinear_z(z)   # warning, not error
  # As double once: with integers RcppArmadillo would convert z and status
  # at every evaluation of the kernels.
  storage.mode(z) <- "double"
  status <- as.numeric(status)

  # Censoring partition (Fine and Gelber, 2001): no stratification by default
  if (!is.null(cens_group)) {
    if (length(cens_group) != length(t))
      stop(sprintf("cens_group must have one value per observation (%d)",
                   length(t)), call. = FALSE)
    if (anyNA(cens_group))
      .qast_stop_at("cens_group contains NA", which(is.na(cens_group)))
  }
  cgroup <- fg_cens_partition(z, cens_group)
  .fg_check_cens_at_L(t, status, cgroup, L)
  
  # --- Precomputation (constants in the optimization) -----------------------
  lt <- h1$linkfun(t)
  lq <- h2$linkfun(qt)
  hL <- h1$linkfun(L)
  n <- length(t)
  
  # --- Helpers ---------------------------------------------------------------
  run_optimizer <- function(par, fn) {
    switch(method,
           "uobyqa"       = uobyqa(par = par, fn = fn),
           "mads"         = mads(par = par, fn = fn,
                                 control = list(trace = FALSE)),
           "Nelder-Mead"  = nmk(par = par, fn = fn),
           "subplex"      = subplex(par = par, fn = fn),
           "BB"           = BBoptim(par = par, fn = fn,
                                    control = list(trace = FALSE),
                                    quiet = TRUE),
           "SANN"         = optim(par = par, fn = fn, method = "SANN"),
           "bobyqa" = ,
           "newuoa" = ,
           "hjkb"   = ,
           "BFGS"   = optimr(par = par, fn = fn, method = method,
                             control = list(trace = FALSE)),
           stop("unknown method: ", method)
    )
  }
  
  extract_conv <- function(res) {
    if (method == "uobyqa") return(res$ierr)
    cc <- res$convergence
    if (method == "mads") cc <- if (cc < 1e-06) 0 else -1
    cc
  }
  
  init_coef <- function(y, zmat) {
    b <- tryCatch(lm(y ~ zmat)$coef[-1],
                  error = function(e) rep(0, ncol(zmat)))
    if (anyNA(b)) b[is.na(b)] <- runif(sum(is.na(b)))
    b
  }
  
  # ---------------------------------------------------------------------------
  # Stage 1: estimation of beta (same equation as the proposed model)
  # ---------------------------------------------------------------------------
  U1min <- run_optimizer(init_coef(lt, z),
                         function(beta) U1_obj_cpp(beta, lt, z,
                                                   status, norm_min))
  stopifnot("The beta estimation algorithm did not converge" =
              extract_conv(U1min) == 0)
  betahat <- U1min$par
  
  # ---------------------------------------------------------------------------
  # Stage 2: estimation of theta
  #
  #   stage2_env(): builds the objective of a (sub)sample
  #     - identity: affine form alpha + W theta, O(p^2) per evaluation
  #     - piecewise: sum over the valid pairs without storing them, O(n log^2 n)
  #   solve_theta(): closed form (identity + L2 + fast_path) or optimizer
  # ---------------------------------------------------------------------------
  km_surv <- prep_km_by_column(z, t, status, group = cgroup)
  
  stage2_env <- function(lt_s, lq_s, z_s, status_s, km_s, beta = betahat) {
    if (phi == "identity") {
      pairs <- make_fg_pairs(lt_s, lq_s, z_s, status_s,
                             beta, hL, km_s, phi)
      aff <- make_fg_affine(pairs, lq_s, z_s)
      list(
        fn = function(th) eval_fg_obj_affine(aff, th, norm_min),
        affine = aff,
        can_fast = TRUE
      )
    } else {
      # Without explicit pairs: active members and a divide-and-conquer
      # sum, O(n) memory and O(n log^2 n) per evaluation. Weights outside
      # the safe range (see fg_pairs_core.cpp): explicit pairs.
      lean <- make_fg_pw_lean_cpp(lt_s, lq_s, z_s, status_s, beta, hL,
                                  km_s)
      fn <- if (isTRUE(attr(lean, "meta")$ok)) {
        function(th) eval_fg_obj_pw_lean_cpp(lean, lq_s, z_s, th, norm_min)
      } else {
        pairs <- make_fg_pairs(lt_s, lq_s, z_s, status_s,
                               beta, hL, km_s, phi)
        function(th) eval_fg_obj(pairs, lq_s, z_s, th, "piecewise", norm_min)
      }
      list(fn = fn, affine = NULL, can_fast = FALSE)
    }
  }
  
  # Closed form when possible; rank deficiency or a large residual
  # gradient -> warning and the optimizer on the same affine objective.
  solve_theta <- function(s2, par0, must_converge = TRUE, replicate = FALSE) {
    if (fast_path && s2$can_fast && norm_min == "L2") {
      rcond <- 1e-10
      sol <- solve_fg_theta_affine(s2$affine, rcond = rcond)
      if (isTRUE(sol$full_rank) && is.finite(sol$grad_norm) &&
          sol$grad_norm < 1e-6) {
        return(list(par = sol$par, converged = TRUE))
      }
      # If the whole z was already flagged as collinear, this warning repeats
      # the same cause; it is issued only when the rank loss is new (e.g. a
      # degenerate jackknife/bootstrap replicate).
      # In the variance replicates, a generic message (without the numbers
      # of each replicate) so that it is reported once
      if (replicate) {
        if (z_full_rank) warning("fast_path: rank-deficient system or large residual ",
                "gradient in variance replicates; theta solved there by ",
                "the optimizer on the same affine objective")
      } else if (z_full_rank) warning(
        "fast_path: rank-deficient system or large residual gradient ",
        sprintf("(rcond = %.2g, full_rank = %s, grad_norm = %.2g). ",
                rcond, isTRUE(sol$full_rank), sol$grad_norm),
        "falling back to the optimizer on the same affine objective; ",
        "theta is not identifiable along the null directions."
      )
    }
    U2min <- run_optimizer(par0, s2$fn)
    converged <- extract_conv(U2min) == 0
    if (must_converge)
      stopifnot("The theta estimation algorithm did not converge" = converged)
    list(par = U2min$par, converged = converged)
  }
  
  stage2 <- stage2_env(lt, lq, z, status, km_surv)
  thetahat <- solve_theta(stage2, init_coef(lq, z))$par
  
  # ---------------------------------------------------------------------------
  # Variance (var_comp = TRUE)
  #
  # Each replicate re-estimates the whole estimator on the same subsample
  # (jackknife) or resample (bootstrap): beta, then the KM censoring
  # weights, the pairs and theta with the beta of the replicate (the
  # theta equation depends on beta through the trimming at h1(L) and the
  # residuals of the survival time). A replicate whose covariates lose
  # rank, or whose optimizer does not converge, gives NA and is discarded
  # (with a warning).
  # ---------------------------------------------------------------------------
  vcov_beta <- vcov_theta <- vcov_diff <- NULL
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
      if (extract_conv(rb) != 0) return(c(na_p, na_p))
      km_s <- prep_km_by_column(z_s, t_s, status_s,
                                group = droplevels(cgroup[rows]))
      s2 <- stage2_env(lt_s, lq_s, z_s, status_s, km_s, beta = rb$par)
      sol <- solve_theta(s2, thetahat, must_converge = FALSE,
                         replicate = TRUE)
      c(rb$par, if (sol$converged) sol$par else na_p)
    }
    # both: joint block, for the covariance of theta - beta
    blocks <- list(beta = seq_len(p), theta = p + seq_len(p),
                   both = seq_len(2L * p))

    if (var_method == "jackknife") {
      vr <- .jackknife_se(n, refit = function(i) refit(-i), blocks = blocks,
                          cov = TRUE)
      .qast_warn_discarded(vr, n, "jackknife")
      replicates <- lapply(vr[c("beta", "theta")], .qast_jack_reps)
    } else {
      vr <- .qast_bootstrap(n, M, refit, blocks)
      .qast_warn_discarded(vr, M, "bootstrap")
      replicates <- lapply(vr[c("beta", "theta")], .qast_boot_reps)
    }
    vcov_beta  <- vr$beta$cov
    vcov_theta <- vr$theta$cov
    dr <- .qast_diff_reps(vr$both, p, var_method)
    replicates$diff <- dr$reps
    vcov_diff <- dr$cov
  }

  .new_qast_fit(
    model = "fg", call = cl, z = z, status = status,
    betahat = betahat, thetahat = thetahat,
    vcov_beta = vcov_beta, vcov_theta = vcov_theta, vcov_diff = vcov_diff,
    var_comp = var_comp, var_method = var_method, M = M,
    lt = lt, lq = lq, t = t, qt = qt, formula_info = formula_info,
    replicates = replicates,
    settings = list(h1 = h1_name, h2 = h2_name, method = method,
                    norm_min = norm_min, L = L, phi = phi,
                    fast_path = fast_path, cens_group = cgroup,
                    cens_group_default = is.null(cens_group))
  )
}

# H_G(L) > 0 in every censoring group (Fine and Gelber, 2001)
.fg_check_cens_at_L <- function(t, status, group, L) {
  bad <- character()
  for (g in levels(group)) {
    i <- group == g
    if (!any(status[i] == 0)) next
    f <- survival::survfit(survival::Surv(t[i], 1 - status[i]) ~ 1)
    k <- findInterval(L, f$time)
    if (k > 0L && f$surv[k] <= 0) bad <- c(bad, g)
  }
  if (length(bad))
    warning(sprintf(paste0("the censoring survival is 0 at L in %d censoring ",
                           "group(s); L must satisfy H_G(L) > 0 (reduce L ",
                           "or a coarser cens_group)"), length(bad)),
            call. = FALSE)
  invisible(bad)
}

