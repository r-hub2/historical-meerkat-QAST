# R/qast_fit.R
#
# Common interface of the models: Qsurv() response for formulas, input
# validation and the S3 class "qast_fit" (print, summary, coef, vcov).

# ---------------------------------------------------------------------------
# Response for formulas
# ---------------------------------------------------------------------------

#' Response for QAST model formulas
#'
#' Bundles the survival time, the quality-adjusted time and the event
#' indicator into a response for the left-hand side of a formula, like
#' \code{survival::Surv()}:
#' \code{semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d, L = 5)}.
#'
#' @param time Observed survival time
#' @param qtime Quality-adjusted survival time
#' @param status Event indicator (1 = event, 0 = censored)
#' @return Matrix with columns \code{time}, \code{qtime} and \code{status}
#'   and class \code{"Qsurv"}
#' @seealso \code{\link{semi_aft_fg}}, \code{\link{semi_aft_prop}}
#' @examples
#' Qsurv(c(2, 5, 3), c(1.5, 4, 3), c(1, 0, 1))
#' @export
Qsurv <- function(time, qtime, status) {
  if (length(time) != length(qtime) || length(time) != length(status))
    stop("Qsurv: time, qtime and status must have the same length",
         call. = FALSE)
  out <- cbind(time = as.numeric(time), qtime = as.numeric(qtime),
               status = as.numeric(status))
  class(out) <- c("Qsurv", class(out))
  out
}

# Converts formula + data into the model inputs. The covariate matrix is
# the model.matrix without the intercept column (the model does not
# estimate it); factors enter as k - 1 indicators. NAs go on and are
# reported by the validation.
.qast_model_data <- function(formula, data, env) {
  if (is.null(data)) data <- env
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!inherits(y, "Qsurv"))
    stop("the left-hand side of the formula must be Qsurv(time, qtime, status)",
         call. = FALSE)
  mm <- stats::model.matrix(attr(mf, "terms"), mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) == 0L)
    stop("the formula has no covariates", call. = FALSE)
  info <- list(terms = stats::delete.response(attr(mf, "terms")),
               xlevels = stats::.getXlevels(attr(mf, "terms"), mf),
               contrasts = attr(mm, "contrasts"))
  attr(mm, "assign") <- NULL
  attr(mm, "contrasts") <- NULL
  list(t = unname(y[, "time"]), qt = unname(y[, "qtime"]),
       status = unname(y[, "status"]), z = mm, info = info)
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

# Error with the first offending positions
.qast_stop_at <- function(what, idx) {
  shown <- utils::head(idx, 5L)
  more <- if (length(idx) > 5L) sprintf(" (and %d more)", length(idx) - 5L) else ""
  stop(what, " at observations ", paste(shown, collapse = ", "), more, ".",
       call. = FALSE)
}

# Validation shared by semi_aft_fg and semi_aft_prop. h1/h2 are make.link
# objects; L is optional (only the Fine and Gelber model uses it).
.validate_qast_data <- function(t, qt, status, z, h1, h2, L = NULL) {
  if (!is.numeric(t) || !is.numeric(qt))
    stop("t and qt must be numeric", call. = FALSE)
  z <- as.matrix(z)
  n <- length(t)
  lens <- c(t = n, qt = length(qt), status = length(status), z = nrow(z))
  if (any(lens != n))
    stop("t, qt, status and z must have the same number of observations (",
         paste(names(lens), lens, sep = " = ", collapse = ", "), ")",
         call. = FALSE)
  if (n < 2L) stop("at least 2 observations are required", call. = FALSE)

  for (nm in c("t", "qt", "status")) {
    v <- get(nm)
    if (anyNA(v)) .qast_stop_at(sprintf("%s contains NA", nm), which(is.na(v)))
  }
  if (anyNA(z)) .qast_stop_at("z contains NA", which(!stats::complete.cases(z)))
  if (!is.numeric(z)) stop("z must be numeric", call. = FALSE)

  if (!all(status %in% c(0, 1)))
    .qast_stop_at("status must be 0 (censored) or 1 (event)",
                  which(!status %in% c(0, 1)))
  if (sum(status) == 0) stop("there are no events (status is all 0)", call. = FALSE)

  if (any(!is.finite(t) | t <= 0))
    .qast_stop_at("t must be finite and positive", which(!is.finite(t) | t <= 0))
  if (any(!is.finite(qt) | qt < 0))
    .qast_stop_at("qt must be finite and non-negative",
                  which(!is.finite(qt) | qt < 0))
  # Utilities in [0, 1]: the adjusted time does not exceed the survival
  # time (relative slack for rounding)
  over <- which(qt > t * (1 + 1e-8))
  if (length(over))
    .qast_stop_at("qt (quality-adjusted time) greater than t; check that the arguments are not swapped",
                  over)

  lt <- h1$linkfun(t)
  lq <- h2$linkfun(qt)
  if (any(!is.finite(lt)))
    .qast_stop_at(sprintf("h1(t) is not finite (link \"%s\")", h1$name),
                  which(!is.finite(lt)))
  if (any(!is.finite(lq)))
    .qast_stop_at(sprintf("h2(qt) is not finite (link \"%s\"); with the log link, qt must be > 0",
                          h2$name), which(!is.finite(lq)))

  if (!is.null(L)) {
    if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L <= 0)
      stop("L must be a positive number", call. = FALSE)
    if (!is.finite(h1$linkfun(L)))
      stop(sprintf("h1(L) is not finite (link \"%s\")", h1$name), call. = FALSE)
  }
  invisible(TRUE)
}

# Methods that do not accept a single coefficient (dfoptim::nmk and mads
# call optimize(); minqa::uobyqa/bobyqa and hjkb do not converge with
# p = 1): with one covariate, use subplex, also derivative-free.
.qast_method_for_p <- function(method, p) {
  if (p == 1L && method %in% c("Nelder-Mead", "uobyqa", "mads", "bobyqa",
                               "hjkb")) {
    message(sprintf(paste0("method = \"%s\" does not support a single coefficient; ",
                           "using \"subplex\"."), method))
    return("subplex")
  }
  method
}

# M is required only for the bootstrap
.validate_var_args <- function(var_comp, var_method, M) {
  if (!is.logical(var_comp) || length(var_comp) != 1L || is.na(var_comp))
    stop("var_comp must be TRUE or FALSE", call. = FALSE)
  if (var_comp && var_method == "bootstrap") {
    if (is.null(M) || !is.numeric(M) || length(M) != 1L || !is.finite(M) ||
        M < 2 || M != round(M))
      stop("with var_method = \"bootstrap\", M (number of replicates) must be ",
           "an integer >= 2", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Intercepts and distribution of the residuals
# ---------------------------------------------------------------------------
# The rank estimating equations are invariant to shifts and do not
# identify the intercept. After the fit, the distribution of the residuals
# of each component is estimated and stored in the qast_fit; the intercept
# is its mean and predict(type = "survival") uses it in S(t | z) = S_res(h(t) - z'b).

# Survival time: Kaplan-Meier curve of the residuals
# e = h1(t) - z'beta (right-continuous steps: S(x) = P(e > x)).
.resid_dist_km <- function(e, status) {
  f <- survival::survfit(survival::Surv(e, status) ~ 1)
  list(time = f$time, surv = f$surv)
}

# Quality-adjusted time: Zhao and Tsiatis distribution of the residuals
# r = h2(qt) - z'theta. Each event gets mass 1 / (n G(t_i)), with G the
# censoring KM on the survival time scale (G = 0 -> mass 0, as in the
# surv_qt_fun kernel); the masses of tied residuals are added up. G is
# the marginal KM by default; the Fine and Gelber fits pass their H_G.
.resid_dist_ipcw <- function(r, t, status, G = km_weights_at(t, status)) {
  w <- ifelse(status == 1 & is.finite(G) & G > 0, 1 / G, 0) / length(r)
  keep <- w > 0
  # Exact values (not through names/character, which round to 15 digits
  # and would shift the jumps)
  v <- r[keep]
  u <- sort(unique(v))
  mass <- as.vector(rowsum(w[keep], match(v, u), reorder = TRUE))
  list(value = u, mass = mass)
}

# Mean of the KM curve (Ding and Nan, 2015), restricted to the largest
# residual: the mass left at the end of the curve (largest residual
# censored) stays at that point. Computed from the jumps, and not by the
# rmean of survfit, which integrates from min(e) without adding it back
# when there are negative residuals.
.km_mean <- function(d) {
  s <- d$surv
  k <- length(s)
  sum(-diff(c(1, s)) * d$time) + s[k] * d$time[k]
}

# IPCW mean (Zhao and Tsiatis, 2000): (1/n) sum status_i r_i / G(t_i)
.ipcw_mean <- function(d) sum(d$mass * d$value)

# S(x) = P(e > x) on the KM curve; 1 before the first jump
.surv_km_at <- function(d, x) {
  idx <- findInterval(x, d$time)
  out <- d$surv[pmax(idx, 1L)]
  out[idx == 0L] <- 1
  out
}

# S(x) = P(r > x) = sum of the masses of the residuals > x, bounded to [0, 1]
# (the total IPCW mass is close to 1, but not exactly 1)
.surv_ipcw_at <- function(d, x) {
  tail_mass <- rev(cumsum(rev(d$mass)))
  idx <- findInterval(x, d$value) + 1L
  out <- c(tail_mass, 0)[idx]
  pmin(pmax(out, 0), 1)
}

.qast_resid_dists <- function(lt, lq, t, status, z, betahat, thetahat,
                              G = km_weights_at(t, status)) {
  list(beta  = .resid_dist_km(lt - drop(z %*% betahat), status),
       theta = .resid_dist_ipcw(lq - drop(z %*% thetahat), t, status, G))
}

# ---------------------------------------------------------------------------
# Class qast_fit
# ---------------------------------------------------------------------------

.new_qast_fit <- function(model, call, z, status, betahat, thetahat,
                          vcov_beta, vcov_theta, var_comp, var_method, M,
                          vcov_diff = NULL,
                          settings, lt, lq, t, qt, formula_info = NULL,
                          replicates = NULL) {
  nms <- colnames(z)
  names(betahat) <- names(thetahat) <- nms
  # Censoring: in the Fine and Gelber model, the KM of the partition of the
  # fit (H_G, the same as the equations and the diagnostic); in the proposed one, the marginal KM
  G <- if (model == "fg" && !is.null(settings$cens_group))
    km_weights_by_group(t, status, settings$cens_group)
  else km_weights_at(t, status)
  rd <- .qast_resid_dists(lt, lq, t, status, z, betahat, thetahat, G)
  fit <- list(beta = betahat, theta = thetahat,
              intercept = c(beta = .km_mean(rd$beta),
                            theta = .ipcw_mean(rd$theta)))
  if (var_comp) {
    dimnames(vcov_beta) <- dimnames(vcov_theta) <- list(nms, nms)
    fit$var_beta  <- stats::setNames(diag(vcov_beta), nms)
    fit$var_theta <- stats::setNames(diag(vcov_theta), nms)
    fit$vcov_beta  <- vcov_beta
    fit$vcov_theta <- vcov_theta
    # theta - beta, from the joint replicates
    if (!is.null(vcov_diff)) {
      dimnames(vcov_diff) <- list(nms, nms)
      fit$var_diff  <- stats::setNames(diag(vcov_diff), nms)
      fit$vcov_diff <- vcov_diff
    }
  }
  fit$model <- model
  fit$call <- call
  fit$n <- length(status)
  fit$n_events <- as.integer(sum(status))
  fit$var_method <- if (var_comp) var_method else NULL
  fit$M <- if (var_comp && var_method == "bootstrap") M else NULL
  fit$settings <- settings
  # For predict(): covariates of the fit, distribution of the residuals,
  # event times (default grid) and, for formulas, terms and levels
  fit$z <- z
  fit$y <- cbind(time = t, qtime = qt, status = status)
  fit$resid_dist <- rd
  fit$event_times <- list(beta  = sort(unique(t[status == 1])),
                          theta = sort(unique(qt[status == 1])))
  fit$formula_info <- formula_info
  # Variance replicates (coefficients + rows used) for
  # predict(se.fit = TRUE)
  fit$replicates <- if (var_comp) replicates else NULL
  class(fit) <- "qast_fit"
  fit
}

# L, phi and censoring partition of the Fine and Gelber model
.qast_print_fg_settings <- function(settings, digits) {
  cat("L = ", format(settings$L, digits = digits), ", phi = ",
      settings$phi, sep = "")
  if (!is.null(settings$cens_group))
    cat(", censoring groups = ", nlevels(settings$cens_group),
        if (isTRUE(settings$cens_group_default)) " (no stratification)",
        sep = "")
  cat("\n")
}

.qast_model_label <- function(model) {
  switch(model,
         fg   = "Fine and Gelber (2001) semiparametric regression model",
         prop = "Bivariate semiparametric regression model (proposed)")
}

#' Methods for QAST fits
#'
#' \code{semi_aft_fg()} and \code{semi_aft_prop()} return an object of class
#' \code{"qast_fit"}: a list with \code{beta} (survival time component),
#' \code{theta} (quality-adjusted survival time component),
#' \code{intercept} and, when \code{var_comp = TRUE}, \code{var_beta},
#' \code{var_theta} and the covariance matrices \code{vcov_beta} and
#' \code{vcov_theta}. It also holds \code{model}, \code{call}, \code{n},
#' \code{n_events}, \code{settings} and, for \code{\link{predict.qast_fit}}
#' and \code{\link{qast_pp}}, the data of the fit (\code{y}: columns
#' \code{time}, \code{qtime} and \code{status}; \code{z}: covariate
#' matrix), the estimated distributions of the residuals
#' (\code{resid_dist}), the distinct event times of each component
#' (\code{event_times}) and, for formula fits, the terms, factor levels
#' and contrasts (\code{formula_info}).
#'
#' @section Intercepts:
#' The rank-based estimating equations are invariant to a location shift,
#' so they do not identify an intercept. After the fit, each intercept is
#' estimated as the mean of the distribution of the residuals of its
#' component, and stored in \code{intercept}, a named vector with elements
#' \code{beta} and \code{theta}:
#' \itemize{
#'   \item \code{beta}: mean of the Kaplan-Meier curve of the residuals
#'     \eqn{h_1(T_i) - Z_i^\top\hat\beta}{h1(T_i) - Z_i'beta} (Ding and Nan,
#'     2015), restricted to the largest residual when it is censored;
#'   \item \code{theta}: inverse probability of censoring weighted mean of
#'     the residuals \eqn{h_2(Q_i) - Z_i^\top\hat\theta}{h2(Q_i) - Z_i'theta}
#'     (Zhao and Tsiatis, 2000), with the Kaplan-Meier censoring
#'     survival \eqn{\hat G}{G} evaluated at the survival time:
#'     \eqn{n^{-1}\sum_i \delta_i r_i / \hat G(T_i)}{(1/n) sum_i delta_i r_i / G(T_i)}.
#'     For \code{semi_aft_prop()} fits, \eqn{\hat G}{G} is the marginal
#'     Kaplan-Meier, and this is the mean of the distribution estimated by
#'     \code{\link{surv_qt_fun}}; for \code{semi_aft_fg()} fits, it is the
#'     Kaplan-Meier within the censoring groups of the fit
#'     (\code{cens_group}), the same weights as the estimating function
#'     (the marginal one without stratification, the default).
#' }
#' @section Quality per unit of survival:
#' \code{coef(fit, "diff")} is \eqn{\hat\theta - \hat\beta}{theta - beta}.
#' With \code{h1 = h2 = "log"}, \eqn{\log\{Q(T)/T\} = (\theta - \beta)^\top Z +
#' \varepsilon_2 - \varepsilon_1}{log(Q/T) = (theta - beta)'Z + e2 - e1}, so
#' \eqn{\exp(\theta_j - \beta_j)}{exp(theta_j - beta_j)} is the change in
#' the quality of life per unit of survival for a unit increase in the
#' covariate \eqn{j}, and also the ratio of the mean quality-adjusted
#' survival per unit of mean survival (Fine and Gelber, 2001, Section 2;
#' the ratio of means extrapolates beyond the follow-up). Its variance
#' uses the joint replicates of the fit, in which beta and theta are
#' estimated on the same subsamples or resamples:
#' \eqn{Var(\hat\theta_j - \hat\beta_j) = Var(\hat\theta_j) +
#' Var(\hat\beta_j) - 2 Cov(\hat\beta_j, \hat\theta_j)}{Var(theta_j -
#' beta_j) = Var(theta_j) + Var(beta_j) - 2 Cov(beta_j, theta_j)}.
#' \code{confint()} gives the normal interval of Fine and Gelber (2001,
#' Section 3) or, with bootstrap replicates, the percentile interval;
#' exponentiate the limits for the ratio. \code{summary()} shows the
#' table and the normal interval of the ratio.
#'
#' The linear predictor of the transformed time is then
#' \code{intercept["beta"] + z \%*\% beta} (and likewise with
#' \code{theta}). The intercepts have no standard errors and are not
#' part of \code{coef()} or \code{vcov()}.
#'
#' @references Fine, J. P. and Gelber, R. D. (2001). Joint regression
#'   analysis of survival and quality-adjusted survival. \emph{Biometrics},
#'   57(2), 376-382. \doi{10.1111/j.0006-341X.2001.00376.x}
#'
#'   Ding, Y. and Nan, B. (2015). Estimating mean survival time:
#'   when is it possible? \emph{Scandinavian Journal of Statistics}, 42(2),
#'   397-413. \doi{10.1111/sjos.12112}
#'
#'   Zhao, H. and Tsiatis, A. A. (2000). Estimating mean quality adjusted
#'   lifetime with censored data. \emph{Sankhya: The Indian Journal of
#'   Statistics, Series B}, 62, 175-188.
#'
#' @param x,object Object of class \code{"qast_fit"}
#' @param digits Significant digits when printing
#' @param component \code{"beta"} (survival time), \code{"theta"}
#'   (quality-adjusted time) or \code{"diff"} (\code{theta - beta}); in
#'   \code{coef()}, \code{"both"} returns a matrix with the columns
#'   \code{beta} and \code{theta}
#' @param parm Coefficients (names or positions); all by default
#' @param level,conf.level Confidence level
#' @param type \code{"wald"} (normal interval, the default) or
#'   \code{"percentile"} (quantiles of the replicates; bootstrap fits only)
#' @param signif.stars If TRUE, flags the p-values with stars
#' @param ... Not used
#' @return \code{coef()}: matrix (p x 2) or vector; \code{vcov()}:
#'   (jackknife or bootstrap) covariance matrix of the component;
#'   \code{confint()}: matrix with the lower and upper limits;
#'   \code{summary()}: a \code{"summary.qast_fit"} object with the
#'   coefficient tables; \code{print()}: the object, invisibly.
#' @name qast_fit
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
#' coef(fit)
#'
#' # Quality per unit of survival: exp(theta - beta), with intervals
#' fit_v <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                      L = quantile(d$t, 0.9), var_comp = TRUE)
#' summary(fit_v)
#' exp(confint(fit_v, component = "diff"))
NULL

#' @rdname qast_fit
#' @export
print.qast_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...) {
  cat("\n", .qast_model_label(x$model), "\n\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("n = ", x$n, ", number of events = ", x$n_events, "\n", sep = "")
  if (x$model == "fg") .qast_print_fg_settings(x$settings, digits)
  cat("\nCoefficients:\n")
  print(coef(x), digits = digits)
  cat("\nIntercepts (mean of the residuals):\n")
  print(x$intercept, digits = digits)
  if (is.null(x$var_beta))
    cat("\n(no variance: refit with var_comp = TRUE for standard errors)\n")
  else
    cat("\nUse summary() for standard errors and tests (",
        x$var_method, ").\n", sep = "")
  invisible(x)
}

#' @rdname qast_fit
#' @export
coef.qast_fit <- function(object,
                          component = c("both", "beta", "theta", "diff"),
                          ...) {
  component <- match.arg(component)
  switch(component,
         both  = cbind(beta = object$beta, theta = object$theta),
         beta  = object$beta,
         theta = object$theta,
         diff  = object$theta - object$beta)
}

#' @rdname qast_fit
#' @export
vcov.qast_fit <- function(object, component = c("beta", "theta", "diff"),
                          ...) {
  component <- match.arg(component)
  if (is.null(object$vcov_beta))
    stop("fit has no variance: refit with var_comp = TRUE", call. = FALSE)
  if (component == "diff" && is.null(object$vcov_diff))
    stop("this fit has no joint replicates for theta - beta; refit it ",
         "with the current version of QAST", call. = FALSE)
  switch(component,
         beta  = object$vcov_beta,
         theta = object$vcov_theta,
         diff  = object$vcov_diff)
}

#' @rdname qast_fit
#' @export
confint.qast_fit <- function(object, parm, level = 0.95,
                             component = c("beta", "theta", "diff"),
                             type = c("wald", "percentile"), ...) {
  component <- match.arg(component)
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1))
    stop("level must be a number in (0, 1)", call. = FALSE)
  est <- coef(object, component)
  V <- vcov(object, component)
  if (missing(parm)) parm <- names(est)
  if (is.numeric(parm)) parm <- names(est)[parm]
  if (anyNA(parm) || !all(parm %in% names(est)))
    stop("unknown coefficients in parm", call. = FALSE)
  a <- (1 - level) / 2
  probs <- c(a, 1 - a)
  if (type == "wald") {
    q <- stats::qnorm(1 - a)
    se <- sqrt(diag(V))[parm]
    ci <- cbind(est[parm] - q * se, est[parm] + q * se)
  } else {
    if (!identical(object$var_method, "bootstrap"))
      stop("type = \"percentile\" requires bootstrap replicates ",
           "(var_method = \"bootstrap\")", call. = FALSE)
    reps <- object$replicates[[component]]$reps
    colnames(reps) <- names(est)
    ci <- t(apply(reps[, parm, drop = FALSE], 2L, stats::quantile,
                  probs = probs, names = FALSE))
  }
  dimnames(ci) <- list(parm, paste(format(100 * probs, trim = TRUE,
                                          scientific = FALSE, digits = 3),
                                   "%"))
  ci
}

.qast_coef_table <- function(est, v) {
  if (is.null(v)) {
    tab <- cbind(est, exp(est))
    colnames(tab) <- c("coef", "exp(coef)")
  } else {
    se <- sqrt(v)
    zval <- est / se
    tab <- cbind(est, exp(est), se, zval, 2 * stats::pnorm(-abs(zval)))
    colnames(tab) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  }
  rownames(tab) <- names(est)
  tab
}

#' @rdname qast_fit
#' @export
summary.qast_fit <- function(object, conf.level = 0.95, ...) {
  diff <- object$theta - object$beta
  ratio <- NULL
  if (!is.null(object$var_diff)) {
    q <- stats::qnorm(1 - (1 - conf.level) / 2)
    se <- sqrt(object$var_diff)
    ratio <- cbind(exp(diff), exp(diff - q * se), exp(diff + q * se))
    dimnames(ratio) <- list(names(diff),
                            c("exp(theta - beta)",
                              sprintf("lower .%02d", round(100 * conf.level)),
                              sprintf("upper .%02d", round(100 * conf.level))))
  }
  out <- list(
    model = object$model, call = object$call, n = object$n,
    n_events = object$n_events, settings = object$settings,
    var_method = object$var_method, M = object$M,
    intercept = object$intercept,
    beta  = .qast_coef_table(object$beta, object$var_beta),
    theta = .qast_coef_table(object$theta, object$var_theta),
    diff  = .qast_coef_table(diff, object$var_diff),
    ratio = ratio, conf.level = conf.level
  )
  class(out) <- "summary.qast_fit"
  out
}

#' @rdname qast_fit
#' @export
print.summary.qast_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {
  rule <- strrep("-", 62)
  cat("\n", rule, "\n", .qast_model_label(x$model), "\n", rule, "\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat("n = ", x$n, ", number of events = ", x$n_events, "\n", sep = "")
  if (x$model == "fg") .qast_print_fg_settings(x$settings, digits)
  if (!is.null(x$var_method))
    cat("Variance: ", x$var_method,
        if (!is.null(x$M)) paste0(" (M = ", x$M, ")"), "\n", sep = "")
  cat("\nSurvival time component\n")
  stats::printCoefmat(x$beta, digits = digits, signif.stars = signif.stars,
                      signif.legend = FALSE)
  cat("Intercept (Kaplan-Meier mean of the residuals): ",
      format(x$intercept[["beta"]], digits = digits), "\n", sep = "")
  cat(rule, "\nQuality-adjusted survival time component\n", sep = "")
  stats::printCoefmat(x$theta, digits = digits, signif.stars = signif.stars,
                      signif.legend = FALSE)
  cat("Intercept (IPCW mean of the residuals): ",
      format(x$intercept[["theta"]], digits = digits), "\n", sep = "")
  if (!is.null(x$diff)) {
    cat(rule, "\nQuality per unit of survival (theta - beta)\n", sep = "")
    stats::printCoefmat(x$diff, digits = digits, signif.stars = signif.stars,
                        signif.legend = FALSE)
    if (!is.null(x$ratio)) {
      cat("\n")
      print(x$ratio, digits = digits)
    }
  }
  if (signif.stars && ncol(x$theta) > 2L)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

# Covariate matrix for predict(): formula -> model.matrix with the terms,
# levels and contrasts of the fit; matrix -> columns by name (or by
# position, if the fit had no names)
.qast_newdata_z <- function(object, newdata) {
  if (is.null(newdata)) return(object$z)
  p <- length(object$beta)
  nms <- names(object$beta)
  info <- object$formula_info
  if (!is.null(info)) {
    if (!is.data.frame(newdata))
      newdata <- as.data.frame(newdata)
    mf <- stats::model.frame(info$terms, newdata, xlev = info$xlevels,
                             na.action = stats::na.pass)
    z <- stats::model.matrix(info$terms, mf, contrasts.arg = info$contrasts)
    z <- z[, colnames(z) != "(Intercept)", drop = FALSE]
  } else {
    z <- as.matrix(newdata)
    if (!is.null(nms) && !is.null(colnames(z))) {
      miss <- setdiff(nms, colnames(z))
      if (length(miss))
        stop("newdata lacks the columns: ", paste(miss, collapse = ", "),
             call. = FALSE)
      z <- z[, nms, drop = FALSE]
    }
  }
  if (ncol(z) != p)
    stop(sprintf("newdata gives %d covariate columns; the fit has %d",
                 ncol(z), p), call. = FALSE)
  if (!is.numeric(z)) stop("newdata covariates must be numeric", call. = FALSE)
  if (anyNA(z)) .qast_stop_at("newdata contains NA", which(!stats::complete.cases(z)))
  z
}

# Evaluation times in the domain of the link of the component
.qast_check_times <- function(object, comp, times) {
  link <- stats::make.link(object$settings[[if (comp == "beta") "h1" else "h2"]])
  if (is.null(times)) times <- object$event_times[[comp]]
  if (!is.numeric(times) || anyNA(times))
    stop("times must be a numeric vector without NA", call. = FALSE)
  ht <- suppressWarnings(link$linkfun(times))   # NaN handled below
  if (any(!is.finite(ht)))
    stop(sprintf("times must be in the domain of the link \"%s\" (with log, > 0)",
                 link$name), call. = FALSE)
  list(times = times, ht = ht)
}

# S(t | z_k) on a grid of times (ht = h(times)), given the distribution of
# the residuals and the coefficients of a component
.qast_surv_matrix <- function(dist, coef, z, ht, comp) {
  lp <- drop(z %*% coef)
  x <- outer(-lp, ht, `+`)          # residual: h(t) - z'b
  surv_at <- if (comp == "beta") .surv_km_at else .surv_ipcw_at
  matrix(surv_at(dist, as.vector(x)), nrow = nrow(x))
}

# Restriction time of type = "rmean": a scalar, or a vector named by
# component (beta, theta); default, the largest event time of the fit
# data on the scale of the component
.qast_tau <- function(object, comp, tau) {
  if (is.null(tau)) return(max(object$event_times[[comp]]))
  if (!is.numeric(tau) || anyNA(tau) || length(tau) < 1L)
    stop("tau must be a positive number (or one per component, named ",
         "beta and theta)", call. = FALSE)
  if (length(tau) > 1L) {
    if (is.null(names(tau)) || !comp %in% names(tau))
      stop("with more than one value, tau must be named by component ",
           "(beta, theta)", call. = FALSE)
    tau <- tau[[comp]]
  }
  if (!is.finite(tau) || tau <= 0)
    stop("tau must be finite and positive", call. = FALSE)
  unname(tau)
}

# Increasing link (the restricted mean maps the residual jumps back to
# the time scale with linkinv)
.qast_link_increasing <- function(object, comp) {
  link <- stats::make.link(object$settings[[if (comp == "beta") "h1" else "h2"]])
  if (!(link$linkfun(2) > link$linkfun(1)))
    stop(sprintf("type = \"rmean\" requires an increasing link (link \"%s\")",
                 link$name), call. = FALSE)
  link
}

# Restricted mean int_0^tau S(t | z_k) dt, exact for the step curve:
# S(t | z) = S_res(h(t) - z'b) jumps at u_j = h^{-1}(v_j + z'b), v_j the
# support of the residual distribution (residuals below h(0) map to 0),
# and is constant between jumps (val: value before the first jump and
# after each one)
.qast_rmean <- function(dist, coef, z, comp, link, tau) {
  if (comp == "beta") {
    v <- dist$time
    val <- c(1, dist$surv)
  } else {
    v <- dist$value
    tail_mass <- rev(cumsum(rev(dist$mass)))
    val <- pmin(pmax(c(tail_mass, 0), 0), 1)
  }
  h0 <- suppressWarnings(link$linkfun(0))
  if (is.na(h0)) h0 <- -Inf
  lp <- drop(z %*% coef)
  vapply(lp, function(a) {
    eta <- v + a
    u <- numeric(length(eta))
    pos <- eta > h0
    u[pos] <- link$linkinv(eta[pos])
    lo <- pmin(c(0, u), tau)
    hi <- pmin(c(u, Inf), tau)
    sum(val * (hi - lo))
  }, numeric(1))
}

# Distribution of the residuals of a component recomputed on the rows
# 'rows' of the data of the fit (jackknife or bootstrap replicate), with
# the same estimators and censoring weights as the fit
.qast_dist_rows <- function(object, comp, rows, coef) {
  y <- object$y[rows, , drop = FALSE]
  z <- object$z[rows, , drop = FALSE]
  t <- y[, "time"]; status <- y[, "status"]
  if (comp == "beta") {
    lt <- stats::make.link(object$settings$h1)$linkfun(t)
    return(.resid_dist_km(lt - drop(z %*% coef), status))
  }
  lq <- stats::make.link(object$settings$h2)$linkfun(y[, "qtime"])
  G <- if (object$model == "fg" && !is.null(object$settings$cens_group))
    km_weights_by_group(t, status, object$settings$cens_group[rows])
  else km_weights_at(t, status)
  .resid_dist_ipcw(lq - drop(z %*% coef), t, status, G)
}

# Rows and coefficients of each replicate of a component
.qast_rep_list <- function(object, comp) {
  rp <- object$replicates[[comp]]
  if (is.null(rp))
    stop(sprintf(paste0("the fit has no variance replicates for %s; refit ",
                        "with var_comp = TRUE"), comp), call. = FALSE)
  n <- object$n
  lapply(seq_len(nrow(rp$reps)), function(k) list(
    coef = rp$reps[k, ],
    rows = if (!is.null(rp$drop)) seq_len(n)[-rp$drop[k]] else rp$index[k, ]))
}

# Standard deviation of the replicates, with the vcov() formula: jackknife
# sqrt((n - 1) / n_valid * sum of squared deviations); bootstrap sd()
.qast_rep_se <- function(object, vals) {
  if (identical(object$var_method, "jackknife")) {
    nv <- ncol(vals)
    sqrt((object$n - 1) / nv * rowSums((vals - rowMeans(vals))^2))
  } else {
    apply(vals, 1L, stats::sd)
  }
}

#' Predictions from QAST fits
#'
#' Linear predictors and survival probabilities for the survival time
#' (\code{beta}) and quality-adjusted survival time (\code{theta})
#' components of a \code{"qast_fit"}.
#'
#' Under the model \eqn{h_1(T) = \alpha + Z^\top\beta + \varepsilon}{h1(T)
#' = alpha + Z'beta + e}, with the error independent of the covariates,
#' the survival function given the covariates is a shift of the residual
#' survival function:
#' \deqn{S(t \mid z) = P\{h_1(T) - Z^\top\beta > h_1(t) - z^\top\beta\} =
#'   S_{e}\{h_1(t) - z^\top\hat\beta\},}{S(t | z) = S_e(h1(t) - z'beta),}
#' with \eqn{e = h_1(T) - Z^\top\hat\beta}{e = h1(T) - Z'beta} (the
#' intercept cancels out). \eqn{S_e}{S_e} is estimated by:
#' \itemize{
#'   \item \code{beta}: the Kaplan-Meier curve of the residuals
#'     \eqn{h_1(T_i) - Z_i^\top\hat\beta}{h1(T_i) - Z_i'beta};
#'   \item \code{theta}: the estimator of Zhao and Tsiatis (1997) applied
#'     to the residuals \eqn{h_2(Q_i) - Z_i^\top\hat\theta}{h2(Q_i) -
#'     Z_i'theta}, which weights each event by the inverse of the
#'     Kaplan-Meier censoring survival at its survival time and corrects
#'     the induced informative censoring of the quality-adjusted time.
#'     Its total mass is close to, but not exactly, 1, so the estimates
#'     are truncated to [0, 1].
#' }
#' The censoring is the Kaplan-Meier: within the censoring groups of the
#' fit for \code{semi_aft_fg()} (\code{cens_group}; marginal without
#' stratification, the default), as in the estimating function; marginal
#' for \code{semi_aft_prop()}, which assumes censoring independent of the
#' covariates (the fit itself uses a Cox model). The survival
#' probabilities are \eqn{P(T > t \mid z)}{P(T > t | z)} (right-continuous
#' steps), the same convention as \code{surv_qt_fun()}.
#'
#' For the survival time, this is the prediction of Ying, Wei and Lin
#' (1992): the Kaplan-Meier curve of the residuals of the rank estimate,
#' evaluated at \eqn{h_1(t) - z^\top\hat\beta}{h1(t) - z'beta}. They use the
#' left-continuous version of the curve, \eqn{P(e \ge x)}{P(e >= x)}; the
#' two versions differ only at the jumps.
#'
#' With \code{se.fit = TRUE} (fits with \code{var_comp = TRUE}), the
#' prediction is recomputed in every jackknife or bootstrap replicate of
#' the fit: the replicate coefficients, and the intercept and the residual
#' distribution re-estimated on the data of the replicate, with the same
#' censoring weights. The standard errors use the same formula as
#' \code{vcov()}: the delete-1 jackknife variance or the bootstrap sample
#' variance. The survival intervals are pointwise: normal
#' (\eqn{\pm z_{1-\alpha/2}}{+/- z} times the standard error) for the
#' jackknife, percentile for the bootstrap, truncated to [0, 1]. The
#' censoring weights are re-estimated in every replicate, as in the
#' variance of the coefficients.
#'
#' \code{type = "rmean"} gives the restricted mean
#' \eqn{\int_0^\tau S(t \mid z)\,dt = E\{\min(T, \tau) \mid z\}}{int_0^tau
#' S(t | z) dt = E[min(T, tau) | z]} (and the same for the quality-adjusted
#' time), computed exactly on the step curve of each subject: the curves
#' are steps between the jumps \eqn{h^{-1}(v_j + z^\top\hat\beta)}{h^-1(v_j
#' + z'beta)}, \eqn{v_j} the residuals. The unrestricted mean is not used
#' because the right tail of the curves is not identified when the
#' largest times are censored. It requires an increasing link (such as
#' the default \code{"log"}).
#'
#' @param object Object of class \code{"qast_fit"}
#' @param newdata Covariates to predict for. For fits from a formula, a
#'   data frame with the covariates of the formula (factors with the
#'   levels of the fit); for fits from a matrix, a matrix or data frame
#'   with the columns of \code{z} (matched by name when the fit had
#'   column names). \code{NULL} (default) uses the covariates of the fit.
#' @param type \code{"lp"} (linear predictor including the intercept,
#'   \code{intercept + z \%*\% coef}, on the scale of the link),
#'   \code{"survival"} (survival probabilities at \code{times}) or
#'   \code{"rmean"} (mean time restricted to \code{tau})
#' @param component \code{"both"}, \code{"beta"} (survival time) or
#'   \code{"theta"} (quality-adjusted time)
#' @param times Times at which to evaluate the survival probabilities, on
#'   the original scale (survival time for \code{beta}, quality-adjusted
#'   time for \code{theta}). \code{NULL} (default) uses the distinct
#'   event times of the fit data of each component.
#' @param tau Restriction time of \code{type = "rmean"}, on the original
#'   scale: a number, or a vector named by component
#'   (\code{c(beta = 5, theta = 3)}). \code{NULL} (default) uses the
#'   largest event time of the fit data of each component.
#' @param se.fit If TRUE, also returns standard errors (and, for
#'   \code{type = "survival"}, pointwise confidence limits) from the
#'   variance replicates of the fit; requires \code{var_comp = TRUE}
#' @param level Confidence level of the survival intervals
#' @param ... Not used
#' @return For \code{type = "lp"}: a vector (one component) or a matrix
#'   with columns \code{beta} and \code{theta}; with \code{se.fit = TRUE},
#'   a list with \code{fit} and \code{se.fit} of the same shape. The same
#'   for \code{type = "rmean"}, with the restriction times in the
#'   attribute \code{"tau"}. For
#'   \code{type = "survival"}: a matrix with one row per subject of
#'   \code{newdata} and one column per time (one component), or a list
#'   with the matrices \code{beta} and \code{theta}; with
#'   \code{se.fit = TRUE}, each matrix is replaced by a list with
#'   \code{surv}, \code{se}, \code{lower} and \code{upper}.
#' @references Ying, Z., Wei, L. J. and Lin, J. S. (1992). Prediction of
#'   survival probability based on a linear regression model.
#'   \emph{Biometrika}, 79(1), 205-209. \doi{10.1093/biomet/79.1.205}
#'
#'   Zhao, H. and Tsiatis, A. A. (1997). A consistent estimator
#'   for the distribution of quality adjusted survival time.
#'   \emph{Biometrika}, 84(2), 339-348. \doi{10.1093/biomet/84.2.339}
#' @seealso \code{\link{qast_fit}} for the intercepts;
#'   \code{\link{surv_qt_fun}}
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(age = rnorm(n), trt = rbinom(n, 1, 0.5))
#' d$t <- rexp(n, rate = exp(-0.3 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#' fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d)
#'
#' new <- data.frame(age = c(0, 0), trt = c(0, 1))
#' predict(fit, new)                       # linear predictors
#' s <- predict(fit, new, type = "survival", times = c(0.5, 1, 2))
#' s$beta                                  # P(T > t | z)
#' s$theta                                 # P(Q > t | z)
#' predict(fit, new, type = "rmean", tau = 2)  # E[min(T, 2) | z], E[min(Q, 2) | z]
#'
#' # Standard errors and pointwise 95% intervals (jackknife replicates)
#' fit_v <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d,
#'                        var_comp = TRUE)
#' predict(fit_v, new, se.fit = TRUE)
#' predict(fit_v, new, type = "survival", component = "theta",
#'         times = c(0.5, 1), se.fit = TRUE)
#' @export
predict.qast_fit <- function(object, newdata = NULL,
                             type = c("lp", "survival", "rmean"),
                             component = c("both", "beta", "theta"),
                             times = NULL, tau = NULL, se.fit = FALSE,
                             level = 0.95, ...) {
  type <- match.arg(type)
  component <- match.arg(component)
  if (is.null(object$resid_dist))
    stop("this fit has no stored residual distribution; refit it with ",
         "the current version of QAST", call. = FALSE)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
    stop("se.fit must be TRUE or FALSE", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1))
    stop("level must be a number in (0, 1)", call. = FALSE)
  if (se.fit && is.null(object$replicates))
    stop("se.fit = TRUE requires the variance replicates: refit with ",
         "var_comp = TRUE", call. = FALSE)
  z <- .qast_newdata_z(object, newdata)
  comps <- if (component == "both") c("beta", "theta") else component
  reps <- if (se.fit)
    lapply(stats::setNames(comps, comps), function(cp) {
      lapply(.qast_rep_list(object, cp), function(r) {
        r$dist <- .qast_dist_rows(object, cp, r$rows, r$coef)
        r
      })
    })

  if (type == "lp") {
    lp_of <- function(cp) object$intercept[[cp]] + drop(z %*% object[[cp]])
    out <- vapply(comps, lp_of, numeric(nrow(z)))
    if (!is.matrix(out)) out <- matrix(out, nrow = nrow(z),
                                       dimnames = list(NULL, comps))
    rownames(out) <- rownames(z)
    if (se.fit) {
      se <- vapply(comps, function(cp) {
        vals <- vapply(reps[[cp]], function(r) {
          a <- if (cp == "beta") .km_mean(r$dist) else .ipcw_mean(r$dist)
          a + drop(z %*% r$coef)
        }, numeric(nrow(z)))
        .qast_rep_se(object, matrix(vals, nrow = nrow(z)))
      }, numeric(nrow(z)))
      if (!is.matrix(se)) se <- matrix(se, nrow = nrow(z),
                                       dimnames = list(NULL, comps))
      rownames(se) <- rownames(z)
      if (length(comps) == 1L) return(list(fit = out[, 1L], se.fit = se[, 1L]))
      return(list(fit = out, se.fit = se))
    }
    return(if (length(comps) == 1L) out[, 1L] else out)
  }

  if (type == "rmean") {
    links <- lapply(stats::setNames(comps, comps),
                    function(cp) .qast_link_increasing(object, cp))
    taus <- vapply(comps, function(cp) .qast_tau(object, cp, tau), numeric(1))
    rm_of <- function(cp, dist, coef)
      .qast_rmean(dist, coef, z, cp, links[[cp]], taus[[cp]])
    out <- vapply(comps, function(cp)
      rm_of(cp, object$resid_dist[[cp]], object[[cp]]), numeric(nrow(z)))
    if (!is.matrix(out)) out <- matrix(out, nrow = nrow(z),
                                       dimnames = list(NULL, comps))
    rownames(out) <- rownames(z)
    attr(out, "tau") <- taus
    if (se.fit) {
      se <- vapply(comps, function(cp) {
        vals <- vapply(reps[[cp]], function(r) rm_of(cp, r$dist, r$coef),
                       numeric(nrow(z)))
        .qast_rep_se(object, matrix(vals, nrow = nrow(z)))
      }, numeric(nrow(z)))
      if (!is.matrix(se)) se <- matrix(se, nrow = nrow(z),
                                       dimnames = list(NULL, comps))
      rownames(se) <- rownames(z)
      if (length(comps) == 1L) {
        f <- out[, 1L]; attr(f, "tau") <- taus
        return(list(fit = f, se.fit = se[, 1L]))
      }
      return(list(fit = out, se.fit = se))
    }
    if (length(comps) == 1L) {
      f <- out[, 1L]; attr(f, "tau") <- taus
      return(f)
    }
    return(out)
  }

  out <- lapply(stats::setNames(comps, comps), function(cp) {
    tm <- .qast_check_times(object, cp, times)
    est <- .qast_surv_matrix(object$resid_dist[[cp]], object[[cp]], z,
                             tm$ht, cp)
    dn <- list(rownames(z), format(tm$times, digits = 6L, trim = TRUE))
    dimnames(est) <- dn
    if (!se.fit) return(est)
    # replicates: n_rep matrices (rows x times) stacked by column
    vals <- vapply(reps[[cp]], function(r)
      as.vector(.qast_surv_matrix(r$dist, r$coef, z, tm$ht, cp)),
      numeric(length(est)))
    vals <- matrix(vals, nrow = length(est))
    se <- matrix(.qast_rep_se(object, vals), nrow = nrow(est),
                 dimnames = dn)
    if (identical(object$var_method, "jackknife")) {
      q <- stats::qnorm(1 - (1 - level) / 2)
      lower <- est - q * se; upper <- est + q * se
    } else {
      a <- (1 - level) / 2
      lower <- matrix(apply(vals, 1L, stats::quantile, probs = a,
                            names = FALSE), nrow = nrow(est))
      upper <- matrix(apply(vals, 1L, stats::quantile, probs = 1 - a,
                            names = FALSE), nrow = nrow(est))
    }
    lower <- pmin(pmax(lower, 0), 1); upper <- pmin(pmax(upper, 0), 1)
    dimnames(lower) <- dimnames(upper) <- dn
    list(surv = est, se = se, lower = lower, upper = upper)
  })
  if (length(comps) == 1L) out[[1L]] else out
}
