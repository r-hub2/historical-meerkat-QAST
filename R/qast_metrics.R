# R/qast_metrics.R
#
# Predictive performance of survival curves for the survival time (T) and
# the quality-adjusted time (Q): Brier score with inverse probability of
# censoring weights (IPCW), integrated Brier score, explained variation
# against the null model and IPCW mean absolute error.
#
# Censoring of Q is induced informative (Q(C) depends on the quality
# path), so the usual weights from the censoring KM of Q are biased. All
# weights here use the censoring KM G on the SURVIVAL time scale, as the
# estimator of Zhao and Tsiatis (1997): each uncensored subject has
# weight status / G(T_i), with expectation 1 when C is independent of
# (T, Q).

# Censoring KM on the survival time scale: G(x) = P(C > x), or the left
# limit G(x-) with left = TRUE
.qast_cens_km <- function(time, status) {
  if (all(status == 1)) return(function(x, left = FALSE) rep(1, length(x)))
  f <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  function(x, left = FALSE) {
    idx <- findInterval(x, f$time, left.open = left)
    out <- f$surv[pmax(idx, 1L)]
    out[idx == 0L] <- 1
    out
  }
}

# 1 / G with 0 where G = 0 (no information left)
.qast_inv <- function(g) ifelse(g > 0, 1 / g, 0)

# Evaluation data: matrix with columns time, status (and qtime)
.qast_check_y <- function(y, n, need_q) {
  y <- as.matrix(unclass(y))
  need <- c("time", "status", if (need_q) "qtime")
  miss <- setdiff(need, colnames(y))
  if (length(miss))
    stop("y must have the columns ", paste(need, collapse = ", "),
         " (e.g. Qsurv(time, qtime, status)); missing: ",
         paste(miss, collapse = ", "), call. = FALSE)
  if (nrow(y) != n)
    stop(sprintf("y has %d rows; the predictions have %d", nrow(y), n),
         call. = FALSE)
  y <- y[, need, drop = FALSE]
  if (anyNA(y)) .qast_stop_at("y contains NA", which(!stats::complete.cases(y)))
  if (!all(y[, "status"] %in% c(0, 1)))
    stop("status must be 0 (censored) or 1 (event)", call. = FALSE)
  if (sum(y[, "status"]) == 0) stop("y has no events", call. = FALSE)
  if (any(y[, "time"] <= 0)) stop("time must be positive", call. = FALSE)
  if (need_q && any(y[, "qtime"] < 0))
    stop("qtime must be non-negative", call. = FALSE)
  y
}

# Trapezoidal integral of f over the grid, divided by its range
.qast_integrate <- function(times, f) {
  k <- length(times)
  if (k == 1L) return(f)
  sum(diff(times) * (f[-1L] + f[-k]) / 2) / (times[k] - times[1L])
}

# Core: surv (n x length(times)), y (time, status[, qtime]), outcome
.qast_metrics_core <- function(surv, times, y, outcome, pred_time, tau) {
  n <- nrow(y)
  time <- y[, "time"]
  status <- y[, "status"]
  obs <- y[, outcome]
  G <- .qast_cens_km(time, status)
  ind <- outer(obs, times, ">")                 # I(Y_i > t)

  if (outcome == "time") {
    # Graf et al. (1999): events up to t weighted by 1/G(T_i-), subjects
    # still at risk by 1/G(t)
    W <- (!ind & status == 1) * .qast_inv(G(time, left = TRUE)) +
      ind * matrix(.qast_inv(G(times)), n, length(times), byrow = TRUE)
    null_surv <- .surv_km_at(.resid_dist_km(time, status), times)
  } else {
    # Only the uncensored contribute, with weight 1/G(T_i) on the survival
    # time scale: whether Q_i > q is known for a censored subject with
    # Q_i(C) > q, but the probability of that is not, so these subjects
    # cannot be reweighted
    wq <- status * .qast_inv(G(time))
    W <- matrix(wq, n, length(times))
    null_surv <- .surv_ipcw_at(.resid_dist_ipcw(obs, time, status, G(time)),
                               times)
  }
  null_mat <- matrix(null_surv, n, length(times), byrow = TRUE)
  brier <- colMeans(W * (ind - surv)^2)
  brier_null <- colMeans(W * (ind - null_mat)^2)
  r2 <- ifelse(brier_null > 0, 1 - brier / brier_null, NA_real_)

  ibs <- .qast_integrate(times, brier)
  ibs_null <- .qast_integrate(times, brier_null)

  # IPCW mean absolute error of min(Y, tau), Hajek weights (only with
  # predicted restricted means)
  w <- status * .qast_inv(G(time))
  mae <- if (is.null(pred_time)) NA_real_ else
    sum(w * abs(pmin(obs, tau) - pred_time)) / sum(w)
  # Total IPCW mass of the uncensored, sum(w) / n: close to 1 when the
  # censoring support covers the survival times (G(T) > 0); clearly below
  # 1, subjects with long survival are never observed and the metrics of
  # the quality-adjusted time are biased at large times
  ipcw_mass <- sum(w) / n

  structure(list(
    outcome = outcome,
    brier = data.frame(time = times, brier = brier, brier_null = brier_null,
                       r2 = r2),
    ibs = ibs, ibs_null = ibs_null,
    r2_ibs = if (ibs_null > 0) 1 - ibs / ibs_null else NA_real_,
    mae = mae, tau = tau, n = n, n_events = as.integer(sum(status)),
    ipcw_mass = ipcw_mass,
    pred_time = pred_time
  ), class = "qast_metrics")
}

#' Predictive performance for survival and quality-adjusted survival
#'
#' Brier score, integrated Brier score, explained variation against the
#' null model and mean absolute error of the predicted curves of the
#' survival time (\eqn{T}) or of the quality-adjusted survival time
#' (\eqn{Q}), with inverse probability of censoring weights (IPCW) that are
#' valid for \eqn{Q}.
#'
#' The censoring of \eqn{Q} is induced informative: the quality-adjusted
#' time accumulated up to the censoring, \eqn{Q(C)}, depends on the
#' quality of life of the subject. The usual metrics, which treat
#' \code{(qtime, status)} as right-censored data and weight by the
#' Kaplan-Meier censoring curve of \eqn{Q}, are therefore biased. Here,
#' as in the estimator of Zhao and Tsiatis (1997), every weight uses the
#' Kaplan-Meier censoring survival \eqn{G} on the \emph{survival time}
#' scale, which assumes censoring independent of \eqn{(T, Q)} and of the
#' covariates:
#' \itemize{
#'   \item Survival time (\code{outcome = "time"}): the Brier score of
#'     Graf et al. (1999),
#'     \deqn{BS(t) = \frac{1}{n}\sum_i \left[\frac{I(T_i \le t,
#'     \delta_i = 1)}{G(T_i-)}\hat S_i(t)^2 + \frac{I(T_i > t)}{G(t)}
#'     \{1 - \hat S_i(t)\}^2\right].}{BS(t) = (1/n) sum_i [I(T_i <= t,
#'     delta_i = 1) S_i(t)^2 / G(T_i-) + I(T_i > t) (1 - S_i(t))^2 / G(t)].}
#'   \item Quality-adjusted time (\code{outcome = "qtime"}): only the
#'     uncensored subjects contribute,
#'     \deqn{BS(q) = \frac{1}{n}\sum_i \frac{\delta_i}{G(T_i)}\{I(Q_i > q) -
#'     \hat S_i(q)\}^2.}{BS(q) = (1/n) sum_i delta_i / G(T_i) (I(Q_i > q) -
#'     S_i(q))^2.}
#'     A censored subject with \eqn{Q_i(C) > q} is known to have
#'     \eqn{Q_i > q}, but the probability of that observation depends on
#'     its quality path, so it cannot be reweighted.
#' }
#' The null model is the curve without covariates: the Kaplan-Meier curve
#' of \eqn{T}, or the Zhao and Tsiatis curve of \eqn{Q}. \code{r2} is
#' \eqn{1 - BS/BS_{null}}{1 - BS / BS_null} at each time and \code{r2_ibs}
#' the same for the integrated scores. The integrated Brier score is the
#' trapezoidal integral over the grid \code{times}, divided by its range.
#'
#' The weights require positivity: every survival time must have a
#' positive probability of being observed, \eqn{G(T) > 0}{G(T) > 0}. When the
#' follow-up ends before the longest survival times, those subjects are
#' never observed; for the survival time the Graf weights still use the
#' subjects at risk, but for the quality-adjusted time the Brier score is
#' then underestimated at large times. \code{ipcw_mass},
#' \eqn{n^{-1}\sum_i \delta_i/G(T_i)}{(1/n) sum_i delta_i / G(T_i)}, has
#' expectation 1 under positivity; values clearly below 1 flag the problem
#' (choose a smaller \code{tau} or grid).
#'
#' The mean absolute error compares the predicted restricted mean with
#' \eqn{\min(Y_i, \tau)}{min(Y_i, tau)}, \eqn{Y = T} or \eqn{Q}, over the
#' uncensored subjects with weights \eqn{\delta_i/G(T_i)}{delta_i / G(T_i)}
#' normalized to sum 1. Restricting at \eqn{\tau} avoids the right tail of
#' the curves, which is not identified when the largest times are
#' censored.
#'
#' @param object A \code{"qast_fit"}, or a matrix of predicted survival
#'   probabilities with one row per subject of \code{y} and one column per
#'   time of \code{times} (any model).
#' @param ... Arguments of the methods.
#' @return An object of class \code{"qast_metrics"}: a list with
#'   \code{brier} (data frame with \code{time}, \code{brier},
#'   \code{brier_null} and \code{r2}), \code{ibs}, \code{ibs_null},
#'   \code{r2_ibs}, \code{mae}, \code{tau}, \code{n}, \code{n_events},
#'   \code{ipcw_mass} (see Details) and \code{pred_time} (predicted
#'   restricted means, \code{NULL} if not given). For a \code{"qast_fit"}
#'   with \code{component = "both"}, a list of two such objects
#'   (\code{beta} and \code{theta}) of class \code{"qast_metrics_list"}.
#' @references Graf, E., Schmoor, C., Sauerbrei, W. and Schumacher, M.
#'   (1999). Assessment and comparison of prognostic classification
#'   schemes for survival data. \emph{Statistics in Medicine}, 18(17-18),
#'   2529-2545.
#'   \doi{10.1002/(SICI)1097-0258(19990915/30)18:17/18<2529::AID-SIM274>3.0.CO;2-5}
#'
#'   Zhao, H. and Tsiatis, A. A. (1997). A consistent estimator for the
#'   distribution of quality adjusted survival time. \emph{Biometrika},
#'   84(2), 339-348. \doi{10.1093/biomet/84.2.339}
#' @seealso \code{\link{predict.qast_fit}} (\code{type = "survival"} and
#'   \code{"rmean"}), \code{\link{surv_qt_fun}}
#' @examples
#' set.seed(1)
#' n <- 120
#' d <- data.frame(age = rnorm(n), trt = rbinom(n, 1, 0.5))
#' d$t <- rexp(n, rate = exp(-0.5 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#' train <- d[1:80, ]
#' test <- d[81:120, ]
#' fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = train)
#'
#' # Out-of-sample performance, survival and quality-adjusted time
#' m <- qast_metrics(fit, newdata = test,
#'                   y = with(test, Qsurv(t, qt, status)))
#' m
#' head(m$theta$brier)
#'
#' # Any model, from its survival matrix: here, the Kaplan-Meier curve
#' # without covariates (the null model, explained variation 0)
#' library(survival)
#' km <- survfit(Surv(t, status) ~ 1, data = train)
#' tt <- sort(unique(test$t[test$status == 1]))
#' s_km <- matrix(summary(km, times = tt, extend = TRUE)$surv,
#'                nrow(test), length(tt), byrow = TRUE)
#' qast_metrics(s_km, y = with(test, Qsurv(t, qt, status)), times = tt)
#' @export
qast_metrics <- function(object, ...) UseMethod("qast_metrics")

#' @rdname qast_metrics
#' @param y Evaluation data: a \code{Qsurv(time, qtime, status)} object or
#'   a matrix with columns \code{time}, \code{status} and, for the
#'   quality-adjusted time, \code{qtime}. For a \code{"qast_fit"}, required
#'   with \code{newdata}; \code{NULL} uses the data of the fit (apparent
#'   performance).
#' @param outcome \code{"time"} (survival time) or \code{"qtime"}
#'   (quality-adjusted time): what the columns of the matrix predict
#' @param times Evaluation times, on the original scale. For a matrix,
#'   the times of its columns (required). For a \code{"qast_fit"},
#'   \code{NULL} (default) uses the distinct event times of \code{y} up to
#'   \code{tau}; a vector applies to both components, a list named
#'   \code{beta} and \code{theta} gives one grid per component.
#' @param tau Restriction time of the mean absolute error (and of the
#'   default grid): a number, or for a \code{"qast_fit"} a vector named by
#'   component. Default: the largest time of \code{times} (matrix) or the
#'   largest event time of \code{y} (\code{"qast_fit"}).
#' @param pred_time Predicted restricted means \eqn{E\{\min(Y, \tau)\}}{E[min(Y,
#'   tau)]}, one per subject, for the mean absolute error. \code{NULL}
#'   (default): no mean absolute error (\code{mae = NA}).
#' @export
qast_metrics.default <- function(object, y, outcome = c("time", "qtime"),
                                 times, tau = NULL, pred_time = NULL, ...) {
  outcome <- match.arg(outcome)
  surv <- as.matrix(object)
  if (!is.numeric(surv)) stop("object must be a numeric matrix of survival ",
                              "probabilities", call. = FALSE)
  if (missing(times) || !is.numeric(times) || anyNA(times) ||
      any(times < 0))
    stop("times (the non-negative times of the columns) are required",
         call. = FALSE)
  if (ncol(surv) != length(times))
    stop(sprintf("object has %d columns; times has %d", ncol(surv),
                 length(times)), call. = FALSE)
  if (anyNA(surv) || any(surv < 0 | surv > 1))
    stop("the survival probabilities must be in [0, 1], without NA",
         call. = FALSE)
  if (anyDuplicated(times)) stop("times must be distinct", call. = FALSE)
  o <- order(times)
  times <- times[o]
  surv <- surv[, o, drop = FALSE]
  y <- .qast_check_y(y, nrow(surv), outcome == "qtime")
  if (is.null(tau)) tau <- max(times)
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0)
    stop("tau must be a positive number", call. = FALSE)
  if (!is.null(pred_time) && (!is.numeric(pred_time) || length(pred_time) != nrow(surv) ||
             anyNA(pred_time))) {
    stop("pred_time must be numeric, one value per subject", call. = FALSE)
  }
  .qast_metrics_core(surv, times, y, outcome, pred_time, tau)
}

#' @rdname qast_metrics
#' @param newdata Covariates of the evaluation subjects, as in
#'   \code{\link{predict.qast_fit}}; \code{NULL} uses the data of the fit
#' @param component \code{"both"}, \code{"beta"} (survival time) or
#'   \code{"theta"} (quality-adjusted time)
#' @export
qast_metrics.qast_fit <- function(object, newdata = NULL, y = NULL,
                                  component = c("both", "beta", "theta"),
                                  times = NULL, tau = NULL, ...) {
  component <- match.arg(component)
  if (is.null(newdata)) {
    if (!is.null(y))
      stop("y without newdata: pass the covariates of y in newdata",
           call. = FALSE)
    y <- object$y
  } else if (is.null(y)) {
    stop("with newdata, y = Qsurv(time, qtime, status) of the same ",
         "subjects is required", call. = FALSE)
  }
  n_new <- if (is.null(newdata)) object$n else
    nrow(.qast_newdata_z(object, newdata))
  y <- .qast_check_y(y, n_new, TRUE)
  comps <- if (component == "both") c("beta", "theta") else component

  out <- lapply(stats::setNames(comps, comps), function(cp) {
    outcome <- if (cp == "beta") "time" else "qtime"
    ev <- sort(unique(y[y[, "status"] == 1, outcome]))
    tau_c <- if (is.null(tau)) max(ev) else {
      tmp <- list(event_times = stats::setNames(list(ev), cp))
      .qast_tau(tmp, cp, tau)
    }
    tm <- if (is.list(times)) times[[cp]] else times
    if (is.null(tm)) tm <- ev[ev <= tau_c]
    if (!length(tm)) stop(sprintf("no evaluation times up to tau for %s", cp),
                          call. = FALSE)
    tm <- sort(unique(tm))
    s <- predict(object, newdata, type = "survival", component = cp,
                 times = tm)
    pt <- predict(object, newdata, type = "rmean", component = cp,
                  tau = tau_c)
    .qast_metrics_core(unname(s), tm, y, outcome, as.vector(pt), tau_c)
  })
  if (length(out) == 1L) return(out[[1L]])
  class(out) <- "qast_metrics_list"
  out
}

#' @rdname qast_metrics
#' @param x Object to print
#' @param digits Significant digits
#' @export
print.qast_metrics <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  lab <- if (x$outcome == "time") "survival time" else
    "quality-adjusted survival time"
  cat(sprintf("Predictive performance: %s (n = %d, %d events)\n", lab, x$n,
              x$n_events))
  b <- x$brier
  cat(sprintf("  Integrated Brier score over [%s, %s]: %s (null model %s)\n",
              format(min(b$time), digits = digits),
              format(max(b$time), digits = digits),
              format(x$ibs, digits = digits),
              format(x$ibs_null, digits = digits)))
  cat(sprintf("  Explained variation (1 - IBS / IBS null): %s\n",
              format(x$r2_ibs, digits = digits)))
  if (!is.na(x$mae))
    cat(sprintf("  Mean absolute error of min(Y, %s): %s\n",
                format(x$tau, digits = digits),
                format(x$mae, digits = digits)))
  cat(sprintf(paste0("  IPCW weights: censoring Kaplan-Meier on the survival ",
                     "time scale (total mass %s)\n"),
              format(x$ipcw_mass, digits = digits)))
  invisible(x)
}

#' @rdname qast_metrics
#' @export
print.qast_metrics_list <- function(x, ...) {
  for (k in seq_along(x)) {
    if (k > 1L) cat("\n")
    print(x[[k]], ...)
  }
  invisible(x)
}
