#' Zhao and Tsiatis estimator of the quality-adjusted survival curve
#'
#' Estimates the survival function of the quality-adjusted survival time
#' with the inverse probability of censoring weighted (IPCW) estimator of
#' Zhao and Tsiatis (1997). Each uncensored subject is weighted by the
#' inverse of the censoring survival at its survival time, which corrects
#' the induced informative censoring on the quality-adjusted scale:
#'
#'   \code{pred_surv[i] = (1/n) * sum_\{j : qalys[j] >= qalys[i]\} status[j] / G(t_j)}
#'
#' where G is the Kaplan-Meier survival of the censoring distribution and
#' \code{qalys = log(qaly) - lin_pred} when a linear predictor is given.
#'
#' @param qaly Quality-adjusted survival times (QALY), on the raw
#'   (positive) scale
#' @param time Observed survival times
#' @param status Event indicator (1 = event, 0 = censored)
#' @param lin_pred Optional linear predictor (\code{X \%*\% theta}). When
#'   given, the transformation \code{log(qaly) - lin_pred} is applied. Must
#'   have the same length as \code{qaly}
#' @param km_fit Optional \code{survfit} object with the censoring KM
#'   curve, precomputed (e.g. reused across calls). If NULL, it is fitted
#'   as \code{Surv(time, 1 - status) ~ 1}
#' @return Numeric vector of length n with the estimated survival
#'   probabilities
#' @references Zhao, H. and Tsiatis, A. A. (1997). A consistent estimator
#'   for the distribution of quality adjusted survival time.
#'   \emph{Biometrika}, 84(2), 339-348. \doi{10.1093/biomet/84.2.339}
#' @seealso \code{\link{semi_aft_fg}}
#' @examples
#' set.seed(1)
#' n <- 50
#' time   <- rexp(n)
#' status <- rbinom(n, 1, 0.8)
#' qaly   <- time * runif(n, 0.5, 1)
#'
#' # Curve on the raw scale
#' head(surv_qt_fun(qaly, time, status))
#'
#' # With a linear predictor (e.g. X \%*\% theta from a fit)
#' lp <- rnorm(n, sd = 0.1)
#' head(surv_qt_fun(qaly, time, status, lin_pred = lp))
#' @export
surv_qt_fun <- function(qaly, time, status, lin_pred = NULL, km_fit = NULL) {
  
  # --- Validacoes -----------------------------------------------------------
  n <- length(qaly)
  if (length(time) != n) {
    stop("time must have the same length as qaly")
  }
  if (length(status) != n) {
    stop("status must have the same length as qaly")
  }
  if (!is.null(lin_pred) && length(lin_pred) != n) {
    stop("lin_pred must have the same length as qaly")
  }
  if (!is.null(lin_pred) && any(qaly <= 0)) {
    stop("qaly must be strictly positive when lin_pred is given ",
         "(log undefined)")
  }
  
  # --- Sobrevivencia da censura (Kaplan-Meier) -------------------------------
  # survfit em R; o kernel C++ avalia a step function por busca binaria
  if (is.null(km_fit)) {
    km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  }
  
  # --- Chamada ao kernel C++ -------------------------------------------------
  # A transformacao log(qaly) - lin_pred acontece dentro do kernel.
  surv_qt_with_km_cpp(
    qaly     = qaly,
    time     = time,
    status   = status,
    km_time  = km_fit$time,
    km_surv  = km_fit$surv,
    lin_pred = lin_pred
  )
}