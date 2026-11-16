#' Censoring weights via a Cox model
#'
#' Survival probabilities of the censoring distribution from
#' \code{survival::coxph}, used as IPCW weights in the proposed model.
#'
#' @param X Covariate matrix (n x p)
#' @param t Observed times
#' @param status Event indicator (1 = event, 0 = censored)
#' @return Numeric vector of length n with the censoring survival
#'   probabilities
#' @importFrom survival coxph
#' @importFrom stats predict
#' @keywords internal
HG <- function(X, t, status) {
  cox_fit <- coxph(Surv(t, 1 - status) ~ X)
  cox_pred <- exp(-predict(object = cox_fit, type = 'expected'))
  return(cox_pred)
}