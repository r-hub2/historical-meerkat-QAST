#' Kaplan-Meier censoring survival at each subject's observed time
#'
#' The censoring KM curve is estimated by \code{survival::survfit} and the
#' step function is evaluated at each \code{t[i]} with \code{findInterval}
#' (O(n log m) per group).
#'
#' Conventions (checked by the tests):
#'   - empty group -> numeric(0)
#'   - group without events -> rep(1, n)
#'   - fit without valid times/survival -> rep(1, n)
#'   - t before the first jump -> 1
#'   - t beyond the last jump -> last value of the curve
#'
#' @param t Observed times of the group
#' @param status Event indicator (1 = event)
#' @param use_legacy If TRUE, evaluates with \code{which()} (reference
#'   implementation for the A/B tests)
#' @return Vector of length \code{length(t)} with \code{S_cens(t_i)}
#' @keywords internal
km_weights_at <- function(t, status, use_legacy = FALSE) {
  if (use_legacy) return(.km_weights_at_legacy(t, status))
  .km_weights_at_fast(t, status)
}

#' Censoring KM evaluation: survfit + findInterval
#' @keywords internal
.km_weights_at_fast <- function(t, status) {
  n <- length(t)
  if (n == 0L) return(numeric(0))
  
  # Group without events: weights 1
  if (sum(status) == 0L) return(rep(1, n))
  
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  if (length(fit$time) == 0L || all(is.na(fit$surv))) {
    return(rep(1, n))
  }
  
  # findInterval(x, v) = number of elements of v that are <= x
  # (left.open = FALSE, rightmost.closed = FALSE) == max(which(v <= x))
  idx <- findInterval(t, fit$time)
  
  # idx == 0 (before the first jump): S = 1. Otherwise: surv[idx].
  w <- fit$surv[pmax(idx, 1L)]
  w[idx == 0L] <- 1
  w
}

#' Reference censoring KM evaluation with which(), O(n m)
#'
#' Used only by the A/B tests of km_weights_at; keep unchanged.
#' @keywords internal
.km_weights_at_legacy <- function(t, status) {
  n <- length(t)
  if (n == 0L) return(numeric(0))
  
  # Group without events: weights 1
  if (sum(status) == 0L) return(rep(1, n))
  
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  if (length(fit$time) == 0L || all(is.na(fit$surv))) {
    return(rep(1, n))
  }
  
  times <- fit$time
  surv <- fit$surv
  
  vapply(t, function(ti) {
    idx <- which(times <= ti)
    if (length(idx) == 0L) {
      1
    } else {
      surv[max(idx)]
    }
  }, numeric(1))
}

#' Censoring partition of the Fine and Gelber model
#'
#' Fine and Gelber (2001) assume that the censoring depends on a discrete
#' partition G of the covariates, chosen by the analyst. As in their
#' application, the default is no stratification (a single group);
#' \code{group} sets another partition.
#'
#' @param z Binary covariate matrix (n x p)
#' @param group NULL (a single group) or a vector of length n
#' @return Factor of length n
#' @keywords internal
fg_cens_partition <- function(z, group = NULL) {
  if (is.null(group)) return(factor(rep("all", nrow(z))))
  factor(group)
}

#' Censoring KM survival of each subject in its censoring group
#'
#' \eqn{\hat H_{G_i}(X_i)}{H_G(X_i)}: the censoring Kaplan-Meier curve of
#' the group of subject i, evaluated at its observed time.
#'
#' @param t Observed times
#' @param status Event indicator
#' @param group Factor with the censoring groups
#' @param use_legacy Passed to km_weights_at (reference for the tests)
#' @return Numeric vector of length n
#' @keywords internal
km_weights_by_group <- function(t, status, group, use_legacy = FALSE) {
  H <- numeric(length(t))
  for (g in levels(group)) {
    i <- which(group == g)
    if (length(i)) H[i] <- km_weights_at(t[i], status[i],
                                         use_legacy = use_legacy)
  }
  H
}

#' Censoring weights of the Fine and Gelber pairs, by column
#'
#' Each subject has one censoring weight \eqn{\hat H_{G_i}(X_i)}{H_G(X_i)}
#' from the censoring partition (Fine and Gelber, 2001), the same in every
#' component of the estimating function. The kernels read it split by the
#' groups of each binary column: \code{w0} (rows with \code{z[, j] == 0})
#' and \code{w1} (rows with \code{z[, j] == 1}), in the original order.
#' Recomputed in every variance replicate.
#'
#' @param z Binary covariate matrix (n x p)
#' @param t Observed times
#' @param status Event indicator
#' @param use_legacy Passed to km_weights_at (reference for the tests)
#' @param group Censoring partition: NULL (default; a single group, no
#'   stratification) or a vector of length n
#' @return List of length p; element j holds w0 and w1
#' @keywords internal
prep_km_by_column <- function(z, t, status, use_legacy = FALSE,
                              group = NULL) {
  H <- km_weights_by_group(t, status, fg_cens_partition(z, group),
                           use_legacy = use_legacy)
  lapply(seq_len(ncol(z)), function(j) {
    x <- z[, j]
    list(w0 = H[x == 0], w1 = H[x == 1])
  })
}
