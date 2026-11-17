# R/qast_pp.R
#
# Lack-of-fit diagnostic by P-P plots of the residuals (Fine and Gelber,
# 2001, Section 4): under the model, the distribution of the residuals
# does not depend on the covariates, so the distributions estimated in two
# covariate groups must form a line from (0, 0) to (1, 1).

#' P-P plots of the residuals by covariate group
#'
#' Graphical lack-of-fit diagnostic of Fine and Gelber (2001, Section 4).
#' Under the model, the distribution of the residuals does not depend on
#' the covariates, so the distributions estimated in two covariate groups
#' should give a P-P plot close to the line from (0, 0) to (1, 1).
#'
#' With residuals \eqn{\hat\epsilon_{1j} = h_1(X_j) -
#' \hat\beta^\top Z_j}{e1_j = h1(X_j) - beta'Z_j} and
#' \eqn{\hat\epsilon_{2j} = h_2(Q_j) - \hat\theta^\top Z_j}{e2_j = h2(Q_j) -
#' theta'Z_j}, for each group \eqn{s}:
#' \itemize{
#'   \item \code{beta} (survival time): the distribution of
#'     \eqn{\epsilon_1}{e1} is estimated by the Kaplan-Meier estimator
#'     applied to the pairs \eqn{(\hat\epsilon_{1j}, \Delta_j)}{(e1_j,
#'     Delta_j)} of the group.
#'   \item \code{theta} (quality-adjusted time): the inverse probability of
#'     censoring weighted distribution of \eqn{\hat\epsilon_{2j}}{e2_j}
#'     among the events of the group, normalized by the total weight of
#'     the group:
#'     \deqn{\hat F_s(t) = \frac{\sum_j \Delta_j \hat H(X_j)^{-1}
#'       I\{A_j, Z_j \in S_s, \hat\epsilon_{2j} < t\}}{\sum_j \Delta_j
#'       \hat H(X_j)^{-1} I\{A_j, Z_j \in S_s\}}.}{F_s(t) = sum_j
#'       Delta_j / H(X_j) I(A_j, j in group s, e2_j < t) / sum_j
#'       Delta_j / H(X_j) I(A_j, j in group s).}
#'     For \code{semi_aft_fg()} fits, as in Fine and Gelber (2001),
#'     \eqn{A_j}{A_j} is the trimming
#'     \eqn{\hat\epsilon_{1j} < h_1(L) - \max_k \hat\beta^\top Z_k}{e1_j <
#'     h1(L) - max_k beta'Z_k}, so \eqn{\hat F_s}{F_s} estimates the
#'     conditional distribution
#'     \eqn{\Pr\{\epsilon_2 < t \mid \epsilon_1 < h_1(L) -
#'     \max_Z(\beta^\top Z)\}}{Pr(e2 < t | e1 < h1(L) - max beta'Z)}
#'     required by the estimating function of \eqn{\theta}{theta}, and
#'     \eqn{\hat H = \hat H_{G_j}(X_j)}{H = H_G(X_j)} is the Kaplan-Meier
#'     censoring survival in the censoring group of the fit
#'     (\code{cens_group} of \code{\link{semi_aft_fg}}; the marginal
#'     Kaplan-Meier without stratification, the default). For
#'     \code{semi_aft_prop()} fits, which have no truncation time, there
#'     is no trimming (\eqn{A_j}{A_j} always holds) and
#'     \eqn{\hat H}{H} are the censoring weights of the fit: the Cox model
#'     of the censoring given the covariates when \code{cens_weight = TRUE},
#'     the Kaplan-Meier within the group otherwise.
#' }
#' Fine and Gelber (2001) recommend checking the model for the survival
#' time first: the quality-adjusted diagnostic assumes it is correct.
#' Each group is compared with the reference group \code{ref} on a common
#' grid of residual values (the event residuals of the two groups; for
#' \code{beta}, up to the smallest of the two largest residuals, beyond
#' which a Kaplan-Meier curve is undefined). As in the article, the
#' quality-adjusted distribution uses the strict indicator
#' \eqn{\hat\epsilon_{2j} < t}{e2_j < t}, so at each grid value it is the
#' left limit; the last point of the grid closes the curve at (1, 1). Continuous covariates must be
#' discretized to form the groups (e.g. with \code{cut()}).
#'
#' @param fit Object of class \code{"qast_fit"}
#' @param group Grouping of the observations of the fit: a vector (factor,
#'   character, numeric or logical) of length \code{fit$n}
#' @param ref Reference level; \code{NULL} (default) uses the first level
#' @param x Object of class \code{"qast_pp"}
#' @param component \code{"both"}, \code{"beta"} or \code{"theta"}
#' @param ... Passed to \code{plot()} (\code{plot} method); not used
#'   (\code{print} method)
#' @return Object of class \code{"qast_pp"}: a list with \code{beta} and
#'   \code{theta}, each a list of data frames (one per group compared with
#'   \code{ref}) with the residual grid \code{t} and the estimated
#'   distributions \code{ref} (x axis) and \code{group} (y axis); plus
#'   \code{levels}, \code{ref}, \code{model}, \code{trim} (the trimming
#'   bound, \code{NULL} without trimming) and \code{n_used} (events per
#'   group in the \code{theta} distribution).
#' @references Fine, J. P. and Gelber, R. D. (2001). Joint regression
#'   analysis of survival and quality-adjusted survival. \emph{Biometrics},
#'   57(2), 376-382. \doi{10.1111/j.0006-341X.2001.00376.x}
#' @seealso \code{\link{semi_aft_fg}}, \code{\link{semi_aft_prop}}
#' @examples
#' set.seed(2)
#' n <- 200
#' d <- data.frame(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
#' d$t <- rexp(n, rate = exp(-0.4 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#' fit <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                    L = quantile(d$t, 0.9))
#' pp <- qast_pp(fit, group = d$trt)
#' pp
#' plot(pp)
#' @export
qast_pp <- function(fit, group, ref = NULL) {
  if (!inherits(fit, "qast_fit"))
    stop("fit must be a qast_fit object", call. = FALSE)
  if (is.null(fit$y))
    stop("this fit does not store its data; refit it with the current ",
         "version of QAST", call. = FALSE)
  n <- fit$n
  if (length(group) != n)
    stop(sprintf("group must have one value per observation of the fit (%d)",
                 n), call. = FALSE)
  if (anyNA(group)) .qast_stop_at("group contains NA", which(is.na(group)))
  g <- factor(group)
  lev <- levels(g)
  if (length(lev) < 2L)
    stop("group must have at least 2 levels", call. = FALSE)
  if (is.null(ref)) ref <- lev[1L]
  ref <- as.character(ref)
  if (!ref %in% lev)
    stop("ref must be one of the levels of group: ",
         paste(lev, collapse = ", "), call. = FALSE)

  t <- fit$y[, "time"]; qt <- fit$y[, "qtime"]
  status <- fit$y[, "status"]
  z <- fit$z
  h1 <- stats::make.link(fit$settings$h1)
  h2 <- stats::make.link(fit$settings$h2)
  e1 <- h1$linkfun(t) - drop(z %*% fit$beta)
  e2 <- h2$linkfun(qt) - drop(z %*% fit$theta)

  # Fine and Gelber trimming: e1 < h1(L) - max_k beta'Z_k
  trim <- NULL
  keep <- rep(TRUE, n)
  if (fit$model == "fg") {
    trim <- h1$linkfun(fit$settings$L) - max(drop(z %*% fit$beta))
    keep <- e1 < trim
  }

  # Censoring survival at X_j
  H <- numeric(n)
  if (fit$model == "prop" && isTRUE(fit$settings$cens_weight)) {
    H <- HG(X = z, t = t, status = status)
  } else if (fit$model == "fg" && !is.null(fit$settings$cens_group)) {
    # H_{G_j}(X_j) from the censoring partition of the fit, as in the article
    H <- km_weights_by_group(t, status, fit$settings$cens_group)
  } else {
    for (l in lev) {
      i <- g == l
      H[i] <- km_weights_at(t[i], status[i])
    }
  }
  w <- ifelse(status == 1 & keep & is.finite(H) & H > 0, 1 / H, 0)

  for (l in lev) {
    i <- g == l
    if (!any(status[i] == 1))
      stop(sprintf("group \"%s\" has no events", l), call. = FALSE)
    if (sum(w[i]) == 0)
      stop(sprintf(paste0("group \"%s\" has no events left for the ",
                          "quality-adjusted distribution%s"), l,
                   if (!is.null(trim)) " after the trimming (increase L?)"
                   else ""), call. = FALSE)
  }

  # Distributions per group (distribution functions)
  cdf1 <- lapply(stats::setNames(lev, lev), function(l) {
    i <- g == l
    d <- .resid_dist_km(e1[i], status[i])
    list(fun = function(x) 1 - .surv_km_at(d, x),
         jumps = sort(unique(e1[i & status == 1])), max = max(e1[i]))
  })
  cdf2 <- lapply(stats::setNames(lev, lev), function(l) {
    i <- g == l & w > 0
    v <- e2[i]; wi <- w[i] / sum(w[i])
    o <- order(v); v <- v[o]; cw <- cumsum(wi[o])
    # Pr(e2 < t), with the strict indicator of the article: sum of v < x
    list(fun = function(x) {
      k <- findInterval(x, v, left.open = TRUE)
      c(0, cw)[k + 1L]
    }, jumps = sort(unique(v)))
  })

  others <- setdiff(lev, ref)
  pp1 <- lapply(stats::setNames(others, others), function(l) {
    tt <- sort(unique(c(cdf1[[ref]]$jumps, cdf1[[l]]$jumps)))
    tt <- tt[tt <= min(cdf1[[ref]]$max, cdf1[[l]]$max)]
    data.frame(t = c(-Inf, tt), ref = c(0, cdf1[[ref]]$fun(tt)),
               group = c(0, cdf1[[l]]$fun(tt)))
  })
  pp2 <- lapply(stats::setNames(others, others), function(l) {
    # With <, the value at a jump is the left limit; Inf closes at (1, 1)
    tt <- c(sort(unique(c(cdf2[[ref]]$jumps, cdf2[[l]]$jumps))), Inf)
    data.frame(t = tt, ref = cdf2[[ref]]$fun(tt),
               group = cdf2[[l]]$fun(tt))
  })

  out <- list(beta = pp1, theta = pp2, levels = lev, ref = ref,
              model = fit$model, trim = trim,
              n_used = vapply(lev, function(l) sum(g == l & w > 0),
                              integer(1)))
  class(out) <- "qast_pp"
  out
}

#' @rdname qast_pp
#' @export
print.qast_pp <- function(x, ...) {
  cat("\nP-P plots of the residuals (Fine and Gelber, 2001)\n")
  cat("Groups: ", paste(x$levels, collapse = ", "), " (reference: ", x$ref,
      ")\n", sep = "")
  if (!is.null(x$trim))
    cat("Quality-adjusted residuals restricted to e1 < h1(L) - max(beta'Z) = ",
        format(x$trim, digits = 4), "\n", sep = "")
  cat("Events in the quality-adjusted distribution: ",
      paste(names(x$n_used), x$n_used, sep = " = ", collapse = ", "), "\n",
      sep = "")
  dev <- function(lst) vapply(lst, function(d) max(abs(d$group - d$ref)), 0)
  cat("Largest vertical distance to the diagonal:\n")
  print(rbind(beta = dev(x$beta), theta = dev(x$theta)), digits = 3)
  cat("Use plot() to draw them.\n")
  invisible(x)
}

#' @rdname qast_pp
#' @export
plot.qast_pp <- function(x, component = c("both", "beta", "theta"), ...) {
  component <- match.arg(component)
  comps <- if (component == "both") c("beta", "theta") else component
  if (length(comps) == 2L) {
    op <- graphics::par(mfrow = c(1L, 2L))
    on.exit(graphics::par(op))
  }
  titles <- c(beta = "(a) Survival",
              theta = if (x$model == "fg")
                "(b) Conditional quality-adjusted survival"
              else "(b) Quality-adjusted survival")
  for (cp in comps) {
    lst <- x[[cp]]
    graphics::plot(c(0, 1), c(0, 1), type = "n", main = titles[[cp]],
                   xlab = paste0("Estimated probabilities - ", x$ref),
                   ylab = if (length(lst) == 1L)
                     paste0("Estimated probabilities - ", names(lst))
                   else "Estimated probabilities - group", ...)
    graphics::abline(0, 1, col = "grey60", lwd = 2)
    for (k in seq_along(lst))
      graphics::lines(lst[[k]]$ref, lst[[k]]$group, type = "s", lty = k)
    if (length(lst) > 1L)
      graphics::legend("bottomright", names(lst), lty = seq_along(lst),
                       bty = "n")
  }
  invisible(x)
}
