# R/utils_jackknife.R
#
# Delete-1 jackknife, with the same standard-error formula as
# jackknifeR::jackknife(d = 1, weight = FALSE) (checked in the tests):
#
#   jk        = n x k matrix of the estimates without observation i
#   valid     = complete.cases(jk)          (rows with NA dropped)
#   theta_dot = colMeans(jk[valid, ])
#   se        = sqrt( (n - 1) / n_valid * colSums((jk - theta_dot)^2) )
#
# Execution: fork (parallel::mclapply); several parameter blocks (e.g.
# beta and theta) in a single pass, with NA removal per block; warnings
# of the workers passed to the main session, once per distinct message.
# The bootstrap (.qast_bootstrap) uses the same runner.
#
# Workers: getOption("QAST.cores"), otherwise the PHYSICAL cores
# (detectCores(logical = FALSE)), at most 2 under R CMD check
# (_R_CHECK_LIMIT_CORES_). On Windows (no fork): sequential.

#' Number of workers for the variance refits
#' @keywords internal
.qast_workers <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  opt <- getOption("QAST.cores")
  nc <- if (!is.null(opt)) {
    as.integer(opt)
  } else {
    suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  if (is.na(nc) || nc < 1L) nc <- 1L
  # R CMD check --as-cran (and CRAN) limit to 2 simultaneous processes;
  # parallel aborts above that when _R_CHECK_LIMIT_CORES_ is set
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(toupper(lim), "FALSE")) nc <- min(nc, 2L)
  as.integer(nc)
}

# Runs refit(task) for every task (fork when cores > 1) and returns the
# matrix of estimates, one row per task. Warnings of the replicates are
# reissued once per distinct message; an error in any replicate is
# re-raised in the main session (a degenerate replicate must return NA
# instead of failing).
.qast_run_reps <- function(tasks, refit, cores = .qast_workers()) {

  one <- function(k) {
    warns <- character()
    val <- withCallingHandlers(
      refit(k),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = as.numeric(val), warnings = warns)
  }

  res <- if (cores > 1L) {
    # mclapply's own warnings ("scheduled core ... error") are redundant:
    # the original error is re-raised below. The warnings of the statistic
    # were already captured inside one().
    suppressWarnings(
      parallel::mclapply(tasks, one, mc.cores = cores, mc.preschedule = TRUE)
    )
  } else {
    lapply(tasks, one)
  }

  # Error in any replicate: re-raised in the main session
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"),
                logical(1))
  if (any(bad)) {
    r <- res[[which(bad)[1]]]
    if (is.null(r)) stop("variance replicates: a worker ended without a result")
    cond <- attr(r, "condition")
    if (inherits(cond, "condition")) stop(cond)
    stop(as.character(r))
  }

  for (msg in unique(unlist(lapply(res, `[[`, "warnings"))))
    warning(msg, call. = FALSE)

  do.call(rbind, lapply(res, `[[`, "value"))
}

# Warning when replicates were discarded (NA: rank-deficient covariates in
# the resample or an optimizer that did not converge)
.qast_warn_discarded <- function(blk, total, what) {
  nd <- max(vapply(blk, function(b) b$n_drop, numeric(1)))
  if (nd > 0)
    warning(sprintf(paste0("%d of %d %s replicates discarded (covariates ",
                           "without full rank in the resample or no ",
                           "convergence of the optimizer)"),
                    nd, total, what), call. = FALSE)
  invisible(nd)
}

#' Delete-1 jackknife standard errors, by parameter block
#'
#' @param n Number of observations
#' @param refit Function of i (the removed index) returning the vector of
#'   estimates on the sample without observation i (NA for a degenerate
#'   replicate, which is discarded)
#' @param blocks Named list of column indices of the vector returned by
#'   refit (e.g. \code{list(beta = 1:p, theta = p + 1:p)})
#' @param cores Number of workers (fork)
#' @param cov If TRUE, returns \code{list(se, cov, reps, drop, n_drop)} per
#'   block, with the jackknife covariance matrix
#'   \code{(n - 1) / n_valid * crossprod(jk - theta_dot)} (whose diagonal is
#'   \code{se^2}), the valid replicates (rows of \code{reps}), the
#'   observation removed in each of them (\code{drop}) and the number of
#'   discarded replicates (\code{n_drop})
#' @return Named list (as blocks) of standard-error vectors, or of
#'   \code{list(se, cov, reps, drop, n_drop)} when \code{cov = TRUE}
#' @keywords internal
.jackknife_se <- function(n, refit, blocks, cores = .qast_workers(),
                          cov = FALSE) {

  jk <- .qast_run_reps(seq_len(n), refit, cores)

  lapply(blocks, function(cols) {
    b <- jk[, cols, drop = FALSE]
    ok <- stats::complete.cases(b)
    b <- b[ok, , drop = FALSE]
    nv <- nrow(b)
    theta_dot <- colMeans(b)
    dev <- sweep(b, 2, theta_dot, "-")
    se <- sqrt(((n - 1) / nv) * colSums(dev^2))
    if (!cov) return(se)
    list(se = se, cov = ((n - 1) / nv) * crossprod(dev),
         reps = unname(b), drop = which(ok), n_drop = sum(!ok))
  })
}

# Nonparametric bootstrap (resampling of whole subjects), by parameter
# block. The resamples are drawn as in boot::boot(sim = "ordinary"):
# matrix(sample.int(n, n * M, replace = TRUE), M), so a given seed gives
# the same resamples. Each replicate refits ALL the components on the same
# resample (joint distribution of the estimates); replicates with NA are
# discarded per block. Per block: list(se, cov, reps, index, n_drop), with
# cov = cov() of the valid replicates and index = their resampled rows.
.qast_bootstrap <- function(n, M, refit, blocks, cores = .qast_workers()) {
  index <- matrix(sample.int(n, n * M, replace = TRUE), M)
  vals <- .qast_run_reps(seq_len(M), function(b) refit(index[b, ]), cores)
  lapply(blocks, function(cols) {
    b <- vals[, cols, drop = FALSE]
    ok <- stats::complete.cases(b)
    b <- b[ok, , drop = FALSE]
    V <- if (nrow(b) > 1L) stats::cov(b)
         else matrix(NA_real_, length(cols), length(cols))
    list(se = sqrt(diag(V)), cov = V, reps = unname(b),
         index = index[ok, , drop = FALSE], n_drop = sum(!ok))
  })
}

# Rank deficiency of a resample (bootstrap) or subsample (jackknife):
# covariates with lower rank than in the full data (e.g. a binary
# covariate that became constant) do not identify the coefficients.
.qast_rank_deficient <- function(z_sub, rank_full) {
  qr(z_sub)$rank < rank_full
}

# Replicates stored for predict(): coefficients and rows of each replicate
.qast_boot_reps <- function(blk) list(reps = blk$reps, index = blk$index)

# Jackknife replicates (output of .jackknife_se with cov = TRUE)
.qast_jack_reps <- function(blk) list(reps = blk$reps, drop = blk$drop)

# Difference theta - beta from the joint block (columns beta, then theta)
# of the same replicates: replicates of the difference (with the rows of
# each replicate, as the other components) and its covariance
# A V A', A = [-I, I] (the jackknife and bootstrap covariances are
# quadratic forms of the replicates, so this equals the covariance of the
# differences)
.qast_diff_reps <- function(blk, p, var_method) {
  r <- blk$reps
  d <- r[, p + seq_len(p), drop = FALSE] - r[, seq_len(p), drop = FALSE]
  A <- cbind(-diag(p), diag(p))
  reps <- if (var_method == "jackknife") list(reps = d, drop = blk$drop)
          else list(reps = d, index = blk$index)
  list(reps = reps, cov = A %*% blk$cov %*% t(A))
}
