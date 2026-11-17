# tests/testthat/test-fg-equivalence-extended.R
#
# Numerical equivalence of the FG kernel (U2_parts_fg_cpp) with the R
# reference. Extends test-fg-equivalence.R (grid of theta x phi x norm) with:
#
#   1. grid of beta   - the precomputed pairs depend on beta
#   2. empty group / constant column / group without events (path rep(1, n))
#   3. heavy censoring with the KM survival reaching zero (denom filter)
#   4. exact tie at the boundary of the h1L indicator (strict inequality)
#   5. smaller L (quantile 0.5) to exercise the truncation on more pairs
#   6. smoke: n = 200, p = 1
#
# All blocks compare U2_parts_fg_cpp against U2_parts_fg_ref (pure R
# reference) and the aggregation of the norm, tolerance 1e-10.

# --------------------------------------------------------------------------
# Helper: full comparison (parts + aggregated obj) for one point
# --------------------------------------------------------------------------
expect_fg_equivalent <- function(theta, lt, lq, z, status, beta, hL,
                                 km_surv, phi, norm_min, tol = 1e-10) {
  parts_cpp <- QAST:::U2_parts_fg_cpp(theta, lt, lq, z, status,
                                      beta, hL, km_surv, phi)
  parts_ref <- QAST:::U2_parts_fg_ref(theta, lt, lq, z, status,
                                      beta, hL, km_surv, phi)
  expect_equal(parts_cpp, parts_ref, tolerance = tol)
  
  obj_cpp <- QAST:::U2_obj_fg_cpp(theta, lt, lq, z, status,
                                  beta, hL, km_surv, phi, norm_min)
  if (norm_min == "L1") {
    obj_ref <- sum(abs(parts_ref)) / sqrt(length(parts_ref))
  } else {
    obj_ref <- sqrt(sum(parts_ref^2)) / sqrt(length(parts_ref))
  }
  expect_equal(obj_cpp, obj_ref, tolerance = tol)
}

# Goes through the grids of phi/norm/theta for a fixed configuration
compare_over_grids <- function(lt, lq, z, status, beta, hL, km_surv,
                               theta_grid) {
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      for (i in seq_len(nrow(theta_grid))) {
        expect_fg_equivalent(theta_grid[i, ], lt, lq, z, status,
                             beta, hL, km_surv, phi, norm_min)
      }
    }
  }
}

# --------------------------------------------------------------------------
# 1. Grid of beta x theta (base case: "clean" random data, L at q90)
# --------------------------------------------------------------------------
test_that("FG: C++ kernel equals the reference on a grid of beta and theta", {
  set.seed(42)
  
  n <- 40
  z <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.4))
  colnames(z) <- c("trt", "sex")
  
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  
  lt <- log(t); lq <- log(qt)
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  # beta = 0 checks the degenerate case z'beta = 0 in the h1L indicator
  beta_grid <- rbind(
    c(0.0,  0.0),
    c(0.3, -0.2),
    c(0.8,  0.5),   # strong effects: push more pairs outside L
    c(-0.5, 0.7)    # mixed signs
  )
  
  theta_grid <- rbind(
    c(0.0, 0.0),
    c(0.2, -0.1),
    c(-0.4, 0.5),
    c(0.7, 0.7)
  )
  
  for (b in seq_len(nrow(beta_grid))) {
    beta <- beta_grid[b, ]
    
    # L at the 0.9 quantile (default of the original test)
    hL <- log(stats::quantile(t, 0.9))
    compare_over_grids(lt, lq, z, status, beta, hL, km_surv, theta_grid)
    
    # Smaller L: the truncation cuts many more pairs
    hL_mid <- log(stats::quantile(t, 0.5))
    compare_over_grids(lt, lq, z, status, beta, hL_mid, km_surv, theta_grid)
  }
})

# --------------------------------------------------------------------------
# 2. Constant column / empty group / group without events
#    Covers: parts(j) = 0 (nk = 0 or nl = 0) and the rep(1, n) path
#    of km_weights_at (group with status == 0 for all).
# --------------------------------------------------------------------------
test_that("FG: constant column, empty group and group without events", {
  set.seed(7)
  
  n <- 60
  z <- cbind(rbinom(n, 1, 0.5),   # normal
             rep(0, n),           # constant 0  -> group 1 empty
             rep(1, n))           # constant 1  -> group 0 empty
  colnames(z) <- c("trt", "cte0", "cte1")
  
  t <- rexp(n)
  # Group without events: all subjects with trt = 1 are censored.
  # (exercises the rep(1, n) path of km_weights_at inside prep_km_by_column)
  status <- ifelse(z[, 1] == 1, 0L, rbinom(n, 1, 0.8))
  
  qt <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  hL <- log(stats::quantile(t, 0.9))
  
  beta  <- c(0.3, 0.0, 0.0)
  theta <- c(0.2, 0.0, 0.0)
  
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      expect_fg_equivalent(theta, lt, lq, z, status, beta, hL,
                           km_surv, phi, norm_min)
    }
  }
  
  # The constant columns must have an identically null component
  parts <- QAST:::U2_parts_fg_cpp(theta, lt, lq, z, status,
                                  beta, hL, km_surv, "identity")
  expect_identical(unname(parts[2]), 0.0)
  expect_identical(unname(parts[3]), 0.0)
})

# --------------------------------------------------------------------------
# 3. Heavy censoring: the KM survival gets close to / reaches zero
#    (exercises the filter !(denom > 0) || !isfinite(denom))
# --------------------------------------------------------------------------
test_that("FG: heavy censoring with near-zero KM weights", {
  set.seed(123)
  
  n <- 120
  z <- cbind(rbinom(n, 1, 0.6), rbinom(n, 1, 0.3))
  colnames(z) <- c("trt", "sex")
  
  t <- rexp(n, rate = exp(0.3 * z[, 1]))
  
  # 70% censoring + 1 subject censored after the last event: ensures a
  # tail with a very low S(t)
  status <- rbinom(n, 1, 0.30)
  status[n] <- 0L
  t[n] <- max(t) * 10
  
  qt <- pmin(t, t * runif(n, 0.3, 1))
  lt <- log(t); lq <- log(qt)
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  # Sanity: the KM weights of the last subject must be valid (not NaN/Inf)
  w0 <- km_surv[[1]]$w0
  expect_true(all(is.finite(w0)))
  
  hL <- log(stats::quantile(t, 0.9))
  compare_over_grids(lt, lq, z, status, c(0.3, -0.2), hL, km_surv,
                     rbind(c(0.0, 0.0), c(0.2, -0.1), c(0.5, 0.3)))
})

# --------------------------------------------------------------------------
# 4. EXACT tie at the boundary of the h1L indicator
#
#   The indicator uses a STRICT inequality:
#     h1L = 1 if  hL - max(z'beta) > max(e1k, e1l)   else 0
#
#   Exact ties are built with beta = 0 (z'beta identically 0) and
#   DUPLICATED times equal to L: for every pair whose two members have
#   t == L, max(e1k, e1l) == hL bit for bit (same call to log()), and
#   the strict comparison must give 0. Pairs with max(e1) just below
#   must give 1. It is the case most sensitive to floating-point
#   reordering: the exact ties built here must stay identical to the
#   reference.
# --------------------------------------------------------------------------
test_that("FG: exact tie at the h1L boundary (strict inequality)", {
  set.seed(99)
  
  n <- 50
  z <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.4))
  colnames(z) <- c("trt", "sex")
  
  t <- rexp(n)
  status <- rbinom(n, 1, 0.8)
  
  # Duplicates times at the value of L: 4 subjects with t exactly equal to L
  L  <- stats::quantile(t, 0.7, names = FALSE)
  idx_dup <- c(1, 2, 3, 4)
  t[idx_dup] <- L
  
  # beta = 0 => z'beta == 0 exactly; the indicator becomes log(L) > max(e1)
  beta <- c(0.0, 0.0)
  
  qt <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  hL <- log(L)  # same log() on the same double => bit for bit equal to the lt of the duplicates
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  theta_grid <- rbind(c(0.0, 0.0), c(0.2, -0.1))
  
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      for (i in seq_len(nrow(theta_grid))) {
        expect_fg_equivalent(theta_grid[i, ], lt, lq, z, status,
                             beta, hL, km_surv, phi, norm_min)
      }
    }
  }
})

# --------------------------------------------------------------------------
# 6. Smoke: p = 1, n = 200
# --------------------------------------------------------------------------
test_that("FG: smoke p = 1, n = 200", {
  set.seed(2024)
  
  n <- 200
  z <- matrix(rbinom(n, 1, 0.45), n, 1)
  colnames(z) <- "x1"
  
  t      <- rexp(n, rate = exp(0.2 * z[, 1]))
  status <- rbinom(n, 1, 0.75)
  qt     <- pmin(t, t * runif(n, 0.5, 1))
  
  lt <- log(t); lq <- log(qt)
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  hL <- log(stats::quantile(t, 0.9))
  
  compare_over_grids(lt, lq, z, status, c(0.2), hL, km_surv,
                     rbind(c(0.0), c(0.3), c(-0.2)))
})
