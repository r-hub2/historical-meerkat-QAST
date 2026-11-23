# tests/testthat/test-km-findinterval-ab.R
#
# The oracle is the evaluation by which() (km_weights_at_legacy, in
# helper-references.R); the
# findInterval path must reproduce it in every edge case. The estimator
# is survfit in both; what is checked is the evaluation of the step
# function.

km_compare <- function(t, status, tol = 1e-12) {
  expect_equal(QAST:::km_weights_at(t, status),
               km_weights_at_legacy(t, status),
               tolerance = tol)
  invisible(TRUE)
}

test_that("KM findInterval: base case with moderate censoring", {
  set.seed(42)
  n <- 150
  km_compare(rexp(n), rbinom(n, 1, 0.8))
})

test_that("KM findInterval: many ties (duplicated times)", {
  set.seed(1)
  n <- 200
  t <- sample(1:10, n, replace = TRUE) * 1.0
  km_compare(t, rbinom(n, 1, 0.7))
})

test_that("KM findInterval: exact ties between failure and censoring", {
  # The tie convention is survfit's; both paths query the SAME fit,
  # only the interpolation is checked
  t      <- c(1, 1, 1, 2, 2, 3, 3, 3, 3)
  status <- c(0, 1, 1, 0, 1, 0, 0, 1, 1)
  km_compare(t, status)
})

test_that("KM findInterval: heavy censoring with a long tail", {
  set.seed(123)
  n <- 120
  t <- rexp(n)
  status <- rbinom(n, 1, 0.30)
  t[n] <- max(t) * 10
  status[n] <- 0L
  km_compare(t, status)
  
  # Subject beyond the tail: BOTH paths return the LAST value of the
  # curve (never NaN)
  w_new <- QAST:::km_weights_at(t, status)
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  expect_identical(w_new[n], fit$surv[length(fit$surv)])
})

test_that("KM findInterval: group without events (rep(1, n))", {
  set.seed(7)
  t <- rexp(50)
  expect_identical(QAST:::km_weights_at(t, rep(0L, 50)), rep(1, 50))
})

test_that("KM findInterval: group without censoring (S == 1)", {
  set.seed(9)
  t <- rexp(80)
  status <- rep(1L, 80)
  km_compare(t, status)
  expect_true(all(QAST:::km_weights_at(t, status) == 1))
})

test_that("KM findInterval: empty group and before the first jump", {
  expect_identical(QAST:::km_weights_at(numeric(0), integer(0)),
                   numeric(0))
  
  t <- c(1, 10, 20, 30)
  status <- c(1, 0, 1, 1)
  expect_identical(QAST:::km_weights_at(t, status)[1], 1)
})

test_that("KM findInterval: names/attributes of t do not propagate", {
  t <- c(a = 1.5, b = 2.5, c = 0.5)
  status <- c(1, 0, 1)
  w <- QAST:::km_weights_at(t, status)
  expect_equal(length(w), 3L)
  expect_null(names(w))
})

test_that("prep_km_by_column: one weight per subject, from the given partition", {
  set.seed(99)
  n <- 100
  z <- cbind(a = rbinom(n, 1, 0.5), b = rbinom(n, 1, 0.3))
  t <- rexp(n)
  status <- rbinom(n, 1, 0.75)
  g <- paste(z[, 1], z[, 2])
  km <- QAST:::prep_km_by_column(z, t, status, group = g)
  km_leg <- prep_km_by_column_legacy(z, t, status, group = g)

  H <- numeric(n)
  for (l in unique(g)) {
    i <- g == l
    H[i] <- QAST:::km_weights_at(t[i], status[i])
  }
  for (j in 1:2) {
    expect_equal(km[[j]]$w0, H[z[, j] == 0], tolerance = 1e-12)
    expect_equal(km[[j]]$w1, H[z[, j] == 1], tolerance = 1e-12)
    expect_equal(km[[j]]$w0, km_leg[[j]]$w0, tolerance = 1e-12)
  }
  expect_length(km, 2L)
})

test_that("prep_km_by_column: partition given by group", {
  set.seed(5)
  n <- 120
  z <- cbind(rbinom(n, 1, .5), rbinom(n, 1, .3), rbinom(n, 1, .6))
  t <- rexp(n); status <- rbinom(n, 1, .75)
  km1 <- QAST:::prep_km_by_column(z, t, status)   # default: no stratification
  H1 <- QAST:::km_weights_at(t, status)
  for (j in 1:3) expect_equal(km1[[j]]$w1, H1[z[, j] == 1], tolerance = 0)
  kmz <- QAST:::prep_km_by_column(z, t, status, group = z[, 2])
  Hz <- QAST:::km_weights_by_group(t, status, factor(z[, 2]))
  expect_equal(kmz[[3]]$w0, Hz[z[, 3] == 0], tolerance = 0)
  expect_identical(nlevels(QAST:::fg_cens_partition(z)), 1L)
})

test_that("prep_km_by_column: golden chain - weights feed the FG kernel", {
  # Smoke: the whole FG pipeline (prep_km -> pairs -> affine) gives an
  # identical objective with KM weights by findInterval and by which()
  set.seed(2024)
  n <- 120
  z <- cbind(rbinom(n, 1, 0.45), rbinom(n, 1, 0.35))
  colnames(z) <- c("trt", "sex")
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  hL <- log(quantile(t, 0.9))
  beta <- c(0.3, -0.2)
  
  pairs_fast <- QAST:::make_fg_pairs(lt, lq, z, status, beta, hL,
                                     QAST:::prep_km_by_column(z, t, status),
                                     "identity")
  pairs_leg <- QAST:::make_fg_pairs(lt, lq, z, status, beta, hL,
                                    prep_km_by_column_legacy(z, t, status),
                                    "identity")
  aff_fast <- QAST:::make_fg_affine(pairs_fast, lq, z)
  aff_leg  <- QAST:::make_fg_affine(pairs_leg, lq, z)
  
  for (th in list(c(0, 0), c(0.2, -0.1))) {
    expect_equal(QAST:::eval_fg_obj_affine(aff_fast, unlist(th), "L2"),
                 QAST:::eval_fg_obj_affine(aff_leg,  unlist(th), "L2"),
                 tolerance = 1e-12)
  }
})


