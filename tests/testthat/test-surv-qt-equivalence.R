# tests/testthat/test-surv-qt-equivalence.R
#
# Tests the equivalence between the original pure R implementation
# (kept here as the oracle) and the optimized C++ kernel.
#
# The oracle faithfully reproduces the original script, including the
# edge behaviors that the kernel replicates on purpose:
#   - summary(km_fit, times = time) returns NA for times beyond the
#     last event -> the na.rm = TRUE of the original drops the
#     contribution -> the kernel treats the weight as 0;
#   - ties in qalys get the same probability.
# The comparison is > (survival P(Q > q), right-continuous, the same
# convention as predict and the KM); the original script used >=.

# Oracle: original pure R implementation, with the comparison >
surv_qt_ref <- function(QALY, object, pred_lin_qt = NULL) {
  time   <- object[, 1]
  status <- object[, 2]
  if (is.null(pred_lin_qt)) {
    qalys <- QALY
  } else {
    qalys <- log(QALY) - pred_lin_qt
  }
  km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  cens_surv <- summary(object = km_fit, times = time)$surv
  pred_surv <- sapply(1:length(qalys), function(i) {
    sum(x = ((status * (qalys > qalys[i])) / cens_surv), na.rm = TRUE) /
      length(qalys)
  })
  pred_surv
}

gera_dados_survqt <- function(n = 50, seed = 42) {
  set.seed(seed)
  z1     <- rbinom(n, 1, 0.5)
  t      <- rexp(n, rate = exp(0.3 * z1))
  cens   <- rexp(n, rate = 0.2)
  time   <- pmin(t, cens)
  status <- as.integer(t <= cens)
  qaly   <- time * runif(n, 0.4, 1)
  list(time = time, status = status, qaly = qaly, z = z1)
}

test_that("C++ kernel equals the R oracle (without linear predictor)", {
  d <- gera_dados_survqt()
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("C++ kernel equals the R oracle (with linear predictor)", {
  d <- gera_dados_survqt(n = 60, seed = 7)
  
  lin_pred <- d$z * 0.2
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status),
                     pred_lin_qt = lin_pred)
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status,
                            lin_pred = lin_pred)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("Ties in qaly get the same probability", {
  d <- gera_dados_survqt(n = 30, seed = 99)
  # Forces ties: rounds the QALYs to a coarse grid
  d$qaly <- round(d$qaly, 1)
  
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  
  expect_equal(cpp, ref, tolerance = 1e-10)
  
  # Direct check of the invariant: same input, same output
  if (any(duplicated(d$qaly))) {
    tab <- tapply(cpp, d$qaly, function(x) length(unique(round(x, 12))))
    expect_true(all(tab == 1))
  }
})

test_that("Times beyond the last KM event have zero weight (NA in the original)", {
  d <- gera_dados_survqt(n = 20, seed = 5)
  
  # Forces a very large time: beyond the largest time of the censoring KM.
  # The oracle's summary() would return NA and na.rm = TRUE would drop
  # the contribution; the kernel must replicate that behavior.
  d$time[1] <- max(d$time) * 10
  d$qaly[1] <- d$qaly[1] * 10
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("Wrapper validates the inputs", {
  d <- gera_dados_survqt(n = 10)
  
  expect_error(
    QAST:::surv_qt_fun(qaly = d$qaly, time = d$time[-1], status = d$status),
    "same length as qaly"
  )
  expect_error(
    QAST:::surv_qt_fun(qaly = -d$qaly, time = d$time, status = d$status,
                       lin_pred = rep(0, 10)),
    "strictly positive"
  )
})

test_that("An external km_fit is reused correctly", {
  d <- gera_dados_survqt(n = 40, seed = 3)
  
  km_fit <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
  
  cpp_interno <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time,
                                    status = d$status)
  cpp_externo <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time,
                                    status = d$status, km_fit = km_fit)
  
  expect_equal(cpp_interno, cpp_externo, tolerance = 1e-12)
})