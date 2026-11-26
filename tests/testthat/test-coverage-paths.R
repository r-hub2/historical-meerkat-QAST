# Paths that the other tests do not exercise: optimizers, bootstrap of the
# proposed model, cens_weight = FALSE and the choice of the number of workers.

mk_cov_data <- function(n = 60, seed = 21) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
  t <- rexp(n, rate = exp(-0.4 * z[, "trt"]))
  status <- rbinom(n, 1, 0.8)
  qt <- t * runif(n, 0.5, 1)
  x <- cbind(age = rnorm(n), trt = z[, "trt"])
  list(t = t, qt = qt, status = status, z = z, x = x, L = quantile(t, 0.9))
}

test_that("all accepted optimizers give finite estimates", {
  d <- mk_cov_data()
  methods <- c("uobyqa", "mads", "Nelder-Mead", "subplex", "BB", "SANN",
               "bobyqa", "newuoa", "hjkb", "BFGS")
  for (m in methods) {
    set.seed(1)
    fit_p <- suppressWarnings(semi_aft_prop(d$t, d$qt, d$status, d$x,
                                            method = m))
    expect_true(all(is.finite(coef(fit_p))), info = m)
    expect_identical(fit_p$settings$method, m)
    # FG with fast_path = FALSE: theta also goes through the optimizer
    set.seed(1)
    fit_f <- suppressWarnings(semi_aft_fg(d$t, d$qt, d$status, d$z,
                                          L = d$L, method = m,
                                          fast_path = FALSE))
    expect_true(all(is.finite(coef(fit_f))), info = m)
  }
})

test_that("semi_aft_prop: bootstrap and cens_weight = FALSE", {
  # larger n: in small resamples the coxph of the IPCW weights may warn
  # about separation in a binary covariate (behavior of the data)
  d <- mk_cov_data(n = 150)
  set.seed(5)
  fit_b <- semi_aft_prop(d$t, d$qt, d$status, d$x, var_comp = TRUE,
                         var_method = "bootstrap", M = 15)
  expect_identical(fit_b$var_method, "bootstrap")
  expect_identical(fit_b$M, 15)
  for (comp in c("beta", "theta")) {
    V <- vcov(fit_b, comp)
    expect_true(isSymmetric(V))
    expect_true(all(diag(V) > 0))
  }
  expect_output(print(summary(fit_b)), "Variance: bootstrap \\(M = 15\\)")

  fit_nw <- semi_aft_prop(d$t, d$qt, d$status, d$x, cens_weight = FALSE)
  expect_false(fit_nw$settings$cens_weight)
  expect_true(all(is.finite(coef(fit_nw))))
  # without weights, theta uses the unweighted equation: it differs from the IPCW fit
  fit_w <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_false(isTRUE(all.equal(fit_nw$theta, fit_w$theta)))
})

test_that("semi_aft_fg: piecewise with L1 and jackknife", {
  d <- mk_cov_data()
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, phi = "piecewise",
                     norm_min = "L1", var_comp = TRUE)
  expect_identical(fit$settings$phi, "piecewise")
  expect_true(all(is.finite(fit$var_theta)))
  expect_output(print(fit), "phi = piecewise")
})

test_that(".qast_workers honors the QAST.cores option and the check limit", {
  skip_on_os("windows")
  old_opt <- options(QAST.cores = 3L)
  old_env <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit({
    options(old_opt)
    if (is.na(old_env)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_env)
  })
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  expect_identical(QAST:::.qast_workers(), 3L)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_identical(QAST:::.qast_workers(), 2L)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  expect_identical(QAST:::.qast_workers(), 3L)
  options(QAST.cores = 0L)
  expect_identical(QAST:::.qast_workers(), 1L)
})

test_that("C++ kernels do not touch the R random number generator", {
  # With an RNGScope in the wrappers, every call reloaded .Random.seed and
  # optim(method = "SANN") restarted its random walk at each evaluation.
  d <- mk_cov_data(n = 80)
  lt <- log(d$t); st <- as.numeric(d$status); z <- d$x * 1
  visited <- function(fn) {
    pts <- list()
    g <- function(b) { pts[[length(pts) + 1L]] <<- b; fn(b) }
    set.seed(7)
    stats::optim(c(0, 0), g, method = "SANN", control = list(maxit = 200))
    do.call(rbind, pts)
  }
  target <- function(b) sum((b - c(0.5, -0.2))^2)
  with_cpp <- function(b) {
    QAST:::U1_obj_cpp(b, lt, z, st, "L2")
    target(b)
  }
  expect_identical(visited(with_cpp), visited(target))
  # and a call leaves the stream where it was
  set.seed(3); a <- runif(1)
  set.seed(3); QAST:::U1_obj_cpp(c(0, 0), lt, z, st, "L2"); b <- runif(1)
  expect_identical(a, b)
})

test_that("bobyqa uses NLopt, not minqa", {
  d <- mk_cov_data(n = 80)
  testthat::local_mocked_bindings(
    bobyqa = function(...) stop("minqa::bobyqa called"), .package = "minqa"
  )
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x, method = "bobyqa")
  expect_true(all(is.finite(c(fit$beta, fit$theta))))
  ref <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_equal(fit$beta, ref$beta, tolerance = 0.1)
  fg <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, method = "bobyqa")
  expect_true(all(is.finite(c(fg$beta, fg$theta))))
})
