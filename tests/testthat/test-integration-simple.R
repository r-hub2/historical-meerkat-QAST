# tests/testthat/test-integration-simple.R

test_that("betahat and thetahat reproduce reference values of the original (seed=42)", {
  skip_on_cran()

  # Same data generator as the original. The intercept (column of 1s)
  # enters ONLY in the generation of t; the model does not estimate an
  # intercept (rank equations, invariant to shifts), so it leaves z in the fit.
  set.seed(42)
  n <- 300
  z_gen <- cbind(1, matrix(rnorm(n * 2), n, 2))
  t <- exp(as.vector(z_gen %*% c(0.5, -0.3, 0.2)) + rnorm(n, sd = 0.5))
  qt <- t * runif(n, 0.5, 1.0)
  status <- rbinom(n, 1, 0.85)
  z <- z_gen[, -1]

  fit <- semi_aft_prop(t, qt, status, z,
                       h1 = "log", h2 = "log",
                       method = "Nelder-Mead",
                       norm_min = "L2",
                       cens_weight = TRUE,
                       M = 10, var_comp = FALSE,
                       var_method = "jackknife")

  # References of the original script on the SAME data (set.seed(42), n = 300),
  # slopes only. The original also optimized the flat direction of the
  # intercept (value from runif); without it, the Nelder-Mead path over
  # the step objective changes and the end point moves ~2e-4 within the
  # set of minima. Hence: (1) slopes within 1e-3 and (2) objective at the
  # estimated point not worse than at the reference.
  beta_ref  <- c(-0.320483557958109, 0.177030351648878)
  theta_ref <- c(-0.341859946407396, 0.182101852818395)

  expect_equal(unname(fit$beta),  beta_ref,  tolerance = 1e-3)
  expect_equal(unname(fit$theta), theta_ref, tolerance = 1e-3)

  lt <- log(t); lq <- log(qt)
  HG_surv <- QAST:::HG(X = z, t = t, status = status)
  o1 <- function(b) QAST:::U1_obj_cpp(b, lt, z, status, "L2")
  o2 <- function(th) QAST:::U2_obj_cpp(th, lq, z, status, HG_surv, "L2")
  expect_lte(o1(unname(fit$beta)),  o1(beta_ref)  + 1e-12)
  expect_lte(o2(unname(fit$theta)), o2(theta_ref) + 1e-12)
})

test_that("constant column (intercept) in z is rejected in both models", {
  set.seed(1)
  n <- 40
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)

  expect_error(
    semi_aft_prop(t, qt, st, cbind(1, rnorm(n)), "log", "log",
                  "Nelder-Mead", "L2", M = 0, var_comp = FALSE,
                  var_method = "jackknife"),
    "intercept")

  z_fg <- cbind(zb, const = 0)
  expect_error(
    semi_aft_fg(t, qt, st, z_fg, "log", "log", "Nelder-Mead", "L2",
                L = max(t), phi = "identity", M = 0, var_comp = FALSE,
                var_method = "jackknife"),
    "const")

  # Names of the offending columns appear in the message
  expect_error(QAST:::validate_no_constant_z(cbind(a = 1:3, b = 5)), "b")
  expect_silent(QAST:::validate_no_constant_z(zb))
})

test_that("collinear columns give a warning (not an error) in both models", {
  set.seed(3)
  n <- 50
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)

  # duplicated (FG) and linear combination (proposed)
  z_dup <- cbind(zb, trt2 = zb[, "trt"])
  x1 <- rnorm(n); x2 <- rnorm(n)
  z_lin <- cbind(x1 = x1, x2 = x2, x3 = x1 + 2 * x2)

  expect_warning(
    fit <- semi_aft_fg(t, qt, st, z_dup, "log", "log", "Nelder-Mead", "L2",
                       L = max(t), phi = "identity", M = 0, var_comp = FALSE,
                       var_method = "jackknife"),
    "full rank")
  expect_length(fit$beta, 3)

  expect_warning(
    semi_aft_prop(t, qt, st, z_lin, "log", "log", "Nelder-Mead", "L2",
                  M = 0, var_comp = FALSE, var_method = "jackknife"),
    "x3")

  expect_silent(QAST:::warn_collinear_z(zb))
})

test_that("collinear z: a single warning (without repeating the fast_path one)", {
  set.seed(3)
  n <- 50
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)
  w <- character()
  withCallingHandlers(
    semi_aft_fg(t, qt, st, cbind(zb, trt2 = zb[, "trt"]), "log", "log",
                "Nelder-Mead", "L2", L = max(t), phi = "identity", M = 0,
                var_comp = FALSE, var_method = "jackknife"),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "full rank")
})
