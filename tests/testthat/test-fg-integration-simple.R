gera_dados_fg <- function(n = 60, seed = 123) {
  set.seed(seed)
  z1 <- rbinom(n, 1, 0.5)
  z2 <- rbinom(n, 1, 0.4)
  z  <- cbind(z1, z2)
  colnames(z) <- c("trt", "sex")
  
  t      <- rexp(n, rate = exp(0.3 * z1))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  
  list(t = t, qt = qt, status = status, z = z)
}

test_that("semi_aft_fg fits without variance and returns the right structure", {
  d <- gera_dados_fg()
  
  fit <- semi_aft_fg(
    t = d$t, qt = d$qt, status = d$status, z = d$z,
    h1 = "log", h2 = "log",
    method = "Nelder-Mead", norm_min = "L2",
    L = stats::quantile(d$t, 0.9), phi = "identity",
    M = 0, var_comp = FALSE, var_method = "jackknife"
  )
  
  expect_type(fit, "list")
  expect_true(all(c("beta", "theta") %in% names(fit)))
  expect_false("var_beta" %in% names(fit))
  expect_length(fit$beta, 2)
  expect_length(fit$theta, 2)
  expect_equal(names(fit$beta), c("trt", "sex"))
  expect_equal(names(fit$theta), c("trt", "sex"))
  expect_true(all(is.finite(fit$beta)))
  expect_true(all(is.finite(fit$theta)))
})

test_that("semi_aft_fg rejects non-binary covariates", {
  d <- gera_dados_fg()
  z_cont <- cbind(d$z, age = runif(nrow(d$z), 40, 80))
  
  expect_error(
    semi_aft_fg(
      t = d$t, qt = d$qt, status = d$status, z = z_cont,
      h1 = "log", h2 = "log",
      method = "Nelder-Mead", norm_min = "L2",
      L = stats::quantile(d$t, 0.9), phi = "identity",
      M = 0, var_comp = FALSE, var_method = "jackknife"
    ),
    "binary \\(0/1\\) covariates"
  )
})

test_that("semi_aft_fg with bootstrap returns variances", {
  d <- gera_dados_fg(n = 30)
  
  fit <- semi_aft_fg(
    t = d$t, qt = d$qt, status = d$status, z = d$z,
    h1 = "log", h2 = "log",
    method = "Nelder-Mead", norm_min = "L2",
    L = stats::quantile(d$t, 0.9), phi = "identity",
    M = 5, var_comp = TRUE, var_method = "bootstrap"
  )
  
  expect_true(all(c("var_beta", "var_theta") %in% names(fit)))
  expect_length(fit$var_beta, 2)
  expect_length(fit$var_theta, 2)
  expect_true(all(is.finite(fit$var_beta)))
  expect_true(all(fit$var_beta >= 0))
  expect_true(all(fit$var_theta >= 0))
})