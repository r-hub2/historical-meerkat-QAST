# Caminhos que os demais testes nao exercitam: otimizadores, bootstrap do
# modelo proposto, cens_weight = FALSE e a escolha do numero de workers.

mk_cov_data <- function(n = 60, seed = 21) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
  t <- rexp(n, rate = exp(-0.4 * z[, "trt"]))
  status <- rbinom(n, 1, 0.8)
  qt <- t * runif(n, 0.5, 1)
  x <- cbind(age = rnorm(n), trt = z[, "trt"])
  list(t = t, qt = qt, status = status, z = z, x = x, L = quantile(t, 0.9))
}

test_that("todos os otimizadores aceitos produzem estimativas finitas", {
  d <- mk_cov_data()
  methods <- c("uobyqa", "mads", "Nelder-Mead", "subplex", "BB", "SANN",
               "bobyqa", "newuoa", "hjkb", "BFGS")
  for (m in methods) {
    set.seed(1)
    fit_p <- suppressWarnings(semi_aft_prop(d$t, d$qt, d$status, d$x,
                                            method = m))
    expect_true(all(is.finite(coef(fit_p))), info = m)
    expect_identical(fit_p$settings$method, m)
    # FG com fast_path = FALSE: theta tambem passa pelo otimizador
    set.seed(1)
    fit_f <- suppressWarnings(semi_aft_fg(d$t, d$qt, d$status, d$z,
                                          L = d$L, method = m,
                                          fast_path = FALSE))
    expect_true(all(is.finite(coef(fit_f))), info = m)
  }
})

test_that("semi_aft_prop: bootstrap e cens_weight = FALSE", {
  # n maior: em reamostras pequenas o coxph dos pesos IPCW pode avisar
  # sobre separacao numa covariavel binaria (comportamento dos dados)
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
  # sem pesos, theta usa a equacao nao ponderada: difere do ajuste IPCW
  fit_w <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_false(isTRUE(all.equal(fit_nw$theta, fit_w$theta)))
})

test_that("semi_aft_fg: piecewise com L1 e jackknife", {
  d <- mk_cov_data()
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, phi = "piecewise",
                     norm_min = "L1", var_comp = TRUE)
  expect_identical(fit$settings$phi, "piecewise")
  expect_true(all(is.finite(fit$var_theta)))
  expect_output(print(fit), "phi = piecewise")
})

test_that(".qast_workers respeita a opcao QAST.cores e o limite do check", {
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
