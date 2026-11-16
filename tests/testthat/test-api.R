# Interface dos modelos: validacao das entradas, padroes, formula com
# Qsurv() e a classe qast_fit (print, summary, coef, vcov).

mk_api_data <- function(n = 60, seed = 11) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
  t <- rexp(n, rate = exp(-0.4 * z[, "trt"]))
  status <- rbinom(n, 1, 0.8)
  qt <- t * runif(n, 0.5, 1)
  x <- cbind(idade = rnorm(n), trt = z[, "trt"])
  list(t = t, qt = qt, status = status, z = z, x = x, L = quantile(t, 0.9))
}

test_that("validacao: entradas invalidas param com mensagem clara", {
  d <- mk_api_data()
  fg <- function(t = d$t, qt = d$qt, status = d$status, z = d$z, ...)
    semi_aft_fg(t, qt, status, z, L = d$L, ...)
  prop <- function(t = d$t, qt = d$qt, status = d$status, z = d$x, ...)
    semi_aft_prop(t, qt, status, z, ...)

  t_na <- d$t; t_na[3] <- NA
  expect_error(fg(t = t_na), "t contains NA at observations 3")
  expect_error(prop(t = t_na), "t contains NA")

  t_neg <- d$t; t_neg[2] <- -1
  expect_error(fg(t = t_neg), "t must be finite and positive at observations 2")

  qt_big <- d$qt; qt_big[5] <- d$t[5] * 3
  expect_error(fg(qt = qt_big), "qt .* greater than t")
  expect_error(prop(qt = qt_big), "greater than t")

  qt_zero <- d$qt; qt_zero[4] <- 0
  expect_error(fg(qt = qt_zero), "h2\\(qt\\) is not finite")

  st_bad <- d$status; st_bad[1] <- 2
  expect_error(fg(status = st_bad), "status must be 0")
  expect_error(fg(status = rep(0, length(d$t))), "no events")

  expect_error(fg(t = d$t[-1]), "same number of observations")
  z_na <- d$z; z_na[7, 1] <- NA
  expect_error(fg(z = z_na), "z contains NA at observations 7")

  expect_error(semi_aft_fg(d$t, d$qt, d$status, d$z), "L .* is required")
  expect_error(semi_aft_fg(d$t, d$qt, d$status, d$z, L = -1), "L must be")

  expect_error(fg(var_comp = TRUE, var_method = "bootstrap"),
               "M .* integer >= 2")
  expect_error(fg(var_method = "outro"), "should be one of")
  expect_error(prop(cens_weight = FALSE, var_comp = TRUE),
               "requires the IPCW weights")
})

test_that("padroes: so os dados (e L no FG) sao obrigatorios", {
  d <- mk_api_data()
  full <- semi_aft_fg(d$t, d$qt, d$status, d$z, "log", "log",
                      "Nelder-Mead", "L2", L = d$L, phi = "identity")
  short <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L)
  expect_identical(coef(short), coef(full))

  p_full <- semi_aft_prop(d$t, d$qt, d$status, d$x, "log", "log",
                          "Nelder-Mead", "L2")
  p_short <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_identical(coef(p_short), coef(p_full))
})

test_that("formula com Qsurv() equivale a interface de matrizes", {
  d <- mk_api_data()
  df <- data.frame(tt = d$t, q = d$qt, ev = d$status, d$z, d$x[, "idade",
                                                               drop = FALSE])
  f_mat <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L)
  f_frm <- semi_aft_fg(Qsurv(tt, q, ev) ~ trt + sex, data = df, L = d$L)
  expect_identical(coef(f_frm), coef(f_mat))

  p_mat <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  p_frm <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + trt, data = df)
  expect_identical(coef(p_frm), coef(p_mat))

  # Fator de 3 niveis entra como 2 indicadoras (sem intercepto)
  df$g <- factor(sample(c("a", "b", "c"), nrow(df), TRUE))
  p_fac <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + g, data = df)
  expect_named(p_fac$beta, c("idade", "gb", "gc"))

  expect_error(semi_aft_prop(tt ~ idade, data = df), "Qsurv")
  df$tt[2] <- NA
  expect_error(semi_aft_prop(Qsurv(tt, q, ev) ~ idade, data = df),
               "t contains NA")
  expect_error(Qsurv(1:3, 1:2, c(1, 0, 1)), "same length")
})

test_that("qast_fit: print, summary, coef e vcov", {
  d <- mk_api_data()
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE)

  expect_s3_class(fit, "qast_fit")
  expect_identical(fit$model, "fg")
  expect_identical(fit$n, length(d$t))
  expect_identical(fit$n_events, sum(d$status))

  cf <- coef(fit)
  expect_identical(dim(cf), c(2L, 2L))
  expect_identical(colnames(cf), c("beta", "theta"))
  expect_identical(coef(fit, "beta"), fit$beta)
  expect_identical(coef(fit, "theta"), fit$theta)

  for (comp in c("beta", "theta")) {
    V <- vcov(fit, comp)
    expect_identical(dim(V), c(2L, 2L))
    expect_true(isSymmetric(V))
    expect_equal(diag(V), fit[[paste0("var_", comp)]], tolerance = 1e-12)
  }

  s <- summary(fit)
  expect_s3_class(s, "summary.qast_fit")
  expect_identical(colnames(s$beta),
                   c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
  expect_equal(unname(s$theta[, "se(coef)"]), unname(sqrt(fit$var_theta)))

  expect_output(print(fit), "Fine and Gelber")
  expect_output(print(s), "Quality-adjusted survival time component")

  # Sem variancia: vcov falha com mensagem; summary so com coef/exp(coef)
  fit0 <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_error(vcov(fit0), "var_comp = TRUE")
  expect_identical(colnames(summary(fit0)$beta), c("coef", "exp(coef)"))
  expect_output(print(fit0), "no variance")
})

test_that("vcov do bootstrap e a covariancia das replicas", {
  d <- mk_api_data(n = 50)
  set.seed(3)
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                     var_method = "bootstrap", M = 20)
  expect_identical(fit$var_method, "bootstrap")
  expect_identical(fit$M, 20)
  expect_equal(diag(vcov(fit, "theta")), fit$var_theta, tolerance = 1e-12)
})

test_that("uma unica covariavel: metodos sem suporte a 1D usam subplex", {
  d <- mk_api_data()
  z1 <- d$z[, "trt", drop = FALSE]
  expect_message(fit <- semi_aft_prop(d$t, d$qt, d$status, z1),
                 "using \"subplex\"")
  expect_identical(fit$settings$method, "subplex")
  expect_true(all(is.finite(coef(fit))))
  expect_message(semi_aft_fg(d$t, d$qt, d$status, z1, L = d$L), "subplex")
  # metodo que aceita 1D: sem troca e sem mensagem
  expect_silent(fit_b <- semi_aft_prop(d$t, d$qt, d$status, z1,
                                       method = "BFGS"))
  expect_identical(fit_b$settings$method, "BFGS")
})

test_that("interceptos: media KM (beta) e media IPCW (theta) dos residuos", {
  d <- mk_api_data(n = 80, seed = 21)
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_named(fit$intercept, c("beta", "theta"))

  # Referencia independente: media pelos saltos do KM dos residuos e
  # media IPCW com o KM da censura avaliado em t
  e <- log(d$t) - drop(d$x %*% fit$beta)
  km <- survival::survfit(survival::Surv(e, d$status) ~ 1)
  s <- km$surv
  a_t <- sum(-diff(c(1, s)) * km$time) + s[length(s)] * km$time[length(s)]
  r <- log(d$qt) - drop(d$x %*% fit$theta)
  G <- summary(survival::survfit(survival::Surv(d$t, 1 - d$status) ~ 1),
               times = d$t)$surv
  # censurado com G = 0 da 0/0; como no script da tese, entra com peso 0
  a_q <- sum(d$status * r / G, na.rm = TRUE) / length(r)
  expect_equal(unname(fit$intercept), c(a_t, a_q), tolerance = 1e-12)

  expect_output(print(fit), "Intercepts")
  expect_output(print(summary(fit)), "IPCW mean of the residuals")
})

test_that("interceptos sem censura: medias simples, inclusive com residuos negativos", {
  d <- mk_api_data(n = 60, seed = 5)
  st <- rep(1, length(d$t))
  t <- d$t / 100; qt <- d$qt / 100       # log-tempos negativos
  fit <- semi_aft_fg(t, qt, st, d$z, L = max(t))
  e <- log(t) - drop(d$z %*% fit$beta)
  r <- log(qt) - drop(d$z %*% fit$theta)
  expect_true(all(e < 0))
  expect_equal(fit$intercept[["beta"]], mean(e), tolerance = 1e-12)
  expect_equal(fit$intercept[["theta"]], mean(r), tolerance = 1e-12)
})
