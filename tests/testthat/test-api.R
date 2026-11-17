# Model interface: input validation, defaults, formula with Qsurv() and
# the qast_fit class (print, summary, coef, vcov).

mk_api_data <- function(n = 60, seed = 11) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
  t <- rexp(n, rate = exp(-0.4 * z[, "trt"]))
  status <- rbinom(n, 1, 0.8)
  qt <- t * runif(n, 0.5, 1)
  x <- cbind(idade = rnorm(n), trt = z[, "trt"])
  list(t = t, qt = qt, status = status, z = z, x = x, L = quantile(t, 0.9))
}

test_that("validation: invalid inputs stop with a clear message", {
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

test_that("defaults: only the data (and L in FG) are required", {
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

test_that("formula with Qsurv() equals the matrix interface", {
  d <- mk_api_data()
  df <- data.frame(tt = d$t, q = d$qt, ev = d$status, d$z, d$x[, "idade",
                                                               drop = FALSE])
  f_mat <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L)
  f_frm <- semi_aft_fg(Qsurv(tt, q, ev) ~ trt + sex, data = df, L = d$L)
  expect_identical(coef(f_frm), coef(f_mat))

  p_mat <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  p_frm <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + trt, data = df)
  expect_identical(coef(p_frm), coef(p_mat))

  # A 3-level factor enters as 2 indicators (no intercept)
  df$g <- factor(sample(c("a", "b", "c"), nrow(df), TRUE))
  p_fac <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + g, data = df)
  expect_named(p_fac$beta, c("idade", "gb", "gc"))

  expect_error(semi_aft_prop(tt ~ idade, data = df), "Qsurv")
  df$tt[2] <- NA
  expect_error(semi_aft_prop(Qsurv(tt, q, ev) ~ idade, data = df),
               "t contains NA")
  expect_error(Qsurv(1:3, 1:2, c(1, 0, 1)), "same length")
})

test_that("qast_fit: print, summary, coef and vcov", {
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

  # Without variance: vcov fails with a message; summary only with coef/exp(coef)
  fit0 <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_error(vcov(fit0), "var_comp = TRUE")
  expect_identical(colnames(summary(fit0)$beta), c("coef", "exp(coef)"))
  expect_output(print(fit0), "no variance")
})

test_that("bootstrap vcov is the covariance of the replicates", {
  d <- mk_api_data(n = 50)
  set.seed(3)
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                     var_method = "bootstrap", M = 20)
  expect_identical(fit$var_method, "bootstrap")
  expect_identical(fit$M, 20)
  expect_equal(diag(vcov(fit, "theta")), fit$var_theta, tolerance = 1e-12)
})

test_that("single covariate: methods without 1D support use subplex", {
  d <- mk_api_data()
  z1 <- d$z[, "trt", drop = FALSE]
  expect_message(fit <- semi_aft_prop(d$t, d$qt, d$status, z1),
                 "using \"subplex\"")
  expect_identical(fit$settings$method, "subplex")
  expect_true(all(is.finite(coef(fit))))
  expect_message(semi_aft_fg(d$t, d$qt, d$status, z1, L = d$L), "subplex")
  # method that accepts 1D: no switch and no message
  expect_silent(fit_b <- semi_aft_prop(d$t, d$qt, d$status, z1,
                                       method = "BFGS"))
  expect_identical(fit_b$settings$method, "BFGS")
})

test_that("intercepts: KM mean (beta) and IPCW mean (theta) of the residuals", {
  d <- mk_api_data(n = 80, seed = 21)
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_named(fit$intercept, c("beta", "theta"))

  # Independent reference: mean from the jumps of the KM of the residuals
  # and IPCW mean with the censoring KM evaluated at t
  e <- log(d$t) - drop(d$x %*% fit$beta)
  km <- survival::survfit(survival::Surv(e, d$status) ~ 1)
  s <- km$surv
  a_t <- sum(-diff(c(1, s)) * km$time) + s[length(s)] * km$time[length(s)]
  r <- log(d$qt) - drop(d$x %*% fit$theta)
  G <- summary(survival::survfit(survival::Surv(d$t, 1 - d$status) ~ 1),
               times = d$t)$surv
  # a censored subject with G = 0 gives 0/0; as in the thesis script, it gets weight 0
  a_q <- sum(d$status * r / G, na.rm = TRUE) / length(r)
  expect_equal(unname(fit$intercept), c(a_t, a_q), tolerance = 1e-12)

  expect_output(print(fit), "Intercepts")
  expect_output(print(summary(fit)), "IPCW mean of the residuals")
})

test_that("intercepts without censoring: plain means, also with negative residuals", {
  d <- mk_api_data(n = 60, seed = 5)
  st <- rep(1, length(d$t))
  t <- d$t / 100; qt <- d$qt / 100       # negative log-times
  fit <- semi_aft_fg(t, qt, st, d$z, L = max(t))
  e <- log(t) - drop(d$z %*% fit$beta)
  r <- log(qt) - drop(d$z %*% fit$theta)
  expect_true(all(e < 0))
  expect_equal(fit$intercept[["beta"]], mean(e), tolerance = 1e-12)
  expect_equal(fit$intercept[["theta"]], mean(r), tolerance = 1e-12)
})

test_that("predict lp: intercept + z'coef, formula and matrix", {
  d <- mk_api_data(n = 70, seed = 31)
  df <- data.frame(tt = d$t, q = d$qt, ev = d$status, d$x)
  fit <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + trt, data = df)
  lp <- predict(fit)
  expect_identical(dim(lp), c(70L, 2L))
  expect_equal(unname(lp[, "beta"]),
               unname(fit$intercept[["beta"]] + drop(d$x %*% fit$beta)))
  new <- data.frame(idade = c(-1, 1), trt = c(1, 0))
  expect_equal(unname(predict(fit, new, component = "theta")),
               unname(fit$intercept[["theta"]] +
                        drop(as.matrix(new) %*% fit$theta)))

  # Matrix: columns matched by name, in any order
  fm <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_equal(predict(fm, d$x[1:3, c("trt", "idade")], component = "beta"),
               predict(fm, d$x[1:3, ], component = "beta"))
  expect_error(predict(fm, d$x[, "trt", drop = FALSE]), "lacks the columns")
})

test_that("predict survival: S(t | z) from the distribution of the residuals", {
  d <- mk_api_data(n = 90, seed = 41)
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x)

  # beta: KM of the residuals evaluated at each residual (thesis script)
  e <- log(d$t) - drop(d$x %*% fit$beta)
  km <- survival::survfit(survival::Surv(e, d$status) ~ 1)
  ref_t <- summary(km, times = e)$surv
  S <- predict(fit, type = "survival", component = "beta", times = d$t)
  expect_equal(diag(S), ref_t, tolerance = 1e-12)

  # theta: Zhao and Tsiatis on the residuals, P(r > x), in [0, 1]
  r <- log(d$qt) - drop(d$x %*% fit$theta)
  G <- summary(survival::survfit(survival::Surv(d$t, 1 - d$status) ~ 1),
               times = d$t)$surv
  ok <- d$status == 1 & G > 0
  tg <- c(0.2, 0.7, 1.5)
  Sq <- predict(fit, d$x[1:4, ], type = "survival", component = "theta",
                times = tg)
  for (k in 1:4) {
    x <- log(tg) - sum(d$x[k, ] * fit$theta)
    ref <- vapply(x, function(xx) sum((r > xx)[ok] / G[ok]) / length(r), 0)
    expect_equal(unname(Sq[k, ]), pmin(ref, 1), tolerance = 1e-12)
  }

  # Non-increasing curves in t and within [0, 1]; default grid = events
  both <- predict(fit, d$x[1:5, ], type = "survival")
  for (m in both) {
    expect_true(all(m >= 0 & m <= 1))
    expect_true(all(apply(m, 1, function(v) all(diff(v) <= 1e-12))))
  }
  expect_identical(ncol(both$beta), length(unique(d$t[d$status == 1])))
  expect_error(predict(fit, type = "survival", times = -1), "domain")
})

test_that("predict with a factor: levels of the fit and error for a new level", {
  d <- mk_api_data(n = 80, seed = 51)
  df <- data.frame(tt = d$t, q = d$qt, ev = d$status, idade = d$x[, "idade"],
                   g = factor(rep(c("a", "b", "c"), length.out = 80)))
  fit <- semi_aft_prop(Qsurv(tt, q, ev) ~ idade + g, data = df)
  new <- data.frame(idade = 0, g = c("c", "a"))
  lp <- predict(fit, new, component = "beta")
  expect_equal(unname(lp), unname(fit$intercept[["beta"]] +
                                    c(fit$beta[["gc"]], 0)))
  expect_error(predict(fit, data.frame(idade = 0, g = "z")))
})

test_that("cens_group: validation, default partition and warning for H_G(L) = 0", {
  d <- mk_api_data(n = 80, seed = 61)
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L)
  expect_s3_class(fit$settings$cens_group, "factor")
  expect_identical(nlevels(fit$settings$cens_group), 1L)
  expect_true(fit$settings$cens_group_default)
  expect_output(print(fit), "censoring groups = 1 \\(no stratification\\)")
  # equal to passing a constant
  f1 <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L,
                    cens_group = rep(1, 80))
  expect_equal(coef(f1), coef(fit))
  expect_output(print(f1), "censoring groups = 1\\n")

  fj <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L,
                    cens_group = paste(d$z[, 1], d$z[, 2]))
  expect_identical(nlevels(fj$settings$cens_group), nrow(unique(d$z)))

  # the FG intercept uses the same H_G as the equations
  r <- log(d$qt) - drop(d$z %*% fj$theta)
  H <- QAST:::km_weights_by_group(d$t, d$status, fj$settings$cens_group)
  expect_equal(fj$intercept[["theta"]],
               sum(ifelse(d$status == 1, r / H, 0)) / length(r),
               tolerance = 1e-12)

  expect_error(semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L,
                           cens_group = 1:3), "one value per observation")
  g <- d$z[, "trt"]; g[4] <- NA
  expect_error(semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L,
                           cens_group = g), "cens_group contains NA")

  # L beyond the last censoring of a group with the censoring at 0
  st <- d$status; t <- d$t
  i <- which(d$z[, "trt"] == 1)
  last <- i[which.max(t[i])]; st[last] <- 0
  expect_warning(semi_aft_fg(t, d$qt, st, d$z, L = max(t) * 2,
                             cens_group = d$z[, "trt"]), "H_G\\(L\\) > 0")
})

test_that("predict theta == surv_qt_fun at the own residuals (convention >)", {
  d <- mk_api_data(n = 120, seed = 71)
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  s_fun <- surv_qt_fun(d$qt, d$t, d$status,
                       lin_pred = drop(d$x %*% fit$theta))
  S <- predict(fit, type = "survival", component = "theta", times = d$qt)
  expect_equal(diag(S), pmin(s_fun, 1), tolerance = 1e-12)
  # distribution values stored without rounding
  r <- log(d$qt) - drop(d$x %*% fit$theta)
  expect_true(all(fit$resid_dist$theta$value %in% r))
})

test_that("predict se.fit: jackknife and bootstrap replicates", {
  d <- mk_api_data(n = 70, seed = 81)
  fit <- semi_aft_prop(d$t, d$qt, d$status, d$x, var_comp = TRUE)
  expect_error(predict(semi_aft_prop(d$t, d$qt, d$status, d$x),
                       se.fit = TRUE), "var_comp = TRUE")
  expect_error(predict(fit, se.fit = NA), "TRUE or FALSE")
  expect_error(predict(fit, type = "survival", se.fit = TRUE, level = 2),
               "level")

  zn <- d$x[1:4, ]
  p <- predict(fit, zn, se.fit = TRUE)
  expect_named(p, c("fit", "se.fit"))
  expect_identical(dim(p$se.fit), c(4L, 2L))
  expect_equal(p$fit, predict(fit, zn))

  # without the intercept, the standard error of the replicates is sqrt(z V z')
  reps <- QAST:::.qast_rep_list(fit, "theta")
  vals <- vapply(reps, function(r) drop(zn %*% r$coef), numeric(4))
  expect_equal(QAST:::.qast_rep_se(fit, vals),
               sqrt(diag(zn %*% vcov(fit, "theta") %*% t(zn))),
               tolerance = 1e-12)
  # replicate with all rows and the coefficient of the fit == fit
  expect_identical(QAST:::.qast_dist_rows(fit, "theta", seq_len(70),
                                          fit$theta),
                   fit$resid_dist$theta)

  s <- predict(fit, zn, type = "survival", times = c(0.3, 1), se.fit = TRUE)
  for (cp in c("beta", "theta")) {
    expect_named(s[[cp]], c("surv", "se", "lower", "upper"))
    expect_equal(s[[cp]]$surv,
                 predict(fit, zn, type = "survival", times = c(0.3, 1))[[cp]])
    expect_true(all(s[[cp]]$lower <= s[[cp]]$surv + 1e-12 &
                      s[[cp]]$surv <= s[[cp]]$upper + 1e-12))
    expect_true(all(s[[cp]]$lower >= 0 & s[[cp]]$upper <= 1))
  }

  # bootstrap: regenerated indices == resamples of each replicate
  set.seed(4)
  fb <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                    var_method = "bootstrap", M = 8)
  ix <- fb$replicates$beta$index
  expect_identical(dim(ix), c(8L, 70L))
  b1 <- dfoptim::nmk(unname(fb$beta), function(b)
    QAST:::U1_obj_cpp(b, log(d$t[ix[1, ]]), d$z[ix[1, ], ] * 1,
                      as.numeric(d$status[ix[1, ]]), "L2"))$par
  expect_equal(b1, fb$replicates$beta$reps[1, ], tolerance = 1e-10)
  pb <- predict(fb, d$z[1:2, ], type = "survival", component = "theta",
                times = 1, se.fit = TRUE)
  expect_true(all(pb$se >= 0))
})
