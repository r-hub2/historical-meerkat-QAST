# P-P plots of the residuals (Fine and Gelber, 2001, Section 4)

mk_pp_data <- function(n = 240, seed = 3) {
  set.seed(seed)
  d <- data.frame(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4),
                  age = rnorm(n))
  d$t <- exp(0.4 * d$trt - 0.2 * d$age + rnorm(n))
  d$status <- rbinom(n, 1, .75)
  d$qt <- d$t * runif(n, .4, 1)
  d
}

test_that("FG: KM per group (beta) and the article formula with trimming (theta)", {
  d <- mk_pp_data()
  L <- unname(quantile(d$t, .8))
  fit <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d, L = L)
  pp <- qast_pp(fit, d$sex)
  expect_s3_class(pp, "qast_pp")
  z <- as.matrix(d[, c("trt", "sex")])
  e1 <- log(d$t) - drop(z %*% fit$beta)
  e2 <- log(d$qt) - drop(z %*% fit$theta)

  km <- survival::survfit(survival::Surv(e1, d$status) ~ d$sex)
  tt <- pp$beta[["1"]]$t[-1]
  s <- summary(km, times = tt, extend = TRUE)
  expect_equal(pp$beta[["1"]]$ref[-1], 1 - s$surv[s$strata == "d$sex=0"],
               tolerance = 1e-12)
  expect_equal(pp$beta[["1"]]$group[-1], 1 - s$surv[s$strata == "d$sex=1"],
               tolerance = 1e-12)

  trim <- log(L) - max(z %*% fit$beta)
  expect_equal(pp$trim, trim)
  # H_{G_j}(X_j): censoring KM in the partition of the fit (default: no
  # stratification, as in the application of the article)
  H <- summary(survival::survfit(
    survival::Surv(d$t, 1 - d$status) ~ 1), times = d$t)$surv
  paper <- function(s, x) {
    i <- d$sex == s
    A <- d$status == 1 & e1 < trim & i
    vapply(x, function(v) sum((1 / H)[A & e2 < v]) / sum((1 / H)[A]), 0)
  }
  tt2 <- pp$theta[["1"]]$t
  expect_identical(tail(tt2, 1), Inf)
  expect_equal(pp$theta[["1"]]$ref, paper(0, tt2), tolerance = 1e-12)
  expect_equal(pp$theta[["1"]]$group, paper(1, tt2), tolerance = 1e-12)

  # Stratified: H_G of the cens_group of the fit
  fs <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d, L = L,
                    cens_group = d$trt)
  ps <- qast_pp(fs, d$sex)
  Hs <- QAST:::km_weights_by_group(d$t, d$status, factor(d$trt))
  e1s <- log(d$t) - drop(z %*% fs$beta)
  e2s <- log(d$qt) - drop(z %*% fs$theta)
  As <- d$status == 1 & e1s < ps$trim & d$sex == 1
  x <- ps$theta[["1"]]$t
  expect_equal(ps$theta[["1"]]$group,
               vapply(x, function(v) sum((1 / Hs)[As & e2s < v]) /
                        sum((1 / Hs)[As]), 0), tolerance = 1e-12)
  expect_equal(tail(pp$theta[["1"]]$ref, 1), 1)
  expect_equal(unname(pp$n_used),
               c(sum(d$status == 1 & e1 < trim & d$sex == 0),
                 sum(d$status == 1 & e1 < trim & d$sex == 1)))
})

test_that("proposed: no trimming and the censoring weights of the fit (Cox)", {
  d <- mk_pp_data(seed = 8)
  fit <- semi_aft_prop(Qsurv(t, qt, status) ~ trt + age, data = d)
  grp <- cut(d$age, c(-Inf, 0, Inf), labels = c("young", "old"))
  pp <- qast_pp(fit, grp, ref = "old")
  expect_null(pp$trim)
  expect_identical(pp$ref, "old")
  expect_named(pp$theta, "young")

  z <- as.matrix(d[, c("trt", "age")])
  e2 <- log(d$qt) - drop(z %*% fit$theta)
  H <- QAST:::HG(X = z, t = d$t, status = d$status)
  i <- grp == "young" & d$status == 1
  tt <- pp$theta$young$t
  ref <- vapply(tt, function(v) sum((1 / H)[i & e2 < v]) / sum((1 / H)[i]), 0)
  expect_equal(pp$theta$young$group, ref, tolerance = 1e-12)
})

test_that("several groups, print, plot and validation", {
  d <- mk_pp_data(seed = 5)
  fit <- semi_aft_prop(Qsurv(t, qt, status) ~ trt + age, data = d)
  g3 <- cut(d$age, quantile(d$age, 0:3 / 3), include.lowest = TRUE)
  pp <- qast_pp(fit, g3)
  expect_length(pp$beta, 2L)
  expect_output(print(pp), "Fine and Gelber")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(pp))
  expect_silent(plot(pp, component = "theta"))

  expect_error(qast_pp(fit, d$trt[-1]), "one value per observation")
  expect_error(qast_pp(fit, rep(1, nrow(d))), "at least 2 levels")
  expect_error(qast_pp(fit, d$trt, ref = "9"), "ref must be one of")
  gna <- d$trt; gna[2] <- NA
  expect_error(qast_pp(fit, gna), "group contains NA")
})
