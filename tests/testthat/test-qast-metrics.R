# predict(type = "rmean") and qast_metrics(): exact restricted means,
# IPCW Brier scores on the survival time scale and the null models.

mk_met_data <- function(n = 120, seed = 41) {
  set.seed(seed)
  d <- data.frame(age = rnorm(n), trt = rbinom(n, 1, 0.5))
  d$t <- rexp(n, rate = exp(-0.5 * d$trt))
  d$status <- rbinom(n, 1, 0.75)
  d$qt <- d$t * runif(n, 0.4, 1)
  d
}

test_that("rmean is the exact integral of the predicted step curves", {
  d <- mk_met_data()
  fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d)
  new <- d[1:4, ]
  tau <- 1.5
  for (cp in c("beta", "theta")) {
    rm <- predict(fit, new, type = "rmean", component = cp, tau = tau)
    expect_identical(attr(rm, "tau"), stats::setNames(tau, cp))
    # midpoint rule on a fine grid: error O(h)
    h <- tau / 20000
    s <- predict(fit, new, type = "survival", component = cp,
                 times = seq(h / 2, tau, by = h))
    expect_equal(as.vector(rm), unname(rowSums(s) * h), tolerance = 1e-4)
  }
  # tau before the first jump of every curve: S = 1, mean = tau
  tiny <- 1e-6
  expect_equal(as.vector(predict(fit, new, type = "rmean",
                                 component = "beta", tau = tiny)),
               rep(tiny, 4))
  # increasing in tau; default tau = largest event time of each component
  r <- sapply(c(0.5, 1, 2), function(ta)
    predict(fit, new[1, ], type = "rmean", component = "beta", tau = ta))
  expect_true(all(diff(r) > 0))
  both <- predict(fit, new, type = "rmean")
  expect_identical(colnames(both), c("beta", "theta"))
  expect_equal(attr(both, "tau"),
               c(beta = max(d$t[d$status == 1]),
                 theta = max(d$qt[d$status == 1])))
  expect_equal(predict(fit, new, type = "rmean", tau = c(theta = 1, beta = 2)),
               cbind(beta = predict(fit, new, type = "rmean",
                                    component = "beta", tau = 2),
                     theta = predict(fit, new, type = "rmean",
                                     component = "theta", tau = 1)),
               ignore_attr = TRUE)
  expect_error(predict(fit, new, type = "rmean", tau = -1), "positive")
  expect_error(predict(fit, new, type = "rmean", tau = c(1, 2)), "named")
})

test_that("rmean with standard errors from the replicates", {
  d <- mk_met_data(n = 70)
  fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d,
                       var_comp = TRUE)
  p <- predict(fit, d[1:3, ], type = "rmean", tau = 1, se.fit = TRUE)
  expect_named(p, c("fit", "se.fit"))
  expect_identical(dim(p$se.fit), c(3L, 2L))
  expect_true(all(p$se.fit > 0))
  expect_equal(p$fit, predict(fit, d[1:3, ], type = "rmean", tau = 1))
})

test_that("rmean of the Fine and Gelber model", {
  d <- mk_met_data()
  z <- cbind(trt = d$trt, old = as.numeric(d$age > 0))
  fit <- semi_aft_fg(d$t, d$qt, d$status, z, L = quantile(d$t, 0.9))
  rm <- predict(fit, z[1:3, ], type = "rmean", tau = 1)
  expect_true(all(rm > 0 & rm <= 1))
})

test_that("without censoring the metrics are the plain averages", {
  d <- mk_met_data()
  d$status <- 1
  y <- with(d, Qsurv(t, qt, status))
  set.seed(2)
  tt <- c(0.2, 0.5, 1, 2)
  S <- matrix(runif(nrow(d) * 4), nrow(d))
  for (oc in c("time", "qtime")) {
    obs <- if (oc == "time") d$t else d$qt
    m <- qast_metrics(S, y = y, outcome = oc, times = tt, tau = 2,
                      pred_time = rep(0.7, nrow(d)))
    expect_equal(m$brier$brier, colMeans((outer(obs, tt, ">") - S)^2))
    expect_equal(m$mae, mean(abs(pmin(obs, 2) - 0.7)))
    expect_equal(m$ipcw_mass, 1)
  }
})

test_that("the null model reproduces the reference scores", {
  d <- mk_met_data()
  y <- with(d, Qsurv(t, qt, status))
  n <- nrow(d)
  # survival time: Kaplan-Meier predictions
  km <- survival::survfit(survival::Surv(t, status) ~ 1, data = d)
  tt <- sort(unique(d$t[d$status == 1]))
  S <- matrix(summary(km, times = tt)$surv, n, length(tt), byrow = TRUE)
  m <- qast_metrics(S, y = y, times = tt)
  expect_equal(m$brier$brier, m$brier$brier_null)
  expect_equal(m$ibs, m$ibs_null)
  expect_false(any(c("r2", "r2_ibs") %in% c(names(m), names(m$brier))))
  # without predicted means, no mean absolute error
  expect_true(is.na(m$mae))
  expect_null(m$pred_time)
  expect_false(any(grepl("absolute", capture.output(print(m)))))
  # quality-adjusted time: Zhao and Tsiatis curve
  tq <- sort(unique(d$qt[d$status == 1]))
  zt <- QAST:::.surv_ipcw_at(QAST:::.resid_dist_ipcw(d$qt, d$t, d$status),
                              tq)
  Sq <- matrix(zt, n, length(tq), byrow = TRUE)
  mq <- qast_metrics(Sq, y = y, outcome = "qtime", times = tq)
  expect_equal(mq$brier$brier, mq$brier$brier_null)
})

test_that("Brier of the survival time follows Graf et al. (1999)", {
  d <- mk_met_data()
  y <- with(d, Qsurv(t, qt, status))
  set.seed(3)
  tt <- c(0.3, 0.8, 1.5)
  S <- matrix(runif(nrow(d) * 3), nrow(d))
  m <- qast_metrics(S, y = y, times = tt)
  G <- survival::survfit(survival::Surv(t, 1 - status) ~ 1, data = d)
  Gat <- function(x) stats::stepfun(G$time, c(1, G$surv))(x)
  Gleft <- function(x) sapply(x, function(u) {
    k <- G$time < u
    if (any(k)) G$surv[max(which(k))] else 1
  })
  ref <- sapply(seq_along(tt), function(j) {
    s <- S[, j]; u <- tt[j]
    died <- d$t <= u & d$status == 1
    mean(ifelse(died, s^2 / Gleft(d$t), 0) +
           ifelse(d$t > u, (1 - s)^2 / Gat(u), 0))
  })
  expect_equal(m$brier$brier, ref)
  # quality-adjusted time: uncensored only, weight 1 / G(T_i)
  mq <- qast_metrics(S, y = y, outcome = "qtime", times = tt)
  w <- d$status / Gat(d$t)
  refq <- colMeans(w * (outer(d$qt, tt, ">") - S)^2)
  expect_equal(mq$brier$brier, refq)
})

test_that("qast_fit method == matrix method with the predictions", {
  d <- mk_met_data(n = 150)
  tr <- d[1:100, ]; te <- d[101:150, ]
  fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = tr)
  y <- with(te, Qsurv(t, qt, status))
  m <- qast_metrics(fit, newdata = te, y = y)
  expect_s3_class(m, "qast_metrics_list")
  expect_named(m, c("beta", "theta"))
  tq <- sort(unique(te$qt[te$status == 1]))
  S <- predict(fit, te, type = "survival", component = "theta", times = tq)
  pt <- predict(fit, te, type = "rmean", component = "theta", tau = max(tq))
  ref <- qast_metrics(unname(S), y = y, outcome = "qtime", times = tq,
                      pred_time = as.vector(pt))
  expect_equal(m$theta$brier, ref$brier)
  expect_equal(m$theta$mae, ref$mae)
  expect_output(print(m), "quality-adjusted survival time")
  # in-sample (apparent) performance and a single component
  m1 <- qast_metrics(fit, component = "beta", tau = 2)
  expect_s3_class(m1, "qast_metrics")
  expect_identical(m1$tau, 2)
  expect_true(all(m1$brier$time <= 2))

  expect_error(qast_metrics(fit, newdata = te), "Qsurv")
  expect_error(qast_metrics(fit, y = y), "newdata")
  expect_error(qast_metrics(fit, newdata = te, y = y[1:10, ]), "rows")
  expect_error(qast_metrics(matrix(0.5, 50, 2), y = y, times = 1),
               "columns")
  expect_error(qast_metrics(matrix(2, 50, 1), y = y, times = 1), "\\[0, 1\\]")
})
