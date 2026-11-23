# Variance replicates (jackknife and bootstrap): every replicate refits the
# whole estimator on the same subsample/resample, with the beta of the
# replicate in the theta equation of the Fine and Gelber model, and
# degenerate replicates are discarded.

mk_rep_data <- function(n = 60, seed = 31) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
  t <- rexp(n, rate = exp(-0.4 * z[, "trt"]))
  status <- rbinom(n, 1, 0.8)
  qt <- t * runif(n, 0.5, 1)
  x <- cbind(age = rnorm(n), trt = z[, "trt"])
  list(t = t, qt = qt, status = status, z = z, x = x, L = quantile(t, 0.9))
}

# theta of the Fine and Gelber model on the rows 'rows' with a given beta
# (closed form, phi = "identity", norm_min = "L2")
fg_theta_rows <- function(d, rows, beta) {
  z <- d$z[rows, , drop = FALSE] * 1
  t <- d$t[rows]; st <- as.numeric(d$status[rows])
  km <- QAST:::prep_km_by_column(z, t, st,
                                 group = QAST:::fg_cens_partition(z))
  pairs <- QAST:::make_fg_pairs(log(t), log(d$qt[rows]), z, st, beta,
                                log(d$L), km, "identity")
  aff <- QAST:::make_fg_affine(pairs, log(d$qt[rows]), z)
  QAST:::solve_fg_theta_affine(aff, rcond = 1e-10)$par
}

test_that("bootstrap: beta and theta share the resamples", {
  d <- mk_rep_data()
  set.seed(8)
  fit <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                     var_method = "bootstrap", M = 12)
  expect_identical(fit$replicates$beta$index, fit$replicates$theta$index)
  # same resamples as boot::boot(sim = "ordinary") with the same seed
  set.seed(8)
  expect_identical(fit$replicates$beta$index,
                   matrix(sample.int(60, 60 * 12, replace = TRUE), 12))
  expect_equal(unname(vcov(fit, "theta")), cov(fit$replicates$theta$reps),
               tolerance = 1e-12)
})

test_that("Fine and Gelber: theta of each replicate uses its own beta", {
  d <- mk_rep_data()
  set.seed(9)
  fb <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                    var_method = "bootstrap", M = 6)
  for (k in 1:3) {
    rows <- fb$replicates$theta$index[k, ]
    b_k <- fb$replicates$beta$reps[k, ]
    expect_equal(fb$replicates$theta$reps[k, ], fg_theta_rows(d, rows, b_k),
                 tolerance = 1e-10)
    # differs from the theta with the beta of the full sample
    expect_false(isTRUE(all.equal(fb$replicates$theta$reps[k, ],
                                  fg_theta_rows(d, rows, fb$beta))))
  }

  fj <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE)
  i <- fj$replicates$theta$drop[5]
  expect_equal(fj$replicates$theta$reps[5, ],
               fg_theta_rows(d, seq_len(60)[-i], fj$replicates$beta$reps[5, ]),
               tolerance = 1e-10)
})

collect_warnings <- function(expr) {
  w <- character()
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("replicates whose covariates lose rank are discarded", {
  d <- mk_rep_data(n = 40)
  # jackknife, proposed model: a single subject with trt = 1, so removing
  # it leaves a constant column
  d$x[, "trt"] <- 0
  d$x[7, "trt"] <- 1
  expect_warning(
    fp <- semi_aft_prop(d$t, d$qt, d$status, d$x, var_comp = TRUE),
    "1 of 40 jackknife replicates discarded"
  )
  expect_identical(nrow(fp$replicates$beta$reps), 39L)
  expect_false(7L %in% fp$replicates$beta$drop)
  expect_false(7L %in% fp$replicates$theta$drop)
  expect_true(all(is.finite(fp$var_beta)) && all(is.finite(fp$var_theta)))

  # bootstrap, both models: resamples without the subject are discarded
  set.seed(2)
  r <- collect_warnings(semi_aft_prop(d$t, d$qt, d$status, d$x,
                                      var_comp = TRUE,
                                      var_method = "bootstrap", M = 30))
  kept <- apply(r$value$replicates$beta$index, 1, function(i) 7L %in% i)
  expect_true(all(kept))
  nd <- 30L - nrow(r$value$replicates$beta$reps)
  expect_gt(nd, 0L)
  expect_true(any(grepl(sprintf("^%d of 30 bootstrap replicates discarded",
                                nd), r$warnings)))

  # Fine and Gelber: three subjects with sex = 1
  d$z[, "sex"] <- 0
  d$z[c(3, 7, 11), "sex"] <- 1
  set.seed(3)
  r <- collect_warnings(semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L,
                                    var_comp = TRUE,
                                    var_method = "bootstrap", M = 60))
  ix <- r$value$replicates$theta$index
  expect_true(all(apply(ix, 1, function(i) any(i %in% c(3, 7, 11)))))
  expect_true(any(grepl("of 60 bootstrap replicates discarded", r$warnings)))
  # warnings of the replicates are reported once each
  expect_identical(anyDuplicated(r$warnings), 0L)
})

test_that(".qast_bootstrap: blocks with NA are discarded separately", {
  set.seed(1)
  vr <- QAST:::.qast_bootstrap(5, 10, function(rows)
    c(mean(rows), if (any(rows == 1L)) NA else sum(rows)),
    blocks = list(a = 1, b = 2), cores = 1L)
  expect_identical(vr$a$n_drop, 0L)
  expect_identical(nrow(vr$b$reps) + vr$b$n_drop, 10L)
  expect_false(any(vr$b$index == 1L))
  expect_equal(vr$a$cov[1, 1], var(rowMeans(vr$a$index)))
})

test_that("theta - beta: joint covariance, confint and summary", {
  d <- mk_rep_data(n = 70)
  set.seed(5)
  fb <- semi_aft_fg(d$t, d$qt, d$status, d$z, L = d$L, var_comp = TRUE,
                    var_method = "bootstrap", M = 30)
  fj <- semi_aft_prop(d$t, d$qt, d$status, d$x, var_comp = TRUE)
  for (f in list(fb, fj)) {
    dd <- f$replicates$diff$reps
    expect_equal(coef(f, "diff"), f$theta - f$beta)
    # covariance of the differences of the same replicates
    ref <- if (f$var_method == "bootstrap") cov(dd) else {
      (f$n - 1) / nrow(dd) * crossprod(sweep(dd, 2, colMeans(dd)))
    }
    expect_equal(unname(vcov(f, "diff")), ref, tolerance = 1e-12)
    # = Var(theta) + Var(beta) - 2 Cov(beta, theta), from the joint rows
    rb <- f$replicates$beta$reps; rt <- f$replicates$theta$reps
    if (f$var_method == "bootstrap") {
      expect_equal(unname(f$var_diff),
                   apply(rt, 2, var) + apply(rb, 2, var) -
                     2 * diag(cov(rb, rt)), tolerance = 1e-12)
    }
    # Wald interval
    ci <- confint(f, component = "diff", level = 0.9)
    se <- sqrt(f$var_diff)
    expect_equal(unname(ci[, 1]), unname(coef(f, "diff") - qnorm(0.95) * se))
    expect_identical(colnames(ci), c("5 %", "95 %"))
    s <- summary(f)
    expect_equal(unname(s$ratio[, 1]), unname(exp(coef(f, "diff"))))
    expect_equal(unname(s$ratio[, 2]),
                 unname(exp(confint(f, component = "diff")[, 1])))
    expect_output(print(s), "Quality per unit of survival")
  }
  # percentile interval: bootstrap only
  pc <- confint(fb, "trt", component = "diff", type = "percentile")
  expect_equal(unname(pc[1, ]),
               unname(quantile(fb$replicates$diff$reps[, 1], c(0.025, 0.975))))
  expect_error(confint(fj, component = "diff", type = "percentile"),
               "bootstrap")
  expect_error(confint(fj, "nope"), "unknown")
  # beta and theta intervals use their own variances
  expect_equal(unname(confint(fj, component = "beta")[, 2]),
               unname(fj$beta + qnorm(0.975) * sqrt(fj$var_beta)))
  # without variance
  f0 <- semi_aft_prop(d$t, d$qt, d$status, d$x)
  expect_error(confint(f0), "var_comp")
  expect_null(summary(f0)$ratio)
})
