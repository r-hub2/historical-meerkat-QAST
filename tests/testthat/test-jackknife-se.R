# .jackknife_se reproduces jackknifeR::jackknife(d = 1) — same
# statistic, same SE formula, NA removal per block.

mk_jk_case <- function(n = 30, seed = 7) {
  set.seed(seed)
  z <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n, exp(.3 * z[, 1])); st <- rbinom(n, 1, .8)
  list(n = n, z = z, lt = log(t), st = st)
}

refit_case <- function(cas, keep) {
  f <- function(b) QAST:::U1_obj_cpp(b, cas$lt[keep],
                                     cas$z[keep, , drop = FALSE],
                                     cas$st[keep], "L2")
  dfoptim::nmk(c(0, 0), f)$par
}

test_that(".jackknife_se == jackknifeR::jackknife (sequential and fork)", {
  skip_if_not_installed("jackknifeR")
  cas <- mk_jk_case()
  dat <- data.frame(id = seq_len(cas$n), pad = 0)
  ref <- suppressWarnings(jackknifeR::jackknife(
    statistic = function(d) refit_case(cas, d$id),
    data = dat, d = 1, numCores = 1))$se

  se1 <- QAST:::.jackknife_se(cas$n, function(i) refit_case(cas, -i),
                              blocks = list(b = 1:2), cores = 1L)$b
  expect_equal(se1, ref, tolerance = 1e-12)

  skip_on_os("windows")
  se2 <- QAST:::.jackknife_se(cas$n, function(i) refit_case(cas, -i),
                              blocks = list(b = 1:2), cores = 2L)$b
  expect_identical(se2, se1)
})

test_that(".jackknife_se: independent blocks (NA removed per block)", {
  # Block 'a' with NA in replicate 3 does not affect block 'b'
  vals <- function(i) c(i^2, if (i == 3) NA else sqrt(i), i)
  se <- QAST:::.jackknife_se(6, vals, blocks = list(a = 1, b = 2:3),
                             cores = 1L)
  m <- t(sapply(1:6, vals))
  jk_se <- function(x) { x <- x[stats::complete.cases(x), , drop = FALSE]
    sqrt((6 - 1) / nrow(x) * colSums(sweep(x, 2, colMeans(x))^2)) }
  expect_equal(se$a, jk_se(m[, 1, drop = FALSE]))
  expect_equal(se$b, jk_se(m[, 2:3, drop = FALSE]))
})

test_that(".jackknife_se: errors and warnings of the workers reach the session", {
  expect_error(
    QAST:::.jackknife_se(4, function(i) if (i == 2) stop("falhou") else i,
                         blocks = list(a = 1), cores = 1L),
    "falhou")
  expect_warning(
    QAST:::.jackknife_se(4, function(i) { warning("aviso-x"); i },
                         blocks = list(a = 1), cores = 1L),
    "aviso-x")
  skip_on_os("windows")
  expect_error(
    QAST:::.jackknife_se(4, function(i) if (i == 2) stop("falhou") else i,
                         blocks = list(a = 1), cores = 2L),
    "falhou")
})
