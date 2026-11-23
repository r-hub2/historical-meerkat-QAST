library(testthat)

mk_pairs_case <- function(seed = 11, n = 60, p = 2, tie_at_L = FALSE) {
  set.seed(seed)
  z <- cbind(rbinom(n, 1, 0.45), rbinom(n, 1, 0.35))
  colnames(z) <- c("trt", "sex")
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  L      <- stats::quantile(t, 0.9)
  if (tie_at_L) {
    # Exact tie: beta = 0 => z'beta = 0 exactly;
    # times duplicated exactly at L; log() of the same double.
    t[1:4] <- L
    status[1:4] <- 1L
  }
  list(t = t, qt = qt, status = status, z = z,
       lt = log(t), lq = log(qt), hL = log(L))
}

pairs_theta_grid <- rbind(c(0.0, 0.0), c(0.2, -0.1),
                          c(-0.4, 0.5), c(0.7, 0.7))
pairs_beta_grid  <- rbind(c(0.0, 0.0), c(0.3, -0.2), c(-0.5, 0.8))


test_that("Pairs (identity): equal the legacy kernel and the reference on a grid", {
  skip_if_not(exists("U2_parts_fg_ref", envir = asNamespace("QAST")),
              "referencia R indisponivel")
  
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      
      cas   <- mk_pairs_case(seed = 10 + bi)
      beta  <- pairs_beta_grid[bi, ]
      theta <- pairs_theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity")
      
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, "identity")
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "identity")
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
      
      # Triangulation: reference kernel vs pure R reference
      parts_ref <- U2_parts_fg_ref(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "identity")
      expect_equal(parts_leg, parts_ref, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_new <- QAST:::eval_fg_obj_cpp(pairs, cas$lq, cas$z,
                                          theta, "identity", nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq,
                                        cas$z, cas$status, beta,
                                        cas$hL, km, "identity", nm)
        expect_equal(obj_new, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Pairs (piecewise): equal the legacy kernel on a grid of beta x theta", {
  # piecewise phi freezes the h1L INDICATOR per pair at construction, but
  # the sign of e2[pa]-e2[pb] is re-evaluated at each theta — here the
  # pairs path must reproduce the reference kernel (same summation order,
  # same division by n2 at the end), hence a tolerance of 1e-10.
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      
      cas   <- mk_pairs_case(seed = 50 + bi)
      beta  <- pairs_beta_grid[bi, ]
      theta <- pairs_theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise")
      
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, "piecewise")
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "piecewise")
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_new <- QAST:::eval_fg_obj_cpp(pairs, cas$lq, cas$z,
                                          theta, "piecewise", nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq,
                                        cas$z, cas$status, beta,
                                        cas$hL, km, "piecewise", nm)
        expect_equal(obj_new, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Pairs: a constant column gives an empty group and parts identically 0", {
  cas <- mk_pairs_case()
  cas$z[, 2] <- 0L   # group 1 empty in the sex column
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  
  for (ph in c("identity", "piecewise")) {
    pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                      cas$status, beta, cas$hL,
                                      km, ph)
    col <- pairs[[2]]
    expect_length(col$kb, 0)
    # Ra aligned with ka (zeros), not empty
    expect_length(col$Ra, length(col$ka))
    expect_true(all(col$Ra == 0))
    expect_length(col$pc, 0)
    
    th <- c(0.2, -0.1)
    parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z, th, ph)
    expect_identical(parts_new[2], 0.0)
    
    parts_leg <- QAST:::U2_parts_fg_cpp(th, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, ph)
    expect_equal(parts_new, parts_leg, tolerance = 1e-10)
  }
})

test_that("Pairs: exact-tie gate in h1L preserved at construction", {
  # The h1L indicator uses a STRICT inequality resolved AT CONSTRUCTION
  # with the same arithmetic as the reference kernel (dot per member — the
  # value is deterministically identical to the dot per pair). This test
  # fails if any floating-point reordering moves a pair across the gate:
  # the parts diverge by a DISCRETE JUMP, not gradually — hence the tight
  # tolerance, without a relative scale.
  cas <- mk_pairs_case(tie_at_L = TRUE)
  beta <- c(0, 0)
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  for (ph in c("identity", "piecewise")) {
    pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                      cas$status, beta, cas$hL,
                                      km, ph)
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      theta <- pairs_theta_grid[ti, ]
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, ph)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, ph)
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
    }
  }
})

test_that("Pairs: invalid arguments stop with a clear message", {
  cas <- mk_pairs_case()
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  
  # invalid phi
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, c(0, 0), cas$hL, km, "banana"),
               "phi")
  
  # beta of the wrong length
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, c(0.1), cas$hL, km, "identity"),
               "beta")
  
  # object without meta evaluated
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    c(0, 0), cas$hL, km, "identity")
  attr(pairs, "meta") <- NULL
  expect_error(QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                        c(0, 0), "identity"),
               "meta")
})
test_that("PERF-1: identity stores no pairs; Ra/Cb bit-identical to piecewise", {
  cas <- mk_pairs_case(seed = 21)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    beta <- pairs_beta_grid[bi, ]
    # algo = "loop": same summation order as piecewise => bit-identical.
    # algo = "auto" (sweep): same sum, another order => 1e-13.
    p_id <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "identity", "loop")
    p_cq <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "identity")
    p_pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "piecewise")
    expect_false(attr(p_id, "meta")$has_pairs)
    expect_true(attr(p_pw, "meta")$has_pairs)
    for (j in seq_along(p_id)) {
      expect_identical(p_id[[j]]$Ra, p_pw[[j]]$Ra)
      expect_identical(p_id[[j]]$Cb, p_pw[[j]]$Cb)
      expect_equal(p_cq[[j]]$Ra, p_pw[[j]]$Ra, tolerance = 1e-13)
      expect_equal(p_cq[[j]]$Cb, p_pw[[j]]$Cb, tolerance = 1e-13)
      expect_length(p_id[[j]]$pc, 0)
      expect_type(p_pw[[j]]$pa, "integer")
    }
    # An identity object refuses the piecewise evaluation (never a silent zero)
    expect_error(QAST:::eval_fg_parts_cpp(p_id, cas$lq, cas$z, c(0, 0),
                                          "piecewise"), "identity")
  }
})

test_that("pairs: all valid cross pairs (Fine and Gelber, 2001)", {
  cas <- mk_pairs_case(seed = 31)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.1, -0.1)
  pr  <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                  beta, cas$hL, km, "piecewise")
  e1 <- cas$lt - drop(cas$z %*% beta); bz <- drop(cas$z %*% beta)
  H <- QAST:::km_weights_by_group(cas$t, cas$status,
                                  QAST:::fg_cens_partition(cas$z))
  for (j in seq_along(pr)) {
    g0 <- which(cas$z[, j] == 0); g1 <- which(cas$z[, j] == 1)
    all_pairs <- expand.grid(a = g0, b = g1)
    ok <- with(all_pairs,
      cas$status[a] == 1 & cas$status[b] == 1 &
      (cas$hL - pmax(bz[a], bz[b])) > pmax(e1[a], e1[b]))
    expected <- all_pairs[ok, ]
    got <- data.frame(a = pr[[j]]$pa + 1L, b = pr[[j]]$pb + 1L)
    expect_identical(nrow(got), nrow(expected))
    key <- function(d) sort(paste(d$a, d$b))
    expect_identical(key(got), key(expected))
    expect_equal(sort(pr[[j]]$pc),
                 sort(1 / (H[expected$a] * H[expected$b])), tolerance = 1e-12)
  }
})

test_that("theta does not depend on the order of the rows", {
  cas <- mk_pairs_case(seed = 32, n = 90)
  L <- exp(cas$hL)
  fit0 <- semi_aft_fg(cas$t, cas$qt, cas$status, cas$z, L = L)
  for (s in 1:3) {
    set.seed(s); o <- sample(length(cas$t))
    fit <- semi_aft_fg(cas$t[o], cas$qt[o], cas$status[o],
                       cas$z[o, , drop = FALSE], L = L)
    expect_equal(fit$theta, fit0$theta, tolerance = 1e-10)
    fpw <- semi_aft_fg(cas$t[o], cas$qt[o], cas$status[o],
                       cas$z[o, , drop = FALSE], L = L, phi = "piecewise")
    if (s == 1) th_pw <- fpw$theta
    expect_equal(fpw$theta, th_pw, tolerance = 1e-6)
  }
})

test_that("semi_aft_fg: piecewise and fast_path paths run (regression B1/B2)", {
  cas <- mk_pairs_case(seed = 41, n = 80)
  args <- list(t = cas$t, qt = cas$qt, status = cas$status, z = cas$z,
               h1 = "log", h2 = "log", method = "Nelder-Mead",
               norm_min = "L2", L = exp(cas$hL), M = 0,
               var_comp = FALSE, var_method = "jackknife")

  # default piecewise == direct optimization over the reference kernel
  # (same objective, same optimizer, same starting value)
  fit_pw <- do.call(semi_aft_fg, c(args, phi = "piecewise"))
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  b   <- unname(fit_pw$beta)
  f_leg <- function(th) QAST:::U2_obj_fg_cpp(th, cas$lt, cas$lq, cas$z,
                                             cas$status, b, cas$hL, km,
                                             "piecewise", "L2")
  th0 <- unname(stats::lm(cas$lq ~ cas$z)$coef[-1])
  th_leg <- dfoptim::nmk(th0, f_leg)$par
  expect_equal(unname(fit_pw$theta), th_leg, tolerance = 1e-8)

  # the closed form is not worse than the optimizer on the same affine objective
  fit_fast <- do.call(semi_aft_fg, c(args, phi = "identity", fast_path = TRUE))
  fit_slow <- do.call(semi_aft_fg, c(args, phi = "identity", fast_path = FALSE))
  expect_length(fit_fast$theta, 2)
  expect_true(all(is.finite(fit_fast$theta)))

  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  pr  <- QAST:::make_fg_pairs(cas$lt, cas$lq, cas$z, cas$status,
                              unname(fit_fast$beta), cas$hL, km, "identity")
  aff <- QAST:::make_fg_affine(pr, cas$lq, cas$z)
  obj <- function(th) QAST:::eval_fg_obj_affine(aff, unname(th), "L2")
  expect_lte(obj(fit_fast$theta), obj(fit_slow$theta) + 1e-10)
})

# ----------------------------------------------------------------------------
# Marginals Ra/Cb in O(N log N) == O(nk*nl) loop
# ----------------------------------------------------------------------------

cdq_vs_loop <- function(lt, lq, z, st, beta, hL, km) {
  a <- QAST:::make_fg_pairs_cpp(lt, lq, z, st, beta, hL, km, "identity", "auto")
  b <- QAST:::make_fg_pairs_cpp(lt, lq, z, st, beta, hL, km, "identity", "loop")
  for (j in seq_along(b)) {
    for (f in c("Ra", "Cb")) {
      # same support (identical set of pairs) and same values
      expect_identical(a[[j]][[f]] == 0, b[[j]][[f]] == 0)
      expect_equal(a[[j]][[f]], b[[j]][[f]], tolerance = 1e-13)
    }
  }
  invisible(attr(a, "meta")$cdq_cols)
}

test_that("PERF-3: the sweep reproduces the loop (ties, cut at L, beta = 0)", {
  for (sd in 1:30) {
    set.seed(sd)
    n <- c(7, 40, 150, 400)[sd %% 4 + 1]; p <- sd %% 3 + 1
    z <- matrix(rbinom(n * p, 1, 0.5), n, p)
    z[1, ] <- 0; z[2, ] <- 1                     # no empty group
    t  <- if (sd %% 3 == 0) sample(1:5, n, TRUE) + 0 else rexp(n)
    st <- rbinom(n, 1, 0.75)
    L  <- if (sd %% 2 == 0) sort(t)[ceiling(0.7 * n)] else max(t)
    beta <- if (sd %% 5 == 0) rep(0, p) else round(rnorm(p), 1)
    km <- QAST:::prep_km_by_column(z, t, st)
    used <- cdq_vs_loop(log(t), log(t * 0.8), z, st, beta, log(L), km)
    expect_true(all(used))
  }
  # exact tie at the L cut
  cas <- mk_pairs_case(tie_at_L = TRUE)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, c(0, 0), cas$hL, km)
  expect_true(all(used))
})

test_that("PERF-3: fallback to the loop with non-finite values or weights out of range", {
  cas <- mk_pairs_case(seed = 5)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.2, -0.1)

  lt_inf <- cas$lt; lt_inf[3] <- -Inf             # e.g. t = 0
  used <- cdq_vs_loop(lt_inf, cas$lq, cas$z, cas$status, beta, cas$hL, km)
  expect_false(any(used))

  # pathological weights on members WITH an event (censored ones form no
  # pairs and the sweep does not even read their weight — it stays exact)
  ev0 <- which(cas$status[cas$z[, 1] == 0] == 1)[1]   # group 0, column 1
  ev1 <- which(cas$status[cas$z[, 2] == 1] == 1)[1]   # group 1, column 2
  km_tiny <- km; km_tiny[[1]]$w0[ev0] <- 1e-200   # the product could overflow
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_tiny)
  expect_identical(used, c(FALSE, TRUE))

  km_neg <- km; km_neg[[2]]$w1[ev1] <- -0.5       # negative weight
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_neg)
  expect_identical(used, c(TRUE, FALSE))

  km_zero <- km; km_zero[[1]]$w0[ev0] <- 0        # zero weight: the sweep handles it
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_zero)
  expect_true(all(used))
})

test_that("PERF-3: invalid algo and piecewise (always the loop)", {
  cas <- mk_pairs_case()
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                        c(0, 0), cas$hL, km, "identity", "x"),
               "algo")
  pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                 c(0, 0), cas$hL, km, "piecewise")
  expect_false(any(attr(pw, "meta")$cdq_cols))
})

test_that("fast_path = TRUE is the default of semi_aft_fg", {
  expect_true(formals(semi_aft_fg)$fast_path)
  cas <- mk_pairs_case(seed = 43, n = 80)
  args <- list(t = cas$t, qt = cas$qt, status = cas$status, z = cas$z,
               h1 = "log", h2 = "log", method = "Nelder-Mead",
               norm_min = "L2", L = exp(cas$hL), M = 0, phi = "identity",
               var_comp = FALSE, var_method = "jackknife")
  expect_identical(do.call(semi_aft_fg, args)$theta,
                   do.call(semi_aft_fg, c(args, fast_path = TRUE))$theta)
})

test_that("piecewise without pairs == explicit pairs (including ties in e2)", {
  for (cfg in list(c(40, 2, 1), c(300, 2, 2), c(500, 3, 3))) {
    set.seed(cfg[3])
    n <- cfg[1]; p <- cfg[2]
    z <- matrix(rbinom(n * p, 1, .45), n, p) * 1
    t <- rexp(n, exp(0.3 * z[, 1])); st <- as.numeric(rbinom(n, 1, .8))
    qt <- t * runif(n, .4, 1); qt[1:10] <- t[1:10] * 0.5
    lt <- log(t); lq <- log(qt); hL <- log(quantile(t, .9))
    km <- QAST:::prep_km_by_column(z, t, st)
    beta <- seq(0.2, -0.2, length.out = p)
    pr <- QAST:::make_fg_pairs_cpp(lt, lq, z, st, beta, hL, km, "piecewise")
    ln <- QAST:::make_fg_pw_lean_cpp(lt, lq, z, st, beta, hL, km)
    expect_true(attr(ln, "meta")$ok)
    for (th in list(rep(0, p), beta, rep(0.3, p))) {
      expect_equal(QAST:::eval_fg_parts_pw_lean_cpp(ln, lq, z, th),
                   QAST:::eval_fg_parts_cpp(pr, lq, z, th, "piecewise"),
                   tolerance = 1e-12)
      expect_equal(QAST:::eval_fg_obj_pw_lean_cpp(ln, lq, z, th, "L1"),
                   QAST:::eval_fg_obj_cpp(pr, lq, z, th, "piecewise", "L1"),
                   tolerance = 1e-12)
    }
  }
})
