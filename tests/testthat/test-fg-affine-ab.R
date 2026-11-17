# Fixes relative to the previous version:
#   - mk_fg_case now HONORS the argument p (it was ignored; it always
#     built 2 columns — a harness bug that broke the fast path);
#   - generic column names "z1..zp" (p != 2 has no "trt"/"sex");
#   - the fast path uses a beta consistent with p = 1;
#   - simplified exact-tie gate (without the obscure line with
#     tie_safe <- TRUE of the previous version).
# ----------------------------------------------------------------------------

library(testthat)

mk_fg_case <- function(seed = 42, n = 60, p = 2, tie_at_L = FALSE) {
  stopifnot(is.integer(as.integer(p)), p >= 1)
  set.seed(seed)
  # FIX B: p honored — p independent binary columns
  z <- matrix(rbinom(n * p, 1, 0.45), nrow = n, ncol = p)
  colnames(z) <- paste0("z", seq_len(p))
  
  # "True" survival coefficient only in the first column (effect on
  # z1), to keep the case identifiable even with p > 2
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  L      <- stats::quantile(t, 0.9)
  
  if (tie_at_L) {
    # Exact tie: beta = 0 => z'beta = 0 exactly;
    # times duplicated exactly at L; log() of the same double.
    t[1:min(4, n)]        <- L
    status[1:min(4, n)]   <- 1L
  }
  list(t = t, qt = qt, status = status, z = z,
       lt = log(t), lq = log(qt), hL = log(L))
}

# Grid of beta consistent with an arbitrary p: the first element varies,
# the others are zero (for p = 1 it becomes c(0), c(0.3), c(-0.5)).
mk_fg_beta_grid <- function(p) {
  base <- rep(0, p)
  rbind(base,
        replace(base, 1, 0.3),
        replace(base, 1, -0.5))
}

fg_theta_grid <- function(p) {
  base <- rep(0, p)
  g <- rbind(base)
  for (k in seq_len(p)) {
    g <- rbind(g, replace(base, k, 0.25),
               replace(base, k, -0.4))
  }
  g
}


test_that("Affine: equals the legacy kernel and the R reference on a grid of beta x theta", {
  skip_if_not(exists("U2_parts_fg_ref", envir = asNamespace("QAST")),
              "referencia R indisponivel")
  
  p <- 2
  beta_grid  <- mk_fg_beta_grid(p)
  theta_grid <- fg_theta_grid(p)
  
  for (bi in seq_len(nrow(beta_grid))) {
    for (ti in seq_len(nrow(theta_grid))) {
      
      cas <- mk_fg_case(seed = 40 + bi, n = 60, p = p)
      beta  <- beta_grid[bi, ]
      theta <- theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity")
      aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
      
      parts_aff <- QAST:::eval_fg_parts_affine_cpp(aff, theta)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")
      expect_equal(parts_aff, parts_leg, tolerance = 1e-10)
      
      # Triangulation with the pure R reference
      parts_ref <- QAST:::U2_parts_fg_ref(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")
      expect_equal(parts_leg, parts_ref, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_aff <- QAST:::eval_fg_obj_affine_cpp(aff, theta, nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity", nm)
        expect_equal(obj_aff, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Compiled piecewise: equals the legacy kernel on a grid of beta x theta", {
  p <- 2
  beta_grid  <- mk_fg_beta_grid(p)
  theta_grid <- fg_theta_grid(p)
  
  for (bi in seq_len(nrow(beta_grid))) {
    for (ti in seq_len(nrow(theta_grid))) {
      
      cas <- mk_fg_case(seed = 80 + bi, n = 60, p = p)
      beta  <- beta_grid[bi, ]
      theta <- theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise")
      pw <- QAST:::compile_fg_pw_cpp(pairs)
      
      parts_pw  <- QAST:::eval_fg_parts_pw_cpp(pw, cas$lq, cas$z, theta)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "piecewise")
      expect_equal(parts_pw, parts_leg, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_pw  <- QAST:::eval_fg_obj_pw_cpp(pw, cas$lq, cas$z, theta, nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise", nm)
        expect_equal(obj_pw, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Constant column: identically null component in every path", {
  # Degenerate column: Ra/Cb are zeros aligned with ka/kb;
  # make_fg_affine_cpp must accept the column and the component is
  # identically 0.
  cas <- mk_fg_case()
  cas$z[, 2] <- 0L   # constant column 2 => group 1 empty
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  # Schema invariants in the degenerate column
  expect_length(pairs[[2]]$kb, 0L)
  expect_length(pairs[[2]]$Ra, length(pairs[[2]]$ka))
  expect_true(all(pairs[[2]]$Ra == 0))
  
  aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)   # must not fail
  # The identity object stores no pairs; the piecewise one is built apart
  pairs_pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                       beta, cas$hL, km, "piecewise")
  pw  <- QAST:::compile_fg_pw_cpp(pairs_pw)
  
  th <- c(0.2, -0.1)
  expect_identical(QAST:::eval_fg_parts_affine_cpp(aff, th)[2], 0.0)
  expect_identical(QAST:::eval_fg_parts_pw_cpp(pw, cas$lq, cas$z, th)[2], 0.0)
  expect_identical(QAST:::U2_parts_fg_cpp(th, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")[2], 0.0)
})

test_that("Gate: exact tie in h1L preserved (affine)", {
  # Exact tie: beta = 0 => z'beta = 0 exactly; times duplicated exactly
  # at L; the strict inequality of h1L replicates bit for bit because
  # the freezing happens AT THE CONSTRUCTION of the pairs. It fails by a
  # DISCRETE JUMP if any floating-point reordering moves a pair across
  # the gate.
  cas <- mk_fg_case(tie_at_L = TRUE)
  beta <- rep(0, ncol(cas$z))
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
  
  theta <- c(0.1, 0.05)
  parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                      cas$z, cas$status, beta, cas$hL,
                                      km, "identity")
  parts_aff <- QAST:::eval_fg_parts_affine_cpp(aff, theta)
  expect_equal(parts_aff, parts_leg, tolerance = 1e-10)
})

test_that("Fast path: least squares reproduce the optimizer (L2/identity)", {
  skip_if_not_installed("minqa")
  
  # FIX C: uobyqa with a starting par EQUAL TO ZERO breaks in minqa's
  # commonArgs: the default rhobeg = min(0.95, 0.2*max(abs(par))) collapses to 0
  # => rhoend = 1e-6*rhobeg = 0 => stopifnot(0 < rhoend) fails.
  # Two defenses:
  #   (a) p = 2 — the uobyqa documentation asks for length(par) > 1;
  #   (b) NON-zero start + explicit rhobeg: immune to a zero start.
  cas <- mk_fg_case(seed = 777, n = 200, p = 2)
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  stopifnot(ncol(cas$z) == 2L, length(beta) == 2L)
  
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  aff   <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
  
  sol <- QAST:::solve_fg_theta_affine_cpp(aff)
  
  obj_fun <- function(th) QAST:::eval_fg_obj_affine_cpp(aff, th, "L2")
  opt <- minqa::uobyqa(c(0.2, -0.1), obj_fun,
                       control = list(rhobeg = 0.5, rhoend = 1e-8))
  
  # Non-inferiority: the closed-form solution cannot be WORSE than the optimizer
  expect_lte(obj_fun(sol$par), obj_fun(opt$par) + 1e-8)
  
  # Gradient ~ 0 at full rank (semi_aft_fg requires < 1e-6)
  expect_lt(sol$grad_norm, 1e-6)
  expect_true(sol$full_rank)
})
