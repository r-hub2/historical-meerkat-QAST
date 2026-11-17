# The radix sort of aft_core.cpp must produce EXACTLY the permutation
# of arma::stable_sort_index(e, "descend") and of order(-e): decreasing,
# ties in the original index order. U1/U2 depend only on that
# permutation.

test_that("radix sort == arma::stable_sort_index == order(-e)", {
  gens <- list(
    function(n) rnorm(n),
    function(n) round(rnorm(n), 1),                            # ties
    function(n) sample(c(-1, -0, 0, 1, 2), n, TRUE),           # -0 / +0
    function(n) sample(c(rnorm(max(n, 4)), Inf, -Inf, 5e-324, -5e-324), n),
    function(n) rep(3.5, n),                                   # all tied
    function(n) sample(c(.Machine$double.xmax, -.Machine$double.xmax,
                         1e-300, -1e-300, 0), n, TRUE),
    function(n) cumsum(rexp(n)),
    function(n) -cumsum(rexp(n)))
  for (g in gens) for (n in c(1, 2, 3, 17, 257, 2000)) {
    set.seed(n)
    e <- g(n)
    r <- QAST:::aft_order_desc_cpp(e, "radix")
    expect_identical(r, QAST:::aft_order_desc_cpp(e, "arma"))
    expect_identical(r, order(-e))
  }
})

test_that("radix sort: NaN keeps Armadillo's error", {
  expect_error(QAST:::aft_order_desc_cpp(c(1, NaN, 2)), "NaN")
})
