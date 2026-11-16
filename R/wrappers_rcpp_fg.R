#' Theta estimating objective of the Fine and Gelber (2001) model
#'
#' U-statistic over pairs between the groups defined by each binary
#' covariate, with censoring KM weights per group (reverse curve
#' \code{Surv(t, 1 - status)}, estimated by survfit), truncation at L and
#' odd function phi ("identity" or "piecewise"), evaluated by the
#' O(nk * nl) loop of src/fg_core.cpp. It is the reference implementation
#' for the A/B tests and benchmarks; semi_aft_fg uses the pair object
#' (\code{make_fg_pairs}).
#'
#' Argument order: the wrapper takes \code{(par, lq, lt, ...)} and passes
#' them to the C++ kernel as \code{(par, lt, lq, ...)}.
#'
#' @param par Numeric vector of theta coefficients
#' @param lq Transformed qt (precomputed)
#' @param lt Transformed t (precomputed)
#' @param z Binary covariate matrix (n x p)
#' @param status Event indicator (1 = event, 0 = censored)
#' @param beta Survival coefficients (fixed at this stage)
#' @param hL Constant h1(L), precomputed
#' @param km_surv List of length p with w0 and w1 per column
#' @param phi Odd function: "identity" or "piecewise"
#' @param norm_min Minimization norm ("L1" or "L2")
#' @return Scalar value of the objective function
#' @seealso \code{\link{semi_aft_fg}}, \code{\link{prep_km_by_column}}
#' @keywords internal
U2_fun_fg_cpp <- function(par, lq, lt, z, status, beta, hL, km_surv,
                          phi, norm_min) {
  U2_obj_fg_cpp(par, lt, lq, z, status, beta, hL, km_surv, phi, norm_min)
}
