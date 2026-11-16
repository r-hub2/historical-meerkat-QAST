#' Exact affine form of the U-statistic components (Fine and Gelber)
#'
#' Extracts, from the pair object (\code{make_fg_pairs}), the decomposition
#'
#'   parts(theta) = alpha + W * theta
#'
#' valid for phi = "identity". The extraction runs once per fit; each later
#' evaluation costs O(p^2), independent of n.
#'
#' Matches the reference kernel to 1e-10 (sums reordered).
#'
#' @param pairs Pair object (make_fg_pairs)
#' @param lq \code{h2$linkfun(qt)}, precomputed
#' @param z Binary covariate matrix (n x p)
#' @return List with alpha (vector p), W (matrix p x p) and the propagated
#'   meta
#' @keywords internal
make_fg_affine <- function(pairs, lq, z) {
  make_fg_affine_cpp(pairs, lq, z)
}

#' Norm (L1/L2) of the affine components, scaled by sqrt(p)
#' @keywords internal
eval_fg_obj_affine <- function(aff, theta, norm_min) {
  eval_fg_obj_affine_cpp(aff, theta, norm_min)
}

#' Closed-form least-squares solution (L2 fast path)
#'
#' Minimizes ||alpha + W theta||^2 through the normal equations. Rank
#' deficiency (degenerate columns, partial non-identifiability) is detected
#' by an extreme singular-value ratio < rcond or a large residual; it then
#' falls back to the minimum-norm solution with full_rank = FALSE
#' (semi_aft_fg issues the warning).
#'
#' @param aff Affine object (make_fg_affine)
#' @param rcond Singular-value ratio threshold for rank detection (1e-10)
#' @return List: par (vector without dim), full_rank (logical), grad_norm
#'   (check; ~0 at full rank)
#' @keywords internal
solve_fg_theta_affine <- function(aff, rcond = 1e-10) {
  solve_fg_theta_affine_cpp(aff, rcond)
}

