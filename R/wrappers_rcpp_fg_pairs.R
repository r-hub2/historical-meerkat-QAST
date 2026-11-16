#' Precomputed pair object of the U-statistic (Fine and Gelber)
#'
#' Computes, once per fit, the constants of the theta stage: the truncation
#' indicator h1L, the KM weights, the event filters and the denominators.
#' Evaluating the objective then costs O(nk + nl) (identity, via the Ra/Cb
#' marginals) or O(nnz) (piecewise). With phi = "identity", Ra/Cb are
#' computed in O(N log^2 N); see src/fg_pairs_core.cpp.
#'
#' @param lt \code{h1$linkfun(t)}, precomputed
#' @param lq \code{h2$linkfun(qt)}, precomputed
#' @param z BINARY covariate matrix (n x p); check it with
#'   validate_discrete_z before calling
#' @param status Event indicator (1 = event)
#' @param beta Survival coefficients (fixed at this stage)
#' @param hL Constant h1(L), precomputed
#' @param km_surv List of length p with w0/w1 per column (output of
#'   prep_km_by_column, in the original column order of z)
#' @param phi "identity" or "piecewise"
#' @return List per column with ka/kb/Ra/Cb/pa/pb/pc plus attr "meta"
#'   (p, n2, has_pairs, cdq_cols). See the schema in the header of
#'   fg_pairs_core.cpp.
#' @keywords internal
make_fg_pairs <- function(lt, lq, z, status, beta, hL, km_surv, phi) {
  make_fg_pairs_cpp(lt, lq, z, status, beta, hL, km_surv, phi)
}

#' Norm (L1/L2) of the components, scaled by sqrt(p)
#' @keywords internal
eval_fg_obj <- function(pairs, lq, z, theta, phi, norm_min) {
  eval_fg_obj_cpp(pairs, lq, z, theta, phi, norm_min)
}