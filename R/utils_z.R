#' Checks that the covariates are binary (Fine and Gelber model)
#'
#' The Fine and Gelber (2001) model requires binary (0/1) covariates,
#' since its U-statistic runs over pairs between the groups defined by each
#' column of z.
#'
#' @param z Covariate matrix (n x p)
#' @return The matrix z, with column names ensured
#' @keywords internal
validate_discrete_z <- function(z) {
  z <- as.matrix(z)
  if (is.null(colnames(z))) {
    colnames(z) <- paste0("x", seq_len(ncol(z)))
  }
  bin_cols <- apply(z, 2, function(col) all(col %in% c(0, 1)))
  if (!all(bin_cols)) {
    invalid <- colnames(z)[!bin_cols]
    stop(
      "The Fine and Gelber model requires binary (0/1) covariates. ",
      "Non-binary columns: ", paste(invalid, collapse = ", "), ". ",
      "Encode multi-level factors as indicators (e.g. model.matrix(~ f)[, -1], without the intercept), ",
      "or use the formula interface."
    )
  }
  z
}

#' Rejects constant columns in z (including an intercept)
#'
#' The estimating equations for beta and theta are rank-based (they depend
#' on the ordering of the residuals \code{lt - z \%*\% beta}) and invariant to a
#' common shift: a constant column adds the same value to every residual
#' and changes no comparison. Its coefficient is not identifiable (the
#' objective is flat along it), so the models do not estimate an intercept.
#' In the Fine and Gelber model, a constant binary column also leaves one of
#' the groups empty and the corresponding U-statistic component null.
#'
#' @param z Covariate matrix (n x p)
#' @return The matrix z (invisibly)
#' @keywords internal
validate_no_constant_z <- function(z) {
  z <- as.matrix(z)
  if (is.null(colnames(z))) colnames(z) <- paste0("x", seq_len(ncol(z)))
  const_cols <- apply(z, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0L || all(v == v[1L])
  })
  if (any(const_cols)) {
    stop(
      "z has constant column(s): ",
      paste(colnames(z)[const_cols], collapse = ", "), ". ",
      "The model does not estimate an intercept (the estimating equations are ",
      "invariant to shifts), so the coefficient of a constant column is ",
      "not identifiable. Remove these columns from z.",
      call. = FALSE
    )
  }
  invisible(z)
}

#' Warns about linearly dependent columns in z
#'
#' If z does not have full column rank (e.g. duplicated columns, or a column
#' equal to a linear combination of others), the coefficients along those
#' directions are not identifiable: the estimating equations are flat along
#' the null space of z. The fit is not stopped; the starting values of those
#' columns come from \code{runif()} (lm returns NA) and the individual
#' estimates depend on the seed. Only identifiable combinations are
#' interpretable.
#'
#' @param z Covariate matrix (n x p), without constant columns
#' @param tol Rank tolerance (as in \code{qr()}, default 1e-7)
#' @return TRUE (invisibly) if z has full rank; FALSE if the warning was
#'   issued
#' @keywords internal
warn_collinear_z <- function(z, tol = 1e-7) {
  z <- as.matrix(z)
  if (is.null(colnames(z))) colnames(z) <- paste0("x", seq_len(ncol(z)))
  q <- qr(z, tol = tol)
  if (q$rank < ncol(z)) {
    dep <- colnames(z)[q$pivot[seq.int(q$rank + 1L, ncol(z))]]
    warning(
      "z is not of full rank (rank ", q$rank, " of ", ncol(z),
      " columns): linearly dependent column(s): ",
      paste(dep, collapse = ", "), ". The coefficients along these ",
      "directions are not identifiable and the individual estimates depend ",
      "on the (random) starting values. Consider removing these columns.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
