#' Fine and Gelber U-statistic components (reference in R)
#'
#' Reference implementation for the R vs C++ equivalence tests. Not for
#' production use.
#'
#' @param theta Quality-adjusted coefficients
#' @param lt Failure times transformed by h1
#' @param lq Quality-adjusted times transformed by h2
#' @param z Binary covariate matrix (n x p)
#' @param status Event indicator (1 = event)
#' @param beta Survival coefficients (fixed)
#' @param hL Constant h1(L), precomputed
#' @param km_surv List of length p with w0 and w1 per column
#' @param phi Odd function: "identity" or "piecewise"
#' @return Vector of length p with the component of each covariate
#' @keywords internal
U2_parts_fg_ref <- function(theta, lt, lq, z, status, beta, hL, km_surv,
                            phi) {
  n <- nrow(z)
  e1 <- as.vector(lt - z %*% beta)
  e2 <- as.vector(lq - z %*% theta)
  
  phi_fun <- switch(phi,
                    identity  = function(x) x,
                    piecewise = function(x) as.integer((x >= 0) - (x <= 0))
  )
  
  sapply(seq_len(ncol(z)), function(j) {
    X <- z[, j]
    
    e1k <- e1[X == 0]
    e1l <- e1[X == 1]
    e2k <- e2[X == 0]
    e2l <- e2[X == 1]
    dk  <- status[X == 0]
    dl  <- status[X == 1]
    
    zk <- z[X == 0, , drop = FALSE]
    zl <- z[X == 1, , drop = FALSE]
    
    wk <- km_surv[[j]]$w0
    wl <- km_surv[[j]]$w1
    
    acc <- 0
    
    for (a in seq_along(e1k)) {
      for (b in seq_along(e1l)) {
        
        # All cross pairs (a in group 0, b in group 1)
        
        # Events in both
        if (dk[a] == 0 || dl[b] == 0) next
        
        denom <- wk[a] * wl[b]
        if (!(denom > 0) || !is.finite(denom)) next
        
        h1L <- as.numeric(
          hL - max(zk[a, ] %*% beta, zl[b, ] %*% beta) > max(e1k[a], e1l[b])
        )
        
        acc <- acc - h1L * phi_fun(e2k[a] - e2l[b]) / denom
      }
    }
    acc / n^2
  })
}