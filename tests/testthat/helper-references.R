# tests/testthat/helper-references.R
#
# Reference implementations used only as oracles by the tests (and by the
# scripts in benchmarks/, which source this file). Not part of the
# package: the package code uses the C++ kernels and the findInterval
# evaluation of the Kaplan-Meier curve, which the tests compare to these.

# Fine and Gelber U-statistic components, by a direct double loop in R
# (all cross pairs of each binary column, censoring weights km_surv,
# trimming at hL); vector of length p
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

# Censoring KM evaluated at each t by which(), O(n m): oracle of
# km_weights_at (findInterval)
km_weights_at_legacy <- function(t, status) {
  n <- length(t)
  if (n == 0L) return(numeric(0))
  
  # Group without events: weights 1
  if (sum(status) == 0L) return(rep(1, n))
  
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  if (length(fit$time) == 0L || all(is.na(fit$surv))) {
    return(rep(1, n))
  }
  
  times <- fit$time
  surv <- fit$surv
  
  vapply(t, function(ti) {
    idx <- which(times <= ti)
    if (length(idx) == 0L) {
      1
    } else {
      surv[max(idx)]
    }
  }, numeric(1))
}

# prep_km_by_column with the which() evaluation
prep_km_by_column_legacy <- function(z, t, status, group = NULL) {
  g <- QAST:::fg_cens_partition(z, group)
  H <- numeric(length(t))
  for (l in levels(g)) {
    i <- which(g == l)
    if (length(i)) H[i] <- km_weights_at_legacy(t[i], status[i])
  }
  lapply(seq_len(ncol(z)), function(j) {
    x <- z[, j]
    list(w0 = H[x == 0], w1 = H[x == 1])
  })
}
