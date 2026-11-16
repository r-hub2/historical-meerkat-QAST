#' Kaplan-Meier censoring survival at each subject's observed time
#'
#' The censoring KM curve is estimated by \code{survival::survfit} and the
#' step function is evaluated at each \code{t[i]} with \code{findInterval}
#' (O(n log m) per group).
#'
#' Conventions (checked by the tests):
#'   - empty group -> numeric(0)
#'   - group without events -> rep(1, n)
#'   - fit without valid times/survival -> rep(1, n)
#'   - t before the first jump -> 1
#'   - t beyond the last jump -> last value of the curve
#'
#' @param t Observed times of the group
#' @param status Event indicator (1 = event)
#' @param use_legacy If TRUE, evaluates with \code{which()} (reference
#'   implementation for the A/B tests)
#' @return Vector of length \code{length(t)} with \code{S_cens(t_i)}
#' @keywords internal
km_weights_at <- function(t, status, use_legacy = FALSE) {
  if (use_legacy) return(.km_weights_at_legacy(t, status))
  .km_weights_at_fast(t, status)
}

#' Censoring KM evaluation: survfit + findInterval
#' @keywords internal
.km_weights_at_fast <- function(t, status) {
  n <- length(t)
  if (n == 0L) return(numeric(0))
  
  # Grupo sem eventos: pesos 1
  if (sum(status) == 0L) return(rep(1, n))
  
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  if (length(fit$time) == 0L || all(is.na(fit$surv))) {
    return(rep(1, n))
  }
  
  # findInterval(x, v) = numero de elementos de v que sao <= x
  # (left.open = FALSE, rightmost.closed = FALSE) == max(which(v <= x))
  idx <- findInterval(t, fit$time)
  
  # idx == 0 (antes do primeiro salto): S = 1. Nos demais: surv[idx].
  w <- fit$surv[pmax(idx, 1L)]
  w[idx == 0L] <- 1
  w
}

#' Reference censoring KM evaluation with which(), O(n m)
#'
#' Used only by the A/B tests of km_weights_at; keep unchanged.
#' @keywords internal
.km_weights_at_legacy <- function(t, status) {
  n <- length(t)
  if (n == 0L) return(numeric(0))
  
  # Grupo sem eventos: pesos 1
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

#' KM weights for each binary covariate column
#'
#' For each column j of z, the censoring KM weights in the groups
#' \code{z[, j] == 0} (w0) and \code{z[, j] == 1} (w1). Identical columns
#' define the same groups and therefore the same weights: with
#' \code{dedupe_cols = TRUE}, survfit runs once per distinct column.
#' Recomputed in every variance replicate.
#'
#' @param z Binary covariate matrix (n x p)
#' @param t Observed times
#' @param status Event indicator
#' @param use_legacy Passed to km_weights_at (reference for the tests)
#' @param dedupe_cols If TRUE (default), computes the weights once per
#'   distinct column (safe by construction: identical columns give
#'   identical groups and weights)
#' @return List of length p; element j holds w0 and w1
#' @keywords internal
prep_km_by_column <- function(z, t, status, use_legacy = FALSE,
                              dedupe_cols = TRUE) {
  p <- ncol(z)
  
  if (!dedupe_cols || p == 1L) {
    out <- lapply(seq_len(p), function(j) {
      x <- z[, j]
      list(
        w0 = km_weights_at(t[x == 0], status[x == 0],
                           use_legacy = use_legacy),
        w1 = km_weights_at(t[x == 1], status[x == 1],
                           use_legacy = use_legacy)
      )
    })
    return(out)
  }
  
  # key[j] = indice da PRIMEIRA coluna identica a coluna j (j se unica).
  # Impressao digital numerica (soma ponderada com pesos fixos, sem RNG)
  # para achar candidatas em O(n*p); a igualdade e CONFIRMADA com
  # identical(), entao colisoes de impressao digital nunca fundem colunas
  # diferentes.
  fp <- as.vector(crossprod(z, sqrt(seq_len(nrow(z)) + 1)))
  key <- seq_len(p)
  for (j in seq_len(p)[-1L]) {
    cand <- which(fp[seq_len(j - 1L)] == fp[j])
    for (k in cand) {
      if (identical(z[, k], z[, j])) { key[j] <- key[k]; break }
    }
  }
  
  uniq <- sort(unique(key))
  
  # Ajusta cada grupo-distinto UMA vez
  cache <- lapply(uniq, function(k) {
    x <- z[, k]
    list(
      w0 = km_weights_at(t[x == 0], status[x == 0],
                         use_legacy = use_legacy),
      w1 = km_weights_at(t[x == 1], status[x == 1],
                         use_legacy = use_legacy)
    )
  })
  names(cache) <- as.character(uniq)
  
  # Replica a saida para todas as colunas com a mesma chave
  lapply(key, function(k) cache[[as.character(k)]])
}