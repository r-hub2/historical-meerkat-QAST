# R/utils_jackknife.R
#
# Jackknife delete-1, com a mesma formula de erro-padrao de
# jackknifeR::jackknife(d = 1, weight = FALSE) (verificado nos testes):
#
#   jk        = matriz n x k das estimativas sem a observacao i
#   validas   = complete.cases(jk)          (linhas com NA descartadas)
#   theta_dot = colMeans(jk[validas, ])
#   se        = sqrt( (n - 1) / n_validas * colSums((jk - theta_dot)^2) )
#
# Execucao: fork (parallel::mclapply); varios blocos de parametros (ex.:
# beta e theta) numa so passada, com remocao de NA por bloco; warnings dos
# workers repassados a sessao principal, uma vez por mensagem distinta.
#
# Workers: getOption("QAST.cores"), senao nucleos FISICOS
# (detectCores(logical = FALSE)), no maximo 2 sob R CMD check
# (_R_CHECK_LIMIT_CORES_). No Windows (sem fork): sequencial.

#' Number of workers for the variance refits
#' @keywords internal
.qast_workers <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  opt <- getOption("QAST.cores")
  nc <- if (!is.null(opt)) {
    as.integer(opt)
  } else {
    suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  if (is.na(nc) || nc < 1L) nc <- 1L
  # R CMD check --as-cran (e o CRAN) limitam a 2 processos simultaneos;
  # o parallel aborta acima disso quando _R_CHECK_LIMIT_CORES_ esta ativo
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(toupper(lim), "FALSE")) nc <- min(nc, 2L)
  as.integer(nc)
}

#' Delete-1 jackknife standard errors, by parameter block
#'
#' @param n Number of observations
#' @param refit Function of i (the removed index) returning the vector of
#'   estimates on the sample without observation i
#' @param blocks Named list of column indices of the vector returned by
#'   refit (e.g. \code{list(beta = 1:p, theta = p + 1:p)})
#' @param cores Number of workers (fork)
#' @param cov If TRUE, returns \code{list(se, cov)} per block, with the
#'   jackknife covariance matrix
#'   \code{(n - 1) / n_valid * crossprod(jk - theta_dot)} (whose diagonal is
#'   \code{se^2})
#' @return Named list (as blocks) of standard-error vectors, or of
#'   \code{list(se, cov)} when \code{cov = TRUE}
#' @keywords internal
.jackknife_se <- function(n, refit, blocks, cores = .qast_workers(),
                          cov = FALSE) {

  one <- function(i) {
    warns <- character()
    val <- withCallingHandlers(
      refit(i),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = as.numeric(val), warnings = warns)
  }

  idx <- seq_len(n)
  res <- if (cores > 1L) {
    # warnings proprios do mclapply ("scheduled core ... error") sao
    # redundantes: o erro original e re-emitido abaixo. Os warnings da
    # estatistica ja foram capturados dentro de one().
    suppressWarnings(
      parallel::mclapply(idx, one, mc.cores = cores, mc.preschedule = TRUE)
    )
  } else {
    lapply(idx, one)
  }

  # Erro em qualquer replica: re-emitido na sessao principal
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"),
                logical(1))
  if (any(bad)) {
    r <- res[[which(bad)[1]]]
    if (is.null(r)) stop("jackknife: a worker ended without a result")
    cond <- attr(r, "condition")
    if (inherits(cond, "condition")) stop(cond)
    stop(as.character(r))
  }

  for (msg in unique(unlist(lapply(res, `[[`, "warnings"))))
    warning(msg, call. = FALSE)

  jk <- do.call(rbind, lapply(res, `[[`, "value"))

  lapply(blocks, function(cols) {
    b <- jk[, cols, drop = FALSE]
    b <- b[stats::complete.cases(b), , drop = FALSE]
    nv <- nrow(b)
    theta_dot <- colMeans(b)
    dev <- sweep(b, 2, theta_dot, "-")
    se <- sqrt(((n - 1) / nv) * colSums(dev^2))
    if (!cov) return(se)
    list(se = se, cov = ((n - 1) / nv) * crossprod(dev))
  })
}
