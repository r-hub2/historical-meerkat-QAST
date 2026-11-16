#' Bivariate semiparametric AFT regression model (Fine and Gelber, 2001)
#'
#' Fits the bivariate semiparametric accelerated failure time model of Fine
#' and Gelber (2001) for the joint analysis of survival time and
#' quality-adjusted survival time. The model works on pairs of subjects
#' between the groups defined by each binary covariate, with Kaplan-Meier
#' censoring weights per group and truncation at \code{L}.
#'
#' @param t Observed survival times, or a formula
#'   \code{Qsurv(t, qt, status) ~ covariates} (see \code{\link{Qsurv}}); with
#'   a formula, \code{qt}, \code{status} and \code{z} come from it and from
#'   \code{data}, and factors enter as indicators (without an intercept)
#' @param qt Quality-adjusted survival times (\code{0 <= qt <= t})
#' @param status Event indicator (1 = event, 0 = censored)
#' @param z Binary covariate matrix (n x p), WITHOUT an intercept. Binary
#'   covariates are central to this model; multi-level factors must be
#'   coded as indicators (the formula interface does it). Constant columns
#'   are rejected, since the model does not estimate an intercept
#' @param h1 Link function for \code{t} (default \code{"log"})
#' @param h2 Link function for \code{qt} (default \code{"log"})
#' @param method Optimizer: "Nelder-Mead" (default), "uobyqa", "mads",
#'   "subplex", "BB", "SANN", "bobyqa", "newuoa", "hjkb", "BFGS". With a
#'   single covariate, methods that do not support one coefficient
#'   (Nelder-Mead, uobyqa, mads, bobyqa, hjkb) are replaced by "subplex",
#'   with a message
#' @param norm_min Minimization norm, \code{"L2"} (default) or \code{"L1"}
#' @param L Truncation time for the quality-adjusted parameters (required,
#'   positive)
#' @param phi Odd function: \code{"identity"} (default) or
#'   \code{"piecewise"}
#' @param M Number of bootstrap replicates; required only with
#'   \code{var_method = "bootstrap"}
#' @param var_comp If TRUE, estimates the variances (default FALSE)
#' @param var_method \code{"jackknife"} (default) or \code{"bootstrap"}
#' @param fast_path If TRUE (default) and \code{phi = "identity"},
#'   \code{norm_min = "L2"}, theta is estimated in closed form: least
#'   squares on the affine representation, which is the EXACT minimum of the
#'   L2 objective (the optimizer only approximates it). In that case
#'   \code{method} is used for beta only; theta does not go through the
#'   optimizer, also in the variance replicates. If the system is
#'   rank-deficient or the residual gradient is large: a warning and a
#'   fallback to the optimizer. With other phi/norm_min combinations, or
#'   with \code{fast_path = FALSE}, theta is always estimated by
#'   \code{method}
#' @param data Data frame with the formula variables, when \code{t} is a
#'   formula
#'
#' @return Object of class \code{"qast_fit"} (see \code{\link{qast_fit}}),
#'   with \code{beta}, \code{theta}, \code{intercept} and, if
#'   \code{var_comp = TRUE},
#'   \code{var_beta}, \code{var_theta}, \code{vcov_beta} and
#'   \code{vcov_theta}
#' @references Fine, J. P. and Gelber, R. D. (2001). Joint regression
#'   analysis of survival and quality-adjusted survival. \emph{Biometrics},
#'   57(2), 376-382.
#' @seealso \code{\link{semi_aft_prop}}; \code{\link{qast_fit}} for
#'   \code{print}, \code{summary}, \code{coef} and \code{vcov}
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(trt = rbinom(n, 1, 0.5), sex = rbinom(n, 1, 0.4))
#' d$t <- rexp(n, rate = exp(-0.4 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#'
#' fit <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                    L = quantile(d$t, 0.9))
#' fit
#'
#' \donttest{
#' # With jackknife standard errors
#' summary(semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
#'                     L = quantile(d$t, 0.9), var_comp = TRUE))
#' }
#' @export
#' @importFrom stats make.link lm pnorm printCoefmat runif optim cov coef vcov
#' @importFrom minqa uobyqa
#' @importFrom subplex subplex
#' @importFrom dfoptim mads nmk
#' @importFrom optimx optimr
#' @importFrom BB BBoptim
semi_aft_fg <- function(t, qt, status, z, h1 = "log", h2 = "log",
                        method = "Nelder-Mead", norm_min = "L2",
                        L, phi = "identity", M = NULL, var_comp = FALSE,
                        var_method = "jackknife", fast_path = TRUE,
                        data = NULL) {
  
  cl <- match.call()
  formula_info <- NULL
  
  # --- Formula --------------------------------------------------------------
  if (inherits(t, "formula")) {
    names(cl)[names(cl) == "t"] <- ""   # Call: mostra a formula sem "t ="
    md <- .qast_model_data(t, data, parent.frame())
    t <- md$t; qt <- md$qt; status <- md$status; z <- md$z
    formula_info <- md$info
  }
  if (missing(L)) stop("L (truncation time) is required", call. = FALSE)
  L <- unname(L)
  
  # --- Validacao -------------------------------------------------------------
  method <- match.arg(method,
                      c("uobyqa", "mads", "Nelder-Mead", "subplex", "BB", "SANN",
                        "bobyqa", "newuoa", "hjkb", "BFGS"))
  phi <- match.arg(phi, c("identity", "piecewise"))
  norm_min <- match.arg(norm_min, c("L1", "L2"))
  var_method <- match.arg(var_method, c("jackknife", "bootstrap"))
  .validate_var_args(var_comp, var_method, M)
  
  h1_name <- h1; h2_name <- h2
  h1 <- make.link(h1)
  h2 <- make.link(h2)
  .validate_qast_data(t, qt, status, z, h1, h2, L = L)
  method <- .qast_method_for_p(method, NCOL(z))
  
  # Covariaveis discretas: premissa central do modelo de Fine e Gelber
  z <- validate_discrete_z(z)
  validate_no_constant_z(z)   # sem intercepto / grupo vazio
  z_full_rank <- warn_collinear_z(z)   # aviso, nao erro
  # Em double uma vez: com inteiros o RcppArmadillo converteria z e status
  # a cada avaliacao dos kernels.
  storage.mode(z) <- "double"
  status <- as.numeric(status)
  
  # --- Pre-computacao (constantes na otimizacao) -----------------------------
  lt <- h1$linkfun(t)
  lq <- h2$linkfun(qt)
  hL <- h1$linkfun(L)
  n <- length(t)
  
  # --- Helpers ---------------------------------------------------------------
  run_optimizer <- function(par, fn) {
    switch(method,
           "uobyqa"       = uobyqa(par = par, fn = fn),
           "mads"         = mads(par = par, fn = fn,
                                 control = list(trace = FALSE)),
           "Nelder-Mead"  = nmk(par = par, fn = fn),
           "subplex"      = subplex(par = par, fn = fn),
           "BB"           = BBoptim(par = par, fn = fn,
                                    control = list(trace = FALSE),
                                    quiet = TRUE),
           "SANN"         = optim(par = par, fn = fn, method = "SANN"),
           "bobyqa" = ,
           "newuoa" = ,
           "hjkb"   = ,
           "BFGS"   = optimr(par = par, fn = fn, method = method,
                             control = list(trace = FALSE)),
           stop("unknown method: ", method)
    )
  }
  
  extract_conv <- function(res) {
    if (method == "uobyqa") return(res$ierr)
    cc <- res$convergence
    if (method == "mads") cc <- if (cc < 1e-06) 0 else -1
    cc
  }
  
  init_coef <- function(y, zmat) {
    b <- tryCatch(lm(y ~ zmat)$coef[-1],
                  error = function(e) rep(0, ncol(zmat)))
    if (anyNA(b)) b[is.na(b)] <- runif(sum(is.na(b)))
    b
  }
  
  # ---------------------------------------------------------------------------
  # Etapa 1: Estimacao de Beta (mesma equacao do modelo proposto)
  # ---------------------------------------------------------------------------
  U1min <- run_optimizer(init_coef(lt, z),
                         function(beta) U1_obj_cpp(beta, lt, z,
                                                   status, norm_min))
  stopifnot("The beta estimation algorithm did not converge" =
              extract_conv(U1min) == 0)
  betahat <- U1min$par
  
  # ---------------------------------------------------------------------------
  # Variancia de Beta (var_comp = TRUE)
  # ---------------------------------------------------------------------------
  vcov_beta <- NULL
  if (var_comp) {
    
    refit_beta <- function(t_sub, z_sub, status_sub) {
      lt_sub <- h1$linkfun(t_sub)
      run_optimizer(betahat,
                    function(b) U1_obj_cpp(b, lt_sub, z_sub,
                                           status_sub, norm_min))$par
    }
    
    # O jackknife de beta e feito junto com o de theta (abaixo): o refit
    # de theta usa o betahat da amostra completa, entao as estimativas por
    # replica sao as mesmas de passadas separadas.

    if (var_method == "bootstrap") {
      dat <- data.frame(t, status, z)
      boot_out <- boot::boot(
        data = dat, R = M,
        parallel = "multicore",
        ncpus = .qast_workers(),
        statistic = function(data, ind) {
          refit_beta(data[ind, 1],
                     as.matrix(data[ind, -(1:2)]),
                     data[ind, 2])
        }
      )
      vcov_beta <- cov(boot_out$t)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Etapa 2: Estimacao de Theta
  #
  #   stage2_env(): monta o objetivo de uma (sub)amostra
  #     - identity: forma afim alpha + W theta, O(p^2) por avaliacao
  #     - piecewise: soma sobre os pares validos, O(nnz) por avaliacao
  #   solve_theta(): forma fechada (identity + L2 + fast_path) ou otimizador
  # ---------------------------------------------------------------------------
  km_surv <- prep_km_by_column(z, t, status)
  
  stage2_env <- function(lt_s, lq_s, z_s, status_s, km_s) {
    if (phi == "identity") {
      pairs <- make_fg_pairs(lt_s, lq_s, z_s, status_s,
                             betahat, hL, km_s, phi)
      aff <- make_fg_affine(pairs, lq_s, z_s)
      list(
        fn = function(th) eval_fg_obj_affine(aff, th, norm_min),
        affine = aff,
        can_fast = TRUE
      )
    } else {
      pairs <- make_fg_pairs(lt_s, lq_s, z_s, status_s,
                             betahat, hL, km_s, phi)
      list(
        fn = function(th) eval_fg_obj(pairs, lq_s, z_s, th, "piecewise",
                                      norm_min),
        affine = NULL,
        can_fast = FALSE
      )
    }
  }
  
  # Forma fechada quando possivel; posto deficiente ou gradiente residual
  # alto -> warning e otimizador sobre o mesmo objetivo afim.
  solve_theta <- function(s2, par0) {
    if (fast_path && s2$can_fast && norm_min == "L2") {
      rcond <- 1e-10
      sol <- solve_fg_theta_affine(s2$affine, rcond = rcond)
      if (isTRUE(sol$full_rank) && is.finite(sol$grad_norm) &&
          sol$grad_norm < 1e-6) {
        return(list(par = sol$par, converged = TRUE))
      }
      # Se z inteiro ja foi sinalizado como colinear, este aviso repete a
      # mesma causa; so e emitido quando a perda de posto e nova (ex.:
      # uma replica de jackknife/bootstrap degenerada).
      if (z_full_rank) warning(
        "fast_path: rank-deficient system or large residual gradient ",
        sprintf("(rcond = %.2g, full_rank = %s, grad_norm = %.2g). ",
                rcond, isTRUE(sol$full_rank), sol$grad_norm),
        "falling back to the optimizer on the same affine objective; ",
        "theta is not identifiable along the null directions."
      )
    }
    U2min <- run_optimizer(par0, s2$fn)
    stopifnot("The theta estimation algorithm did not converge" =
                extract_conv(U2min) == 0)
    list(par = U2min$par, converged = TRUE)
  }
  
  stage2 <- stage2_env(lt, lq, z, status, km_surv)
  thetahat <- solve_theta(stage2, init_coef(lq, z))$par
  
  # ---------------------------------------------------------------------------
  # Variancia de Theta (var_comp = TRUE)
  #
  # Por replica: pesos KM, pares e objetivo reconstruidos na subamostra e
  # resolvidos pelo mesmo caminho da estimativa pontual.
  # ---------------------------------------------------------------------------
  vcov_theta <- NULL
  if (var_comp) {
    
    refit_theta <- function(t_sub, qt_sub, status_sub, z_sub) {
      lq_sub <- h2$linkfun(qt_sub)
      lt_sub <- h1$linkfun(t_sub)
      km_surv_sub <- prep_km_by_column(z_sub, t_sub, status_sub)
      s2 <- stage2_env(lt_sub, lq_sub, z_sub, status_sub, km_surv_sub)
      solve_theta(s2, thetahat)$par
    }
    
    if (var_method == "jackknife") {
      # beta e theta numa so passada (ver R/utils_jackknife.R)
      p <- ncol(z)
      jk <- .jackknife_se(
        n,
        refit = function(i) {
          z_i <- z[-i, , drop = FALSE]
          c(refit_beta(t[-i], z_i, status[-i]),
            refit_theta(t[-i], qt[-i], status[-i], z_i))
        },
        blocks = list(beta = seq_len(p), theta = p + seq_len(p)),
        cov = TRUE
      )
      vcov_beta  <- jk$beta$cov
      vcov_theta <- jk$theta$cov
    }

    if (var_method == "bootstrap") {
      dat <- data.frame(t, qt, status, z)
      boot_out <- boot::boot(
        data = dat, R = M,
        parallel = "multicore",
        ncpus = .qast_workers(),
        statistic = function(data, ind) {
          refit_theta(data[ind, 1], data[ind, 2], data[ind, 3],
                      as.matrix(data[ind, -(1:3)]))
        }
      )
      vcov_theta <- cov(boot_out$t)
    }
  }
  
  .new_qast_fit(
    model = "fg", call = cl, z = z, status = status,
    betahat = betahat, thetahat = thetahat,
    vcov_beta = vcov_beta, vcov_theta = vcov_theta,
    var_comp = var_comp, var_method = var_method, M = M,
    lt = lt, lq = lq, t = t, qt = qt, formula_info = formula_info,
    settings = list(h1 = h1_name, h2 = h2_name, method = method,
                    norm_min = norm_min, L = L, phi = phi,
                    fast_path = fast_path)
  )
}
