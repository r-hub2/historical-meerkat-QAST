# R/qast_fit.R
#
# Interface comum dos modelos: resposta Qsurv() para formulas, validacao
# das entradas e a classe S3 "qast_fit" (print, summary, coef, vcov).

# ---------------------------------------------------------------------------
# Resposta para formulas
# ---------------------------------------------------------------------------

#' Response for QAST model formulas
#'
#' Bundles the survival time, the quality-adjusted time and the event
#' indicator into a response for the left-hand side of a formula, like
#' \code{survival::Surv()}:
#' \code{semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d, L = 5)}.
#'
#' @param time Observed survival time
#' @param qtime Quality-adjusted survival time
#' @param status Event indicator (1 = event, 0 = censored)
#' @return Matrix with columns \code{time}, \code{qtime} and \code{status}
#'   and class \code{"Qsurv"}
#' @seealso \code{\link{semi_aft_fg}}, \code{\link{semi_aft_prop}}
#' @examples
#' Qsurv(c(2, 5, 3), c(1.5, 4, 3), c(1, 0, 1))
#' @export
Qsurv <- function(time, qtime, status) {
  if (length(time) != length(qtime) || length(time) != length(status))
    stop("Qsurv: time, qtime and status must have the same length",
         call. = FALSE)
  out <- cbind(time = as.numeric(time), qtime = as.numeric(qtime),
               status = as.numeric(status))
  class(out) <- c("Qsurv", class(out))
  out
}

# Converte formula + data nas entradas dos modelos. A matriz de covariaveis
# e o model.matrix sem a coluna do intercepto (o modelo nao o estima);
# fatores entram como k - 1 indicadoras. NAs seguem adiante e sao
# reportados pela validacao.
.qast_model_data <- function(formula, data, env) {
  if (is.null(data)) data <- env
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!inherits(y, "Qsurv"))
    stop("the left-hand side of the formula must be Qsurv(time, qtime, status)",
         call. = FALSE)
  mm <- stats::model.matrix(attr(mf, "terms"), mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) == 0L)
    stop("the formula has no covariates", call. = FALSE)
  info <- list(terms = stats::delete.response(attr(mf, "terms")),
               xlevels = stats::.getXlevels(attr(mf, "terms"), mf),
               contrasts = attr(mm, "contrasts"))
  attr(mm, "assign") <- NULL
  attr(mm, "contrasts") <- NULL
  list(t = unname(y[, "time"]), qt = unname(y[, "qtime"]),
       status = unname(y[, "status"]), z = mm, info = info)
}

# ---------------------------------------------------------------------------
# Validacao das entradas
# ---------------------------------------------------------------------------

# Erro com as primeiras posicoes problematicas
.qast_stop_at <- function(what, idx) {
  shown <- utils::head(idx, 5L)
  more <- if (length(idx) > 5L) sprintf(" (and %d more)", length(idx) - 5L) else ""
  stop(what, " at observations ", paste(shown, collapse = ", "), more, ".",
       call. = FALSE)
}

# Validacao comum a semi_aft_fg e semi_aft_prop. h1/h2 sao os objetos de
# make.link; L e opcional (so o modelo de Fine e Gelber o usa).
.validate_qast_data <- function(t, qt, status, z, h1, h2, L = NULL) {
  if (!is.numeric(t) || !is.numeric(qt))
    stop("t and qt must be numeric", call. = FALSE)
  z <- as.matrix(z)
  n <- length(t)
  lens <- c(t = n, qt = length(qt), status = length(status), z = nrow(z))
  if (any(lens != n))
    stop("t, qt, status and z must have the same number of observations (",
         paste(names(lens), lens, sep = " = ", collapse = ", "), ")",
         call. = FALSE)
  if (n < 2L) stop("at least 2 observations are required", call. = FALSE)

  for (nm in c("t", "qt", "status")) {
    v <- get(nm)
    if (anyNA(v)) .qast_stop_at(sprintf("%s contains NA", nm), which(is.na(v)))
  }
  if (anyNA(z)) .qast_stop_at("z contains NA", which(!stats::complete.cases(z)))
  if (!is.numeric(z)) stop("z must be numeric", call. = FALSE)

  if (!all(status %in% c(0, 1)))
    .qast_stop_at("status must be 0 (censored) or 1 (event)",
                  which(!status %in% c(0, 1)))
  if (sum(status) == 0) stop("there are no events (status is all 0)", call. = FALSE)

  if (any(!is.finite(t) | t <= 0))
    .qast_stop_at("t must be finite and positive", which(!is.finite(t) | t <= 0))
  if (any(!is.finite(qt) | qt < 0))
    .qast_stop_at("qt must be finite and non-negative",
                  which(!is.finite(qt) | qt < 0))
  # Utilidades em [0, 1]: o tempo ajustado nao excede o de sobrevivencia
  # (folga relativa para arredondamento)
  over <- which(qt > t * (1 + 1e-8))
  if (length(over))
    .qast_stop_at("qt (quality-adjusted time) greater than t; check that the arguments are not swapped",
                  over)

  lt <- h1$linkfun(t)
  lq <- h2$linkfun(qt)
  if (any(!is.finite(lt)))
    .qast_stop_at(sprintf("h1(t) is not finite (link \"%s\")", h1$name),
                  which(!is.finite(lt)))
  if (any(!is.finite(lq)))
    .qast_stop_at(sprintf("h2(qt) is not finite (link \"%s\"); with the log link, qt must be > 0",
                          h2$name), which(!is.finite(lq)))

  if (!is.null(L)) {
    if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L <= 0)
      stop("L must be a positive number", call. = FALSE)
    if (!is.finite(h1$linkfun(L)))
      stop(sprintf("h1(L) is not finite (link \"%s\")", h1$name), call. = FALSE)
  }
  invisible(TRUE)
}

# Metodos que nao aceitam um unico coeficiente (dfoptim::nmk e mads pedem
# optimize(); minqa::uobyqa/bobyqa e hjkb nao convergem com p = 1): com
# uma covariavel, usa subplex, tambem sem derivadas.
.qast_method_for_p <- function(method, p) {
  if (p == 1L && method %in% c("Nelder-Mead", "uobyqa", "mads", "bobyqa",
                               "hjkb")) {
    message(sprintf(paste0("method = \"%s\" does not support a single coefficient; ",
                           "using \"subplex\"."), method))
    return("subplex")
  }
  method
}

# M e obrigatorio so no bootstrap
.validate_var_args <- function(var_comp, var_method, M) {
  if (!is.logical(var_comp) || length(var_comp) != 1L || is.na(var_comp))
    stop("var_comp must be TRUE or FALSE", call. = FALSE)
  if (var_comp && var_method == "bootstrap") {
    if (is.null(M) || !is.numeric(M) || length(M) != 1L || !is.finite(M) ||
        M < 2 || M != round(M))
      stop("with var_method = \"bootstrap\", M (number of replicates) must be ",
           "an integer >= 2", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Interceptos e distribuicao dos residuos
# ---------------------------------------------------------------------------
# As equacoes de posto sao invariantes a deslocamento e nao identificam o
# intercepto. Depois do ajuste, a distribuicao dos residuos de cada
# componente e estimada e guardada no qast_fit; o intercepto e a media
# dela e predict(type = "survival") a usa em S(t | z) = S_res(h(t) - z'b).

# Tempo de sobrevivencia: curva de Kaplan-Meier dos residuos
# e = h1(t) - z'beta (degraus continuos a direita: S(x) = P(e > x)).
.resid_dist_km <- function(e, status) {
  f <- survival::survfit(survival::Surv(e, status) ~ 1)
  list(time = f$time, surv = f$surv)
}

# Tempo ajustado por qualidade: distribuicao de Zhao e Tsiatis dos
# residuos r = h2(qt) - z'theta. Cada evento recebe massa
# 1 / (n G(t_i)), com G o KM da censura na escala do tempo de
# sobrevivencia (G = 0 -> massa 0, como no kernel de surv_qt_fun);
# massas de residuos empatados sao somadas.
.resid_dist_ipcw <- function(r, t, status) {
  G <- km_weights_at(t, status)
  w <- ifelse(status == 1 & is.finite(G) & G > 0, 1 / G, 0) / length(r)
  keep <- w > 0
  mass <- tapply(w[keep], r[keep], sum)
  list(value = as.numeric(names(mass)), mass = unname(as.vector(mass)))
}

# Media da curva KM (Ding e Nan, 2015), restrita ao maior residuo: a massa
# que sobra no fim da curva (maior residuo censurado) fica nesse ponto.
# Calculada pelos saltos, e nao pelo rmean do survfit, que integra a partir
# de min(e) sem soma-lo de volta quando ha residuos negativos.
.km_mean <- function(d) {
  s <- d$surv
  k <- length(s)
  sum(-diff(c(1, s)) * d$time) + s[k] * d$time[k]
}

# Media IPCW (Zhao e Tsiatis, 2000): (1/n) sum status_i r_i / G(t_i)
.ipcw_mean <- function(d) sum(d$mass * d$value)

# S(x) = P(e > x) na curva KM; 1 antes do primeiro salto
.surv_km_at <- function(d, x) {
  idx <- findInterval(x, d$time)
  out <- d$surv[pmax(idx, 1L)]
  out[idx == 0L] <- 1
  out
}

# S(x) = P(r > x) = soma das massas dos residuos > x, limitada a [0, 1]
# (a massa total do IPCW e proxima de 1, mas nao exatamente 1)
.surv_ipcw_at <- function(d, x) {
  tail_mass <- rev(cumsum(rev(d$mass)))
  idx <- findInterval(x, d$value) + 1L
  out <- c(tail_mass, 0)[idx]
  pmin(pmax(out, 0), 1)
}

.qast_resid_dists <- function(lt, lq, t, status, z, betahat, thetahat) {
  list(beta  = .resid_dist_km(lt - drop(z %*% betahat), status),
       theta = .resid_dist_ipcw(lq - drop(z %*% thetahat), t, status))
}

# ---------------------------------------------------------------------------
# Classe qast_fit
# ---------------------------------------------------------------------------

.new_qast_fit <- function(model, call, z, status, betahat, thetahat,
                          vcov_beta, vcov_theta, var_comp, var_method, M,
                          settings, lt, lq, t, qt, formula_info = NULL) {
  nms <- colnames(z)
  names(betahat) <- names(thetahat) <- nms
  rd <- .qast_resid_dists(lt, lq, t, status, z, betahat, thetahat)
  fit <- list(beta = betahat, theta = thetahat,
              intercept = c(beta = .km_mean(rd$beta),
                            theta = .ipcw_mean(rd$theta)))
  if (var_comp) {
    dimnames(vcov_beta) <- dimnames(vcov_theta) <- list(nms, nms)
    fit$var_beta  <- stats::setNames(diag(vcov_beta), nms)
    fit$var_theta <- stats::setNames(diag(vcov_theta), nms)
    fit$vcov_beta  <- vcov_beta
    fit$vcov_theta <- vcov_theta
  }
  fit$model <- model
  fit$call <- call
  fit$n <- length(status)
  fit$n_events <- as.integer(sum(status))
  fit$var_method <- if (var_comp) var_method else NULL
  fit$M <- if (var_comp && var_method == "bootstrap") M else NULL
  fit$settings <- settings
  # Para predict(): covariaveis do ajuste, distribuicao dos residuos,
  # tempos de evento (grade padrao) e, nas formulas, termos e niveis
  fit$z <- z
  fit$y <- cbind(time = t, qtime = qt, status = status)
  fit$resid_dist <- rd
  fit$event_times <- list(beta  = sort(unique(t[status == 1])),
                          theta = sort(unique(qt[status == 1])))
  fit$formula_info <- formula_info
  class(fit) <- "qast_fit"
  fit
}

.qast_model_label <- function(model) {
  switch(model,
         fg   = "Fine and Gelber (2001) semiparametric regression model",
         prop = "Bivariate semiparametric regression model (proposed)")
}

#' Methods for QAST fits
#'
#' \code{semi_aft_fg()} and \code{semi_aft_prop()} return an object of class
#' \code{"qast_fit"}: a list with \code{beta} (survival time component),
#' \code{theta} (quality-adjusted survival time component),
#' \code{intercept} and, when \code{var_comp = TRUE}, \code{var_beta},
#' \code{var_theta} and the covariance matrices \code{vcov_beta} and
#' \code{vcov_theta}. It also holds \code{model}, \code{call}, \code{n},
#' \code{n_events}, \code{settings} and, for \code{\link{predict.qast_fit}}
#' and \code{\link{qast_pp}}, the data of the fit (\code{y}: columns
#' \code{time}, \code{qtime} and \code{status}; \code{z}: covariate
#' matrix), the estimated distributions of the residuals
#' (\code{resid_dist}), the distinct event times of each component
#' (\code{event_times}) and, for formula fits, the terms, factor levels
#' and contrasts (\code{formula_info}).
#'
#' @section Intercepts:
#' The rank-based estimating equations are invariant to a location shift,
#' so they do not identify an intercept. After the fit, each intercept is
#' estimated as the mean of the distribution of the residuals of its
#' component, and stored in \code{intercept}, a named vector with elements
#' \code{beta} and \code{theta}:
#' \itemize{
#'   \item \code{beta}: mean of the Kaplan-Meier curve of the residuals
#'     \eqn{h_1(T_i) - Z_i^\top\hat\beta}{h1(T_i) - Z_i'beta} (Ding and Nan,
#'     2015), restricted to the largest residual when it is censored;
#'   \item \code{theta}: inverse probability of censoring weighted mean of
#'     the residuals \eqn{h_2(Q_i) - Z_i^\top\hat\theta}{h2(Q_i) - Z_i'theta}
#'     (Zhao and Tsiatis, 2000), with the Kaplan-Meier censoring
#'     survival \eqn{\hat G}{G} evaluated at the survival time:
#'     \eqn{n^{-1}\sum_i \delta_i r_i / \hat G(T_i)}{(1/n) sum_i delta_i r_i / G(T_i)}.
#'     This is the mean of the distribution estimated by
#'     \code{\link{surv_qt_fun}}.
#' }
#' The linear predictor of the transformed time is then
#' \code{intercept["beta"] + z \%*\% beta} (and likewise with
#' \code{theta}). The intercepts have no standard errors and are not
#' part of \code{coef()} or \code{vcov()}.
#'
#' @references Ding, Y. and Nan, B. (2015). Estimating mean survival time:
#'   when is it possible? \emph{Scandinavian Journal of Statistics}, 42(2),
#'   397-413. \doi{10.1111/sjos.12112}
#'
#'   Zhao, H. and Tsiatis, A. A. (2000). Estimating mean quality adjusted
#'   lifetime with censored data. \emph{Sankhya: The Indian Journal of
#'   Statistics, Series B}, 62, 175-188.
#'
#' @param x,object Object of class \code{"qast_fit"}
#' @param digits Significant digits when printing
#' @param component \code{"beta"} (survival time) or \code{"theta"}
#'   (quality-adjusted time); in \code{coef()}, \code{"both"} returns a
#'   matrix with the two columns
#' @param signif.stars If TRUE, flags the p-values with stars
#' @param ... Not used
#' @return \code{coef()}: matrix (p x 2) or vector; \code{vcov()}:
#'   (jackknife or bootstrap) covariance matrix of the component;
#'   \code{summary()}: a \code{"summary.qast_fit"} object with the
#'   coefficient tables; \code{print()}: the object, invisibly.
#' @name qast_fit
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
#' coef(fit)
NULL

#' @rdname qast_fit
#' @export
print.qast_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...) {
  cat("\n", .qast_model_label(x$model), "\n\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("n = ", x$n, ", number of events = ", x$n_events, "\n", sep = "")
  if (x$model == "fg")
    cat("L = ", format(x$settings$L, digits = digits), ", phi = ",
        x$settings$phi, "\n", sep = "")
  cat("\nCoefficients:\n")
  print(coef(x), digits = digits)
  cat("\nIntercepts (mean of the residuals):\n")
  print(x$intercept, digits = digits)
  if (is.null(x$var_beta))
    cat("\n(no variance: refit with var_comp = TRUE for standard errors)\n")
  else
    cat("\nUse summary() for standard errors and tests (",
        x$var_method, ").\n", sep = "")
  invisible(x)
}

#' @rdname qast_fit
#' @export
coef.qast_fit <- function(object, component = c("both", "beta", "theta"),
                          ...) {
  component <- match.arg(component)
  switch(component,
         both  = cbind(beta = object$beta, theta = object$theta),
         beta  = object$beta,
         theta = object$theta)
}

#' @rdname qast_fit
#' @export
vcov.qast_fit <- function(object, component = c("beta", "theta"), ...) {
  component <- match.arg(component)
  if (is.null(object$vcov_beta))
    stop("fit has no variance: refit with var_comp = TRUE", call. = FALSE)
  if (component == "beta") object$vcov_beta else object$vcov_theta
}

.qast_coef_table <- function(est, v) {
  if (is.null(v)) {
    tab <- cbind(est, exp(est))
    colnames(tab) <- c("coef", "exp(coef)")
  } else {
    se <- sqrt(v)
    zval <- est / se
    tab <- cbind(est, exp(est), se, zval, 2 * stats::pnorm(-abs(zval)))
    colnames(tab) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  }
  rownames(tab) <- names(est)
  tab
}

#' @rdname qast_fit
#' @export
summary.qast_fit <- function(object, ...) {
  out <- list(
    model = object$model, call = object$call, n = object$n,
    n_events = object$n_events, settings = object$settings,
    var_method = object$var_method, M = object$M,
    intercept = object$intercept,
    beta  = .qast_coef_table(object$beta, object$var_beta),
    theta = .qast_coef_table(object$theta, object$var_theta)
  )
  class(out) <- "summary.qast_fit"
  out
}

#' @rdname qast_fit
#' @export
print.summary.qast_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {
  rule <- strrep("-", 62)
  cat("\n", rule, "\n", .qast_model_label(x$model), "\n", rule, "\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat("n = ", x$n, ", number of events = ", x$n_events, "\n", sep = "")
  if (x$model == "fg")
    cat("L = ", format(x$settings$L, digits = digits), ", phi = ",
        x$settings$phi, "\n", sep = "")
  if (!is.null(x$var_method))
    cat("Variance: ", x$var_method,
        if (!is.null(x$M)) paste0(" (M = ", x$M, ")"), "\n", sep = "")
  cat("\nSurvival time component\n")
  stats::printCoefmat(x$beta, digits = digits, signif.stars = signif.stars,
                      signif.legend = FALSE)
  cat("Intercept (Kaplan-Meier mean of the residuals): ",
      format(x$intercept[["beta"]], digits = digits), "\n", sep = "")
  cat(rule, "\nQuality-adjusted survival time component\n", sep = "")
  stats::printCoefmat(x$theta, digits = digits, signif.stars = signif.stars,
                      signif.legend = FALSE)
  cat("Intercept (IPCW mean of the residuals): ",
      format(x$intercept[["theta"]], digits = digits), "\n", sep = "")
  if (signif.stars && ncol(x$theta) > 2L)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

# Matriz de covariaveis para predict(): formula -> model.matrix com os
# termos, niveis e contrastes do ajuste; matriz -> colunas pelo nome (ou
# pela posicao, se o ajuste nao tinha nomes)
.qast_newdata_z <- function(object, newdata) {
  if (is.null(newdata)) return(object$z)
  p <- length(object$beta)
  nms <- names(object$beta)
  info <- object$formula_info
  if (!is.null(info)) {
    if (!is.data.frame(newdata))
      newdata <- as.data.frame(newdata)
    mf <- stats::model.frame(info$terms, newdata, xlev = info$xlevels,
                             na.action = stats::na.pass)
    z <- stats::model.matrix(info$terms, mf, contrasts.arg = info$contrasts)
    z <- z[, colnames(z) != "(Intercept)", drop = FALSE]
  } else {
    z <- as.matrix(newdata)
    if (!is.null(nms) && !is.null(colnames(z))) {
      miss <- setdiff(nms, colnames(z))
      if (length(miss))
        stop("newdata lacks the columns: ", paste(miss, collapse = ", "),
             call. = FALSE)
      z <- z[, nms, drop = FALSE]
    }
  }
  if (ncol(z) != p)
    stop(sprintf("newdata gives %d covariate columns; the fit has %d",
                 ncol(z), p), call. = FALSE)
  if (!is.numeric(z)) stop("newdata covariates must be numeric", call. = FALSE)
  if (anyNA(z)) .qast_stop_at("newdata contains NA", which(!stats::complete.cases(z)))
  z
}

# S(t | z_k) numa grade de tempos, para um componente
.qast_predict_surv <- function(object, z, comp, times) {
  link <- stats::make.link(object$settings[[if (comp == "beta") "h1" else "h2"]])
  if (is.null(times)) times <- object$event_times[[comp]]
  if (!is.numeric(times) || anyNA(times))
    stop("times must be a numeric vector without NA", call. = FALSE)
  ht <- link$linkfun(times)
  if (any(!is.finite(ht)))
    stop(sprintf("times must be in the domain of the link \"%s\" (with log, > 0)",
                 link$name), call. = FALSE)
  lp <- drop(z %*% object[[comp]])
  x <- outer(-lp, ht, `+`)          # residuo: h(t) - z'b
  surv_at <- if (comp == "beta") .surv_km_at else .surv_ipcw_at
  out <- matrix(surv_at(object$resid_dist[[comp]], as.vector(x)),
                nrow = nrow(x))
  dimnames(out) <- list(rownames(z), format(times, digits = 6L,
                                            trim = TRUE))
  out
}

#' Predictions from QAST fits
#'
#' Linear predictors and survival probabilities for the survival time
#' (\code{beta}) and quality-adjusted survival time (\code{theta})
#' components of a \code{"qast_fit"}.
#'
#' Under the model \eqn{h_1(T) = \alpha + Z^\top\beta + \varepsilon}{h1(T)
#' = alpha + Z'beta + e}, with the error independent of the covariates,
#' the survival function given the covariates is a shift of the residual
#' survival function:
#' \deqn{S(t \mid z) = P\{h_1(T) - Z^\top\beta > h_1(t) - z^\top\beta\} =
#'   S_{e}\{h_1(t) - z^\top\hat\beta\},}{S(t | z) = S_e(h1(t) - z'beta),}
#' with \eqn{e = h_1(T) - Z^\top\hat\beta}{e = h1(T) - Z'beta} (the
#' intercept cancels out). \eqn{S_e}{S_e} is estimated by:
#' \itemize{
#'   \item \code{beta}: the Kaplan-Meier curve of the residuals
#'     \eqn{h_1(T_i) - Z_i^\top\hat\beta}{h1(T_i) - Z_i'beta};
#'   \item \code{theta}: the estimator of Zhao and Tsiatis (1997) applied
#'     to the residuals \eqn{h_2(Q_i) - Z_i^\top\hat\theta}{h2(Q_i) -
#'     Z_i'theta}, which weights each event by the inverse of the
#'     Kaplan-Meier censoring survival at its survival time and corrects
#'     the induced informative censoring of the quality-adjusted time.
#'     Its total mass is close to, but not exactly, 1, so the estimates
#'     are truncated to [0, 1].
#' }
#' Both use the Kaplan-Meier (marginal) censoring distribution, which
#' assumes censoring independent of the covariates. The survival
#' probabilities are \eqn{P(T > t \mid z)}{P(T > t | z)} (right-continuous
#' steps). \code{surv_qt_fun()} evaluates each subject at its own
#' residual with \eqn{\ge}{>=}, i.e. \eqn{P(r \ge r_i)}{P(r >= r_i)}.
#'
#' @param object Object of class \code{"qast_fit"}
#' @param newdata Covariates to predict for. For fits from a formula, a
#'   data frame with the covariates of the formula (factors with the
#'   levels of the fit); for fits from a matrix, a matrix or data frame
#'   with the columns of \code{z} (matched by name when the fit had
#'   column names). \code{NULL} (default) uses the covariates of the fit.
#' @param type \code{"lp"} (linear predictor including the intercept,
#'   \code{intercept + z \%*\% coef}, on the scale of the link) or
#'   \code{"survival"} (survival probabilities at \code{times})
#' @param component \code{"both"}, \code{"beta"} (survival time) or
#'   \code{"theta"} (quality-adjusted time)
#' @param times Times at which to evaluate the survival probabilities, on
#'   the original scale (survival time for \code{beta}, quality-adjusted
#'   time for \code{theta}). \code{NULL} (default) uses the distinct
#'   event times of the fit data of each component.
#' @param ... Not used
#' @return For \code{type = "lp"}: a vector (one component) or a matrix
#'   with columns \code{beta} and \code{theta}. For
#'   \code{type = "survival"}: a matrix with one row per subject of
#'   \code{newdata} and one column per time (one component), or a list
#'   with the matrices \code{beta} and \code{theta}.
#' @references Zhao, H. and Tsiatis, A. A. (1997). A consistent estimator
#'   for the distribution of quality adjusted survival time.
#'   \emph{Biometrika}, 84(2), 339-348. \doi{10.1093/biomet/84.2.339}
#' @seealso \code{\link{qast_fit}} for the intercepts;
#'   \code{\link{surv_qt_fun}}
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(age = rnorm(n), trt = rbinom(n, 1, 0.5))
#' d$t <- rexp(n, rate = exp(-0.3 * d$trt))
#' d$status <- rbinom(n, 1, 0.8)
#' d$qt <- d$t * runif(n, 0.5, 1)
#' fit <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d)
#'
#' new <- data.frame(age = c(0, 0), trt = c(0, 1))
#' predict(fit, new)                       # linear predictors
#' s <- predict(fit, new, type = "survival", times = c(0.5, 1, 2))
#' s$beta                                  # P(T > t | z)
#' s$theta                                 # P(Q > t | z)
#' @export
predict.qast_fit <- function(object, newdata = NULL,
                             type = c("lp", "survival"),
                             component = c("both", "beta", "theta"),
                             times = NULL, ...) {
  type <- match.arg(type)
  component <- match.arg(component)
  if (is.null(object$resid_dist))
    stop("this fit has no stored residual distribution; refit it with ",
         "the current version of QAST", call. = FALSE)
  z <- .qast_newdata_z(object, newdata)
  comps <- if (component == "both") c("beta", "theta") else component

  if (type == "lp") {
    out <- vapply(comps, function(cp)
      object$intercept[[cp]] + drop(z %*% object[[cp]]),
      numeric(nrow(z)))
    if (!is.matrix(out)) out <- matrix(out, nrow = nrow(z),
                                       dimnames = list(NULL, comps))
    rownames(out) <- rownames(z)
    return(if (length(comps) == 1L) out[, 1L] else out)
  }

  out <- lapply(stats::setNames(comps, comps), function(cp)
    .qast_predict_surv(object, z, cp, times))
  if (length(comps) == 1L) out[[1L]] else out
}
