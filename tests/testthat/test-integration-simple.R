# tests/testthat/test-integration-simple.R

test_that("betahat e thetahat reproduzem valores de referencia do original (seed=42)", {
  skip_on_cran()

  # Mesmo gerador de dados do original. O intercepto (coluna de 1) entra
  # SO na geracao de t; o modelo nao estima intercepto (equacoes de
  # posto, invariantes a deslocamento), entao ele sai de z no ajuste.
  set.seed(42)
  n <- 300
  z_gen <- cbind(1, matrix(rnorm(n * 2), n, 2))
  t <- exp(as.vector(z_gen %*% c(0.5, -0.3, 0.2)) + rnorm(n, sd = 0.5))
  qt <- t * runif(n, 0.5, 1.0)
  status <- rbinom(n, 1, 0.85)
  z <- z_gen[, -1]

  fit <- semi_aft_prop(t, qt, status, z,
                       h1 = "log", h2 = "log",
                       method = "Nelder-Mead",
                       norm_min = "L2",
                       cens_weight = TRUE,
                       M = 10, var_comp = FALSE,
                       var_method = "jackknife")

  # Referencias do script original nos MESMOS dados (set.seed(42), n = 300),
  # apenas as inclinacoes. O original otimizava tambem a direcao plana do
  # intercepto (valor vindo de runif); sem ela, a trajetoria do
  # Nelder-Mead sobre o objetivo escalonado muda e o ponto final se move
  # ~2e-4 dentro do conjunto de minimos. Por isso: (1) inclinacoes em
  # 1e-3 e (2) objetivo no ponto estimado nao pior que na referencia.
  beta_ref  <- c(-0.320483557958109, 0.177030351648878)
  theta_ref <- c(-0.341859946407396, 0.182101852818395)

  expect_equal(unname(fit$beta),  beta_ref,  tolerance = 1e-3)
  expect_equal(unname(fit$theta), theta_ref, tolerance = 1e-3)

  lt <- log(t); lq <- log(qt)
  HG_surv <- QAST:::HG(X = z, t = t, status = status)
  o1 <- function(b) QAST:::U1_obj_cpp(b, lt, z, status, "L2")
  o2 <- function(th) QAST:::U2_obj_cpp(th, lq, z, status, HG_surv, "L2")
  expect_lte(o1(unname(fit$beta)),  o1(beta_ref)  + 1e-12)
  expect_lte(o2(unname(fit$theta)), o2(theta_ref) + 1e-12)
})

test_that("coluna constante (intercepto) em z e rejeitada nos dois modelos", {
  set.seed(1)
  n <- 40
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)

  expect_error(
    semi_aft_prop(t, qt, st, cbind(1, rnorm(n)), "log", "log",
                  "Nelder-Mead", "L2", M = 0, var_comp = FALSE,
                  var_method = "jackknife"),
    "intercept")

  z_fg <- cbind(zb, const = 0)
  expect_error(
    semi_aft_fg(t, qt, st, z_fg, "log", "log", "Nelder-Mead", "L2",
                L = max(t), phi = "identity", M = 0, var_comp = FALSE,
                var_method = "jackknife"),
    "const")

  # Nomes das colunas ofensoras aparecem na mensagem
  expect_error(QAST:::validate_no_constant_z(cbind(a = 1:3, b = 5)), "b")
  expect_silent(QAST:::validate_no_constant_z(zb))
})

test_that("colunas colineares geram warning (nao erro) nos dois modelos", {
  set.seed(3)
  n <- 50
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)

  # duplicada (FG) e combinacao linear (proposto)
  z_dup <- cbind(zb, trt2 = zb[, "trt"])
  x1 <- rnorm(n); x2 <- rnorm(n)
  z_lin <- cbind(x1 = x1, x2 = x2, x3 = x1 + 2 * x2)

  expect_warning(
    fit <- semi_aft_fg(t, qt, st, z_dup, "log", "log", "Nelder-Mead", "L2",
                       L = max(t), phi = "identity", M = 0, var_comp = FALSE,
                       var_method = "jackknife"),
    "full rank")
  expect_length(fit$beta, 3)

  expect_warning(
    semi_aft_prop(t, qt, st, z_lin, "log", "log", "Nelder-Mead", "L2",
                  M = 0, var_comp = FALSE, var_method = "jackknife"),
    "x3")

  expect_silent(QAST:::warn_collinear_z(zb))
})

test_that("z colinear: um unico aviso (sem repetir o do fast_path)", {
  set.seed(3)
  n <- 50
  zb <- cbind(trt = rbinom(n, 1, .5), sex = rbinom(n, 1, .4))
  t <- rexp(n); qt <- t * runif(n, .5, 1); st <- rbinom(n, 1, .8)
  w <- character()
  withCallingHandlers(
    semi_aft_fg(t, qt, st, cbind(zb, trt2 = zb[, "trt"]), "log", "log",
                "Nelder-Mead", "L2", L = max(t), phi = "identity", M = 0,
                var_comp = FALSE, var_method = "jackknife"),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "full rank")
})
