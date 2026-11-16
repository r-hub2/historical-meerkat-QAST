# tests/testthat/test-km-findinterval-ab.R
#
# O oraculo e a avaliacao por which() (.km_weights_at_legacy); o caminho
# com findInterval deve reproduzi-lo em todos os casos de borda. O
# estimador e o survfit nos dois; o que se valida e a avaliacao da step
# function.

km_compare <- function(t, status, tol = 1e-12) {
  expect_equal(QAST:::km_weights_at(t, status),
               QAST:::km_weights_at(t, status, use_legacy = TRUE),
               tolerance = tol)
  invisible(TRUE)
}

test_that("KM findInterval: cenario base com censura moderada", {
  set.seed(42)
  n <- 150
  km_compare(rexp(n), rbinom(n, 1, 0.8))
})

test_that("KM findInterval: empates abundantes (tempos duplicados)", {
  set.seed(1)
  n <- 200
  t <- sample(1:10, n, replace = TRUE) * 1.0
  km_compare(t, rbinom(n, 1, 0.7))
})

test_that("KM findInterval: empates exatos entre falha e censura", {
  # A convencao de empates e do survfit; ambos os caminhos
  # consultam o MESMO fit, apenas a interpolacao e validada
  t      <- c(1, 1, 1, 2, 2, 3, 3, 3, 3)
  status <- c(0, 1, 1, 0, 1, 0, 0, 1, 1)
  km_compare(t, status)
})

test_that("KM findInterval: censura pesada com cauda longa", {
  set.seed(123)
  n <- 120
  t <- rexp(n)
  status <- rbinom(n, 1, 0.30)
  t[n] <- max(t) * 10
  status[n] <- 0L
  km_compare(t, status)
  
  # Sujeito alem da cauda: AMBOS os caminhos retornam o ULTIMO valor da
  # curva (nunca NaN)
  w_new <- QAST:::km_weights_at(t, status)
  fit <- survival::survfit(survival::Surv(t, 1 - status) ~ 1)
  expect_identical(w_new[n], fit$surv[length(fit$surv)])
})

test_that("KM findInterval: grupo sem eventos (rep(1, n))", {
  set.seed(7)
  t <- rexp(50)
  expect_identical(QAST:::km_weights_at(t, rep(0L, 50)), rep(1, 50))
})

test_that("KM findInterval: grupo sem censuras (S == 1)", {
  set.seed(9)
  t <- rexp(80)
  status <- rep(1L, 80)
  km_compare(t, status)
  expect_true(all(QAST:::km_weights_at(t, status) == 1))
})

test_that("KM findInterval: grupo vazio e antes do primeiro salto", {
  expect_identical(QAST:::km_weights_at(numeric(0), integer(0)),
                   numeric(0))
  
  t <- c(1, 10, 20, 30)
  status <- c(1, 0, 1, 1)
  expect_identical(QAST:::km_weights_at(t, status)[1], 1)
})

test_that("KM findInterval: t com nomes/atributos nao propagam lixo", {
  t <- c(a = 1.5, b = 2.5, c = 0.5)
  status <- c(1, 0, 1)
  w <- QAST:::km_weights_at(t, status)
  expect_equal(length(w), 3L)
  expect_null(names(w))
})

test_that("prep_km_by_column: dedupe de colunas identicas e consistente", {
  set.seed(99)
  n <- 100
  z <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.3))
  z <- cbind(z, z[, 1])           # terceira coluna duplicada da primeira
  colnames(z) <- c("a", "b", "a_dup")
  t <- rexp(n)
  status <- rbinom(n, 1, 0.75)
  
  km_dd  <- QAST:::prep_km_by_column(z, t, status, dedupe_cols = TRUE)
  km_nod <- QAST:::prep_km_by_column(z, t, status, dedupe_cols = FALSE)
  km_leg <- QAST:::prep_km_by_column(z, t, status, use_legacy = TRUE)
  
  for (j in 1:3) {
    expect_equal(km_dd[[j]]$w0, km_leg[[j]]$w0, tolerance = 1e-12)
    expect_equal(km_dd[[j]]$w1, km_leg[[j]]$w1, tolerance = 1e-12)
    expect_equal(km_nod[[j]]$w0, km_leg[[j]]$w0, tolerance = 1e-12)
  }
  
  # Colunas duplicadas produzem pesos identicos com e sem dedupe
  expect_equal(km_dd[[1]]$w0, km_dd[[3]]$w0, tolerance = 0)
  
  # E a chave de dedupe realmente colapsa (sem NA): saida tem comprimento p
  expect_equal(length(km_dd), 3L)
})

test_that("prep_km_by_column: golden chain - pesos alimentam kernel FG", {
  # Fumaca: o pipeline FG inteiro (prep_km -> pares -> afim) produz
  # objetivo identico com pesos KM por findInterval e por which()
  set.seed(2024)
  n <- 120
  z <- cbind(rbinom(n, 1, 0.45), rbinom(n, 1, 0.35))
  colnames(z) <- c("trt", "sex")
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  hL <- log(quantile(t, 0.9))
  beta <- c(0.3, -0.2)
  
  pairs_fast <- QAST:::make_fg_pairs(lt, lq, z, status, beta, hL,
                                     QAST:::prep_km_by_column(z, t, status),
                                     "identity")
  pairs_leg <- QAST:::make_fg_pairs(lt, lq, z, status, beta, hL,
                                    QAST:::prep_km_by_column(z, t, status,
                                                             use_legacy = TRUE),
                                    "identity")
  aff_fast <- QAST:::make_fg_affine(pairs_fast, lq, z)
  aff_leg  <- QAST:::make_fg_affine(pairs_leg, lq, z)
  
  for (th in list(c(0, 0), c(0.2, -0.1))) {
    expect_equal(QAST:::eval_fg_obj_affine(aff_fast, unlist(th), "L2"),
                 QAST:::eval_fg_obj_affine(aff_leg,  unlist(th), "L2"),
                 tolerance = 1e-12)
  }
})


test_that("prep_km_by_column: dedupe por impressao digital == sem dedupe", {
  set.seed(5)
  n <- 120
  base <- cbind(rbinom(n, 1, .5), rbinom(n, 1, .3), rbinom(n, 1, .6))
  # duplicatas intercaladas, inclusive da ultima e complementares (1 - x)
  z <- cbind(base[, 1], base[, 2], base[, 1], 1 - base[, 2], base[, 3],
             base[, 2], base[, 3], 1 - base[, 2])
  t <- rexp(n); status <- rbinom(n, 1, .75)
  km_dd  <- QAST:::prep_km_by_column(z, t, status, dedupe_cols = TRUE)
  km_nod <- QAST:::prep_km_by_column(z, t, status, dedupe_cols = FALSE)
  expect_identical(km_dd, km_nod)
})
