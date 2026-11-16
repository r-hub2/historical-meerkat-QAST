# tests/testthat/test-fg-equivalence-extended.R
#
# Equivalencia numerica do kernel FG (U2_parts_fg_cpp) com a referencia
# em R. Estende test-fg-equivalence.R (grade de theta x phi x norm) com:
#
#   1. grade de beta   - os pares pre-computados dependem de beta
#   2. grupo vazio / coluna constante / grupo sem eventos (path rep(1, n))
#   3. censura pesada com sobrevivencia KM chegando a zero (filtro denom)
#   4. empate exato na fronteira do indicador h1L (desigualdade estrita)
#   5. L menor (quantil 0.5) para exercitar a truncagem em mais pares
#   6. fumaca: n = 200, p = 1
#
# Todos os blocos comparam U2_parts_fg_cpp contra U2_parts_fg_ref
# (referencia R puro) e a agregacao da norma, tolerancia 1e-10.

# --------------------------------------------------------------------------
# Helper: comparacao completa (parts + obj agregado) para um ponto
# --------------------------------------------------------------------------
expect_fg_equivalent <- function(theta, lt, lq, z, status, beta, hL,
                                 km_surv, phi, norm_min, tol = 1e-10) {
  parts_cpp <- QAST:::U2_parts_fg_cpp(theta, lt, lq, z, status,
                                      beta, hL, km_surv, phi)
  parts_ref <- QAST:::U2_parts_fg_ref(theta, lt, lq, z, status,
                                      beta, hL, km_surv, phi)
  expect_equal(parts_cpp, parts_ref, tolerance = tol)
  
  obj_cpp <- QAST:::U2_obj_fg_cpp(theta, lt, lq, z, status,
                                  beta, hL, km_surv, phi, norm_min)
  if (norm_min == "L1") {
    obj_ref <- sum(abs(parts_ref)) / sqrt(length(parts_ref))
  } else {
    obj_ref <- sqrt(sum(parts_ref^2)) / sqrt(length(parts_ref))
  }
  expect_equal(obj_cpp, obj_ref, tolerance = tol)
}

# Percorre as grades de phi/norm/theta para uma configuracao fixa
compare_over_grids <- function(lt, lq, z, status, beta, hL, km_surv,
                               theta_grid) {
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      for (i in seq_len(nrow(theta_grid))) {
        expect_fg_equivalent(theta_grid[i, ], lt, lq, z, status,
                             beta, hL, km_surv, phi, norm_min)
      }
    }
  }
}

# --------------------------------------------------------------------------
# 1. Grade de beta x theta (caso base: aleatorio "limpo", L no q90)
# --------------------------------------------------------------------------
test_that("FG: kernel C++ equivale a referencia em grade de beta e theta", {
  set.seed(42)
  
  n <- 40
  z <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.4))
  colnames(z) <- c("trt", "sex")
  
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  
  lt <- log(t); lq <- log(qt)
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  # beta = 0 verifica o caso degenerado z'beta = 0 no indicador h1L
  beta_grid <- rbind(
    c(0.0,  0.0),
    c(0.3, -0.2),
    c(0.8,  0.5),   # efeitos fortes: empurra mais pares para fora de L
    c(-0.5, 0.7)    # sinais mistos
  )
  
  theta_grid <- rbind(
    c(0.0, 0.0),
    c(0.2, -0.1),
    c(-0.4, 0.5),
    c(0.7, 0.7)
  )
  
  for (b in seq_len(nrow(beta_grid))) {
    beta <- beta_grid[b, ]
    
    # L no quantil 0.9 (padrao do teste original)
    hL <- log(stats::quantile(t, 0.9))
    compare_over_grids(lt, lq, z, status, beta, hL, km_surv, theta_grid)
    
    # L menor: truncagem corta muitos mais pares
    hL_mid <- log(stats::quantile(t, 0.5))
    compare_over_grids(lt, lq, z, status, beta, hL_mid, km_surv, theta_grid)
  }
})

# --------------------------------------------------------------------------
# 2. Coluna constante / grupo vazio / grupo sem eventos
#    Coberturas: parts(j) = 0 (nk = 0 ou nl = 0) e o caminho rep(1, n)
#    de km_weights_at (grupo com status == 0 para todos).
# --------------------------------------------------------------------------
test_that("FG: coluna constante, grupo vazio e grupo sem eventos", {
  set.seed(7)
  
  n <- 60
  z <- cbind(rbinom(n, 1, 0.5),   # normal
             rep(0, n),           # constante 0  -> grupo 1 vazio
             rep(1, n))           # constante 1  -> grupo 0 vazio
  colnames(z) <- c("trt", "cte0", "cte1")
  
  t <- rexp(n)
  # Grupo sem eventos: todos os sujeitos com trt = 1 sao censurados.
  # (exerce o caminho rep(1, n) de km_weights_at dentro de prep_km_by_column)
  status <- ifelse(z[, 1] == 1, 0L, rbinom(n, 1, 0.8))
  
  qt <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  hL <- log(stats::quantile(t, 0.9))
  
  beta  <- c(0.3, 0.0, 0.0)
  theta <- c(0.2, 0.0, 0.0)
  
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      expect_fg_equivalent(theta, lt, lq, z, status, beta, hL,
                           km_surv, phi, norm_min)
    }
  }
  
  # As colunas constantes devem ter componente identicamente nula
  parts <- QAST:::U2_parts_fg_cpp(theta, lt, lq, z, status,
                                  beta, hL, km_surv, "identity")
  expect_identical(unname(parts[2]), 0.0)
  expect_identical(unname(parts[3]), 0.0)
})

# --------------------------------------------------------------------------
# 3. Censura pesada: sobrevivencia KM chega perto de / atinge zero
#    (exerce o filtro !(denom > 0) || !isfinite(denom))
# --------------------------------------------------------------------------
test_that("FG: censura pesada com pesos KM quase nulos", {
  set.seed(123)
  
  n <- 120
  z <- cbind(rbinom(n, 1, 0.6), rbinom(n, 1, 0.3))
  colnames(z) <- c("trt", "sex")
  
  t <- rexp(n, rate = exp(0.3 * z[, 1]))
  
  # 70% de censura + 1 sujeito censurado num tempo posterior ao ultimo
  # evento: garante cauda com S(t) muito baixa
  status <- rbinom(n, 1, 0.30)
  status[n] <- 0L
  t[n] <- max(t) * 10
  
  qt <- pmin(t, t * runif(n, 0.3, 1))
  lt <- log(t); lq <- log(qt)
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  # Sanity: pesos KH do ultimo sujeito devem ser validos (nao NaN/Inf)
  w0 <- km_surv[[1]]$w0
  expect_true(all(is.finite(w0)))
  
  hL <- log(stats::quantile(t, 0.9))
  compare_over_grids(lt, lq, z, status, c(0.3, -0.2), hL, km_surv,
                     rbind(c(0.0, 0.0), c(0.2, -0.1), c(0.5, 0.3)))
})

# --------------------------------------------------------------------------
# 4. Empate EXATO na fronteira do indicador h1L
#
#   O indicador usa desigualdade ESTRITA:
#     h1L = 1 se  hL - max(z'beta) > max(e1k, e1l)   senao 0
#
#   Construimos empates exatos com beta = 0 (z'beta identicamente 0) e
#   tempos DUPLICADOS iguais a L: para todo par cujos dois membros tenham
#   t == L, max(e1k, e1l) == hL bit a bit (mesma chamada a log()), e a
#   comparacao estrita deve dar 0. Pares com max(e1) logo abaixo devem
#   dar 1. E o caso mais sensivel a reordenacao de ponto flutuante: os
#   empates exatos construidos aqui devem permanecer identicos a
#   referencia.
# --------------------------------------------------------------------------
test_that("FG: empate exato na fronteira de h1L (desigualdade estrita)", {
  set.seed(99)
  
  n <- 50
  z <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.4))
  colnames(z) <- c("trt", "sex")
  
  t <- rexp(n)
  status <- rbinom(n, 1, 0.8)
  
  # Duplica tempos no valor de L: 4 sujeitos com t exatamente igual a L
  L  <- stats::quantile(t, 0.7, names = FALSE)
  idx_dup <- c(1, 2, 3, 4)
  t[idx_dup] <- L
  
  # beta = 0 => z'beta == 0 exatamente; o indicador vira log(L) > max(e1)
  beta <- c(0.0, 0.0)
  
  qt <- pmin(t, t * runif(n, 0.4, 1))
  lt <- log(t); lq <- log(qt)
  hL <- log(L)  # mesma operacao log() aplicada ao mesmo double => bit a bit igual a lt dos duplicados
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  theta_grid <- rbind(c(0.0, 0.0), c(0.2, -0.1))
  
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      for (i in seq_len(nrow(theta_grid))) {
        expect_fg_equivalent(theta_grid[i, ], lt, lq, z, status,
                             beta, hL, km_surv, phi, norm_min)
      }
    }
  }
})

# --------------------------------------------------------------------------
# 6. Fumaca: p = 1, n = 200
# --------------------------------------------------------------------------
test_that("FG: fumaca p = 1, n = 200", {
  set.seed(2024)
  
  n <- 200
  z <- matrix(rbinom(n, 1, 0.45), n, 1)
  colnames(z) <- "x1"
  
  t      <- rexp(n, rate = exp(0.2 * z[, 1]))
  status <- rbinom(n, 1, 0.75)
  qt     <- pmin(t, t * runif(n, 0.5, 1))
  
  lt <- log(t); lq <- log(qt)
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  hL <- log(stats::quantile(t, 0.9))
  
  compare_over_grids(lt, lq, z, status, c(0.2), hL, km_surv,
                     rbind(c(0.0), c(0.3), c(-0.2)))
})
