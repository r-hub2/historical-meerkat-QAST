# tests/testthat/test-surv-qt-equivalence.R
#
# Testa a equivalencia entre a implementacao original em R puro
# (preservada aqui como oracle) e o kernel C++ otimizado.
#
# O oracle reproduz fielmente o script original, inclusive os
# comportamentos de borda que o kernel replica propositalmente:
#   - summary(km_fit, times = time) retorna NA para tempos alem do
#     ultimo evento -> o na.rm = TRUE do original descarta a
#     contribuicao -> o kernel trata o peso como 0;
#   - empates (ties) em qalys recebem a mesma probabilidade
#     (comparacao >= do original).

# Oracle: copia fiel da implementacao original em R puro
surv_qt_ref <- function(QALY, object, pred_lin_qt = NULL) {
  time   <- object[, 1]
  status <- object[, 2]
  if (is.null(pred_lin_qt)) {
    qalys <- QALY
  } else {
    qalys <- log(QALY) - pred_lin_qt
  }
  km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  cens_surv <- summary(object = km_fit, times = time)$surv
  pred_surv <- sapply(1:length(qalys), function(i) {
    sum(x = ((status * (qalys >= qalys[i])) / cens_surv), na.rm = TRUE) /
      length(qalys)
  })
  pred_surv
}

gera_dados_survqt <- function(n = 50, seed = 42) {
  set.seed(seed)
  z1     <- rbinom(n, 1, 0.5)
  t      <- rexp(n, rate = exp(0.3 * z1))
  cens   <- rexp(n, rate = 0.2)
  time   <- pmin(t, cens)
  status <- as.integer(t <= cens)
  qaly   <- time * runif(n, 0.4, 1)
  list(time = time, status = status, qaly = qaly, z = z1)
}

test_that("Kernel C++ equivale ao oracle em R (sem preditor linear)", {
  d <- gera_dados_survqt()
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("Kernel C++ equivale ao oracle em R (com preditor linear)", {
  d <- gera_dados_survqt(n = 60, seed = 7)
  
  lin_pred <- d$z * 0.2
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status),
                     pred_lin_qt = lin_pred)
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status,
                            lin_pred = lin_pred)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("Empates em qaly recebem a mesma probabilidade", {
  d <- gera_dados_survqt(n = 30, seed = 99)
  # Forca empates: arredonda os QALYs para uma grade grosseira
  d$qaly <- round(d$qaly, 1)
  
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  
  expect_equal(cpp, ref, tolerance = 1e-10)
  
  # Verificacao direta do invariante: mesma entrada, mesma saida
  if (any(duplicated(d$qaly))) {
    tab <- tapply(cpp, d$qaly, function(x) length(unique(round(x, 12))))
    expect_true(all(tab == 1))
  }
})

test_that("Tempos alem do ultimo evento do KM tem peso zero (NA do original)", {
  d <- gera_dados_survqt(n = 20, seed = 5)
  
  # Forca um tempo muito grande: alem do maior tempo do KM da censura.
  # O summary() do oracle retornaria NA e o na.rm = TRUE descartaria
  # a contribuicao; o kernel deve replicar esse comportamento.
  d$time[1] <- max(d$time) * 10
  d$qaly[1] <- d$qaly[1] * 10
  
  ref <- surv_qt_ref(QALY = d$qaly, object = cbind(d$time, d$status))
  cpp <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time, status = d$status)
  
  expect_equal(cpp, ref, tolerance = 1e-10)
})

test_that("Wrapper valida entradas", {
  d <- gera_dados_survqt(n = 10)
  
  expect_error(
    QAST:::surv_qt_fun(qaly = d$qaly, time = d$time[-1], status = d$status),
    "same length as qaly"
  )
  expect_error(
    QAST:::surv_qt_fun(qaly = -d$qaly, time = d$time, status = d$status,
                       lin_pred = rep(0, 10)),
    "strictly positive"
  )
})

test_that("km_fit externo e reutilizado corretamente", {
  d <- gera_dados_survqt(n = 40, seed = 3)
  
  km_fit <- survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1)
  
  cpp_interno <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time,
                                    status = d$status)
  cpp_externo <- QAST:::surv_qt_fun(qaly = d$qaly, time = d$time,
                                    status = d$status, km_fit = km_fit)
  
  expect_equal(cpp_interno, cpp_externo, tolerance = 1e-12)
})