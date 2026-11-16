library(testthat)

mk_pairs_case <- function(seed = 11, n = 60, p = 2, tie_at_L = FALSE) {
  set.seed(seed)
  z <- cbind(rbinom(n, 1, 0.45), rbinom(n, 1, 0.35))
  colnames(z) <- c("trt", "sex")
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  L      <- stats::quantile(t, 0.9)
  if (tie_at_L) {
    # Empate exato: beta = 0 => z'beta = 0 exato;
    # tempos duplicados exatamente em L; log() do mesmo double.
    t[1:4] <- L
    status[1:4] <- 1L
  }
  list(t = t, qt = qt, status = status, z = z,
       lt = log(t), lq = log(qt), hL = log(L))
}

pairs_theta_grid <- rbind(c(0.0, 0.0), c(0.2, -0.1),
                          c(-0.4, 0.5), c(0.7, 0.7))
pairs_beta_grid  <- rbind(c(0.0, 0.0), c(0.3, -0.2), c(-0.5, 0.8))


test_that("Pares (identity): equivale ao legado e a referencia em grade", {
  skip_if_not(exists("U2_parts_fg_ref", envir = asNamespace("QAST")),
              "referencia R indisponivel")
  
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      
      cas   <- mk_pairs_case(seed = 10 + bi)
      beta  <- pairs_beta_grid[bi, ]
      theta <- pairs_theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity")
      
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, "identity")
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "identity")
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
      
      # Triangulacao: kernel de referencia vs referencia R pura
      parts_ref <- QAST:::U2_parts_fg_ref(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "identity")
      expect_equal(parts_leg, parts_ref, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_new <- QAST:::eval_fg_obj_cpp(pairs, cas$lq, cas$z,
                                          theta, "identity", nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq,
                                        cas$z, cas$status, beta,
                                        cas$hL, km, "identity", nm)
        expect_equal(obj_new, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Pares (piecewise): equivale ao legado em grade de beta x theta", {
  # phi piecewise congela o INDICADOR h1L por par na construcao, mas o
  # sinal de e2[pa]-e2[pb] e reavaliado a cada theta — aqui o caminho
  # por pares deve reproduzir o kernel de referencia (mesma ordem de soma,
  # mesma divisao por n2 no fim), logo tolerancia de 1e-10.
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      
      cas   <- mk_pairs_case(seed = 50 + bi)
      beta  <- pairs_beta_grid[bi, ]
      theta <- pairs_theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise")
      
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, "piecewise")
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, "piecewise")
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_new <- QAST:::eval_fg_obj_cpp(pairs, cas$lq, cas$z,
                                          theta, "piecewise", nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq,
                                        cas$z, cas$status, beta,
                                        cas$hL, km, "piecewise", nm)
        expect_equal(obj_new, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Pares: coluna constante gera grupo vazio e parts identicamente 0", {
  cas <- mk_pairs_case()
  cas$z[, 2] <- 0L   # grupo 1 vazio na coluna sex
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  
  for (ph in c("identity", "piecewise")) {
    pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                      cas$status, beta, cas$hL,
                                      km, ph)
    col <- pairs[[2]]
    expect_length(col$kb, 0)
    # Ra alinhado a ka (zeros), nao vazio
    expect_length(col$Ra, length(col$ka))
    expect_true(all(col$Ra == 0))
    expect_length(col$pc, 0)
    
    th <- c(0.2, -0.1)
    parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z, th, ph)
    expect_identical(parts_new[2], 0.0)
    
    parts_leg <- QAST:::U2_parts_fg_cpp(th, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, ph)
    expect_equal(parts_new, parts_leg, tolerance = 1e-10)
  }
})

test_that("Pares: gate de empate exato em h1L preservado na construcao", {
  # O indicador h1L usa desigualdade ESTRETA resolvida NA CONSTRUCAO
  # com a mesma aritmetica do kernel de referencia (dot por membro — o
  # valor e deterministicamente identico ao dot por par). Este teste
  # falha se qualquer reordenacao de floating point mover um par de
  # lado do gate: os parts divergem por um SALTO DISCRETO, nao gradual
  # — por isso a tolerancia apertada, sem escala relativa.
  cas <- mk_pairs_case(tie_at_L = TRUE)
  beta <- c(0, 0)
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  for (ph in c("identity", "piecewise")) {
    pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                      cas$status, beta, cas$hL,
                                      km, ph)
    for (ti in seq_len(nrow(pairs_theta_grid))) {
      theta <- pairs_theta_grid[ti, ]
      parts_new <- QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                            theta, ph)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                          cas$z, cas$status, beta,
                                          cas$hL, km, ph)
      expect_equal(parts_new, parts_leg, tolerance = 1e-10)
    }
  }
})

test_that("Pares: argumentos invalidos param com mensagem clara", {
  cas <- mk_pairs_case()
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  
  # phi invalido
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, c(0, 0), cas$hL, km, "banana"),
               "phi")
  
  # beta de comprimento errado
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, c(0.1), cas$hL, km, "identity"),
               "beta")
  
  # objeto sem meta avaliado
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    c(0, 0), cas$hL, km, "identity")
  attr(pairs, "meta") <- NULL
  expect_error(QAST:::eval_fg_parts_cpp(pairs, cas$lq, cas$z,
                                        c(0, 0), "identity"),
               "meta")
})
test_that("PERF-1: identity nao guarda pares; Ra/Cb bit-identicos ao piecewise", {
  cas <- mk_pairs_case(seed = 21)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  for (bi in seq_len(nrow(pairs_beta_grid))) {
    beta <- pairs_beta_grid[bi, ]
    # algo = "loop": mesma ordem de soma do piecewise => bit-identico.
    # algo = "auto" (CDQ): mesma soma, outra ordem => 1e-13.
    p_id <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "identity", "loop")
    p_cq <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "identity")
    p_pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                     beta, cas$hL, km, "piecewise")
    expect_false(attr(p_id, "meta")$has_pairs)
    expect_true(attr(p_pw, "meta")$has_pairs)
    for (j in seq_along(p_id)) {
      expect_identical(p_id[[j]]$Ra, p_pw[[j]]$Ra)
      expect_identical(p_id[[j]]$Cb, p_pw[[j]]$Cb)
      expect_equal(p_cq[[j]]$Ra, p_pw[[j]]$Ra, tolerance = 1e-13)
      expect_equal(p_cq[[j]]$Cb, p_pw[[j]]$Cb, tolerance = 1e-13)
      expect_length(p_id[[j]]$pc, 0)
      expect_type(p_pw[[j]]$pa, "integer")
    }
    # Objeto identity recusa avaliacao piecewise (nunca zero silencioso)
    expect_error(QAST:::eval_fg_parts_cpp(p_id, cas$lq, cas$z, c(0, 0),
                                          "piecewise"), "identity")
    expect_error(QAST:::compile_fg_pw_cpp(p_id), "identity")
  }
})

test_that("FIX-4: pares respeitam a < b (posicional), como o legado", {
  cas <- mk_pairs_case(seed = 31)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  pr  <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                  c(0.1, -0.1), cas$hL, km, "piecewise")
  for (j in seq_along(pr)) {
    pos_a <- match(pr[[j]]$pa, pr[[j]]$ka)
    pos_b <- match(pr[[j]]$pb, pr[[j]]$kb)
    expect_true(all(pos_a < pos_b))
  }
})

test_that("semi_aft_fg: caminhos piecewise e fast_path rodam (regressao B1/B2)", {
  cas <- mk_pairs_case(seed = 41, n = 80)
  args <- list(t = cas$t, qt = cas$qt, status = cas$status, z = cas$z,
               h1 = "log", h2 = "log", method = "Nelder-Mead",
               norm_min = "L2", L = exp(cas$hL), M = 0,
               var_comp = FALSE, var_method = "jackknife")

  # piecewise padrao == otimizacao direta sobre o kernel de referencia
  # (mesmo objetivo, mesmo otimizador, mesmo valor inicial)
  fit_pw <- do.call(semi_aft_fg, c(args, phi = "piecewise"))
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  b   <- unname(fit_pw$beta)
  f_leg <- function(th) QAST:::U2_fun_fg_cpp(th, cas$lq, cas$lt, cas$z,
                                             cas$status, b, cas$hL, km,
                                             "piecewise", "L2")
  th0 <- unname(stats::lm(cas$lq ~ cas$z)$coef[-1])
  th_leg <- dfoptim::nmk(th0, f_leg)$par
  expect_equal(unname(fit_pw$theta), th_leg, tolerance = 1e-8)

  # forma fechada nao e pior que o otimizador no mesmo objetivo afim
  fit_fast <- do.call(semi_aft_fg, c(args, phi = "identity", fast_path = TRUE))
  fit_slow <- do.call(semi_aft_fg, c(args, phi = "identity", fast_path = FALSE))
  expect_length(fit_fast$theta, 2)
  expect_true(all(is.finite(fit_fast$theta)))

  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  pr  <- QAST:::make_fg_pairs(cas$lt, cas$lq, cas$z, cas$status,
                              unname(fit_fast$beta), cas$hL, km, "identity")
  aff <- QAST:::make_fg_affine(pr, cas$lq, cas$z)
  obj <- function(th) QAST:::eval_fg_obj_affine(aff, unname(th), "L2")
  expect_lte(obj(fit_fast$theta), obj(fit_slow$theta) + 1e-10)
})

# ----------------------------------------------------------------------------
# Marginais Ra/Cb em O(N log^2 N) == laco O(nk*nl)
# ----------------------------------------------------------------------------

cdq_vs_loop <- function(lt, lq, z, st, beta, hL, km) {
  a <- QAST:::make_fg_pairs_cpp(lt, lq, z, st, beta, hL, km, "identity", "auto")
  b <- QAST:::make_fg_pairs_cpp(lt, lq, z, st, beta, hL, km, "identity", "loop")
  for (j in seq_along(b)) {
    for (f in c("Ra", "Cb")) {
      # mesmo suporte (conjunto de pares identico) e mesmos valores
      expect_identical(a[[j]][[f]] == 0, b[[j]][[f]] == 0)
      expect_equal(a[[j]][[f]], b[[j]][[f]], tolerance = 1e-13)
    }
  }
  invisible(attr(a, "meta")$cdq_cols)
}

test_that("PERF-3: CDQ reproduz o laco (empates, corte em L, beta = 0)", {
  for (sd in 1:30) {
    set.seed(sd)
    n <- c(7, 40, 150, 400)[sd %% 4 + 1]; p <- sd %% 3 + 1
    z <- matrix(rbinom(n * p, 1, 0.5), n, p)
    z[1, ] <- 0; z[2, ] <- 1                     # nenhum grupo vazio
    t  <- if (sd %% 3 == 0) sample(1:5, n, TRUE) + 0 else rexp(n)
    st <- rbinom(n, 1, 0.75)
    L  <- if (sd %% 2 == 0) sort(t)[ceiling(0.7 * n)] else max(t)
    beta <- if (sd %% 5 == 0) rep(0, p) else round(rnorm(p), 1)
    km <- QAST:::prep_km_by_column(z, t, st)
    used <- cdq_vs_loop(log(t), log(t * 0.8), z, st, beta, log(L), km)
    expect_true(all(used))
  }
  # caso de empate exato no corte de L
  cas <- mk_pairs_case(tie_at_L = TRUE)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, c(0, 0), cas$hL, km)
  expect_true(all(used))
})

test_that("PERF-3: fallback para o laco com nao finitos ou pesos fora da faixa", {
  cas <- mk_pairs_case(seed = 5)
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.2, -0.1)

  lt_inf <- cas$lt; lt_inf[3] <- -Inf             # ex.: t = 0
  used <- cdq_vs_loop(lt_inf, cas$lq, cas$z, cas$status, beta, cas$hL, km)
  expect_false(any(used))

  # pesos patologicos em membros COM evento (censurados nao formam pares
  # e o CDQ nem consulta o peso deles — continua exato)
  ev0 <- which(cas$status[cas$z[, 1] == 0] == 1)[1]   # grupo 0, coluna 1
  ev1 <- which(cas$status[cas$z[, 2] == 1] == 1)[1]   # grupo 1, coluna 2
  km_tiny <- km; km_tiny[[1]]$w0[ev0] <- 1e-200   # produto poderia estourar
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_tiny)
  expect_identical(used, c(FALSE, TRUE))

  km_neg <- km; km_neg[[2]]$w1[ev1] <- -0.5       # peso negativo
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_neg)
  expect_identical(used, c(TRUE, FALSE))

  km_zero <- km; km_zero[[1]]$w0[ev0] <- 0        # peso zero: CDQ trata
  used <- cdq_vs_loop(cas$lt, cas$lq, cas$z, cas$status, beta, cas$hL, km_zero)
  expect_true(all(used))
})

test_that("PERF-3: algo invalido e piecewise (sempre laco)", {
  cas <- mk_pairs_case()
  km  <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  expect_error(QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                        c(0, 0), cas$hL, km, "identity", "x"),
               "algo")
  pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                 c(0, 0), cas$hL, km, "piecewise")
  expect_false(any(attr(pw, "meta")$cdq_cols))
})

test_that("fast_path = TRUE e o padrao de semi_aft_fg", {
  expect_true(formals(semi_aft_fg)$fast_path)
  cas <- mk_pairs_case(seed = 43, n = 80)
  args <- list(t = cas$t, qt = cas$qt, status = cas$status, z = cas$z,
               h1 = "log", h2 = "log", method = "Nelder-Mead",
               norm_min = "L2", L = exp(cas$hL), M = 0, phi = "identity",
               var_comp = FALSE, var_method = "jackknife")
  expect_identical(do.call(semi_aft_fg, args)$theta,
                   do.call(semi_aft_fg, c(args, fast_path = TRUE))$theta)
})
