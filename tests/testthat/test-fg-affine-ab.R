# Correcoes relativas a versao anterior:
#   - mk_fg_case agora HONRA o argumento p (era ignorado; construia
#     sempre 2 colunas — bug de harness que quebrou o fast path);
#   - nomes de coluna genericos "z1..zp" (p != 2 nao tem "trt"/"sex");
#   - fast path usa beta coerente com p=1;
#   - gate de empate exato simplificado (sem a linha obscura com
#     tie_safe <- TRUE da versao anterior).
# ----------------------------------------------------------------------------

library(testthat)

mk_fg_case <- function(seed = 42, n = 60, p = 2, tie_at_L = FALSE) {
  stopifnot(is.integer(as.integer(p)), p >= 1)
  set.seed(seed)
  # FIX B: p respeitado — p colunas binarias independentes
  z <- matrix(rbinom(n * p, 1, 0.45), nrow = n, ncol = p)
  colnames(z) <- paste0("z", seq_len(p))
  
  # Coeficiente de sobrevivencia "real" so na primeira coluna (efeito
  # em z1), para manter o caso identificavel mesmo com p > 2
  t      <- rexp(n, rate = exp(0.3 * z[, 1]))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  L      <- stats::quantile(t, 0.9)
  
  if (tie_at_L) {
    # Empate exato: beta = 0 => z'beta = 0 exato;
    # tempos duplicados exatamente em L; log() do mesmo double.
    t[1:min(4, n)]        <- L
    status[1:min(4, n)]   <- 1L
  }
  list(t = t, qt = qt, status = status, z = z,
       lt = log(t), lq = log(qt), hL = log(L))
}

# Grade de beta coerente com p arbitrario: primeiro elemento varia,
# restantes zerados (para p=1 vira c(0), c(0.3), c(-0.5)).
mk_fg_beta_grid <- function(p) {
  base <- rep(0, p)
  rbind(base,
        replace(base, 1, 0.3),
        replace(base, 1, -0.5))
}

fg_theta_grid <- function(p) {
  base <- rep(0, p)
  g <- rbind(base)
  for (k in seq_len(p)) {
    g <- rbind(g, replace(base, k, 0.25),
               replace(base, k, -0.4))
  }
  g
}


test_that("Afim: equivale ao legado e a referencia R em grade de beta x theta", {
  skip_if_not(exists("U2_parts_fg_ref", envir = asNamespace("QAST")),
              "referencia R indisponivel")
  
  p <- 2
  beta_grid  <- mk_fg_beta_grid(p)
  theta_grid <- fg_theta_grid(p)
  
  for (bi in seq_len(nrow(beta_grid))) {
    for (ti in seq_len(nrow(theta_grid))) {
      
      cas <- mk_fg_case(seed = 40 + bi, n = 60, p = p)
      beta  <- beta_grid[bi, ]
      theta <- theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity")
      aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
      
      parts_aff <- QAST:::eval_fg_parts_affine_cpp(aff, theta)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")
      expect_equal(parts_aff, parts_leg, tolerance = 1e-10)
      
      # Triangulacao com a referencia R pura
      parts_ref <- QAST:::U2_parts_fg_ref(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")
      expect_equal(parts_leg, parts_ref, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_aff <- QAST:::eval_fg_obj_affine_cpp(aff, theta, nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "identity", nm)
        expect_equal(obj_aff, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Piecewise compilado: equivale ao legado em grade de beta x theta", {
  p <- 2
  beta_grid  <- mk_fg_beta_grid(p)
  theta_grid <- fg_theta_grid(p)
  
  for (bi in seq_len(nrow(beta_grid))) {
    for (ti in seq_len(nrow(theta_grid))) {
      
      cas <- mk_fg_case(seed = 80 + bi, n = 60, p = p)
      beta  <- beta_grid[bi, ]
      theta <- theta_grid[ti, ]
      
      km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
      pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise")
      pw <- QAST:::compile_fg_pw_cpp(pairs)
      
      parts_pw  <- QAST:::eval_fg_parts_pw_cpp(pw, cas$lq, cas$z, theta)
      parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "piecewise")
      expect_equal(parts_pw, parts_leg, tolerance = 1e-10)
      
      for (nm in c("L1", "L2")) {
        obj_pw  <- QAST:::eval_fg_obj_pw_cpp(pw, cas$lq, cas$z, theta, nm)
        obj_leg <- QAST:::U2_obj_fg_cpp(theta, cas$lt, cas$lq, cas$z,
                                        cas$status, beta, cas$hL,
                                        km, "piecewise", nm)
        expect_equal(obj_pw, obj_leg, tolerance = 1e-10)
      }
    }
  }
})

test_that("Coluna constante: componente identicamente nula em todos os caminhos", {
  # Coluna degenerada: Ra/Cb sao zeros alinhados a ka/kb; o
  # make_fg_affine_cpp deve aceitar a coluna e a componente e
  # identicamente 0.
  cas <- mk_fg_case()
  cas$z[, 2] <- 0L   # coluna 2 constante => grupo 1 vazio
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  # Invariantes do schema na coluna degenerada
  expect_length(pairs[[2]]$kb, 0L)
  expect_length(pairs[[2]]$Ra, length(pairs[[2]]$ka))
  expect_true(all(pairs[[2]]$Ra == 0))
  
  aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)   # nao deve dar erro
  # O objeto identity nao guarda pares; o piecewise e construido a parte
  pairs_pw <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                       beta, cas$hL, km, "piecewise")
  pw  <- QAST:::compile_fg_pw_cpp(pairs_pw)
  
  th <- c(0.2, -0.1)
  expect_identical(QAST:::eval_fg_parts_affine_cpp(aff, th)[2], 0.0)
  expect_identical(QAST:::eval_fg_parts_pw_cpp(pw, cas$lq, cas$z, th)[2], 0.0)
  expect_identical(QAST:::U2_parts_fg_cpp(th, cas$lt, cas$lq, cas$z,
                                          cas$status, beta, cas$hL,
                                          km, "identity")[2], 0.0)
})

test_that("Gate: empate exato em h1L preservado (afim)", {
  # Empate exato: beta = 0 => z'beta = 0 exato; tempos
  # duplicados exatamente em L; a desigualdade estrita do h1L replica
  # bit a bit porque o congelamento acontece NA CONSTRUCAO dos pares.
  # Falha por SALTO DISCRETO se qualquer reordenacao de floating point
  # mover um par de lado do gate.
  cas <- mk_fg_case(tie_at_L = TRUE)
  beta <- rep(0, ncol(cas$z))
  
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  aff <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
  
  theta <- c(0.1, 0.05)
  parts_leg <- QAST:::U2_parts_fg_cpp(theta, cas$lt, cas$lq,
                                      cas$z, cas$status, beta, cas$hL,
                                      km, "identity")
  parts_aff <- QAST:::eval_fg_parts_affine_cpp(aff, theta)
  expect_equal(parts_aff, parts_leg, tolerance = 1e-10)
})

test_that("Fast path: minimos quadrados reproduzem o otimizador (L2/identity)", {
  skip_if_not_installed("minqa")
  
  # FIX C: uobyqa com par inicial IGUAL A ZERO quebra no commonArgs do
  # minqa: rhobeg default = min(0.95, 0.2*max(abs(par))) colapsa a 0
  # => rhoend = 1e-6*rhobeg = 0 => stopifnot(0 < rhoend) falha.
  # Defesas duplas:
  #   (a) p = 2 — a documentacao do uobyqa pede length(par) > 1;
  #   (b) start NAO-nulo + rhobeg explicito: imune a start nulo.
  cas <- mk_fg_case(seed = 777, n = 200, p = 2)
  km <- QAST:::prep_km_by_column(cas$z, cas$t, cas$status)
  beta <- c(0.3, -0.2)
  stopifnot(ncol(cas$z) == 2L, length(beta) == 2L)
  
  pairs <- QAST:::make_fg_pairs_cpp(cas$lt, cas$lq, cas$z, cas$status,
                                    beta, cas$hL, km, "identity")
  aff   <- QAST:::make_fg_affine_cpp(pairs, cas$lq, cas$z)
  
  sol <- QAST:::solve_fg_theta_affine_cpp(aff)
  
  obj_fun <- function(th) QAST:::eval_fg_obj_affine_cpp(aff, th, "L2")
  opt <- minqa::uobyqa(c(0.2, -0.1), obj_fun,
                       control = list(rhobeg = 0.5, rhoend = 1e-8))
  
  # Nao-inferioridade: solucao fechada nao pode ser PIOR que o otimizador
  expect_lte(obj_fun(sol$par), obj_fun(opt$par) + 1e-8)
  
  # Gradiente ~ 0 no posto completo (semi_aft_fg exige < 1e-6)
  expect_lt(sol$grad_norm, 1e-6)
  expect_true(sol$full_rank)
})
