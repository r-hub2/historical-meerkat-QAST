# tests/testthat/test-equivalence.R
#
# Golden reference: as funcoes U1_fun / U2_fun ORIGINAIS do script QAST,
# copiadas sem modificacao (exceto tornar 'norm_min' parametro explicito,
# que no original vinha do escopo envolvente via closure).

# ---- Referencia: etapa beta ------------------------------------------------

Slfun_ref <- function(z, e1, l) {
  sapply(X = e1, u = e1, z = z, l = l,
         FUN = function(X, u, z, l) mean((u >= X) * z^l))
}

U1_fun_ref <- function(par, t, z, status, h1, norm_min) {
  beta <- par
  e1 <- as.vector(h1$linkfun(t) - (z %*% beta))
  
  S0 <- Slfun_ref(e1 = e1, z = 1, l = 0)
  S1 <- apply(X = z, e1 = e1, l = 1, MARGIN = 2, FUN = Slfun_ref)
  
  U1 <- colSums(x = (S0 * (z - S1 / S0) * status))
  if (norm_min == "L1") outU1 <- sum(abs(U1)) / sqrt(length(par))
  if (norm_min == "L2") outU1 <- sqrt(sum(U1^2)) / sqrt(length(par))
  outU1
}

# ---- Referencia: etapa theta (ponderada por IPCW) --------------------------

Slfun_w_ref <- function(z, e2, l, status, HG_surv) {
  if (l == 0) {
    Sl <- sapply(1:nrow(z), function(i) mean((e2 >= e2[i]) * status / HG_surv))
  } else {
    Sl <- t(sapply(1:nrow(z), function(i)
      colMeans((e2 >= e2[i]) * z^l * status / HG_surv)))
  }
  Sl
}

U2_fun_ref <- function(par, qt, z, status, h2, HG_surv, norm_min) {
  theta <- par
  e2 <- as.vector(h2$linkfun(qt) - (z %*% theta))
  
  S0 <- Slfun_w_ref(z = z, e2 = e2, l = 0, status = status, HG_surv = HG_surv)
  S1 <- Slfun_w_ref(z = z, e2 = e2, l = 1, status = status, HG_surv = HG_surv)
  
  mu <- replace(x = S1 / S0, list = (S1 == 0 & S0 == 0), values = 0)
  if (length(par) == 1) mu <- t(mu)
  
  U2 <- colSums(x = (S0 * (z - mu) * status / HG_surv))
  if (norm_min == "L1") outU2 <- sum(abs(U2)) / sqrt(length(par))
  if (norm_min == "L2") outU2 <- sqrt(sum(U2^2)) / sqrt(length(par))
  outU2
}

# ---- Dados sinteticos compartilhados --------------------------------------

make_data <- function(n = 150, p = 3, seed = 20261021) {
  set.seed(seed)
  z <- matrix(rnorm(n * p), n, p)
  storage.mode(z) <- "double"
  colnames(z) <- paste0("x", seq_len(p))
  t  <- exp(as.vector(z %*% seq(0.4, -0.2, length.out = p)) + rnorm(n))
  qt <- t * runif(n, 0.3, 1.0)
  list(z = z, t = t, qt = qt, status = rbinom(n, 1, 0.85))
}

# ---- Testes ---------------------------------------------------------------

test_that("U1_obj_cpp reproduz U1_fun original (L1 e L2, varios beta)", {
  d  <- make_data()
  h1 <- make.link("log")
  lt <- h1$linkfun(d$t)
  
  betas <- list(rep(0, 3), c(0.4, -0.3, 0.2), c(-1, 0.5, 1.2))
  
  for (nm in c("L1", "L2")) {
    for (b in betas) {
      ref <- U1_fun_ref(b, t = d$t, z = d$z, status = d$status,
                        h1 = h1, norm_min = nm)
      got <- U1_obj_cpp(b, lt, d$z, d$status, nm)
      expect_equal(got, ref, tolerance = 1e-9,
                   label = sprintf("U1/%s beta=(%s)", nm, paste(b, collapse = ",")))
    }
  }
})

test_that("U2_obj_cpp reproduz U2_fun original (ponderado, L1 e L2)", {
  d  <- make_data(seed = 999)
  h2 <- make.link("log")
  lq <- h2$linkfun(d$qt)
  
  set.seed(7)
  HG_surv <- pmax(runif(nrow(d$z), 0.05, 1.5), 0.05)  # pesos IPCW positivos
  
  thetas <- list(rep(0, 3), c(0.3, -0.5, 0.8), c(1.1, 0.2, -0.7))
  
  for (nm in c("L1", "L2")) {
    for (th in thetas) {
      ref <- U2_fun_ref(th, qt = d$qt, z = d$z, status = d$status,
                        h2 = h2, HG_surv = HG_surv, norm_min = nm)
      got <- U2_obj_cpp(th, lq, d$z, d$status, HG_surv, nm)
      expect_equal(got, ref, tolerance = 1e-9,
                   label = sprintf("U2/%s theta=(%s)", nm, paste(th, collapse = ",")))
    }
  }
})

test_that("caso p = 1 nao quebra dimensoes", {
  d  <- make_data(n = 80, p = 1, seed = 4242)
  h1 <- make.link("log"); h2 <- make.link("log")
  
  ref1 <- U1_fun_ref(0.5, t = d$t, z = d$z, status = d$status,
                     h1 = h1, norm_min = "L2")
  got1 <- U1_obj_cpp(0.5, h1$linkfun(d$t), d$z, d$status, "L2")
  expect_equal(got1, ref1, tolerance = 1e-9)
  
  HG <- pmax(runif(80, 0.2, 1.0), 0.2)
  ref2 <- U2_fun_ref(0.5, qt = d$qt, z = d$z, status = d$status,
                     h2 = h2, HG_surv = HG, norm_min = "L2")
  got2 <- U2_obj_cpp(0.5, h2$linkfun(d$qt), d$z, d$status, HG, "L2")
  expect_equal(got2, ref2, tolerance = 1e-9)
})

test_that("casos degenerados (status todo zero, pesos nulos) retornam 0", {
  d  <- make_data(n = 60, seed = 31337)
  st0 <- rep(0, nrow(d$z))
  
  expect_equal(U1_obj_cpp(rep(0.2, 3), log(d$t), d$z, st0, "L2"),
               U1_fun_ref(rep(0.2, 3), d$t, d$z, st0, make.link("log"), "L2"),
               tolerance = 1e-12)
  
  HG <- rep(1, nrow(d$z))
  expect_equal(U2_obj_cpp(rep(0.2, 3), log(d$qt), d$z, st0, HG, "L2"),
               U2_fun_ref(rep(0.2, 3), d$qt, d$z, st0, make.link("log"), HG, "L2"),
               tolerance = 1e-12)
})