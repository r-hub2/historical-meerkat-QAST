// src/fg_core.cpp
//
// Kernel de referencia da U-estatistica de Fine e Gelber (2001): laco
// direto O(nk * nl) sobre os pares de cada coluna binaria de z. Usado nos
// testes A/B e benchmarks; semi_aft_fg usa o objeto de pares de
// fg_pairs_core.cpp, que deve reproduzir este kernel.
//
// Para cada coluna j, a U-estatistica soma pares (a, b) com a no grupo
// z[, j] == 0 e b no grupo z[, j] == 1.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace fg_internal {

// phi: funcao impar do estimador de theta
inline double phi_apply(const double x, const std::string& phi) {
  if (phi == "identity") {
    return x;
  }
  // "piecewise": sign(x), com 0 -> 0
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

// Agregacao das componentes pela norma L1 ou L2
inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// Componentes da U-estatistica, uma por coluna binaria de z
arma::vec u2_parts(const arma::vec& theta,
                   const arma::vec& lt,
                   const arma::vec& lq,
                   const arma::mat& z,
                   const arma::vec& status,
                   const arma::vec& beta,
                   const double hL,
                   const Rcpp::List& km_surv,
                   const std::string& phi)
{
  const int n = static_cast<int>(z.n_rows);
  const int p = static_cast<int>(z.n_cols);
  const double n2 = static_cast<double>(n) * static_cast<double>(n);
  
  const arma::vec e1 = lt - z * beta;
  const arma::vec e2 = lq - z * theta;
  
  arma::vec parts(p);
  
  for (int j = 0; j < p; ++j) {
    
    // Grupos definidos pela covariavel discreta da coluna j
    const arma::uvec ik = arma::find(z.col(j) == 0.0);
    const arma::uvec il = arma::find(z.col(j) == 1.0);
    const int nk = static_cast<int>(ik.n_elem);
    const int nl = static_cast<int>(il.n_elem);
    
    if (nk == 0 || nl == 0) {
      parts(j) = 0.0;
      continue;
    }
    
    const arma::vec e1k = e1.elem(ik);
    const arma::vec e1l = e1.elem(il);
    const arma::vec e2k = e2.elem(ik);
    const arma::vec e2l = e2.elem(il);
    const arma::vec dk  = status.elem(ik);
    const arma::vec dl  = status.elem(il);
    
    const arma::mat zk = z.rows(ik);
    const arma::mat zl = z.rows(il);
    
    const Rcpp::List kmj = km_surv[j];
    const arma::vec wk = Rcpp::as<arma::vec>(kmj["w0"]);
    const arma::vec wl = Rcpp::as<arma::vec>(kmj["w1"]);
    
    double acc = 0.0;
    
    for (int a = 0; a < nk; ++a) {
      for (int b = 0; b < nl; ++b) {
        
        // Emparelhamento posicional: a < b (indices dentro dos grupos)
        if (a >= b) continue;
        
        // Pares com evento em ambos (deltak * deltal)
        if (dk(a) == 0.0 || dl(b) == 0.0) continue;
        
        // Peso KM por grupo; descarta denominador nulo ou nao finito
        const double denom = wk(a) * wl(b);
        if (!(denom > 0.0) || !std::isfinite(denom)) continue;
        
        // Indicador de truncacao: h1(L) - max(z_k'beta, z_l'beta)
        //   maior que max(e1k, e1l)
        const double zk_beta = arma::dot(zk.row(a), beta);
        const double zl_beta = arma::dot(zl.row(b), beta);
        const double h1L =
          (hL - std::max(zk_beta, zl_beta) > std::max(e1k(a), e1l(b)))
          ? 1.0 : 0.0;
        
        const double contrib =
        h1L * phi_apply(e2k(a) - e2l(b), phi) / denom;
        
        // Fator (X_k - X_l) com X_k = 0 e X_l = 1
        acc += -contrib;
      }
    }
    parts(j) = acc / n2;
  }
  return parts;
}
}

// U2_parts_fg_cpp
//
// Para cada covariavel binaria de z, calcula a componente da
// U-estatistica de pares entre os grupos 0 e 1.
//
// Argumentos:
//   theta  - coeficientes de qualidade-ajustada
//   lt     - tempos de falha transformados por h1 (pre-computados)
//   lq     - tempos ajustados por qualidade transformados por h2
//   z      - matriz de covariaveis binarias (n x p)
//   status - vetor indicador de evento (1 = evento)
//   beta   - coeficientes de sobrevivencia (fixos nesta etapa)
//   hL     - constante h1(L), pre-computada
//   km_surv- lista de comprimento p com w0 e w1 por coluna
//   phi    - funcao impar: "identity" ou "piecewise"
//
// Retorna vetor numerico de comprimento p (sem atributo dim).
// [[Rcpp::export]]
Rcpp::NumericVector U2_parts_fg_cpp(const arma::vec& theta,
                                    const arma::vec& lt,
                                    const arma::vec& lq,
                                    const arma::mat& z,
                                    const arma::vec& status,
                                    const arma::vec& beta,
                                    const double hL,
                                    const Rcpp::List& km_surv,
                                    const std::string& phi)
{
  const arma::vec parts = fg_internal::u2_parts(theta, lt, lq, z, status,
                                                beta, hL, km_surv, phi);
  
  // NumericVector: um arma::vec exportado chegaria ao R como matriz p x 1
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// U2_obj_fg_cpp
//
// Norma (L1 ou L2) das componentes da U-estatistica, padronizada
// por sqrt(p).
//
// Argumentos: os mesmos de U2_parts_fg_cpp, mais:
//   norm_min - norma de minimizacao ("L1" ou "L2")
//
// Retorna escalar double com o valor da funcao objetivo.
// [[Rcpp::export]]
double U2_obj_fg_cpp(const arma::vec& theta,
                     const arma::vec& lt,
                     const arma::vec& lq,
                     const arma::mat& z,
                     const arma::vec& status,
                     const arma::vec& beta,
                     const double hL,
                     const Rcpp::List& km_surv,
                     const std::string& phi,
                     const std::string& norm_min)
{
  const arma::vec parts = fg_internal::u2_parts(theta, lt, lq, z, status,
                                                beta, hL, km_surv, phi);
  return fg_internal::norm_aggregate(parts, norm_min);
}
