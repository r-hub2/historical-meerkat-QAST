// src/aft_core.cpp
//
// Kernels das equacoes de estimacao do modelo proposto (bivariado
// semiparametrico AFT para tempo de sobrevivencia e tempo ajustado por
// qualidade de vida; tese de doutorado de da Silva, 2026, IME-USP). U1_obj_cpp
// tambem e a equacao de beta do modelo de Fine e Gelber (semi_aft_fg).
//
// Cada avaliacao: residuos e = y - z * par, ordenacao descendente
// ESTAVEL de e (empates na ordem original dos indices), um laco O(n) de
// somas sufixas e um gemv sobre z. Custo O(n log n + n p).
//
// Pesos IPCW (U2): w = status / HG_surv. Com HG_surv = 0 os pesos ficam
// NaN/Inf e o objetivo propaga NaN/Inf (falha visivel); os pesos devem
// ser finitos e >= 0.

#include <RcppArmadillo.h>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace aft_internal {
inline void check_dims(const arma::vec& par, const arma::vec& y,
                       const arma::mat& z, const arma::vec& status,
                       const char* who) {
  if (par.n_elem != z.n_cols)
    Rcpp::stop("%s: length(par) must equal ncol(z)", who);
  if (y.n_elem != z.n_rows || status.n_elem != z.n_rows)
    Rcpp::stop("%s: response/status must have nrow(z) elements", who);
}

// -----------------------------------------------------------------------------
// Ordenacao descendente ESTAVEL por radix sort (LSD, 8 passes de 8 bits):
// a mesma permutacao de arma::stable_sort_index(e, "descend") e de
// order(-e) no R (verificado em test-aft-order.R).
//
// Chave: o padrao de bits do double e mapeado para uint64 com ordem
// crescente igual a ordem numerica (positivos: liga o bit de sinal;
// negativos: inverte todos os bits) e depois invertido (descendente).
// -0.0 e normalizado para +0.0 (x + 0.0), pois a comparacao numerica os
// trata como iguais (empate) e os bits nao. Com NaN, delega ao Armadillo
// (que emite erro). Passes em que todas as chaves tem o mesmo digito sao
// pulados (comum nos bits altos do expoente).
// -----------------------------------------------------------------------------
inline arma::uvec stable_order_desc(const arma::vec& e) {
  const arma::uword n = e.n_elem;
  if (n < 2 || e.has_nan()) return arma::stable_sort_index(e, "descend");

  std::vector<std::uint64_t> key(n), key_tmp(n);
  std::vector<arma::uword>   idx(n), idx_tmp(n);
  for (arma::uword i = 0; i < n; ++i) {
    const double x = e[i] + 0.0;                 // -0.0 -> +0.0
    std::uint64_t u;
    std::memcpy(&u, &x, sizeof u);
    u = (u >> 63) ? ~u : (u | 0x8000000000000000ULL);
    key[i] = ~u;                                 // crescente em key == descendente em e
    idx[i] = i;
  }

  std::size_t count[256];
  for (int pass = 0; pass < 8; ++pass) {
    const int shift = 8 * pass;
    std::fill(count, count + 256, std::size_t(0));
    for (arma::uword i = 0; i < n; ++i) ++count[(key[i] >> shift) & 0xFF];
    if (count[(key[0] >> shift) & 0xFF] == n) continue;   // digito unico
    std::size_t pos = 0;
    for (int d = 0; d < 256; ++d) { const std::size_t c = count[d]; count[d] = pos; pos += c; }
    for (arma::uword i = 0; i < n; ++i) {                 // estavel
      const std::size_t dst = count[(key[i] >> shift) & 0xFF]++;
      key_tmp[dst] = key[i];
      idx_tmp[dst] = idx[i];
    }
    key.swap(key_tmp);
    idx.swap(idx_tmp);
  }
  return arma::uvec(idx);
}
}

// Exposto apenas para os testes da ordenacao
// [[Rcpp::export]]
Rcpp::IntegerVector aft_order_desc_cpp(const arma::vec& e,
                                       const std::string& method = "radix") {
  const arma::uvec o = (method == "arma")
    ? arma::stable_sort_index(e, "descend")
    : aft_internal::stable_order_desc(e);
  Rcpp::IntegerVector out(o.n_elem);
  for (arma::uword i = 0; i < o.n_elem; ++i) out[i] = static_cast<int>(o[i]) + 1;
  return out;
}

// -----------------------------------------------------------------------------
// Equacao de estimacao de beta, norma L1/L2 padronizada por sqrt(p).
//
// Com Zs = z ordenado por e decrescente, d = status ordenado,
// S0[r] = (r + 1)/n e C[r,] = sum_{k<=r} Zs[k,]:
//   U = sum_r d[r] * (S0[r] * Zs[r,] - C[r,] / n)
// Trocando a ordem da soma, C' d = sum_k Zs[k,] * R[k], R[k] = sum_{r>=k} d[r],
// logo U = z' v com v[ord[r]] = d[r] * S0[r] - R[r] / n.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double U1_obj_cpp(const arma::vec& par,
                  const arma::vec& lt,
                  const arma::mat& z,
                  const arma::vec& status,
                  const std::string& norm_min) {
  
  const int n = static_cast<int>(lt.n_elem);
  const int p = static_cast<int>(z.n_cols);
  aft_internal::check_dims(par, lt, z, status, "U1_obj_cpp");

  arma::vec e = lt - z * par;
  arma::uvec ord = aft_internal::stable_order_desc(e);

  // v[ord[r]] = d[r] * S0[r] - R[r] / n
  const double dn = static_cast<double>(n);
  arma::vec v(n);
  double tail = 0.0;
  for (int r = n - 1; r >= 0; --r) {
    const arma::uword i = ord[r];
    const double d = status[i];
    tail += d;
    v[i] = d * (static_cast<double>(r + 1) / dn) - tail / dn;
  }

  const arma::vec U = z.t() * v;
  
  if (norm_min == "L1")
    return accu(abs(U)) / std::sqrt(static_cast<double>(p));
  return std::sqrt(accu(square(U))) / std::sqrt(static_cast<double>(p));
}

// -----------------------------------------------------------------------------
// Equacao de estimacao ponderada (theta) do modelo proposto, norma L1/L2.
//
// Pesos IPCW w = status / HG_surv, ordenados (ws). Com cw[r] = sum_{k<=r} ws[k]
// e Cw[r,] = sum_{k<=r} ws[k] Zs[k,]:
//   U = sum_r ws[r] * ((cw[r] / n) * Zs[r,] - Cw[r,] / n) = z' v,
//   v[ord[r]] = ws[r] * cw[r] / n - ws[r] * Rw[r] / n,  Rw[r] = sum_{k>=r} ws[k]
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double U2_obj_cpp(const arma::vec& par,
                  const arma::vec& lq,
                  const arma::mat& z,
                  const arma::vec& status,
                  const arma::vec& HG_surv,
                  const std::string& norm_min) {
  
  const int n = static_cast<int>(lq.n_elem);
  const int p = static_cast<int>(z.n_cols);
  aft_internal::check_dims(par, lq, z, status, "U2_obj_cpp");
  if (HG_surv.n_elem != lq.n_elem)
    Rcpp::stop("U2_obj_cpp: HG_surv must have length n");

  arma::vec e = lq - z * par;
  arma::uvec ord = aft_internal::stable_order_desc(e);

  const arma::vec w  = status / HG_surv;
  const arma::vec ws = w.elem(ord);

  // Produtos separados (e nao ws * (cw - Rw)) para que pesos NaN/Inf
  // propaguem para o objetivo (ver cabecalho).
  const arma::vec cw = arma::cumsum(ws);
  const double dn = static_cast<double>(n);
  arma::vec v(n);
  double tail = 0.0;
  for (int r = n - 1; r >= 0; --r) {
    tail += ws[r];
    v[ord[r]] = ws[r] * (cw[r] / dn) - (ws[r] * tail) / dn;
  }

  const arma::vec U = z.t() * v;
  
  if (norm_min == "L1")
    return accu(abs(U)) / std::sqrt(static_cast<double>(p));
  return std::sqrt(accu(square(U))) / std::sqrt(static_cast<double>(p));
}
