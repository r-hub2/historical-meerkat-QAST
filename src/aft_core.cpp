// src/aft_core.cpp
//
// Kernels of the estimating equations of the proposed model (bivariate
// semiparametric AFT model for survival time and quality-adjusted
// survival time; doctoral thesis of da Silva, 2026, IME-USP). U1_obj_cpp
// is also the beta equation of the Fine and Gelber model (semi_aft_fg).
//
// Each evaluation: residuals e = y - z * par, STABLE descending sort of
// e (ties in the original index order), an O(n) loop of suffix sums and
// a gemv over z. Cost O(n log n + n p).
//
// IPCW weights (U2): w = status / HG_surv. With HG_surv = 0 the weights
// become NaN/Inf and the objective propagates NaN/Inf (visible failure);
// the weights must be finite and >= 0.

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
// STABLE descending sort by radix sort (LSD, 8 passes of 8 bits): the
// same permutation as arma::stable_sort_index(e, "descend") and as
// order(-e) in R (checked in test-aft-order.R).
//
// Key: the bit pattern of the double is mapped to a uint64 whose
// increasing order is the numeric order (positives: set the sign bit;
// negatives: flip all bits) and then inverted (descending). -0.0 is
// normalized to +0.0 (x + 0.0), since the numeric comparison treats them
// as equal (a tie) and the bits do not. With NaN, defers to Armadillo
// (which raises an error). Passes in which all keys share the same digit
// are skipped (common in the high bits of the exponent).
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
    key[i] = ~u;                                 // increasing in key == descending in e
    idx[i] = i;
  }

  std::size_t count[256];
  for (int pass = 0; pass < 8; ++pass) {
    const int shift = 8 * pass;
    std::fill(count, count + 256, std::size_t(0));
    for (arma::uword i = 0; i < n; ++i) ++count[(key[i] >> shift) & 0xFF];
    if (count[(key[0] >> shift) & 0xFF] == n) continue;   // single digit
    std::size_t pos = 0;
    for (int d = 0; d < 256; ++d) { const std::size_t c = count[d]; count[d] = pos; pos += c; }
    for (arma::uword i = 0; i < n; ++i) {                 // stable
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

// Exported only for the tests of the sort
// [[Rcpp::export(rng = false)]]
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
// Estimating equation of beta, L1/L2 norm standardized by sqrt(p).
//
// With Zs = z sorted by decreasing e, d = sorted status,
// S0[r] = (r + 1)/n e C[r,] = sum_{k<=r} Zs[k,]:
//   U = sum_r d[r] * (S0[r] * Zs[r,] - C[r,] / n)
// Swapping the order of the sum, C' d = sum_k Zs[k,] * R[k], R[k] = sum_{r>=k} d[r],
// so U = z' v with v[ord[r]] = d[r] * S0[r] - R[r] / n.
// -----------------------------------------------------------------------------
// [[Rcpp::export(rng = false)]]
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
// Weighted estimating equation (theta) of the proposed model, L1/L2 norm.
//
// IPCW weights w = status / HG_surv, sorted (ws). With cw[r] = sum_{k<=r} ws[k]
// e Cw[r,] = sum_{k<=r} ws[k] Zs[k,]:
//   U = sum_r ws[r] * ((cw[r] / n) * Zs[r,] - Cw[r,] / n) = z' v,
//   v[ord[r]] = ws[r] * cw[r] / n - ws[r] * Rw[r] / n,  Rw[r] = sum_{k>=r} ws[k]
// -----------------------------------------------------------------------------
// [[Rcpp::export(rng = false)]]
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

  // Separate products (and not ws * (cw - Rw)) so that NaN/Inf weights
  // propagate to the objective (see the header).
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
