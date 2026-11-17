// src/fg_core.cpp
//
// Reference kernel of the U-statistic of Fine and Gelber (2001): a direct
// O(nk * nl) loop over the pairs of each binary column of z. Used in the
// A/B tests and benchmarks; semi_aft_fg uses the pairs object of
// fg_pairs_core.cpp, which must reproduce this kernel.
//
// For each column j, the U-statistic sums the pairs (a, b) with a in the
// group z[, j] == 0 and b in the group z[, j] == 1.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace fg_internal {

// phi: odd function of the theta estimator
inline double phi_apply(const double x, const std::string& phi) {
  if (phi == "identity") {
    return x;
  }
  // "piecewise": sign(x), with 0 -> 0
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

// Aggregation of the components by the L1 or L2 norm
inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// Components of the U-statistic, one per binary column of z
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
    
    // Groups defined by the discrete covariate of column j
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
        
        // All cross pairs (a in group 0, b in group 1): with a binary Z,
        // the Fine and Gelber sum over k < l of Z_kl U_kl is the sum over
        // these pairs, each one once.
        // Pairs with an event in both (deltak * deltal)
        if (dk(a) == 0.0 || dl(b) == 0.0) continue;
        
        // KM weight of each member; drops a null or non-finite denominator
        const double denom = wk(a) * wl(b);
        if (!(denom > 0.0) || !std::isfinite(denom)) continue;
        
        // Truncation indicator: h1(L) - max(z_k'beta, z_l'beta)
        //   greater than max(e1k, e1l)
        const double zk_beta = arma::dot(zk.row(a), beta);
        const double zl_beta = arma::dot(zl.row(b), beta);
        const double h1L =
          (hL - std::max(zk_beta, zl_beta) > std::max(e1k(a), e1l(b)))
          ? 1.0 : 0.0;
        
        const double contrib =
        h1L * phi_apply(e2k(a) - e2l(b), phi) / denom;
        
        // Factor (X_k - X_l) with X_k = 0 and X_l = 1
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
// For each binary covariate of z, computes the component of the
// U-statistic of pairs between the groups 0 and 1.
//
// Arguments:
//   theta  - quality-adjusted coefficients
//   lt     - failure times transformed by h1 (precomputed)
//   lq     - quality-adjusted times transformed by h2
//   z      - binary covariate matrix (n x p)
//   status - event indicator (1 = event)
//   beta   - survival coefficients (fixed in this stage)
//   hL     - constant h1(L), precomputed
//   km_surv- list of length p with w0 and w1 per column
//   phi    - odd function: "identity" or "piecewise"
//
// Returns a numeric vector of length p (without a dim attribute).
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
  
  // NumericVector: an exported arma::vec would reach R as a p x 1 matrix
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// U2_obj_fg_cpp
//
// Norm (L1 or L2) of the components of the U-statistic, standardized
// by sqrt(p).
//
// Arguments: the same as U2_parts_fg_cpp, plus:
//   norm_min - minimization norm ("L1" or "L2")
//
// Returns a double scalar with the value of the objective function.
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
