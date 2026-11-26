// src/surv_qt_core.cpp
//
// IPCW estimator of Zhao and Tsiatis (1997) for the survival curve of
// the quality-adjusted time (QALY). The documentation lives in the
// wrapper R/surv_qt_wrappers.R.
//
// pred_surv[i] = (1/n) sum_{j : qalys[j] > qalys[i]} status[j] / G(t_j),
// the survival P(Q > q) at the QALY of each observation, computed by an
// O(n log n) sort and an O(n) cumulative sum, with ties grouped. G is
// the censoring KM (survfit in R), evaluated at each time by binary
// search.
//
// Conventions:
// - time beyond the last KM time: G = NaN and the weight is taken as 0;
// - time before the first KM time: G = 1;
// - QALY <= 0 with a linear predictor: error (log undefined).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace survqt_internal {

// Evaluation of the KM step function at an arbitrary time.
// Returns: 1.0 before the first event time;
//          NaN after the last time (the weight becomes 0 in surv_qt_cpp);
//          km_surv[k] with km_time[k] <= t < km_time[k+1], otherwise.
inline double km_eval(const double t,
                      const NumericVector& km_time,
                      const NumericVector& km_surv)
{
  const int m = static_cast<int>(km_time.size());
  
  if (m == 0) return 1.0;
  if (t < km_time[0]) return 1.0;
  if (t > km_time[m - 1]) return std::numeric_limits<double>::quiet_NaN();
  
  // Binary search: largest k such that km_time[k] <= t
  int lo = 0, hi = m - 1;
  while (lo < hi) {
    const int mid = (lo + hi + 1) / 2;  // rounds up
    if (km_time[mid] <= t) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  return km_surv[lo];
}

} // namespace survqt_internal

// surv_qt_cpp
//
// Core of the survival estimator for QALY:
//
//   pred_surv[i] = (1/n) * sum_{j : qalys[j] > qalys[i]} (status[j] / cens_surv[j])
//
// Takes cens_surv ALREADY interpolated (use surv_qt_with_km_cpp to
// compute the interpolation from the survfit vectors).
//
// Arguments:
//   qaly       - vector of quality-adjusted times (on the raw scale)
//   time       - observed times (used only to check the length)
//   status     - event indicator (1 = event, 0 = censored)
//   cens_surv  - censoring survival evaluated at each observation
//   lin_pred   - optional: linear predictor (X %*% theta); if given,
//                qalys = log(qaly) - lin_pred
//
// Returns a numeric vector of length n. Internal use (called by
// surv_qt_with_km_cpp); not exported to R.
NumericVector surv_qt_cpp(const NumericVector& qaly,
                          const NumericVector& time,
                          const NumericVector& status,
                          const NumericVector& cens_surv,
                          const Nullable<NumericVector>& lin_pred = R_NilValue)
{
  const int n = static_cast<int>(qaly.size());
  
  if (static_cast<int>(time.size()) != n)
    stop("time must have the same length as qaly");
  if (static_cast<int>(status.size()) != n)
    stop("status must have the same length as qaly");
  if (static_cast<int>(cens_surv.size()) != n)
    stop("cens_surv must have the same length as qaly");
  
  // -----------------------------------------------------------------------
  // Transformation: qalys = qaly (raw) or log(qaly) - lin_pred
  // -----------------------------------------------------------------------
  NumericVector qalys(n);
  
  if (lin_pred.isNotNull()) {
    const NumericVector lp = lin_pred.get();
    if (static_cast<int>(lp.size()) != n)
      stop("lin_pred must have the same length as qaly");
    for (int i = 0; i < n; ++i) {
      if (qaly[i] <= 0.0) {
        stop("qaly[%d] <= 0 with lin_pred given: log undefined. "
               "Fix the data before calling.",
               i + 1);
      }
      qalys[i] = std::log(qaly[i]) - lp[i];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      qalys[i] = qaly[i];
    }
  }
  
  // -----------------------------------------------------------------------
  // Weights: status / cens_surv
  // NaN, zero or non-finite in cens_surv -> weight 0
  // -----------------------------------------------------------------------
  NumericVector weight(n);
  for (int i = 0; i < n; ++i) {
    const double cs = cens_surv[i];
    if (R_IsNaN(cs) || !(cs > 0.0) || !std::isfinite(cs)) {
      weight[i] = 0.0;
    } else {
      weight[i] = status[i] / cs;
    }
  }
  
  // -----------------------------------------------------------------------
  // Estimator by a decreasing sort + cumulative sweep
  //
  // Sorting the indices by decreasing qalys makes the cumulative sum,
  // when the tie group of qalys[i] is reached, hold exactly the weights
  // of the elements with qalys > qalys[i].
  // -----------------------------------------------------------------------
  IntegerVector idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;
  
  std::sort(idx.begin(), idx.end(),
            [&qalys](int a, int b) {
              if (qalys[a] != qalys[b]) return qalys[a] > qalys[b];
              return a < b;  // deterministic tie-break by index
            });
  
  NumericVector pred_surv(n);
  const double dn = static_cast<double>(n);
  
  int j = 0;
  double acc = 0.0;
  while (j < n) {
    const double current = qalys[idx[j]];
    const int start = j;

    // Tie group: all get the sum of the strictly greater ones
    // (comparison >), accumulated before adding the group itself
    while (j < n && qalys[idx[j]] == current) {
      pred_surv[idx[j]] = acc / dn;
      ++j;
    }
    for (int k = start; k < j; ++k) acc += weight[idx[k]];
  }
  
  return pred_surv;
}

// surv_qt_with_km_cpp
//
// Convenience version: interpolates the KM step function (vectors
// km_time/km_surv precomputed by survfit in R) at the observed times
// and delegates to the core surv_qt_cpp.
//
// Arguments:
//   qaly      - vector of quality-adjusted times (raw scale)
//   time      - observed times
//   status    - event indicator (1 = event, 0 = censored)
//   km_time   - times of the censoring KM survfit (sorted)
//   km_surv   - corresponding KM survival probabilities
//   lin_pred  - optional: linear predictor (X %*% theta)
//
// Returns a numeric vector of length n.
// [[Rcpp::export(rng = false)]]
NumericVector surv_qt_with_km_cpp(const NumericVector& qaly,
                                  const NumericVector& time,
                                  const NumericVector& status,
                                  const NumericVector& km_time,
                                  const NumericVector& km_surv,
                                  const Nullable<NumericVector>& lin_pred = R_NilValue)
{
  const int n = static_cast<int>(qaly.size());
  
  if (static_cast<int>(time.size()) != n)
    stop("time must have the same length as qaly");
  if (static_cast<int>(status.size()) != n)
    stop("status must have the same length as qaly");
  if (km_time.size() != km_surv.size())
    stop("km_time and km_surv must have the same length");
  
  // KM interpolation at each observed time
  NumericVector cens_surv(n);
  for (int i = 0; i < n; ++i) {
    cens_surv[i] = survqt_internal::km_eval(time[i], km_time, km_surv);
  }
  
  return surv_qt_cpp(qaly, time, status, cens_surv, lin_pred);
}
