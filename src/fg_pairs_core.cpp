// src/fg_pairs_core.cpp
//
// Precomputed pairs object of the U-statistic of Fine and Gelber. The
// constants of the theta stage (truncation indicator h1L, KM weights,
// event filters) are resolved once per fit; evaluating the objective
// costs O(nk + nl) (identity, via the marginals Ra/Cb) or O(nnz)
// (piecewise). The reference kernel is the loop of fg_core.cpp, which
// this object reproduces (A/B tests).
//
// Valid pair (a, b), with a in group 0 and b in group 1 of column j:
//   - all cross pairs between the groups (with a binary Z, the Fine and
//     Gelber sum over k < l of Z_kl U_kl is the sum over these pairs,
//     each pair once; it does not depend on the order of the rows);
//   - an event in both members;
//   - denom = wk[a] * wl[b] finite and > 0;
//   - h1L: hL - max(zbk[a], zbl[b]) > max(e1k[a], e1l[b]) (strict).
// The weight of the pair is c_ab = 1 / denom; the division by n^2 and
// the sign are applied in the evaluation: parts(j) = -(sum) / n^2.
//
// With phi = "identity", Ra/Cb are computed in O(N log N)
// (marginals_cdq, below) and the pairs are not stored; with "piecewise",
// by the O(nk * nl) loop, which also stores pa/pb/pc.
//
// FIDELITY: zbk/zbl use arma::dot per member, the same expression as the
// reference kernel, so that the strict comparison of h1L is bit for bit
// the same (the exact-tie test depends on it). Do not replace it by
// z * beta: the matrix product reorders sums and may flip ties.
//
// OUTPUT SCHEMA:
//   INVARIANTS (validated by the consumers):
//     len(Ra) == len(ka), len(Cb) == len(kb)   (also in degenerate
//                                                columns: zeros)
//     len(pa) == len(pb) == len(pc)
//   pairs[[j]] : Rcpp::List with
//     "ka" - global integer (0-based) of the rows of group 0 (original order)
//     "kb" - global integer of the rows of group 1 (original order)
//     "Ra" - vec, marginal sum_b c_ab per member of group 0 (aligned with ka)
//     "Cb" - vec, marginal sum_a c_ab per member of group 1 (aligned with kb)
//     "pa" - global integer of the row of member 0 of each valid pair
//     "pb" - global integer of the row of member 1 of each valid pair
//     "pc" - vec, weight c_ab per valid pair
//     (pa/pb/pc empty when has_pairs = FALSE)
//   attr(pairs, "meta") : List with
//     "p"         - int, number of columns of z
//     "n2"        - double, n^2
//     "has_pairs" - logical, TRUE if pa/pb/pc were built
//     "cdq_cols"  - logical per column, TRUE if Ra/Cb came from the
//                   O(N log N) sweep (marginals_cdq; the name is historical)

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <algorithm>
#include <functional>
#include <vector>

using namespace Rcpp;
using namespace arma;

namespace fgp {

// ---------------------------------------------------------------------------
// Meta helpers (attr("meta") returns a proxy: materialize it as a List)
// ---------------------------------------------------------------------------

inline Rcpp::List get_meta(const Rcpp::List& obj) {
  if (!obj.hasAttribute("meta"))
    stop("fg_pairs: object without the 'meta' attribute; use make_fg_pairs_cpp");
  const Rcpp::List meta = obj.attr("meta");   // materialize the proxy ONCE
  return meta;
}

inline int meta_p(const Rcpp::List& meta) {
  if (!meta.containsElementNamed("p"))
    stop("fg_pairs: meta has no 'p' field");
  return Rcpp::as<int>(meta["p"]);
}

inline double meta_n2(const Rcpp::List& meta) {
  if (!meta.containsElementNamed("n2"))
    stop("fg_pairs: meta has no 'n2' field");
  return Rcpp::as<double>(meta["n2"]);
}

// ---------------------------------------------------------------------------
// piecewise phi: sign(x) with 0 -> 0
// ---------------------------------------------------------------------------

inline double phi_pw(const double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

inline void validate_phi(const std::string& phi) {
  if (phi != "identity" && phi != "piecewise")
    stop("fg_pairs: invalid phi (use identity or piecewise): %s",
         phi.c_str());
}

// ---------------------------------------------------------------------------
// Aggregation by the L1/L2 norm, standardized by sqrt(p)
// ---------------------------------------------------------------------------

inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// 0-based indices as R integers (read by pointer in the evaluation)
inline Rcpp::IntegerVector to_int_idx(const arma::uvec& v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

// Object without the has_pairs field: treated as holding pairs
inline bool meta_has_pairs(const Rcpp::List& meta) {
  if (!meta.containsElementNamed("has_pairs")) return true;
  return Rcpp::as<bool>(meta["has_pairs"]);
}

// ===========================================================================
// Marginals Ra/Cb in O(N log N) (phi = "identity")
//
// Valid pair (a, b) in the O(nk*nl) loop:
//   dk[a] != 0, dl[b] != 0,
//   denom = wk[a]*wl[b] > 0 and finite,
//   hL - max(zbk[a], zbl[b]) > max(e1k[a], e1l[b]).
//
// EXACT FLOATING-POINT IDENTITY: fl(hL - .) is monotone non-increasing
// (monotone rounding), hence
//   fl(hL - max(x, y)) == min(fl(hL - x), fl(hL - y)).
// With A_a = hL - zbk[a], B_b = hL - zbl[b], u_a = e1k[a], v_b = e1l[b]:
//   min(A_a, B_b) > max(u_a, v_b)
//   <=> A_a > u_a  and  B_b > v_b        (per-element filters)
//       and A_a > v_b  and  B_b > u_a    (cross conditions)
// All comparisons are on the SAME doubles as the loop: the set of pairs
// is IDENTICAL. Only the order of the sums changes (positive terms, no
// cancellation; ~1e-15 relative).
//
// Weight: c_ab = 1/(wk wl) ~ (1/wk)(1/wl), hence
//   Ra[a] = (1/wk[a]) * sum_{b: v_b < A_a, B_b > u_a} 1/wl[b]
//   Cb[b] = (1/wl[b]) * sum_{a: A_a > v_b, u_a < B_b} 1/wk[a]
// — weighted 2D dominance counts, solved by one sweep with a Fenwick
// tree (marg_cross over the two whole groups).
//
// FALLBACK to the loop (returns false) when the identity or the
// factorization is not guaranteed: non-finite values (std::max with NaN
// is asymmetric) or positive weights outside [1e-150, 1e150] (the
// product could leave the normal range and change the denom > 0 / finite filter).
// ===========================================================================

struct Fenwick {
  std::vector<double> t;
  explicit Fenwick(size_t n) : t(n + 1, 0.0) {}
  void add(size_t i, double v) {
    for (++i; i < t.size(); i += i & (~i + 1)) t[i] += v;
  }
  double sum(size_t k) const {          // sum of the ranks [0, k)
    double s = 0.0;
    for (; k > 0; k -= k & (~k + 1)) s += t[k];
    return s;
  }
};

struct MargData {
  // ACTIVE members (passed the per-element filters), in position order
  std::vector<int>    pa, pb;      // position within the group
  std::vector<double> A, u, wa;    // group 0: A_a, u_a, 1/wk[a]
  std::vector<double> B, v, wb;    // group 1: B_b, v_b, 1/wl[b]
  std::vector<double> SRa, SCb;    // accumulated sums
};

// Contributions of all pairs a in [a0, a1) x b in [b0, b1)
inline void marg_cross(MargData& d, int a0, int a1, int b0, int b1) {
  const int na = a1 - a0, nb = b1 - b0;
  if (na == 0 || nb == 0) return;
  std::vector<int> qa(na), pbv(nb);
  for (int i = 0; i < na; ++i) qa[i] = a0 + i;
  for (int k = 0; k < nb; ++k) pbv[k] = b0 + k;

  // --- Ra: points b (v asc), queries a (A asc); sum of wb with B > u
  {
    std::vector<double> Bd(nb);
    for (int k = 0; k < nb; ++k) Bd[k] = d.B[b0 + k];
    std::sort(Bd.begin(), Bd.end(), std::greater<double>());
    Bd.erase(std::unique(Bd.begin(), Bd.end()), Bd.end());
    std::sort(pbv.begin(), pbv.end(),
              [&d](int x, int y) { return d.v[x] < d.v[y]; });
    std::sort(qa.begin(), qa.end(),
              [&d](int x, int y) { return d.A[x] < d.A[y]; });
    Fenwick fw(Bd.size());
    int p = 0;
    for (int q : qa) {
      while (p < nb && d.v[pbv[p]] < d.A[q]) {          // v_b < A_a
        const int b = pbv[p++];
        const size_t r = std::lower_bound(Bd.begin(), Bd.end(), d.B[b],
                                          std::greater<double>()) - Bd.begin();
        fw.add(r, d.wb[b]);
      }
      // ranks with B > u_a (Bd decreasing)
      const size_t k = std::lower_bound(Bd.begin(), Bd.end(), d.u[q],
                                        std::greater<double>()) - Bd.begin();
      d.SRa[q] += fw.sum(k);
    }
  }

  // --- Cb: points a (A desc), queries b (v desc); sum of wa with u < B
  {
    std::vector<double> Ua(na);
    for (int i = 0; i < na; ++i) Ua[i] = d.u[a0 + i];
    std::sort(Ua.begin(), Ua.end());
    Ua.erase(std::unique(Ua.begin(), Ua.end()), Ua.end());
    std::sort(qa.begin(), qa.end(),
              [&d](int x, int y) { return d.A[x] > d.A[y]; });
    std::sort(pbv.begin(), pbv.end(),
              [&d](int x, int y) { return d.v[x] > d.v[y]; });
    Fenwick fw(Ua.size());
    int p = 0;
    for (int q : pbv) {
      while (p < na && d.A[qa[p]] > d.v[q]) {           // A_a > v_b
        const int a = qa[p++];
        const size_t r = std::lower_bound(Ua.begin(), Ua.end(), d.u[a])
                         - Ua.begin();
        fw.add(r, d.wa[a]);
      }
      // ranks with u < B_b (Ua increasing)
      const size_t k = std::lower_bound(Ua.begin(), Ua.end(), d.B[q])
                       - Ua.begin();
      d.SCb[q] += fw.sum(k);
    }
  }
}

// Fills Ra/Cb; false => use the O(nk*nl) loop (see FALLBACK above)
inline bool marginals_cdq(const double hL,
                          const arma::vec& zbk, const arma::vec& zbl,
                          const arma::vec& e1k, const arma::vec& e1l,
                          const arma::vec& dk,  const arma::vec& dl,
                          const arma::vec& wk,  const arma::vec& wl,
                          arma::vec& Ra, arma::vec& Cb) {
  if (!std::isfinite(hL) || !zbk.is_finite() || !zbl.is_finite() ||
      !e1k.is_finite() || !e1l.is_finite())
    return false;

  const double wmin = 1e-150, wmax = 1e150;
  auto weight_ok = [&](double w, bool& in_range) {
    // same filter as denom > 0 && finite, per element (w in (0, inf)).
    // w < 0: two negative weights would give denom > 0 in the loop — the
    // per-element factorization does not reproduce that, so it forces the fallback.
    if (w < 0.0) { in_range = false; return false; }
    if (!(w > 0.0) || !std::isfinite(w)) return false;
    if (w < wmin || w > wmax) in_range = false;
    return true;
  };

  const int nk = static_cast<int>(zbk.n_elem);
  const int nl = static_cast<int>(zbl.n_elem);
  MargData d;
  bool in_range = true;
  for (int a = 0; a < nk; ++a) {
    if (dk[a] == 0.0 || !weight_ok(wk[a], in_range)) continue;
    const double Aa = hL - zbk[a];
    if (!(Aa > e1k[a])) continue;                       // A_a > u_a
    d.pa.push_back(a); d.A.push_back(Aa); d.u.push_back(e1k[a]);
    d.wa.push_back(1.0 / wk[a]);
  }
  for (int b = 0; b < nl; ++b) {
    if (dl[b] == 0.0 || !weight_ok(wl[b], in_range)) continue;
    const double Bb = hL - zbl[b];
    if (!(Bb > e1l[b])) continue;                       // B_b > v_b
    d.pb.push_back(b); d.B.push_back(Bb); d.v.push_back(e1l[b]);
    d.wb.push_back(1.0 / wl[b]);
  }
  if (!in_range) return false;

  d.SRa.assign(d.pa.size(), 0.0);
  d.SCb.assign(d.pb.size(), 0.0);
  marg_cross(d, 0, static_cast<int>(d.pa.size()),
             0, static_cast<int>(d.pb.size()));

  for (size_t i = 0; i < d.pa.size(); ++i) Ra[d.pa[i]] = d.wa[i] * d.SRa[i];
  for (size_t k = 0; k < d.pb.size(); ++k) Cb[d.pb[k]] = d.wb[k] * d.SCb[k];
  return true;
}

} // namespace fgp

// ===========================================================================
// CONSTRUCTOR OF THE PAIRS OBJECT
// ===========================================================================

// make_fg_pairs_cpp
//
// Builds the object of precomputed pairs per binary column of z.
// The filters of a valid pair (header) are applied here, once.
//
// Arguments:
//   lt      - h1$linkfun(t), precomputed
//   lq      - h2$linkfun(qt), precomputed
//   z       - binary covariate matrix (n x p)
//   status  - event indicator (1 = event)
//   beta    - survival coefficients (fixed in this stage)
//   hL      - constant h1(L), precomputed
//   km_surv - list of length p with w0/w1 per column (original order)
//   phi     - "identity" or "piecewise". With "identity" only Ra/Cb are
//             built; "piecewise" also builds pa/pb/pc
//             (the piecewise object serves both phi in the evaluation)
//   algo    - "auto" (default): identity uses marginals_cdq, O(N log N),
//             with a per-column fallback to the loop; "loop" forces the
//             O(nk*nl) loop (reference for the A/B tests). piecewise always
//             uses the loop (it needs the explicit pairs).
//
// Returns: list per column + attr "meta" (see the schema in the header);
// meta$cdq_cols tells, per column, whether the CDQ path was used.
// [[Rcpp::export(rng = false)]]
Rcpp::List make_fg_pairs_cpp(const arma::vec& lt,
                             const arma::vec& lq,
                             const arma::mat& z,
                             const arma::vec& status,
                             const arma::vec& beta,
                             const double hL,
                             const Rcpp::List& km_surv,
                             const std::string& phi,
                             const std::string& algo = "auto") {

  fgp::validate_phi(phi);
  if (algo != "auto" && algo != "loop")
    stop("make_fg_pairs: invalid algo (use auto or loop): %s",
         algo.c_str());
  const bool store_pairs = (phi == "piecewise");
  const bool try_cdq = !store_pairs && algo == "auto";

  const int n = static_cast<int>(z.n_rows);
  const int p = static_cast<int>(z.n_cols);
  const double n2 = static_cast<double>(n) * static_cast<double>(n);
  
  if (static_cast<int>(lt.n_elem) != n ||
      static_cast<int>(lq.n_elem) != n ||
      static_cast<int>(status.n_elem) != n)
    stop("make_fg_pairs: lt/lq/status must have length n");
  if (static_cast<int>(beta.n_elem) != p)
    stop("make_fg_pairs: beta must have length p = ncol(z)");
  if (static_cast<int>(km_surv.size()) != p)
    stop("make_fg_pairs: km_surv must have length p = ncol(z)");
  
  // Residuals e1 = lt - z * beta (same expression as the reference kernel)
  const arma::vec e1 = lt - z * beta;
  
  Rcpp::List out(p);
  Rcpp::LogicalVector cdq_cols(p, false);

  for (int j = 0; j < p; ++j) {
    
    // Groups defined by the binary covariate of column j
    const arma::uvec ik = arma::find(z.col(j) == 0.0);
    const arma::uvec il = arma::find(z.col(j) == 1.0);
    const int nk = static_cast<int>(ik.n_elem);
    const int nl = static_cast<int>(il.n_elem);
    
    // Degenerate column: no pairs, but with the schema invariants
    // (Ra/Cb of zeros aligned with ka/kb); parts(j) = 0. Use (n, fill::zeros):
    // arma::vec(0) is ambiguous (it would match Col(const char*)).
    if (nk == 0 || nl == 0) {
      out[j] = Rcpp::List::create(
        Rcpp::Named("ka") = fgp::to_int_idx(ik),
        Rcpp::Named("kb") = fgp::to_int_idx(il),
        Rcpp::Named("Ra") = arma::vec(static_cast<arma::uword>(nk),
                    arma::fill::zeros),
        Rcpp::Named("Cb") = arma::vec(static_cast<arma::uword>(nl),
                    arma::fill::zeros),
        Rcpp::Named("pa") = Rcpp::IntegerVector(0),
        Rcpp::Named("pb") = Rcpp::IntegerVector(0),
        Rcpp::Named("pc") = Rcpp::NumericVector(0)
      );
      continue;
    }
    
    const arma::vec e1k = e1.elem(ik);
    const arma::vec e1l = e1.elem(il);
    const arma::vec dk  = status.elem(ik);
    const arma::vec dl  = status.elem(il);
    
    const arma::mat zk = z.rows(ik);
    const arma::mat zl = z.rows(il);
    
    // dot per member, once (see FIDELITY in the header)
    arma::vec zbk(nk), zbl(nl);
    for (int a = 0; a < nk; ++a) zbk[a] = arma::dot(zk.row(a), beta);
    for (int b = 0; b < nl; ++b) zbl[b] = arma::dot(zl.row(b), beta);
    
    // KM weights of the group (original order — aligned with ik/il)
    const Rcpp::List kmj = km_surv[j];
    if (!kmj.containsElementNamed("w0") ||
        !kmj.containsElementNamed("w1"))
        stop("make_fg_pairs: km_surv[[%d]] has no w0/w1 fields", j + 1);
    const arma::vec wk = Rcpp::as<arma::vec>(kmj["w0"]);
    const arma::vec wl = Rcpp::as<arma::vec>(kmj["w1"]);
    if (static_cast<int>(wk.n_elem) != nk ||
        static_cast<int>(wl.n_elem) != nl)
      stop("make_fg_pairs: km_surv[[%d]] w0/w1 do not match the "
             "groups (nk = %d, nl = %d)", j + 1, nk, nl);
    
    // Output accumulators
    arma::vec  Ra(nk, arma::fill::zeros);
    arma::vec  Cb(nl, arma::fill::zeros);
    std::vector<int>    pa_v, pb_v;
    std::vector<double> pc_v;
    if (store_pairs) {
      pa_v.reserve(static_cast<size_t>(nk) * 4);   // conservative guess
      pb_v.reserve(static_cast<size_t>(nk) * 4);
      pc_v.reserve(static_cast<size_t>(nk) * 4);
    }

    // Marginals in O(N log N); false => loop below
    const bool used_cdq = try_cdq &&
      fgp::marginals_cdq(hL, zbk, zbl, e1k, e1l, dk, dl, wk, wl, Ra, Cb);
    cdq_cols[j] = used_cdq;

    // Double loop over all cross pairs; a censored member a is skipped,
    // since none of its pairs is valid.
    for (int a = 0; a < nk && !used_cdq; ++a) {
      if (dk[a] == 0.0) continue;

      for (int b = 0; b < nl; ++b) {
        if (dl[b] == 0.0) continue;
        
        const double denom = wk[a] * wl[b];
        if (!(denom > 0.0) || !std::isfinite(denom)) continue;
        
        // Truncation indicator: strict inequality (see FIDELITY)
        const double h1L_ind =
          (hL - std::max(zbk[a], zbl[b]) > std::max(e1k[a], e1l[b]))
          ? 1.0 : 0.0;
        
        if (h1L_ind == 0.0) continue;   // null weight: not counted in nnz
        
        const double c_ab = h1L_ind / denom;
        
        Ra[a] += c_ab;
        Cb[b] += c_ab;
        if (store_pairs) {
          pa_v.push_back(static_cast<int>(ik[a]));
          pb_v.push_back(static_cast<int>(il[b]));
          pc_v.push_back(c_ab);
        }
      }
    }

    out[j] = Rcpp::List::create(
      Rcpp::Named("ka") = fgp::to_int_idx(ik),
      Rcpp::Named("kb") = fgp::to_int_idx(il),
      Rcpp::Named("Ra") = Ra,
      Rcpp::Named("Cb") = Cb,
      Rcpp::Named("pa") = Rcpp::IntegerVector(pa_v.begin(), pa_v.end()),
      Rcpp::Named("pb") = Rcpp::IntegerVector(pb_v.begin(), pb_v.end()),
      Rcpp::Named("pc") = Rcpp::NumericVector(pc_v.begin(), pc_v.end())
    );
  }

  Rcpp::List meta = Rcpp::List::create(
    Rcpp::Named("p")         = p,
    Rcpp::Named("n2")        = n2,
    Rcpp::Named("has_pairs") = store_pairs,
    Rcpp::Named("cdq_cols")  = cdq_cols
  );
  out.attr("meta") = meta;
  return out;
}

// ===========================================================================
// EVALUATION
// ===========================================================================

// Shared core: components per column, returns arma::vec.
// (parts never goes through NumericVector in the internal path — there
// is no arma::vec(begin, end).)
namespace fgp {

arma::vec eval_parts(const Rcpp::List& pairs,
                     const arma::vec& lq,
                     const arma::mat& z,
                     const arma::vec& theta,
                     const std::string& phi) {
  
  validate_phi(phi);
  
  const Rcpp::List meta = fgp::get_meta(pairs);
  const int    p  = fgp::meta_p(meta);
  const double n2 = fgp::meta_n2(meta);
  
  if (static_cast<int>(z.n_cols) != p)
    stop("eval_fg_parts: z has %d columns, meta['p'] = %d",
         static_cast<int>(z.n_cols), p);
  if (static_cast<int>(lq.n_elem) != static_cast<int>(z.n_rows) ||
      static_cast<int>(theta.n_elem) != p)
    stop("eval_fg_parts: incompatible lq/z/theta dimensions");
  
  // Residuals for ALL rows, once: O(n*p)
  const arma::vec e2 = lq - z * theta;
  
  arma::vec parts(p, arma::fill::zeros);
  
  if (phi == "identity") {
    // O(nk + nl) per column via the marginals:
    //   parts(j) = -(1/n^2) [ sum_a e2k(a) Ra(a) - sum_b e2l(b) Cb(b) ]
    // Exact algebraic decomposition of the sum over pairs; reordered sum
    // (1e-10 from the reference kernel, not bit-exact). Degenerate columns
    // are skipped (parts(j) = 0).
    for (int j = 0; j < p; ++j) {
      const Rcpp::List colj = pairs[j];
      
      const arma::uvec ka = Rcpp::as<arma::uvec>(colj["ka"]);
      const arma::uvec kb = Rcpp::as<arma::uvec>(colj["kb"]);
      if (ka.n_elem == 0 || kb.n_elem == 0) continue;   // degenerate
      
      const arma::vec Ra = Rcpp::as<arma::vec>(colj["Ra"]);
      const arma::vec Cb = Rcpp::as<arma::vec>(colj["Cb"]);
      
      const arma::vec e2k = e2.elem(ka);
      const arma::vec e2l = e2.elem(kb);
      
      parts[j] = -(arma::dot(e2k, Ra) - arma::dot(e2l, Cb)) / n2;
    }
  } else {
    // piecewise: sequential O(nnz) over the frozen pairs
    //   parts(j) = -(1/n^2) sum_k pc(k) * phi(e2[pa(k)] - e2[pb(k)])
    // pa/pb/pc read by pointer (IntegerVector/NumericVector over the
    // SEXP do not copy when the type is already integer/double).
    if (!fgp::meta_has_pairs(meta))
      stop("eval_fg_parts: an object built with phi = 'identity' has no "
             "pairs; rebuild it with make_fg_pairs(..., 'piecewise')");
    const double* e2p = e2.memptr();
    for (int j = 0; j < p; ++j) {
      const Rcpp::List colj = pairs[j];

      const Rcpp::IntegerVector pa = colj["pa"];
      if (pa.size() == 0) continue;   // no valid pairs (or degenerate)

      const Rcpp::IntegerVector pb = colj["pb"];
      const Rcpp::NumericVector pc = colj["pc"];
      if (pa.size() != pb.size() || pa.size() != pc.size())
        stop("eval_fg_parts: column %d has misaligned pa/pb/pc", j + 1);

      const int*    pap = pa.begin();
      const int*    pbp = pb.begin();
      const double* pcp = pc.begin();
      double acc = 0.0;
      const R_xlen_t nnz = pa.size();
      for (R_xlen_t k = 0; k < nnz; ++k) {
        acc -= pcp[k] * fgp::phi_pw(e2p[pap[k]] - e2p[pbp[k]]);
      }
      parts[j] = acc / n2;
    }
  }
  
  return parts;
}

} // namespace fgp

// eval_fg_parts_cpp
//
// Components of the U-statistic from the pairs object.
// Returns a numeric vector of length p without a dim attribute.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector eval_fg_parts_cpp(const Rcpp::List& pairs,
                                      const arma::vec& lq,
                                      const arma::mat& z,
                                      const arma::vec& theta,
                                      const std::string& phi) {
  const arma::vec parts = fgp::eval_parts(pairs, lq, z, theta, phi);
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// eval_fg_obj_cpp
//
// Norm (L1/L2) of the components, standardized by sqrt(p). Interface
// compatible with U2_fun_fg_cpp for the A/B tests.
// [[Rcpp::export(rng = false)]]
double eval_fg_obj_cpp(const Rcpp::List& pairs,
                       const arma::vec& lq,
                       const arma::mat& z,
                       const arma::vec& theta,
                       const std::string& phi,
                       const std::string& norm_min) {
  const arma::vec parts = fgp::eval_parts(pairs, lq, z, theta, phi);
  return fgp::norm_aggregate(parts, norm_min);
}

// ===========================================================================
// PIECEWISE WITHOUT PAIRS (O(n) memory, O(n log^2 n) per evaluation)
//
// With phi = "piecewise", component j is
//   parts(j) = -(1/n^2) sum_{valid a, b} c_ab sign(e2[a] - e2[b]),
// and the pair (a, b) is valid when the members pass the per-element
// filters (event, weight, A_a > u_a, B_b > v_b) and the cross conditions
// v_b < A_a and u_a < B_b (the same decomposition as marginals_cdq), with
// c_ab = (1/wk[a]) (1/wl[b]). Construction stores only the active
// members; the evaluation solves the sum by divide and conquer on
// v_b < A_a, with a sweep on u_a < B_b and a Fenwick tree on the ranks of e2.
//
// The set of pairs is the same as the loop's (exact floating-point
// identity, see marginals_cdq); only the order of the sums changes
// (~1e-15 relative). Same FALLBACK as marginals_cdq: with non-finite
// values or weights outside [1e-150, 1e150], the column is flagged
// (ok = FALSE) and the fit uses the explicit pairs object.
// ===========================================================================

namespace fgp {

struct LeanCol {
  std::vector<int>    ia, ib;      // global rows of the active members
  std::vector<double> A, u, wa;    // group 0
  std::vector<double> B, v, wb;    // group 1
};

// Divide-and-conquer item: key of the condition v_b < A_a
struct PwItem {
  double key;      // A_a (member a) or v_b (member b)
  bool   is_a;
  int    k;        // position in the group vector
};

struct PwEval {
  const LeanCol* c;
  std::vector<double> e2a, e2b;    // e2 of the active members
  std::vector<double> vals;        // distinct e2 values (ranks)
  std::vector<PwItem> items;       // sorted by key (a before b)
  double acc;
};

inline void pw_cross(PwEval& E, int lo, int mid, int hi) {
  const LeanCol& c = *E.c;
  std::vector<int> bs, as;          // b in the left half, a in the right half
  for (int i = lo; i < mid; ++i) if (!E.items[i].is_a) bs.push_back(E.items[i].k);
  for (int i = mid; i < hi; ++i) if (E.items[i].is_a) as.push_back(E.items[i].k);
  if (bs.empty() || as.empty()) return;
  // u_a < B_b: b by decreasing B, a by decreasing u
  std::sort(bs.begin(), bs.end(), [&c](int x, int y) { return c.B[x] > c.B[y]; });
  std::sort(as.begin(), as.end(), [&c](int x, int y) { return c.u[x] > c.u[y]; });
  Fenwick fw(E.vals.size());
  double total = 0.0;
  size_t p = 0;
  std::vector<size_t> added;
  for (int a : as) {
    while (p < bs.size() && c.B[bs[p]] > c.u[a]) {
      const int b = bs[p++];
      const size_t r = std::lower_bound(E.vals.begin(), E.vals.end(),
                                        E.e2b[b]) - E.vals.begin();
      fw.add(r, c.wb[b]);
      total += c.wb[b];
    }
    if (total == 0.0) continue;
    const double x = E.e2a[a];
    const size_t lo_r = std::lower_bound(E.vals.begin(), E.vals.end(), x)
                        - E.vals.begin();               // e2b < x
    const size_t hi_r = std::upper_bound(E.vals.begin(), E.vals.end(), x)
                        - E.vals.begin();               // e2b <= x
    const double less = fw.sum(lo_r);
    const double greater = total - fw.sum(hi_r);
    E.acc += c.wa[a] * (less - greater);
  }
}

inline void pw_brute(PwEval& E, int lo, int hi) {
  const LeanCol& c = *E.c;
  for (int i = lo; i < hi; ++i) {
    if (E.items[i].is_a) continue;              // b comes before a
    const int b = E.items[i].k;
    for (int j = i + 1; j < hi; ++j) {
      if (!E.items[j].is_a) continue;
      const int a = E.items[j].k;
      if (!(c.u[a] < c.B[b])) continue;
      const double d = E.e2a[a] - E.e2b[b];
      if (d > 0.0) E.acc += c.wa[a] * c.wb[b];
      else if (d < 0.0) E.acc -= c.wa[a] * c.wb[b];
    }
  }
}

inline void pw_cdq(PwEval& E, int lo, int hi) {
  if (hi - lo <= 64) { pw_brute(E, lo, hi); return; }
  const int mid = lo + (hi - lo) / 2;
  pw_cross(E, lo, mid, hi);
  pw_cdq(E, lo, mid);
  pw_cdq(E, mid, hi);
}

// Sum of c_ab sign(e2a - e2b) of a column
inline double pw_column_sum(const LeanCol& c, const arma::vec& e2) {
  PwEval E;
  E.c = &c;
  E.acc = 0.0;
  const size_t na = c.ia.size(), nb = c.ib.size();
  if (na == 0 || nb == 0) return 0.0;
  E.e2a.resize(na); E.e2b.resize(nb);
  for (size_t k = 0; k < na; ++k) E.e2a[k] = e2[c.ia[k]];
  for (size_t k = 0; k < nb; ++k) E.e2b[k] = e2[c.ib[k]];
  E.vals = E.e2b;
  std::sort(E.vals.begin(), E.vals.end());
  E.vals.erase(std::unique(E.vals.begin(), E.vals.end()), E.vals.end());
  E.items.reserve(na + nb);
  for (size_t k = 0; k < na; ++k)
    E.items.push_back({c.A[k], true, static_cast<int>(k)});
  for (size_t k = 0; k < nb; ++k)
    E.items.push_back({c.v[k], false, static_cast<int>(k)});
  // increasing key; tie: a before b (b before a => v_b < A_a)
  std::sort(E.items.begin(), E.items.end(),
            [](const PwItem& x, const PwItem& y) {
              if (x.key != y.key) return x.key < y.key;
              return x.is_a && !y.is_a;
            });
  pw_cdq(E, 0, static_cast<int>(E.items.size()));
  return E.acc;
}

// Active members of a column; false => fallback (see the header)
inline bool lean_column(const double hL,
                        const arma::vec& zbk, const arma::vec& zbl,
                        const arma::vec& e1k, const arma::vec& e1l,
                        const arma::vec& dk,  const arma::vec& dl,
                        const arma::vec& wk,  const arma::vec& wl,
                        const arma::uvec& ik, const arma::uvec& il,
                        LeanCol& out) {
  if (!std::isfinite(hL) || !zbk.is_finite() || !zbl.is_finite() ||
      !e1k.is_finite() || !e1l.is_finite())
    return false;
  const double wmin = 1e-150, wmax = 1e150;
  bool in_range = true;
  auto weight_ok = [&](double w) {
    if (w < 0.0) { in_range = false; return false; }
    if (!(w > 0.0) || !std::isfinite(w)) return false;
    if (w < wmin || w > wmax) in_range = false;
    return true;
  };
  for (arma::uword a = 0; a < zbk.n_elem; ++a) {
    if (dk[a] == 0.0 || !weight_ok(wk[a])) continue;
    const double Aa = hL - zbk[a];
    if (!(Aa > e1k[a])) continue;
    out.ia.push_back(static_cast<int>(ik[a]));
    out.A.push_back(Aa); out.u.push_back(e1k[a]); out.wa.push_back(1.0 / wk[a]);
  }
  for (arma::uword b = 0; b < zbl.n_elem; ++b) {
    if (dl[b] == 0.0 || !weight_ok(wl[b])) continue;
    const double Bb = hL - zbl[b];
    if (!(Bb > e1l[b])) continue;
    out.ib.push_back(static_cast<int>(il[b]));
    out.B.push_back(Bb); out.v.push_back(e1l[b]); out.wb.push_back(1.0 / wl[b]);
  }
  return in_range;
}

inline LeanCol lean_from_list(const Rcpp::List& col) {
  LeanCol c;
  c.ia = Rcpp::as<std::vector<int>>(col["ia"]);
  c.ib = Rcpp::as<std::vector<int>>(col["ib"]);
  c.A  = Rcpp::as<std::vector<double>>(col["A"]);
  c.u  = Rcpp::as<std::vector<double>>(col["u"]);
  c.wa = Rcpp::as<std::vector<double>>(col["wa"]);
  c.B  = Rcpp::as<std::vector<double>>(col["B"]);
  c.v  = Rcpp::as<std::vector<double>>(col["v"]);
  c.wb = Rcpp::as<std::vector<double>>(col["wb"]);
  return c;
}

} // namespace fgp

// make_fg_pw_lean_cpp
//
// Piecewise object without pairs: per column, the active members (ia/A/u/wa
// in group 0, ib/B/v/wb in group 1). attr "meta": p, n2 and ok (FALSE
// if some column needs the fallback to the explicit pairs).
// [[Rcpp::export(rng = false)]]
Rcpp::List make_fg_pw_lean_cpp(const arma::vec& lt,
                               const arma::vec& lq,
                               const arma::mat& z,
                               const arma::vec& status,
                               const arma::vec& beta,
                               const double hL,
                               const Rcpp::List& km_surv) {
  const int n = static_cast<int>(z.n_rows);
  const int p = static_cast<int>(z.n_cols);
  if (static_cast<int>(lt.n_elem) != n || static_cast<int>(lq.n_elem) != n ||
      static_cast<int>(status.n_elem) != n)
    stop("make_fg_pw_lean: lt/lq/status must have length n");
  if (static_cast<int>(beta.n_elem) != p ||
      static_cast<int>(km_surv.size()) != p)
    stop("make_fg_pw_lean: beta and km_surv must have length p = ncol(z)");

  const arma::vec e1 = lt - z * beta;
  Rcpp::List out(p);
  bool ok = true;
  for (int j = 0; j < p; ++j) {
    const arma::uvec ik = arma::find(z.col(j) == 0.0);
    const arma::uvec il = arma::find(z.col(j) == 1.0);
    const Rcpp::List kmj = km_surv[j];
    const arma::vec wk = Rcpp::as<arma::vec>(kmj["w0"]);
    const arma::vec wl = Rcpp::as<arma::vec>(kmj["w1"]);
    if (wk.n_elem != ik.n_elem || wl.n_elem != il.n_elem)
      stop("make_fg_pw_lean: km_surv[[%d]] w0/w1 do not match the groups",
           j + 1);
    const arma::mat zk = z.rows(ik), zl = z.rows(il);
    // dot per member, as in the reference kernel (see FIDELITY)
    arma::vec zbk(ik.n_elem), zbl(il.n_elem);
    for (arma::uword a = 0; a < ik.n_elem; ++a) zbk[a] = arma::dot(zk.row(a), beta);
    for (arma::uword b = 0; b < il.n_elem; ++b) zbl[b] = arma::dot(zl.row(b), beta);
    fgp::LeanCol c;
    if (ik.n_elem > 0 && il.n_elem > 0)
      ok = fgp::lean_column(hL, zbk, zbl, e1.elem(ik), e1.elem(il),
                            status.elem(ik), status.elem(il), wk, wl,
                            ik, il, c) && ok;
    out[j] = Rcpp::List::create(
      Rcpp::Named("ia") = Rcpp::wrap(c.ia), Rcpp::Named("A") = Rcpp::wrap(c.A),
      Rcpp::Named("u")  = Rcpp::wrap(c.u),  Rcpp::Named("wa") = Rcpp::wrap(c.wa),
      Rcpp::Named("ib") = Rcpp::wrap(c.ib), Rcpp::Named("B") = Rcpp::wrap(c.B),
      Rcpp::Named("v")  = Rcpp::wrap(c.v),  Rcpp::Named("wb") = Rcpp::wrap(c.wb));
  }
  out.attr("meta") = Rcpp::List::create(
    Rcpp::Named("p") = p,
    Rcpp::Named("n2") = static_cast<double>(n) * static_cast<double>(n),
    Rcpp::Named("ok") = ok);
  return out;
}

// eval_fg_parts_pw_lean_cpp: piecewise components from the object
// without pairs
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector eval_fg_parts_pw_lean_cpp(const Rcpp::List& lean,
                                              const arma::vec& lq,
                                              const arma::mat& z,
                                              const arma::vec& theta) {
  const Rcpp::List meta = fgp::get_meta(lean);
  const int p = fgp::meta_p(meta);
  const double n2 = fgp::meta_n2(meta);
  if (static_cast<int>(z.n_cols) != p ||
      static_cast<int>(theta.n_elem) != p ||
      lq.n_elem != z.n_rows)
    stop("eval_fg_parts_pw_lean: incompatible lq/z/theta dimensions");
  const arma::vec e2 = lq - z * theta;
  Rcpp::NumericVector parts(p);
  for (int j = 0; j < p; ++j) {
    const fgp::LeanCol c = fgp::lean_from_list(lean[j]);
    parts[j] = -fgp::pw_column_sum(c, e2) / n2;
  }
  return parts;
}

// eval_fg_obj_pw_lean_cpp: norm (L1/L2) of the components, by sqrt(p)
// [[Rcpp::export(rng = false)]]
double eval_fg_obj_pw_lean_cpp(const Rcpp::List& lean,
                               const arma::vec& lq,
                               const arma::mat& z,
                               const arma::vec& theta,
                               const std::string& norm_min) {
  const Rcpp::NumericVector parts = eval_fg_parts_pw_lean_cpp(lean, lq, z, theta);
  const arma::vec pv(parts.begin(), parts.size());
  return fgp::norm_aggregate(pv, norm_min);
}
