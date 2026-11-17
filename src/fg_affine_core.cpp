// src/fg_affine_core.cpp
//
// Exact affine form (phi = "identity") and piecewise evaluation over
// flattened vectors for the U-statistic of Fine and Gelber, from the
// pairs object of fg_pairs_core.cpp.
//
// SCHEMA OF THE `pairs` OBJECT (produced by make_fg_pairs_cpp):
//   pairs[[j]] (list per binary covariate column j):
//     "ka"  - integer, global (0-based) indices of the rows of group 0
//     "kb"  - integer, global indices of the rows of group 1
//     "Ra"  - vec, marginal sum_b c_ab per member of group 0 (aligned with ka)
//     "Cb"  - vec, marginal sum_a c_ab per member of group 1 (aligned with kb)
//     "pa"  - integer, global index of the group 0 member, per valid pair
//     "pb"  - integer, global index of the group 1 member, per valid pair
//     "pc"  - vec, weight c_ab = h1L(a,b)/denom(a,b) per valid pair
//   attr(pairs, "meta"): "p" (int), "n2" (double, n^2), "has_pairs",
//     "cdq_cols"
//
// FIDELITY (checked by the A/B tests against fg_core.cpp):
//   - identity: parts_j(theta) = alpha_j + W_j theta is the exact
//     algebraic decomposition of -(1/n^2)[ sum_a e2k(a) Ra_a - sum_b e2l(b) Cb_b ],
//     e2 = lq - z theta; reordered sums (1e-10, not bit-exact).
//   - piecewise: same pair order and arithmetic as the reference kernel;
//     indicator 0 -> 0, > 0 -> 1, < 0 -> -1; division by n^2 at the
//     end.
//   - degenerate columns (empty group): Ra/Cb zeroed, no pairs =>
//     parts(j) = 0.

#include <RcppArmadillo.h>
#include <climits>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace fga {

// attr("meta") returns a proxy without operator[]: materialize it as a List
inline Rcpp::List get_meta(const Rcpp::List& obj) {
  const Rcpp::List meta = obj.attr("meta");   // materialize as a List
  return meta;
}

inline int meta_p(const Rcpp::List& meta) {
  return Rcpp::as<int>(meta["p"]);
}

inline double meta_n2(const Rcpp::List& meta) {
  return Rcpp::as<double>(meta["n2"]);
}

// Column field access with a clear diagnostic if the name diverges from
// the schema of fg_pairs_core.cpp (avoids a silent segfault).
inline Rcpp::List get_col(const Rcpp::List& pairs, int j) {
  if (j < 0 || j >= pairs.size())
    stop("fg_affine: invalid column index: %d", j);
  const Rcpp::List colj = pairs[j];
  return colj;
}

// Aggregation by the L1/L2 norm, standardized by sqrt(p)
inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// piecewise phi: sign(x) with 0 -> 0
inline double phi_pw(const double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

} // namespace fga

// ===========================================================================
// AFFINE PART (phi = "identity")
// ===========================================================================

// make_fg_affine_cpp
//
// Extracts the exact affine form of the components of the U-statistic:
//
//   parts(theta) = alpha + W * theta,   alpha in R^p, W in R^{p x p}
//
// Derivation (exact algebraic decomposition of the sum over pairs):
//
//   parts_j(theta) = -(1/n^2) [ sum_a e2k(a) Ra_a - sum_b e2l(b) Cb_b ]
//     with e2k(a) = lq[ka] - zk(a) * theta,  e2l(b) = lq[kb] - zl(b) * theta
//
//   => alpha_j = -(1/n^2) [ lqk' Ra - lql' Cb ]
//   => W(j,:)  =  (1/n^2) [ Ra' zk - Cb' zl ]
//
// Each later evaluation costs O(p^2), independent of n.
//
// Arguments:
//   pairs - pairs object (schema in the header of the file)
//   lq    - h2$linkfun(qt), precomputed (constant in the optimization)
//   z     - binary covariate matrix (n x p)
//
// Returns List(alpha = vec p, W = mat p x p) with the meta propagated.
// [[Rcpp::export]]
Rcpp::List make_fg_affine_cpp(const Rcpp::List& pairs,
                              const arma::vec& lq,
                              const arma::mat& z) {
  
  const Rcpp::List meta = fga::get_meta(pairs);
  const int    p  = fga::meta_p(meta);
  const double n2 = fga::meta_n2(meta);
  
  const int nz = static_cast<int>(z.n_cols);
  if (nz != p)
    stop("make_fg_affine: z has %d columns, meta['p'] = %d", nz, p);
  if (static_cast<int>(lq.n_elem) != static_cast<int>(z.n_rows))
    stop("make_fg_affine: lq and z have different numbers of rows");
  
  arma::vec alpha(p, arma::fill::zeros);
  arma::mat W(p, p, arma::fill::zeros);
  
  for (int j = 0; j < p; ++j) {
    const Rcpp::List colj = fga::get_col(pairs, j);
    
    const arma::uvec ka = Rcpp::as<arma::uvec>(colj["ka"]);
    const arma::uvec kb = Rcpp::as<arma::uvec>(colj["kb"]);
    const arma::vec  Ra = Rcpp::as<arma::vec>( colj["Ra"]);
    const arma::vec  Cb = Rcpp::as<arma::vec>( colj["Cb"]);
    
    if (ka.n_elem != Ra.n_elem || kb.n_elem != Cb.n_elem)
      stop("make_fg_affine: column %d has misaligned ka/Ra or kb/Cb", j);
    
    // Degenerate column (empty group): alpha = W_j = 0
    if (ka.n_elem == 0 && kb.n_elem == 0) continue;
    
    const arma::vec lqk = lq.elem(ka);
    const arma::vec lql = lq.elem(kb);
    const arma::mat zk  = z.rows(ka);
    const arma::mat zl  = z.rows(kb);
    
    alpha(j) = -(arma::dot(lqk, Ra) - arma::dot(lql, Cb)) / n2;
    
    const arma::rowvec ra_t = arma::conv_to<arma::rowvec>::from(Ra);
    const arma::rowvec cb_t = arma::conv_to<arma::rowvec>::from(Cb);
    W.row(j) = (ra_t * zk - cb_t * zl) / n2;
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("W")     = W
  );
  out.attr("meta") = meta;
  return out;
}

// eval_fg_parts_affine_cpp
//
// Components in the affine form: parts(theta) = alpha + W * theta.
// O(p^2). Returns an R vector without a dim attribute (converted at the boundary).
// [[Rcpp::export]]
Rcpp::NumericVector eval_fg_parts_affine_cpp(const Rcpp::List& aff,
                                             const arma::vec& theta) {
  const arma::vec alpha = Rcpp::as<arma::vec>(aff["alpha"]);
  const arma::mat W     = Rcpp::as<arma::mat>(aff["W"]);
  
  const arma::vec parts = alpha + W * theta;
  
  // NumericVector: an exported arma::vec would reach R as a p x 1 matrix
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// eval_fg_obj_affine_cpp
//
// Norm (L1/L2) of the affine components, standardized by sqrt(p).
// [[Rcpp::export]]
double eval_fg_obj_affine_cpp(const Rcpp::List& aff,
                              const arma::vec& theta,
                              const std::string& norm_min) {
  const arma::vec alpha = Rcpp::as<arma::vec>(aff["alpha"]);
  const arma::mat W     = Rcpp::as<arma::mat>(aff["W"]);
  
  const arma::vec parts = alpha + W * theta;
  
  return fga::norm_aggregate(parts, norm_min);
}

// solve_fg_theta_affine_cpp
//
// Fast path (L2 + identity): minimizes ||alpha + W theta||^2 in closed
// form by the normal equations
//
//     G theta = -g,   G = W'W,  g = W'alpha
//
// Rank deficiency (degenerate columns, partial identifiability): detected
// by an extreme singular value ratio < rcond OR a large residual of the
// direct solution; in these cases the MINIMUM NORM solution via pinv is
// used and full_rank = FALSE — the point estimate remains valid, but the
// caller (semi_aft_fg) warns the user.
//
// The gradient at the solution (2 W'(alpha + W theta)) is returned as a
// check: at full rank it must be ~0 (semi_aft_fg requires < 1e-6).
//
// Returns List:
//   par        - theta coefficient (R vector WITHOUT dim)
//   full_rank  - logical
//   grad_norm  - norm of the gradient of the quadratic objective at the solution
// [[Rcpp::export]]
Rcpp::List solve_fg_theta_affine_cpp(const Rcpp::List& aff,
                                     const double rcond = 1e-10) {
  
  const arma::vec alpha = Rcpp::as<arma::vec>(aff["alpha"]);
  const arma::mat W     = Rcpp::as<arma::mat>(aff["W"]);
  const int p = static_cast<int>(alpha.n_elem);
  
  const arma::mat G = W.t() * W;   // p x p, symmetric positive semidefinite
  const arma::vec g = W.t() * alpha;
  
  // --- Rank detection by the extreme singular value ratio ----------
  bool full_rank = false;
  if (p > 0) {
    arma::mat U, V;
    arma::vec sv;
    if (arma::svd(U, sv, V, G) && sv.n_elem > 0) {
      const double smax = sv(0);
      const double smin = sv(sv.n_elem - 1);
      full_rank = (smax > 0.0) && (smin / smax >= rcond);
    }
  }
  
  // --- Direct solution (full rank) ---------------------------------
  arma::vec theta;
  bool solved = false;
  if (full_rank) {
    const arma::vec cand = arma::solve(G, -g,
                                       arma::solve_opts::likely_sympd);
    if (cand.n_elem == static_cast<arma::uword>(p) && cand.is_finite()) {
      // Residual check: solves even with poor conditioning
      const double res  = arma::norm(G * cand + g, 2);
      const double gsc  = arma::norm(g, 2);
      const double tol  = 1e-8 * (1.0 + gsc);
      if (res <= tol) {
        theta  = cand;
        solved = true;
      }
    }
  }
  
  // --- Fallback: minimum norm (rank deficient) ----------------------
  if (!solved) {
    if (p > 0) {
      theta = -(arma::pinv(G, rcond) * g);
    } else {
      theta.reset();
    }
    full_rank = false;
  }
  
  if (!theta.is_finite()) {
    // Pathological case (G == 0, alpha != 0): theta = 0 and a warning via grad
    theta.zeros(p);
    full_rank = false;
  }
  
  // Gradient of the objective ||alpha + W theta||^2 at the solution: 2 W' parts
  const arma::vec parts    = alpha + W * theta;
  const double   grad_norm = (p > 0) ? arma::norm(W.t() * parts, 2) : 0.0;
  
  return Rcpp::List::create(
    Rcpp::Named("par")       = Rcpp::NumericVector(theta.begin(),
                theta.end()),
                Rcpp::Named("full_rank") = full_rank,
                Rcpp::Named("grad_norm") = grad_norm
  );
}

// ===========================================================================
// PIECEWISE PART OVER FLATTENED VECTORS
// ===========================================================================

// compile_fg_pw_cpp
//
// Flattens the pairs object into global vectors for the sequential
// O(nnz) piecewise evaluation (the event/weight/h1L filters were already
// resolved when the pairs were built).
//
// Output layout (CSR per column):
//   pa, pb  - integer: global rows of the group 0/1 members of each pair
//   pc      - vec: weight c_ab of each pair
//   col_ptr - integer (p+1): column boundaries (col_ptr[j] .. col_ptr[j+1])
//   meta    - propagated
// [[Rcpp::export]]
Rcpp::List compile_fg_pw_cpp(const Rcpp::List& pairs) {
  
  const Rcpp::List meta = fga::get_meta(pairs);
  const int    p  = fga::meta_p(meta);
  const double n2 = fga::meta_n2(meta);
  const double n  = std::sqrt(n2);

  // An object built with phi = "identity" stores no pairs: an explicit
  // error instead of null components.
  if (meta.containsElementNamed("has_pairs") &&
      !Rcpp::as<bool>(meta["has_pairs"]))
    stop("compile_fg_pw: an object built with phi = 'identity' has no "
           "pairs; rebuild it with make_fg_pairs(..., 'piecewise')");

  // Indices as integer, read by pointer in pw_parts; the range [0, n)
  // is validated here, once, and not at each evaluation.
  std::vector<int>    pa_v, pb_v;
  std::vector<double> pc_v;
  std::vector<int>    col_ptr_v;
  col_ptr_v.reserve(static_cast<size_t>(p) + 1);
  col_ptr_v.push_back(0);

  for (int j = 0; j < p; ++j) {
    const Rcpp::List colj = fga::get_col(pairs, j);

    // The coercion to integer also accepts double indices
    const Rcpp::IntegerVector paj = colj["pa"];
    const Rcpp::IntegerVector pbj = colj["pb"];
    const Rcpp::NumericVector pcj = colj["pc"];

    if (paj.size() != pbj.size() || paj.size() != pcj.size())
      stop("compile_fg_pw: column %d has misaligned pa/pb/pc", j);

    for (R_xlen_t k = 0; k < paj.size(); ++k) {
      if (paj[k] < 0 || pbj[k] < 0 || paj[k] >= n || pbj[k] >= n)
        stop("compile_fg_pw: column %d has an index outside [0, n)", j);
    }

    pa_v.insert(pa_v.end(), paj.begin(), paj.end());
    pb_v.insert(pb_v.end(), pbj.begin(), pbj.end());
    pc_v.insert(pc_v.end(), pcj.begin(), pcj.end());

    if (pa_v.size() > static_cast<size_t>(INT_MAX))
      stop("compile_fg_pw: number of pairs exceeds INT_MAX");
    col_ptr_v.push_back(static_cast<int>(pa_v.size()));
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("pa")      = Rcpp::IntegerVector(pa_v.begin(), pa_v.end()),
    Rcpp::Named("pb")      = Rcpp::IntegerVector(pb_v.begin(), pb_v.end()),
    Rcpp::Named("pc")      = Rcpp::NumericVector(pc_v.begin(), pc_v.end()),
    Rcpp::Named("col_ptr") = Rcpp::IntegerVector(col_ptr_v.begin(),
                                                 col_ptr_v.end())
  );
  out.attr("meta") = meta;
  return out;
}

// Shared core: piecewise components as arma::vec
namespace fga {

arma::vec pw_parts(const Rcpp::List& pw,
                   const arma::vec& lq,
                   const arma::mat& z,
                   const arma::vec& theta) {
  
  const Rcpp::List meta = fga::get_meta(pw);
  const int    p  = fga::meta_p(meta);
  const double n2 = fga::meta_n2(meta);
  
  // Copy-free views over the R vectors (IntegerVector/NumericVector over
  // a SEXP of the right type do not allocate)
  const Rcpp::IntegerVector pa      = pw["pa"];
  const Rcpp::IntegerVector pb      = pw["pb"];
  const Rcpp::NumericVector pc      = pw["pc"];
  const Rcpp::IntegerVector col_ptr = pw["col_ptr"];

  if (static_cast<int>(col_ptr.size()) != p + 1)
    stop("eval_fg_parts_pw: col_ptr has length %d, expected %d",
         static_cast<int>(col_ptr.size()), p + 1);
  if (pa.size() != pb.size() || pa.size() != pc.size() ||
      col_ptr[p] != pa.size())
    stop("eval_fg_parts_pw: misaligned pa/pb/pc/col_ptr");
  if (std::sqrt(n2) != static_cast<double>(z.n_rows) ||
      lq.n_elem != z.n_rows || static_cast<int>(theta.n_elem) != p ||
      static_cast<int>(z.n_cols) != p)
    stop("eval_fg_parts_pw: incompatible lq/z/theta dimensions");

  // Residuals e2 = lq - z*theta for ALL rows, once: O(n p)
  const arma::vec e2 = lq - z * theta;
  const double* e2p = e2.memptr();
  const int*    pap = pa.begin();
  const int*    pbp = pb.begin();
  const double* pcp = pc.begin();

  arma::vec parts(p, arma::fill::zeros);

  // Pure sequential sweep: each valid pair contributes
  //   -c_ab * phi(e2[pa] - e2[pb]),  piecewise phi with 0 -> 0
  for (int j = 0; j < p; ++j) {
    const int lo = col_ptr[j];
    const int hi = col_ptr[j + 1];
    double acc = 0.0;
    for (int k = lo; k < hi; ++k) {
      acc -= pcp[k] * fga::phi_pw(e2p[pap[k]] - e2p[pbp[k]]);
    }
    parts[j] = acc / n2;
  }
  return parts;
}

} // namespace fga

// eval_fg_parts_pw_cpp
//
// Piecewise components. Returns an R vector without dim.
// [[Rcpp::export]]
Rcpp::NumericVector eval_fg_parts_pw_cpp(const Rcpp::List& pw,
                                         const arma::vec& lq,
                                         const arma::mat& z,
                                         const arma::vec& theta) {
  const arma::vec parts = fga::pw_parts(pw, lq, z, theta);
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// eval_fg_obj_pw_cpp
//
// Norm (L1/L2) of the piecewise components, standardized by sqrt(p).
// [[Rcpp::export]]
double eval_fg_obj_pw_cpp(const Rcpp::List& pw,
                          const arma::vec& lq,
                          const arma::mat& z,
                          const arma::vec& theta,
                          const std::string& norm_min) {
  const arma::vec parts = fga::pw_parts(pw, lq, z, theta);
  return fga::norm_aggregate(parts, norm_min);
}
