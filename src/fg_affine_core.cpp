// src/fg_affine_core.cpp
//
// Forma afim exata (phi = "identity") e avaliacao piecewise sobre vetores
// achatados para a U-estatistica de Fine e Gelber, a partir do objeto de
// pares de fg_pairs_core.cpp.
//
// SCHEMA DO OBJETO `pairs` (produzido por make_fg_pairs_cpp):
//   pairs[[j]] (lista por coluna j de covariavel binaria):
//     "ka"  - integer, indices globais (0-based) das linhas do grupo 0
//     "kb"  - integer, indices globais das linhas do grupo 1
//     "Ra"  - vec, marginal sum_b c_ab por membro do grupo 0 (alinhada a ka)
//     "Cb"  - vec, marginal sum_a c_ab por membro do grupo 1 (alinhada a kb)
//     "pa"  - integer, indice global do membro do grupo 0, por par valido
//     "pb"  - integer, indice global do membro do grupo 1, por par valido
//     "pc"  - vec, peso c_ab = h1L(a,b)/denom(a,b) por par valido
//   attr(pairs, "meta"): "p" (int), "n2" (double, n^2), "has_pairs",
//     "cdq_cols"
//
// FIDELIDADE (verificada nos testes A/B contra fg_core.cpp):
//   - identity: parts_j(theta) = alpha_j + W_j theta e a decomposicao
//     algebrica exata de -(1/n^2)[ sum_a e2k(a) Ra_a - sum_b e2l(b) Cb_b ],
//     e2 = lq - z theta; somas reordenadas (1e-10, nao bit-exata).
//   - piecewise: mesma ordem de pares e mesma aritmetica do kernel de
//     referencia; indicador 0 -> 0, > 0 -> 1, < 0 -> -1; divisao por n^2
//     no fim.
//   - colunas degeneradas (grupo vazio): Ra/Cb zerados, sem pares =>
//     parts(j) = 0.

#include <RcppArmadillo.h>
#include <climits>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace fga {

// attr("meta") devolve um proxy sem operator[]: materializa em List
inline Rcpp::List get_meta(const Rcpp::List& obj) {
  const Rcpp::List meta = obj.attr("meta");   // materializa em List
  return meta;
}

inline int meta_p(const Rcpp::List& meta) {
  return Rcpp::as<int>(meta["p"]);
}

inline double meta_n2(const Rcpp::List& meta) {
  return Rcpp::as<double>(meta["n2"]);
}

// Acesso a campo de coluna com diagnostico claro se o nome divergir
// do schema do fg_pairs_core.cpp (evita segfault silencioso).
inline Rcpp::List get_col(const Rcpp::List& pairs, int j) {
  if (j < 0 || j >= pairs.size())
    stop("fg_affine: invalid column index: %d", j);
  const Rcpp::List colj = pairs[j];
  return colj;
}

// Agregacao pela norma L1/L2, padronizada por sqrt(p)
inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// phi piecewise: sign(x) com 0 -> 0
inline double phi_pw(const double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

} // namespace fga

// ===========================================================================
// PARTE AFIM (phi = "identity")
// ===========================================================================

// make_fg_affine_cpp
//
// Extrai a forma afim exata das componentes da U-estatistica:
//
//   parts(theta) = alpha + W * theta,   alpha em R^p, W em R^{p x p}
//
// Derivacao (decomposicao algebrica exata da soma de pares):
//
//   parts_j(theta) = -(1/n^2) [ sum_a e2k(a) Ra_a - sum_b e2l(b) Cb_b ]
//     com e2k(a) = lq[ka] - zk(a) * theta,  e2l(b) = lq[kb] - zl(b) * theta
//
//   => alpha_j = -(1/n^2) [ lqk' Ra - lql' Cb ]
//   => W(j,:)  =  (1/n^2) [ Ra' zk - Cb' zl ]
//
// Avaliacao subsequente custa O(p^2), independente de n.
//
// Argumentos:
//   pairs - objeto de pares (schema no cabecalho do arquivo)
//   lq    - h2$linkfun(qt), pre-computado (constante na otimizacao)
//   z     - matriz de covariaveis binarias (n x p)
//
// Retorna List(alpha = vec p, W = mat p x p) com meta propagada.
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
    
    // Coluna degenerada (grupo vazio): alpha = W_j = 0
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
// Componentes na forma afim: parts(theta) = alpha + W * theta.
// O(p^2). Retorna vetor R sem atributo dim (conversao na fronteira).
// [[Rcpp::export]]
Rcpp::NumericVector eval_fg_parts_affine_cpp(const Rcpp::List& aff,
                                             const arma::vec& theta) {
  const arma::vec alpha = Rcpp::as<arma::vec>(aff["alpha"]);
  const arma::mat W     = Rcpp::as<arma::mat>(aff["W"]);
  
  const arma::vec parts = alpha + W * theta;
  
  // NumericVector: um arma::vec exportado chegaria ao R como matriz p x 1
  return Rcpp::NumericVector(parts.begin(), parts.end());
}

// eval_fg_obj_affine_cpp
//
// Norma (L1/L2) das componentes afins, padronizada por sqrt(p).
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
// Caminho rapido (L2 + identity): minimiza ||alpha + W theta||^2 em
// forma fechada via equacoes normais
//
//     G theta = -g,   G = W'W,  g = W'alpha
//
// Posto deficiente (colunas degeneradas, identificabilidade parcial):
// detectado por razao singular extremo < rcond OU residuo alto da
// solucao direta; nesses casos usa-se a solucao de NORMA MINIMA via
// pinv e full_rank = FALSE — a estimativa pontual segue valida, mas o
// caller (semi_aft_fg) avisa o usuario.
//
// O gradiente na solucao (2 W'(alpha + W theta)) e retornado como
// verificacao: no posto completo deve ser ~0 (semi_aft_fg exige < 1e-6).
//
// Retorna List:
//   par        - coeficiente theta (vetor R SEM dim)
//   full_rank  - logical
//   grad_norm  - norma do gradiente do objetivo quadratico na solucao
// [[Rcpp::export]]
Rcpp::List solve_fg_theta_affine_cpp(const Rcpp::List& aff,
                                     const double rcond = 1e-10) {
  
  const arma::vec alpha = Rcpp::as<arma::vec>(aff["alpha"]);
  const arma::mat W     = Rcpp::as<arma::mat>(aff["W"]);
  const int p = static_cast<int>(alpha.n_elem);
  
  const arma::mat G = W.t() * W;   // p x p, simetrica semidefinida positiva
  const arma::vec g = W.t() * alpha;
  
  // --- Deteccao de posto via razao singular extremos -------------
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
  
  // --- Solucao direta (posto completo) -----------------------------
  arma::vec theta;
  bool solved = false;
  if (full_rank) {
    const arma::vec cand = arma::solve(G, -g,
                                       arma::solve_opts::likely_sympd);
    if (cand.n_elem == static_cast<arma::uword>(p) && cand.is_finite()) {
      // Verificacao por residuo: resolve mesmo com condicionamento ruim
      const double res  = arma::norm(G * cand + g, 2);
      const double gsc  = arma::norm(g, 2);
      const double tol  = 1e-8 * (1.0 + gsc);
      if (res <= tol) {
        theta  = cand;
        solved = true;
      }
    }
  }
  
  // --- Fallback: norma minima (posto deficiente) ---------------------
  if (!solved) {
    if (p > 0) {
      theta = -(arma::pinv(G, rcond) * g);
    } else {
      theta.reset();
    }
    full_rank = false;
  }
  
  if (!theta.is_finite()) {
    // Caso patologico (G == 0, alpha != 0): theta = 0 e aviso via grad
    theta.zeros(p);
    full_rank = false;
  }
  
  // Gradiente do objetivo ||alpha + W theta||^2 na solucao: 2 W' parts
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
// PARTE PIECEWISE SOBRE VETORES ACHATADOS
// ===========================================================================

// compile_fg_pw_cpp
//
// Achata o objeto de pares em vetores globais para avaliacao piecewise
// O(nnz) sequencial (os filtros de eventos/pesos/h1L ja foram
// resolvidos na construcao dos pares).
//
// Layout de saida (CSR por coluna):
//   pa, pb  - integer: linhas globais dos membros do grupo 0/1 de cada par
//   pc      - vec: peso c_ab de cada par
//   col_ptr - integer (p+1): fronteiras de coluna (col_ptr[j] .. col_ptr[j+1])
//   meta    - propagada
// [[Rcpp::export]]
Rcpp::List compile_fg_pw_cpp(const Rcpp::List& pairs) {
  
  const Rcpp::List meta = fga::get_meta(pairs);
  const int    p  = fga::meta_p(meta);
  const double n2 = fga::meta_n2(meta);
  const double n  = std::sqrt(n2);

  // Objeto construido com phi = "identity" nao guarda pares: erro
  // explicito em vez de componentes nulas.
  if (meta.containsElementNamed("has_pairs") &&
      !Rcpp::as<bool>(meta["has_pairs"]))
    stop("compile_fg_pw: an object built with phi = 'identity' has no "
           "pairs; rebuild it with make_fg_pairs(..., 'piecewise')");

  // Indices como integer, lidos por ponteiro em pw_parts; a faixa [0, n)
  // e validada aqui, uma vez, e nao a cada avaliacao.
  std::vector<int>    pa_v, pb_v;
  std::vector<double> pc_v;
  std::vector<int>    col_ptr_v;
  col_ptr_v.reserve(static_cast<size_t>(p) + 1);
  col_ptr_v.push_back(0);

  for (int j = 0; j < p; ++j) {
    const Rcpp::List colj = fga::get_col(pairs, j);

    // A coercao para integer aceita tambem indices double
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

// Nucleo compartilhado: componentes piecewise como arma::vec
namespace fga {

arma::vec pw_parts(const Rcpp::List& pw,
                   const arma::vec& lq,
                   const arma::mat& z,
                   const arma::vec& theta) {
  
  const Rcpp::List meta = fga::get_meta(pw);
  const int    p  = fga::meta_p(meta);
  const double n2 = fga::meta_n2(meta);
  
  // Visoes sem copia sobre os vetores do R (IntegerVector/NumericVector
  // sobre SEXP do tipo certo nao alocam)
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

  // Residuos e2 = lq - z*theta para TODAS as linhas, uma vez: O(n p)
  const arma::vec e2 = lq - z * theta;
  const double* e2p = e2.memptr();
  const int*    pap = pa.begin();
  const int*    pbp = pb.begin();
  const double* pcp = pc.begin();

  arma::vec parts(p, arma::fill::zeros);

  // Varredura sequencial pura: cada par valido contribui
  //   -c_ab * phi(e2[pa] - e2[pb]),  phi piecewise com 0 -> 0
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
// Componentes piecewise. Retorna vetor R sem dim.
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
// Norma (L1/L2) das componentes piecewise, padronizada por sqrt(p).
// [[Rcpp::export]]
double eval_fg_obj_pw_cpp(const Rcpp::List& pw,
                          const arma::vec& lq,
                          const arma::mat& z,
                          const arma::vec& theta,
                          const std::string& norm_min) {
  const arma::vec parts = fga::pw_parts(pw, lq, z, theta);
  return fga::norm_aggregate(parts, norm_min);
}
