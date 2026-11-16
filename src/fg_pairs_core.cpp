// src/fg_pairs_core.cpp
//
// Objeto de pares pre-computados da U-estatistica de Fine e Gelber. As
// constantes da etapa de theta (indicador de truncagem h1L, pesos KM,
// filtros de eventos) sao resolvidas uma vez por ajuste; a avaliacao do
// objetivo fica O(nk + nl) (identity, via marginais Ra/Cb) ou O(nnz)
// (piecewise). O kernel de referencia e o laco de fg_core.cpp, que este
// objeto reproduz (testes A/B).
//
// Par (a, b) valido, com a no grupo 0 e b no grupo 1 da coluna j:
//   - emparelhamento posicional a < b (indices dentro dos grupos);
//   - evento nos dois membros;
//   - denom = wk[a] * wl[b] finito e > 0;
//   - h1L: hL - max(zbk[a], zbl[b]) > max(e1k[a], e1l[b]) (estrita).
// O peso do par e c_ab = 1 / denom; a divisao por n^2 e o sinal sao
// aplicados na avaliacao: parts(j) = -(soma) / n^2.
//
// Com phi = "identity", Ra/Cb sao calculados em O(N log^2 N)
// (marginals_cdq, abaixo) e os pares nao sao guardados; com "piecewise",
// pelo laco O(nk * nl), que tambem guarda pa/pb/pc.
//
// FIDELIDADE: zbk/zbl usam arma::dot por membro, a mesma expressao do
// kernel de referencia, para que a comparacao estrita de h1L seja bit a
// bit a mesma (o teste de empate exato depende disso). Nao trocar por
// z * beta: o produto matricial reordena somas e pode inverter empates.
//
// SCHEMA DE SAIDA:
//   INVARIANTES (validadas por quem consome):
//     len(Ra) == len(ka), len(Cb) == len(kb)   (inclusive em colunas
//                                                degeneradas: zeros)
//     len(pa) == len(pb) == len(pc)
//   pairs[[j]] : Rcpp::List com
//     "ka" - integer global (0-based) das linhas do grupo 0 (ordem original)
//     "kb" - integer global das linhas do grupo 1 (ordem original)
//     "Ra" - vec, marginal sum_b c_ab por membro do grupo 0 (alinhada a ka)
//     "Cb" - vec, marginal sum_a c_ab por membro do grupo 1 (alinhada a kb)
//     "pa" - integer global da linha do membro 0 de cada par valido
//     "pb" - integer global da linha do membro 1 de cada par valido
//     "pc" - vec, peso c_ab por par valido
//     (pa/pb/pc vazios quando has_pairs = FALSE)
//   attr(pairs, "meta") : List com
//     "p"         - int, numero de colunas de z
//     "n2"        - double, n^2
//     "has_pairs" - logical, TRUE se pa/pb/pc foram construidos
//     "cdq_cols"  - logical por coluna, TRUE se Ra/Cb vieram do CDQ

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
// Helpers de meta (attr("meta") devolve proxy: materializa em List)
// ---------------------------------------------------------------------------

inline Rcpp::List get_meta(const Rcpp::List& obj) {
  if (!obj.hasAttribute("meta"))
    stop("fg_pairs: object without the 'meta' attribute; use make_fg_pairs_cpp");
  const Rcpp::List meta = obj.attr("meta");   // materializa o proxy UMA vez
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
// phi piecewise: sign(x) com 0 -> 0
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
// Agregacao pela norma L1/L2, padronizada por sqrt(p)
// ---------------------------------------------------------------------------

inline double norm_aggregate(const arma::vec& parts,
                             const std::string& norm_min) {
  const double q = std::sqrt(static_cast<double>(parts.n_elem));
  if (norm_min == "L1") {
    return arma::accu(arma::abs(parts)) / q;
  }
  return std::sqrt(arma::accu(arma::square(parts))) / q;
}

// Indices 0-based como integer do R (lidos por ponteiro na avaliacao)
inline Rcpp::IntegerVector to_int_idx(const arma::uvec& v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}

// Objeto sem o campo has_pairs: tratado como contendo pares
inline bool meta_has_pairs(const Rcpp::List& meta) {
  if (!meta.containsElementNamed("has_pairs")) return true;
  return Rcpp::as<bool>(meta["has_pairs"]);
}

// ===========================================================================
// Marginais Ra/Cb em O(N log^2 N) (phi = "identity")
//
// Par (a, b) valido no laco O(nk*nl):
//   a < b (posicional), dk[a] != 0, dl[b] != 0,
//   denom = wk[a]*wl[b] > 0 e finito,
//   hL - max(zbk[a], zbl[b]) > max(e1k[a], e1l[b]).
//
// IDENTIDADE EXATA EM PONTO FLUTUANTE: fl(hL - .) e monotona nao
// crescente (arredondamento monotono), logo
//   fl(hL - max(x, y)) == min(fl(hL - x), fl(hL - y)).
// Com A_a = hL - zbk[a], B_b = hL - zbl[b], u_a = e1k[a], v_b = e1l[b]:
//   min(A_a, B_b) > max(u_a, v_b)
//   <=> A_a > u_a  e  B_b > v_b          (filtros por elemento)
//       e A_a > v_b  e  B_b > u_a        (condicoes cruzadas)
// Todas as comparacoes sao sobre os MESMOS doubles do laco: o conjunto
// de pares e IDENTICO. So a ordem das somas muda (parcelas positivas,
// sem cancelamento; ~1e-15 relativo).
//
// Peso: c_ab = 1/(wk wl) ~ (1/wk)(1/wl), logo
//   Ra[a] = (1/wk[a]) * sum_{b > a, v_b < A_a, B_b > u_a} 1/wl[b]
//   Cb[b] = (1/wl[b]) * sum_{a < b, A_a > v_b, u_a < B_b} 1/wk[a]
// — contagens ponderadas de dominancia 3D, resolvidas por divisao e
// conquista na posicao (CDQ) + varredura + arvore de Fenwick.
//
// FALLBACK para o laco (retorna false) quando a identidade ou a
// fatoracao nao sao garantidas: valores nao finitos (std::max com NaN e
// assimetrico) ou pesos positivos fora de [1e-150, 1e150] (o produto
// poderia sair da faixa normal e mudar o filtro denom > 0 / finito).
// ===========================================================================

struct Fenwick {
  std::vector<double> t;
  explicit Fenwick(size_t n) : t(n + 1, 0.0) {}
  void add(size_t i, double v) {
    for (++i; i < t.size(); i += i & (~i + 1)) t[i] += v;
  }
  double sum(size_t k) const {          // soma dos ranks [0, k)
    double s = 0.0;
    for (; k > 0; k -= k & (~k + 1)) s += t[k];
    return s;
  }
};

struct MargData {
  // membros ATIVOS (passaram nos filtros por elemento), em ordem de posicao
  std::vector<int>    pa, pb;      // posicao dentro do grupo
  std::vector<double> A, u, wa;    // grupo 0: A_a, u_a, 1/wk[a]
  std::vector<double> B, v, wb;    // grupo 1: B_b, v_b, 1/wl[b]
  std::vector<double> SRa, SCb;    // somas acumuladas
};

// Contribuicoes cruzadas: a em [a0, a1) (metade esquerda das posicoes),
// b em [b0, b1) (metade direita) — todo par cruza com a < b.
inline void marg_cross(MargData& d, int a0, int a1, int b0, int b1) {
  const int na = a1 - a0, nb = b1 - b0;
  if (na == 0 || nb == 0) return;
  std::vector<int> qa(na), pbv(nb);
  for (int i = 0; i < na; ++i) qa[i] = a0 + i;
  for (int k = 0; k < nb; ++k) pbv[k] = b0 + k;

  // --- Ra: pontos b (v asc), consultas a (A asc); soma de wb com B > u
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
      // ranks com B > u_a (Bd decrescente)
      const size_t k = std::lower_bound(Bd.begin(), Bd.end(), d.u[q],
                                        std::greater<double>()) - Bd.begin();
      d.SRa[q] += fw.sum(k);
    }
  }

  // --- Cb: pontos a (A desc), consultas b (v desc); soma de wa com u < B
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
      // ranks com u < B_b (Ua crescente)
      const size_t k = std::lower_bound(Ua.begin(), Ua.end(), d.B[q])
                       - Ua.begin();
      d.SCb[q] += fw.sum(k);
    }
  }
}

// Blocos pequenos: forca bruta com a MESMA decomposicao das condicoes
inline void marg_brute(MargData& d, int a0, int a1, int b0, int b1) {
  for (int i = a0; i < a1; ++i) {
    for (int k = b0; k < b1; ++k) {
      if (d.pb[k] > d.pa[i] && d.A[i] > d.v[k] && d.B[k] > d.u[i]) {
        d.SRa[i] += d.wb[k];
        d.SCb[k] += d.wa[i];
      }
    }
  }
}

// Divisao e conquista sobre o intervalo de posicoes [lo, hi)
inline void marg_cdq(MargData& d, int lo, int hi,
                     int a0, int a1, int b0, int b1) {
  if (a0 == a1 || b0 == b1) return;
  if (hi - lo <= 1 ||
      static_cast<long long>(a1 - a0) * (b1 - b0) <= 1024) {
    marg_brute(d, a0, a1, b0, b1);
    return;
  }
  const int mid = lo + (hi - lo) / 2;
  const int am = static_cast<int>(
    std::lower_bound(d.pa.begin() + a0, d.pa.begin() + a1, mid) - d.pa.begin());
  const int bm = static_cast<int>(
    std::lower_bound(d.pb.begin() + b0, d.pb.begin() + b1, mid) - d.pb.begin());
  marg_cross(d, a0, am, bm, b1);          // a esquerda x b direita
  marg_cdq(d, lo, mid, a0, am, b0, bm);
  marg_cdq(d, mid, hi, am, a1, bm, b1);
}

// Preenche Ra/Cb; false => usar o laco O(nk*nl) (ver FALLBACK acima)
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
    // mesmo filtro que denom > 0 && finito, por elemento (w em (0, inf)).
    // w < 0: dois pesos negativos dariam denom > 0 no laco — a fatoracao
    // por elemento nao reproduz isso, entao forca o fallback.
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
  marg_cdq(d, 0, std::max(nk, nl), 0, static_cast<int>(d.pa.size()),
           0, static_cast<int>(d.pb.size()));

  for (size_t i = 0; i < d.pa.size(); ++i) Ra[d.pa[i]] = d.wa[i] * d.SRa[i];
  for (size_t k = 0; k < d.pb.size(); ++k) Cb[d.pb[k]] = d.wb[k] * d.SCb[k];
  return true;
}

} // namespace fgp

// ===========================================================================
// CONSTRUTOR DO OBJETO DE PARES
// ===========================================================================

// make_fg_pairs_cpp
//
// Constroi o objeto de pares pre-computados por coluna binaria de z.
// Os filtros de par valido (cabecalho) sao aplicados aqui, uma vez.
//
// Argumentos:
//   lt      - h1$linkfun(t), pre-computado
//   lq      - h2$linkfun(qt), pre-computado
//   z       - matriz de covariaveis binarias (n x p)
//   status  - indicador de evento (1 = evento)
//   beta    - coeficientes de sobrevivencia (fixos nesta etapa)
//   hL      - constante h1(L), pre-computada
//   km_surv - lista de comprimento p com w0/w1 por coluna (ordem original)
//   phi     - "identity" ou "piecewise". Com "identity" so Ra/Cb sao
//             construidos; "piecewise" constroi tambem pa/pb/pc
//             (o objeto piecewise serve aos dois phi na avaliacao)
//   algo    - "auto" (padrao): identity usa marginals_cdq, O(N log^2 N),
//             com fallback por coluna para o laco; "loop" forca o laco
//             O(nk*nl) (referencia para testes A/B). piecewise sempre
//             usa o laco (precisa dos pares explicitos).
//
// Retorna: lista por coluna + attr "meta" (ver schema no cabecalho);
// meta$cdq_cols indica, por coluna, se o caminho CDQ foi usado.
// [[Rcpp::export]]
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
  
  // Residuos e1 = lt - z * beta (mesma expressao do kernel de referencia)
  const arma::vec e1 = lt - z * beta;
  
  Rcpp::List out(p);
  Rcpp::LogicalVector cdq_cols(p, false);

  for (int j = 0; j < p; ++j) {
    
    // Grupos definidos pela covariavel binaria da coluna j
    const arma::uvec ik = arma::find(z.col(j) == 0.0);
    const arma::uvec il = arma::find(z.col(j) == 1.0);
    const int nk = static_cast<int>(ik.n_elem);
    const int nl = static_cast<int>(il.n_elem);
    
    // Coluna degenerada: sem pares, mas com as invariantes do schema
    // (Ra/Cb de zeros alinhados a ka/kb); parts(j) = 0. Use (n, fill::zeros):
    // arma::vec(0) e ambiguo (casaria com Col(const char*)).
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
    
    // dot por membro, uma vez (ver FIDELIDADE no cabecalho)
    arma::vec zbk(nk), zbl(nl);
    for (int a = 0; a < nk; ++a) zbk[a] = arma::dot(zk.row(a), beta);
    for (int b = 0; b < nl; ++b) zbl[b] = arma::dot(zl.row(b), beta);
    
    // Pesos KM do grupo (ordem original — alinhados com ik/il)
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
    
    // Acumuladores das saidas
    arma::vec  Ra(nk, arma::fill::zeros);
    arma::vec  Cb(nl, arma::fill::zeros);
    std::vector<int>    pa_v, pb_v;
    std::vector<double> pc_v;
    if (store_pairs) {
      pa_v.reserve(static_cast<size_t>(nk) * 2);   // chute conservador
      pb_v.reserve(static_cast<size_t>(nk) * 2);
      pc_v.reserve(static_cast<size_t>(nk) * 2);
    }

    // Marginais em O(N log^2 N); false => laco abaixo
    const bool used_cdq = try_cdq &&
      fgp::marginals_cdq(hL, zbk, zbl, e1k, e1l, dk, dl, wk, wl, Ra, Cb);
    cdq_cols[j] = used_cdq;

    // Laco duplo: b em [a + 1, nl) (a < b); membro a censurado e pulado,
    // pois nenhum par seu e valido.
    for (int a = 0; a < nk && !used_cdq; ++a) {
      if (dk[a] == 0.0) continue;

      for (int b = a + 1; b < nl; ++b) {
        if (dl[b] == 0.0) continue;
        
        const double denom = wk[a] * wl[b];
        if (!(denom > 0.0) || !std::isfinite(denom)) continue;
        
        // Indicador de truncagem: desigualdade estrita (ver FIDELIDADE)
        const double h1L_ind =
          (hL - std::max(zbk[a], zbl[b]) > std::max(e1k[a], e1l[b]))
          ? 1.0 : 0.0;
        
        if (h1L_ind == 0.0) continue;   // peso nulo: nao entra em nnz
        
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
// AVALIACAO
// ===========================================================================

// Nucleo compartilhado: componentes por coluna, retorna arma::vec.
// (parts nunca passa por NumericVector no caminho interno — sem
// arma::vec(begin, end), que nao existe.)
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
  
  // Residuos para TODAS as linhas, uma vez: O(n*p)
  const arma::vec e2 = lq - z * theta;
  
  arma::vec parts(p, arma::fill::zeros);
  
  if (phi == "identity") {
    // O(nk + nl) por coluna via marginais:
    //   parts(j) = -(1/n^2) [ sum_a e2k(a) Ra(a) - sum_b e2l(b) Cb(b) ]
    // Decomposicao algebrica exata da soma de pares; soma reordenada
    // (1e-10 do kernel de referencia, nao bit-exata). Colunas degeneradas
    // sao puladas (parts(j) = 0).
    for (int j = 0; j < p; ++j) {
      const Rcpp::List colj = pairs[j];
      
      const arma::uvec ka = Rcpp::as<arma::uvec>(colj["ka"]);
      const arma::uvec kb = Rcpp::as<arma::uvec>(colj["kb"]);
      if (ka.n_elem == 0 || kb.n_elem == 0) continue;   // degenerada
      
      const arma::vec Ra = Rcpp::as<arma::vec>(colj["Ra"]);
      const arma::vec Cb = Rcpp::as<arma::vec>(colj["Cb"]);
      
      const arma::vec e2k = e2.elem(ka);
      const arma::vec e2l = e2.elem(kb);
      
      parts[j] = -(arma::dot(e2k, Ra) - arma::dot(e2l, Cb)) / n2;
    }
  } else {
    // piecewise: O(nnz) sequencial sobre os pares congelados
    //   parts(j) = -(1/n^2) sum_k pc(k) * phi(e2[pa(k)] - e2[pb(k)])
    // pa/pb/pc lidos por ponteiro (IntegerVector/NumericVector sobre o
    // SEXP nao copiam quando o tipo ja e integer/double).
    if (!fgp::meta_has_pairs(meta))
      stop("eval_fg_parts: an object built with phi = 'identity' has no "
             "pairs; rebuild it with make_fg_pairs(..., 'piecewise')");
    const double* e2p = e2.memptr();
    for (int j = 0; j < p; ++j) {
      const Rcpp::List colj = pairs[j];

      const Rcpp::IntegerVector pa = colj["pa"];
      if (pa.size() == 0) continue;   // sem pares validos (ou degenerada)

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
// Componentes da U-estatistica a partir do objeto de pares.
// Retorna vetor numerico de comprimento p sem atributo dim.
// [[Rcpp::export]]
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
// Norma (L1/L2) das componentes, padronizada por sqrt(p). Interface
// compativel com U2_fun_fg_cpp para testes A/B.
// [[Rcpp::export]]
double eval_fg_obj_cpp(const Rcpp::List& pairs,
                       const arma::vec& lq,
                       const arma::mat& z,
                       const arma::vec& theta,
                       const std::string& phi,
                       const std::string& norm_min) {
  const arma::vec parts = fgp::eval_parts(pairs, lq, z, theta, phi);
  return fgp::norm_aggregate(parts, norm_min);
}
