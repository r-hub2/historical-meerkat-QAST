// src/surv_qt_core.cpp
//
// Estimador IPCW de Zhao e Tsiatis (1997) para a curva de sobrevivencia do
// tempo ajustado por qualidade de vida (QALY). A documentacao vive no
// wrapper R/surv_qt_wrappers.R.
//
// pred_surv[i] = (1/n) sum_{j : qalys[j] >= qalys[i]} status[j] / G(t_j),
// calculado por uma ordenacao O(n log n) e uma soma acumulada O(n), com
// empates agrupados. G e o KM da censura (survfit em R), avaliado em cada
// tempo por busca binaria.
//
// Convencoes:
// - tempo alem do ultimo tempo do KM: G = NaN e o peso e tratado como 0;
// - tempo antes do primeiro tempo do KM: G = 1;
// - QALY <= 0 com preditor linear: erro (log indefinido).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace survqt_internal {

// Avaliacao da step function KM em um tempo arbitrario.
// Retorna: 1.0 antes do primeiro tempo de evento;
//          NaN depois do ultimo tempo (o peso vira 0 em surv_qt_cpp);
//          km_surv[k] com km_time[k] <= t < km_time[k+1], caso contrario.
inline double km_eval(const double t,
                      const NumericVector& km_time,
                      const NumericVector& km_surv)
{
  const int m = static_cast<int>(km_time.size());
  
  if (m == 0) return 1.0;
  if (t < km_time[0]) return 1.0;
  if (t > km_time[m - 1]) return std::numeric_limits<double>::quiet_NaN();
  
  // Busca binaria: maior k tal que km_time[k] <= t
  int lo = 0, hi = m - 1;
  while (lo < hi) {
    const int mid = (lo + hi + 1) / 2;  // arredonda para cima
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
// Nucleo do estimador de sobrevivencia para QALY:
//
//   pred_surv[i] = (1/n) * sum_{j : qalys[j] >= qalys[i]} (status[j] / cens_surv[j])
//
// Recebe cens_surv JA interpolado (use surv_qt_with_km_cpp para calcular
// a interpolacao a partir dos vetores do survfit).
//
// Argumentos:
//   qaly       - vetor de tempos ajustados por qualidade (na escala crua)
//   time       - tempos observados (usados apenas para validacao de tamanho)
//   status     - indicador de evento (1 = evento, 0 = censura)
//   cens_surv  - sobrevivencia da censura avaliada em cada observacao
//   lin_pred   - opcional: predicao linear (X %*% theta); se fornecida,
//                qalys = log(qaly) - lin_pred
//
// Retorna vetor numerico de comprimento n. Uso interno (chamado por
// surv_qt_with_km_cpp); nao exportado para o R.
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
  // Transformacao: qalys = qaly (crus) ou log(qaly) - lin_pred
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
  // Pesos: status / cens_surv
  // NaN, zero ou nao finito em cens_surv -> peso 0
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
  // Estimador via ordenacao decrescente + sweep acumulado
  //
  // Ordenar os indices por qalys decrescente faz com que, ao processar o
  // elemento i, a soma acumulada contenha exatamente os pesos dos
  // elementos com qalys >= qalys[i] (com ties agrupados).
  // -----------------------------------------------------------------------
  IntegerVector idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;
  
  std::sort(idx.begin(), idx.end(),
            [&qalys](int a, int b) {
              if (qalys[a] != qalys[b]) return qalys[a] > qalys[b];
              return a < b;  // desempate deterministico por indice
            });
  
  NumericVector pred_surv(n);
  const double dn = static_cast<double>(n);
  
  int j = 0;
  double acc = 0.0;
  while (j < n) {
    const double current = qalys[idx[j]];
    
    // Acumula todos os elementos do grupo de empate atual
    while (j < n && qalys[idx[j]] == current) {
      acc += weight[idx[j]];
      ++j;
    }
    
    // Todos os membros do grupo de empate recebem a mesma probabilidade
    // (comparacao >=)
    for (int k = j - 1; k >= 0 && qalys[idx[k]] == current; --k) {
      pred_surv[idx[k]] = acc / dn;
    }
  }
  
  return pred_surv;
}

// surv_qt_with_km_cpp
//
// Versao de conveniencia: interpola a step function KM (vetores
// km_time/km_surv pre-computados via survfit em R) nos tempos
// observados e delega ao nucleo surv_qt_cpp.
//
// Argumentos:
//   qaly      - vetor de tempos ajustados por qualidade (escala crua)
//   time      - tempos observados
//   status    - indicador de evento (1 = evento, 0 = censura)
//   km_time   - tempos do survfit KM da censura (ordenados)
//   km_surv   - sobrevivencias KM correspondentes
//   lin_pred  - opcional: predicao linear (X %*% theta)
//
// Retorna vetor numerico de comprimento n.
// [[Rcpp::export]]
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
  
  // Interpolacao do KM em cada tempo observado
  NumericVector cens_surv(n);
  for (int i = 0; i < n; ++i) {
    cens_surv[i] = survqt_internal::km_eval(time[i], km_time, km_surv);
  }
  
  return surv_qt_cpp(qaly, time, status, cens_surv, lin_pred);
}
