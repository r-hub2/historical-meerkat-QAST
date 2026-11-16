# QAST <img src="man/figures/logo.png" align="right" height="139" alt="QAST hex logo" />

*[English version below](#english)*

O **QAST** implementa a análise de regressão conjunta do **tempo de
sobrevivência** e do **tempo de sobrevivência ajustado por qualidade de vida**
(QALY) com modelos semiparamétricos bivariados de tempo de falha acelerado
(AFT, *accelerated failure time*).

Os parâmetros são estimados por equações de estimação baseadas em postos, com
pesos de probabilidade inversa de censura (IPCW) e estimação da variância por
jackknife ou bootstrap. Os núcleos numéricos são escritos em C++
(Rcpp/RcppArmadillo).

## Modelos

| Função | Modelo | Covariáveis |
|---|---|---|
| `semi_aft_prop()` | Modelo proposto na tese de doutorado de da Silva (2026, IME-USP) | contínuas ou discretas |
| `semi_aft_fg()` | Fine e Gelber (2001) | binárias (0/1) |
| `surv_qt_fun()` | Estimador de Zhao e Tsiatis (1997) da curva de sobrevivência do tempo ajustado por qualidade | — |

Nos dois modelos de regressão, `beta` descreve o tempo de sobrevivência e
`theta` o tempo ajustado por qualidade.

## Instalação

O pacote ainda não está no CRAN. Para instalar a partir do Codeberg (requer um
compilador C++):

```r
# install.packages("remotes")
remotes::install_git("https://codeberg.org/landsilva/QAST")
```

### Compilador C++

O QAST é compilado na instalação, então o R precisa encontrar as
ferramentas de compilação do sistema:

- **Windows:** instale o [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  correspondente à sua versão do R (por exemplo, Rtools 4.4 para o R 4.4.x).
- **macOS:** instale as ferramentas de linha de comando do Xcode no Terminal
  com `xcode-select --install` e o compilador Fortran indicado em
  [mac.r-project.org/tools](https://mac.r-project.org/tools/), exigido pelo
  RcppArmadillo.
- **Linux:** instale o compilador e os cabeçalhos de desenvolvimento do R:
  - Debian/Ubuntu: `sudo apt install build-essential r-base-dev`
  - Fedora: `sudo dnf install gcc-c++ gcc-gfortran R-devel`
  - Arch e derivadas (Manjaro, EndeavourOS):
    `sudo pacman -S --needed base-devel gcc-fortran` (os cabeçalhos do R
    já vêm no pacote `r`)

Para conferir se o R encontra o compilador:

```r
# install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)
```

## Uso

```r
library(QAST)

set.seed(2026)
n <- 200
d <- data.frame(trt = rbinom(n, 1, 0.5), sexo = rbinom(n, 1, 0.4),
                idade = rnorm(n))
t <- rexp(n, rate = exp(-0.4 * d$trt))   # tempo de sobrevivência
c <- rexp(n, rate = 0.3)                 # censura
d$status <- as.integer(t <= c)
d$t <- pmin(t, c)
d$qt <- d$t * runif(n, 0.5, 1)           # tempo ajustado por qualidade
```

A resposta é escrita com `Qsurv(tempo, tempo_ajustado, status)`, no lado
esquerdo da fórmula, como o `Surv()` do pacote `survival`. Fatores entram
como indicadoras, sem intercepto.

### Modelo de Fine e Gelber

```r
fit_fg <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sexo, data = d,
                      L = quantile(d$t, 0.9), var_comp = TRUE)
fit_fg            # coeficientes dos dois componentes
summary(fit_fg)   # erros-padrão, estatísticas z e p-valores
```

As covariáveis devem ser binárias. Com `phi = "identity"` e
`norm_min = "L2"` (os padrões), `theta` é obtido em forma fechada
(`fast_path = TRUE`); `method` é usado para `beta`.

### Modelo proposto

```r
fit_prop <- semi_aft_prop(Qsurv(t, qt, status) ~ idade + trt, data = d,
                          var_comp = TRUE)
coef(fit_prop)            # matriz com as colunas beta e theta
vcov(fit_prop, "theta")   # covariância (jackknife) do componente theta
```

Os dois modelos também aceitam vetores e uma matriz de covariáveis:
`semi_aft_prop(t, qt, status, z)`. Os padrões são `h1 = h2 = "log"`,
`method = "Nelder-Mead"`, `norm_min = "L2"` e `var_comp = FALSE`; com
`var_method = "bootstrap"`, informe o número de réplicas em `M`.

### Curva de sobrevivência do tempo ajustado por qualidade

```r
s_qt <- surv_qt_fun(qaly = d$qt, time = d$t, status = d$status)
```

Os ajustes devolvem um objeto de classe `qast_fit`, com `beta`, `theta` e,
se `var_comp = TRUE`, as variâncias e matrizes de covariância. As entradas
são validadas: `NA`, tempos não positivos, `qt > t`, `status` fora de
{0, 1} e comprimentos diferentes param com uma mensagem que indica as
observações problemáticas.

## Observações

- **Sem intercepto em `z`.** As equações de estimação são invariantes a
  deslocamentos, então o modelo não estima intercepto: colunas constantes são
  rejeitadas com erro. Colunas linearmente dependentes geram um aviso, pois os
  coeficientes nessas direções não são identificáveis.
- **Paralelismo.** O jackknife e o bootstrap usam os núcleos físicos da máquina
  (via fork; no Windows, execução sequencial). Para fixar o número de
  processos: `options(QAST.cores = 4)`.
- **Otimizadores.** `method` aceita `"uobyqa"`, `"mads"`, `"Nelder-Mead"`,
  `"subplex"`, `"BB"`, `"SANN"`, `"bobyqa"`, `"newuoa"`, `"hjkb"` e `"BFGS"`.

## Documentação

- `vignette("QAST")`: introdução ao método e exemplo completo (em inglês);
  requer instalar com `remotes::install_git(..., build_vignettes = TRUE)`.
- `citation("QAST")`: como citar o pacote e o modelo de Fine e Gelber.
- [NEWS.md](NEWS.md): histórico de versões.

## Referências

da Silva, A. (2026). *Regression analysis of quality-adjusted survival time
in cancer patients admitted to intensive care units: a new estimation method*.
Tese de doutorado, Instituto de Matemática e Estatística, Universidade de São
Paulo. Orientação: Gisela Tunes da Silva.

Fine, J. P. e Gelber, R. D. (2001). Joint regression analysis of survival and
quality-adjusted survival. *Biometrics*, 57(2), 376–382.
<https://doi.org/10.1111/j.0006-341X.2001.00376.x>

Zhao, H. e Tsiatis, A. A. (1997). A consistent estimator for the distribution
of quality adjusted survival time. *Biometrika*, 84(2), 339–348.
<https://doi.org/10.1093/biomet/84.2.339>

## Licença

GPL (>= 3). Ver [LICENSE.md](https://codeberg.org/landsilva/QAST/src/branch/main/LICENSE.md).

---

## English

**QAST** implements joint regression analysis of **survival time** and
**quality-adjusted survival time** (QALY) with bivariate semiparametric
accelerated failure time (AFT) models.

Parameters are estimated through rank-based estimating equations, with
inverse probability of censoring weights (IPCW) and variance estimation by
jackknife or bootstrap. The numerical kernels are written in C++
(Rcpp/RcppArmadillo).

### Models

| Function | Model | Covariates |
|---|---|---|
| `semi_aft_prop()` | Model proposed in the doctoral thesis of da Silva (2026, IME-USP) | continuous or discrete |
| `semi_aft_fg()` | Fine and Gelber (2001) | binary (0/1) |
| `surv_qt_fun()` | Zhao and Tsiatis (1997) estimator of the quality-adjusted survival curve | — |

In both regression models, `beta` describes survival time and `theta` the
quality-adjusted survival time.

### Installation

QAST is not on CRAN yet. To install it from Codeberg (a C++ compiler is
required):

```r
# install.packages("remotes")
remotes::install_git("https://codeberg.org/landsilva/QAST")
```

#### C++ compiler

QAST is compiled at installation, so R must find the system build tools:

- **Windows:** install the [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  version that matches your R (e.g. Rtools 4.4 for R 4.4.x).
- **macOS:** install the Xcode command line tools from the Terminal with
  `xcode-select --install`, plus the Fortran compiler listed at
  [mac.r-project.org/tools](https://mac.r-project.org/tools/), which
  RcppArmadillo requires.
- **Linux:** install the compiler and R's development headers:
  - Debian/Ubuntu: `sudo apt install build-essential r-base-dev`
  - Fedora: `sudo dnf install gcc-c++ gcc-gfortran R-devel`
  - Arch and derivatives (Manjaro, EndeavourOS):
    `sudo pacman -S --needed base-devel gcc-fortran` (R's headers ship with
    the `r` package)

To check that R finds the compiler:

```r
# install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)
```

### Usage

The response is written as `Qsurv(time, adjusted_time, status)` on the
left-hand side of a formula, like `survival::Surv()`; factors enter as
indicators, without an intercept. The data `d` are the ones simulated in the
Portuguese section above.

```r
library(QAST)

# Fine and Gelber (2001): binary covariates
fit_fg <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sexo, data = d,
                      L = quantile(d$t, 0.9), var_comp = TRUE)
fit_fg            # coefficients of both components
summary(fit_fg)   # standard errors, z statistics and p-values

# Proposed model: continuous or discrete covariates
fit_prop <- semi_aft_prop(Qsurv(t, qt, status) ~ idade + trt, data = d,
                          var_comp = TRUE)
coef(fit_prop)            # matrix with columns beta and theta
vcov(fit_prop, "theta")   # (jackknife) covariance of the theta component

# Quality-adjusted survival curve (Zhao and Tsiatis, 1997)
s_qt <- surv_qt_fun(qaly = d$qt, time = d$t, status = d$status)
```

Both models also accept vectors and a covariate matrix:
`semi_aft_prop(t, qt, status, z)`. Defaults are `h1 = h2 = "log"`,
`method = "Nelder-Mead"`, `norm_min = "L2"` and `var_comp = FALSE`; with
`var_method = "bootstrap"`, set the number of replicates in `M`.

Fits return an object of class `qast_fit` with `beta`, `theta` and, when
`var_comp = TRUE`, the variances and covariance matrices. Inputs are
validated: `NA`, non-positive times, `qt > t`, `status` outside {0, 1} and
mismatched lengths stop with a message pointing to the offending
observations.

### Notes

- **No intercept in `z`.** The estimating equations are invariant to shifts,
  so the models do not estimate an intercept: constant columns are rejected
  with an error. Linearly dependent columns trigger a warning, since the
  coefficients along those directions are not identifiable.
- **Fine and Gelber covariates** must be binary. With `phi = "identity"` and
  `norm_min = "L2"` (the defaults), `theta` has a closed form
  (`fast_path = TRUE`) and `method` is used for `beta` only.
- **Parallelism.** Jackknife and bootstrap use the physical cores (via fork;
  sequential on Windows). To set the number of processes:
  `options(QAST.cores = 4)`.
- **Optimizers.** `method` accepts `"uobyqa"`, `"mads"`, `"Nelder-Mead"`,
  `"subplex"`, `"BB"`, `"SANN"`, `"bobyqa"`, `"newuoa"`, `"hjkb"` and
  `"BFGS"`.

### Documentation

- `vignette("QAST")`: introduction to the method and a complete example;
  requires installing with `remotes::install_git(..., build_vignettes = TRUE)`.
- `citation("QAST")`: how to cite the package and the Fine and Gelber model.
- [NEWS.md](NEWS.md): version history.

### References

da Silva, A. (2026). *Regression analysis of quality-adjusted survival time
in cancer patients admitted to intensive care units: a new estimation method*.
Doctoral thesis, Institute of Mathematics and Statistics, University of São
Paulo. Advisor: Gisela Tunes da Silva.

Fine, J. P. and Gelber, R. D. (2001). Joint regression analysis of survival
and quality-adjusted survival. *Biometrics*, 57(2), 376–382.
<https://doi.org/10.1111/j.0006-341X.2001.00376.x>

Zhao, H. and Tsiatis, A. A. (1997). A consistent estimator for the
distribution of quality adjusted survival time. *Biometrika*, 84(2), 339–348.
<https://doi.org/10.1093/biomet/84.2.339>

### License

GPL (>= 3). See [LICENSE.md](https://codeberg.org/landsilva/QAST/src/branch/main/LICENSE.md).
