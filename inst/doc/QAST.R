## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(QAST.cores = 2)

## ----data---------------------------------------------------------------------
library(QAST)

set.seed(2026)
n <- 150
d <- data.frame(
  trt = rbinom(n, 1, 0.5),
  sex = rbinom(n, 1, 0.4),
  age = round(rnorm(n, 60, 8))
)
t_event <- rexp(n, rate = exp(-0.4 * d$trt))   # treatment lengthens survival
t_cens  <- rexp(n, rate = 0.3)
d$status <- as.integer(t_event <= t_cens)
d$t  <- pmin(t_event, t_cens)
d$qt <- d$t * runif(n, 0.5, 1)                  # utility between 0.5 and 1
head(d)

## ----fg-----------------------------------------------------------------------
fit_fg <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
                      L = quantile(d$t, 0.9), var_comp = TRUE)
fit_fg

## ----fg-summary---------------------------------------------------------------
summary(fit_fg)

## ----fg-cens------------------------------------------------------------------
fit_fg1 <- semi_aft_fg(Qsurv(t, qt, status) ~ trt + sex, data = d,
                       L = quantile(d$t, 0.9), cens_group = d$trt)
coef(fit_fg1)

## ----prop---------------------------------------------------------------------
fit_prop <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = d,
                          var_comp = TRUE)
coef(fit_prop)
vcov(fit_prop, "theta")

## ----intercept----------------------------------------------------------------
fit_prop$intercept
# Predicted log survival time and log quality-adjusted time
z <- model.matrix(~ age + trt, d)[, -1]
head(cbind(
  log_t  = fit_prop$intercept[["beta"]]  + drop(z %*% fit_prop$beta),
  log_qt = fit_prop$intercept[["theta"]] + drop(z %*% fit_prop$theta)
))

## ----predict, fig.width = 6, fig.height = 4-----------------------------------
new <- data.frame(age = c(0, 0), trt = c(0, 1))
predict(fit_prop, new)
tg <- seq(0.1, 3, length.out = 60)
s <- predict(fit_prop, new, type = "survival", times = tg)
matplot(tg, t(s$theta), type = "s", lty = 1:2, col = 1, ylim = c(0, 1),
        xlab = "Quality-adjusted survival time",
        ylab = "Estimated survival probability")
legend("topright", c("trt = 0", "trt = 1"), lty = 1:2, bty = "n")

## ----rmean--------------------------------------------------------------------
predict(fit_prop, new, type = "rmean", tau = 2)

## ----metrics------------------------------------------------------------------
train <- d[1:100, ]
test <- d[101:150, ]
fit_tr <- semi_aft_prop(Qsurv(t, qt, status) ~ age + trt, data = train)
m <- qast_metrics(fit_tr, newdata = test,
                  y = with(test, Qsurv(t, qt, status)))
m

## ----pp, fig.width = 7, fig.height = 3.6--------------------------------------
pp <- qast_pp(fit_prop, group = cut(d$age, c(-Inf, median(d$age), Inf),
                                    labels = c("younger", "older")))
pp
plot(pp)

## ----surv-qt, fig.width = 6, fig.height = 4-----------------------------------
s_qt <- surv_qt_fun(qaly = d$qt, time = d$t, status = d$status)
ord <- order(d$qt)
plot(d$qt[ord], s_qt[ord], type = "s", ylim = c(0, 1),
     xlab = "Quality-adjusted survival time",
     ylab = "Estimated survival probability")

