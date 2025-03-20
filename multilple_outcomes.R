suppressMessages(library(survey))
n <- 1e5
devtools::load_all("survivl/")
forms <- list(W ~ 1,
              Z ~ Z_l1 + X_l1,
              X ~ X_l1 + Z_l0,
              list(Y ~ W + X_l0,
              D ~ W + X_l0),
              list(Y = list(Z ~ W),
                   D = list(Z ~ W)))
fams <- list(3, 1, 5, c(3,3), c(1))
pars <- list(W = list(beta=0, phi=1/2),
             Z = list(beta=c(0,0.7,0.2), phi=1),
             X = list(beta=c(-0.5,0.25,0.5)),
             Y = list(beta=c(1,-1/10,1/5), phi=1),
             D = list(beta = c(0, 0.5, 1), phi = 1),
             cop = list(Y=list(Z=list(beta=c(0.5, 0.0))),
                        D = list(Z = list(beta = c(-0.25, 0.3)))))
link <- list("log", "identity", "logit",c("identity", "identity"))

set.seed(123)
browser()
dat <- msm_samp(n, T=10, formulas = forms, family = fams, pars = pars,
                link = link)
datl <- surv_to_long(dat)

temp <- ipwtm(
  exposure = X,
  family = "binomial",
  numerator = ~ 1,
  denominator = ~ Z,
  link = "logit",
  id = id,
  timevar = t,
  type = "all",
  data = datl)


datl$wt <- temp$ipw.weights
glmY <- suppressWarnings(coxph(Surv(t, t_stop,Y)~ X, id = id, data = datl,
                               weights = datl$wt, timefix = FALSE))
