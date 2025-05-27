suppressMessages(library(survey))
n <- 1e5

forms <- list(W ~ 1,
              Z ~ Z_l1 + X_l1,
              X ~ X_l1 + Z_l0,
              Y ~ W + X_l0,
              list(Y = list(Z ~ W)))
fams <- list(3, 1, 5, 3, c(1))
pars <- list(W = list(beta=0, phi=1/2),
             Z = list(beta=c(0,0.7,0.2), phi=1),
             X = list(beta=c(-0.5,0.25,0.5)),
             Y = list(beta=c(1,-1/10,1/5), phi=1),
             cop = list(Y=list(Z=list(beta=c(0.5, 0.0))))
)
link <- list("log", "identity", "logit", "inverse")

set.seed(123)
dat <- msm_samp(n, T=5, formulas = forms, family = fams, pars = pars,
                link = link)
datl <- surv_to_long(dat, lag=1)

## keep W, first set, and then Z,X for t=1
vars <- c("W", paste0(c("Z","X","Y","Z","X"), "_", c(0,0,0,1,1)))
dat2 <- dat[vars]

## first test for static covariates
set.seed(123)
dati <- msm_samp(dat=dat2, T=5, formulas = forms, family = fams, pars = pars,
                link = link)
datil <- surv_to_long(dati, lag=1)

dati$T[dati$T == 0.5] <- dat$T[dati$T == 0.5]
# mean(dat$T)
# mean(dati$T)
#
# table(dat$status)
# table(dati$status)
#
# table(ceiling(dat$T-1))
# table(ceiling(dati$T-1))

# msm_samp(T=3, formulas=forms, family=fams, pars=pars, dat=df)

# summary(svyglm(Z ~ Z_l1 + X_l1, design = svydesign(~ 1, data=datl)))
modZ <- suppressWarnings(summary(svyglm(Z ~ Z_l1 + X_l1, design = svydesign(~ 1, data=datil[datil$t >= 1,])))$coefficients)
modXp <- suppressWarnings(svyglm(X ~ X_l1 + Z, family = binomial, design = svydesign(~ 1, data=datil[datil$t >= 1,])))
ps <- predict(modXp, type="response")
wt <- datil$X[datil$t >= 1]/ps + (1-datil$X[datil$t >= 1])/(1-ps)
modX <- summary(modXp)$coefficients
modY <- suppressWarnings(summary(svyglm(I(1-Y) ~ W + X, family = binomial(log), start = c(-1,1/10,-1/5), design = svydesign(~ 1, data=datil[datil$t >= 1,], weights = wt)))$coefficients)

test_that("plasmode simulation works", {
  expect_equal(dat[vars], dati[vars])
  expect_lt(max(abs((modZ[,1] - pars$Z$beta)/modZ[,2])), 2.5)
  expect_lt(max(abs((modX[,1] - pars$X$beta)/modX[,2])), 2.5)
  # expect_lt(sum(((-modY[-1,1] - pars$Y$beta[-1])/modY[-1,2])^2), qchisq(0.99, df=2))  # intercept is wrong!
})
