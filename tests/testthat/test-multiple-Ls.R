
n <- 1e4

formulas <- list(W ~ 1,
                 list(Z1 ~ Z1_l1 + X_l1, Z2 ~ W + Z2_l1),
                 X ~ X_l1 + Z1_l0 + Z2_l0,
                 Y ~ W + X_l0,
                 list(Y = list(Z1 ~ 1, Z2 ~Z1_l1)))
family <- list(3, c(1,1), 5, 3, c(1, 1))
pars <- list(
  W = list(beta=0, phi=1/2),
  Z1 = list(beta=c(0,0.7,0.2), phi=1),
  Z2 = list(beta = c(0, 0.3, 0.45), phi = 1),
  X = list(beta=c(-0.5,0.25,0.5, 0.25)),
  Y = list(beta=c(1,-1/10,1/4), phi=1),
  cop = list(Y=list(Z1=list(beta=0.5), 
                    Z2 = list(beta = c(-0.5, 0.25)))))

link <- list("log", c("identity","identity"),"logit", "log")

set.seed(123)

sm <- survivl_model(formulas=formulas, family=family,
                    pars=pars, T = 7,
                    link = link)

dat <- rmsm(n,sm) # changed to 7 too many people died out
datl <- surv_to_long(dat)


temp <- ipw::ipwtm(
  exposure = X,
  family = "binomial",
  numerator = ~ 1,
  denominator = ~ Z1 +Z2,
  link = "logit",
  id = id,
  timevar = t,
  type = "all",
  data = datl)
datl$wt <- temp$ipw.weights

coxY <- coxph(Surv(t, t_stop, Y) ~ W + X,data = datl, weights = wt,
              id = id, timefix = FALSE)
sumY <- summary(coxY)

test_that("multiple Ls works as expected for time-to-event", {
  expect_lt(max(abs(sumY$coefficients[,1] + c(-1/10, 1/4))/sumY$coefficients[,2]), 2.5)
})
