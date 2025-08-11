#longitudinal-no-survivl
n <- 1e4
formulas <- list(W ~ 1,
                 Z ~ Z_l1 + X_l1,
                 X ~ X_l1 + Z_l0,
                 Y ~ W + X_l0,
                 list(Y ~ Z_l0))
family <- list(3, 1, 5, 3, c(1))
pars <- list(W = list(beta=0, phi=1/2),
             Z = list(beta=c(0,0.7,0.2), phi=1),
             X = list(beta=c(-0.5,0.25,0.5)),
             Y = list(beta=c(0,-0.45, -0.2), phi=1),
             cop = list(Y=list(rho = -0.8))
)
link <- list("log", "identity", "logit", "log")

set.seed(123)

sm <- survivl_model(formulas=formulas, family=family,
                    pars=pars, T = 7,
                    link = link, method = "bootstrap",
                    control = list(bootsims = 2.5e2))

dat <- rmsm(n,sm) # changed to 7 too many people died out
datl <- surv_to_long(dat)

temp <- ipw::ipwtm(
  exposure = X, 
  family = "binomial",
  link = "logit",
  numerator = ~ W,
  denominator = ~Z,
  id = id,
  timevar = t,
  type = "all",
  data = datl
)
datl$wt <- temp$ipw.weights


sumY <- summary(survreg(Surv(t_stop, Y) ~ W + X, 
                        dist = "exponential",
                        data = datl,
                        weights = wt))
test_that("sim_seaman_keogh works correctly for continuous exponential", {
  expect_lt(max(abs(sumY$table[3,1] - 1/5) /sumY$table[3,2]), 2.5)
})
