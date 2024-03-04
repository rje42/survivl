suppressMessages(library(dplyr))

formulas <- list(Sex ~ 1,
                 Z ~ X_l1,
                 X ~ Z_l0,
                 Y ~ X_l0)
family <- list(5,1,5,3)
pars <- list(Sex = list(beta=0),
             Z = list(beta=c(1,0.5), phi=1),
             X = list(beta=c(-0.25,0.5)),
             Y = list(beta=c(0,1), lambda0=0.01))

# dat <- msm_samp(1e2,T=2, formulas, family, pars,
#                 link = list("logit", "identity", "logit", "inverse", "?"))
dat <- surv_samp(1e2,T=2, formulas, family, pars)
dat <- dat %>% mutate(Sex = factor(Sex, labels=c("female", "male")))
datl <- surv_to_long(dat)

test_that("surv_to_long preserves factors", {
  expect_equal(dat$Sex[dat$id], datl$Sex[match(dat$id, datl$id)])
})
