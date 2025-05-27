formulas <- list(list(),
                 Z ~ X_l1,
                 X ~ Z_l0,
                 list(Y ~ X_l0, Cen ~ X_l0),
                 cop ~ 1)
family <- list(integer(0), 1, 5, list(3,3), 1)

pars <- list(Z = list(beta=c(1,0.5), phi=1),
             X = list(beta=c(-0.25,0.5)),
             Y = list(beta=c(1.1,-0.2), phi=1),
             Cen = list(beta=c(1.1,-0.2), phi=1),
             cop = list(beta=1))

set.seed(123)
dat <- msm_samp(1e3, T=5,
                formulas = formulas,
                family = family,
                pars = pars,
                link = list(character(0), "identity", "logit", c("inverse", "inverse")))

test_that("status recorded correctly", {
  expect_true(all(dat$status - 2*rowSums(dat[,paste0("Y_", 0:4)], na.rm=TRUE) - rowSums(dat[,paste0("Cen_", 0:4)], na.rm=TRUE) == 0))
})
