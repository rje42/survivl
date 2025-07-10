suppressWarnings(library(survey))
suppressWarnings(library(survival))
suppressWarnings(library(ipw))
gamma1 <- 0.1
gamma2 <- 0.15
gamma3 <- 0.1
gamma4 <- 0.3
gamma5 <- 0.3

rho_to_beta <- function(rho){
  x <- (rho + 1)/2
  return(log(x/(1-x)))
}


n <- 1e4


dat <- data.frame(X1=rnorm(n,mean = 0, sd = 1),
                  X2 = rbinom(n, size = 1, prob = 0.5))
dat$B1 <- rnorm(qtls$B1, 0.4*dat$X2 -0.2, sd = 1)

dat[["Z1_0"]] <- rnorm(n, 0.2*dat$X1, sd = 1)
dat[["Z2_0"]] <- rbinom(n, 1, expit(-0.2 + 0.4 * dat$X2))
dat[["A_0"]] = rbinom(n, 1, expit(-1 + 0.1 * dat$X1 + 0.15*dat$X2 + 0.1*dat$B1 + 0.3*dat[["Z1_0"]] + 0.3* dat[["Z2_0"]]))


# test where B2 is missing, but other baslines are there. 
formulas <- list(list(B2 ~ X1),
                 list(Z1 ~  B2 + Z1_l1 + A_l1, Z2 ~ B2 + Z2_l1 + A_l1),
                 A ~ X1 + X2 + B1 + Z1_l0 + Z2_l0 + A_l1,
                 Y ~ A_l0,
                 list(Z1 ~ 1, Z2 ~ 1))

family <- list(list(1),
               list(1,5),
               5,
               3,
               list(2,2))
pars <- list(B2 = list(beta = c(-0.2, 0.4), phi = 1),
             Z1 = list(beta=c(0.3,0.4, 0.7,-0.6 ), phi=1),
             Z2 = list(beta = c(-0.2,0.4,1,-0.6)),
             A = list(beta=c(-1, gamma1, gamma2, gamma3, gamma4, gamma5, 1)),
             Y = list(beta = c(0.5,0.2), phi = 1), # negative 1/exp
             cop = list(Y = list(Z1 = list(beta= rho_to_beta(-0.6), par2 = 5),
                                 Z2 = list(beta = rho_to_beta(0.2), par2 = 5))))

sm <- survivl_model(formulas=formulas, family=family,
                    pars=pars, T = 5, dat = dat, qtls = NULL, # no quantiles given
                    link = list(list("identity"),
                    list("identity", "logit"),
                    "logit", "log"))
dat2 <- rmsm(n, sm)

datl <- surv_to_long(dat2)


temp <- ipwtm(
  exposure = A,
  family = "binomial",
  numerator = ~ 1,
  denominator = ~ Z1 +Z2,
  link = "logit",
  id = id,
  timevar = t,
  type = "all",
  data = datl)


datl$wt <- temp$ipw.weights
glmY <- suppressWarnings(coxph(Surv(t, t_stop,Y)~ A, id = id, data = datl,
                               weights = datl$wt, timefix = FALSE))
sumY <- summary(glmY)



test_that("msm_samp() works as expected for time-to-event", {
  # expect_lt(max(abs(sumZ$coefficients[,1] - c(0, 0.7, 0.2))/sumZ$coefficients[,2]), 2.5)
  # expect_lt(max(abs(sumX$coefficients[,1] - c(-0.5, 0.25, 0.5))/sumX$coefficients[,2]), 2.5)
  expect_lt(max(abs(sumY$coefficients[,1] + c(-1/5))/sumY$coefficients[,2]), 2.5)
})

