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



formulas <- list(list(X1 ~ 1, X2 ~ 1,B1 ~ X2, B2 ~ X1),
                 list(Z1 ~  B2 + Z1_l1 + A_l1, Z2 ~ B2 + Z2_l1 + A_l1),
                 A ~ X1 + X2 + B1 + Z1_l0 + Z2_l0 + A_l1,
                 Y ~ A_l0,
                 list(copZ1 ~ 1, copZ2 ~ 1))

family <- list(list(1,5,1,1),
               list(1,5),
               5,
               3,
               list(2,2))
pars <- list(X1 = list(beta = 0, phi = 1),
             X2 = list(beta = 0), #expit(0) = 0.5
             B1 = list(beta = c(-0.2,0.4),phi = 1),
             B2 = list(beta = c(0, 0.2), phi = 1),

             Z1 = list(beta=c(0.3,0.4, 0.7,-0.6 ), phi=1),
             Z2 = list(beta = c(-0.2,0.4,1,-0.6)),
             A = list(beta=c(-1, gamma1, gamma2, gamma3, gamma4, gamma5, 1)),
             Y = list(beta = c(0.5,0.2), phi = 1), # negative 1/exp
             copZ1 = list(beta= rho_to_beta(-0.6), par2 = 5),
             copZ2 = list(beta = rho_to_beta(0.2), par2 = 5))

dat <- msm_samp(n=1e5, T = 10,
                formulas = formulas,
                family = family,
                pars = pars,
                link = list(list("identity", "logit", "identity", "identity"),
                            list("identity", "logit"),
                            "logit", "log"))

# dat <- as.data.table(dat)
datl <- surv_to_long(dat)


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

