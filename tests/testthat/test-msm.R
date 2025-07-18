suppressMessages(library(survey))
n <- 1e4

formulas <- list(W ~ 1,
              Z ~ Z_l1 + X_l1,
              X ~ X_l1 + Z_l0,
              Y ~ W + X_l0,
              list(Y = list(Z ~ W)))
family <- list(3, 1, 5, 3, c(1))
pars <- list(W = list(beta=0, phi=1/2),
             Z = list(beta=c(0,0.7,0.2), phi=1),
             X = list(beta=c(-0.5,0.25,0.5)),
             Y = list(beta=c(1,-1/10,1/5), phi=1),
             cop = list(Y=list(Z=list(beta=c(0.5, 0.0))))
)
link <- list("log", "identity", "logit", "log")

set.seed(123)

sm <- survivl_model(formulas=formulas, family=family,
                    pars=pars, T = 7,
                    link = link)

dat <- rmsm(n,sm) # changed to 7 too many people died out
datl <- surv_to_long(dat, lag=1)

glmZ <- glm(Z ~ Z_l1 + X_l1, data=datl)
sumZ <- summary(glmZ)
glmX <- glm(X ~ X_l1 + Z, data=datl, family=binomial)
sumX <- summary(glmX)
ps <- predict(glmX, type="response")
wt <- datl$X/ps + (1-datl$X)/(1-ps)



coxY <- coxph(Surv(t, t_stop, Y) ~ W + X,data = datl, weights = wt,
              id = id, timefix = FALSE)
sumY <- summary(coxY)



test_that("msm_samp() works as expected", {
  expect_lt(max(abs(sumZ$coefficients[,1] - c(0, 0.7, 0.2))/sumZ$coefficients[,2]), 2.5)
  expect_lt(max(abs(sumX$coefficients[,1] - c(-0.5, 0.25, 0.5))/sumX$coefficients[,2]), 2.5)
  expect_lt(max(abs(sumY$coefficients[,1] + c(-1/10, 1/5))/sumY$coefficients[,3]), 2.5)
})

