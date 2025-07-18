#longitudinal-no-survivl

## Two Examples

formulas <- list(C ~ 1,
                 Z ~ X_l1 + C,
                 X ~ Z_l0 + C,
                 Y ~ X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + X_6 + C,
                 cop ~ 1)
family <- list(5,1,5,3,1)
link <- list("logit", "identity", "logit", "inverse")
do_pars <- c(0.05,0.5, 0.25, 0, 0.34, 0.25, 0.45, 0.2, 0.05) # needs to be of length 1 + 7 + 1 (intercept, X, C)
pars <- list(C = list(beta=0),
             Z = list(beta = c(-1/2,1/2,0.25), phi=0.5),
             X = list(beta = c(0,1/2,1/10)),
             Y = list(beta = do_pars, phi=1),
             cop = list(beta=0.8472979))  # gives correlation 0.4


set.seed(123)
n <- 1e4

surv_model <- survivl_model(T=7, formulas=formulas, family=family, pars=pars, link=link)

dat <- rmsm(n, surv_model)

glm_X0 <- glm(X_0 ~ Z_0, family = binomial(), data = dat)
glm_X1 <- glm(X_1 ~ X_0 + Z_0 + Z_1, family = binomial(), data = dat)
glm_X2 <- glm(X_2 ~ X_0 + X_1 + Z_0 + Z_1 + Z_2, family = binomial(), data = dat)
glm_X3 <- glm(X_3 ~ X_0 + X_1 + X_2 + Z_0 + Z_1 + Z_2 + Z_3, family = binomial(), data = dat)
glm_X4 <- glm(X_4 ~ X_0 + X_1 + X_2 + X_3 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4, family = binomial(), data = dat)
glm_X5 <- glm(X_5 ~ X_0 + X_1 + X_2 + X_3 + X_4 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4 + Z_5, family = binomial(), data = dat)
glm_X6 <- glm(X_6 ~ X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4 + Z_5 + Z_6, family = binomial(), data = dat)

ps_X0 <- predict(glm_X0, type = "response")
ps_X1 <- predict(glm_X1, type = "response")
ps_X2 <- predict(glm_X2, type = "response")
ps_X3 <- predict(glm_X3, type = "response")
ps_X4 <- predict(glm_X4, type = "response")
ps_X5 <- predict(glm_X5, type = "response")
ps_X6 <- predict(glm_X6, type = "response")

# Compute weights depending on actual values of X_t
w0 <- ifelse(dat$X_0 == 1, ps_X0, 1 - ps_X0)
w1 <- ifelse(dat$X_1 == 1, ps_X1, 1 - ps_X1)
w2 <- ifelse(dat$X_2 == 1, ps_X2, 1 - ps_X2)
w3 <- ifelse(dat$X_3 == 1, ps_X3, 1 - ps_X3)
w4 <- ifelse(dat$X_4 == 1, ps_X4, 1 - ps_X4)
w5 <- ifelse(dat$X_5 == 1, ps_X5, 1 - ps_X5)
w6 <- ifelse(dat$X_6 == 1, ps_X6, 1 - ps_X6)

dat$iptw_weights <- 1 / (w0 * w1 * w2 * w3 * w4 * w5 * w6)

msm_fit <- glm(Y ~ X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + X_6 + C, 
               data = dat, weights = iptw_weights, family = Gamma(link = "log"))


sum_fit <- summary(msm_fit)


test_that("rmsm works as expected for longitudinal data", {
  expect_lt(max(abs(sum_fit$coefficients[,1] - do_pars)/sum_fit$coefficients[,2]), 2.5)
})



# Example 2 cumulative link function much less parameters


formulas <- list(C ~ 1,
                 Z ~ X_l1 + C,
                 X ~ Z_l0 + C,
                 Y ~ I(X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + X_6) + C,
                 cop ~ 1)
family <- list(5,1,5,3,1)
link <- list("logit", "identity", "logit", "inverse")
do_pars <- c(0.05,0.45, 0.05)
pars <- list(C = list(beta=0),
             Z = list(beta = c(-1/2,1/2,0.25), phi=0.5),
             X = list(beta = c(0,1/2,1/10)),
             Y = list(beta = do_pars, phi=1), # needs to be of length 1 + 7 + 1 (intercept, X, C)
             cop = list(beta=0.8472979))  # gives correlation 0.4


set.seed(123)
n <- 1e4

surv_model <- survivl_model(T=7, formulas=formulas, family=family, pars=pars, link=link)

dat <- rmsm(n, surv_model)


glm_X0 <- glm(X_0 ~ Z_0, family = binomial(), data = dat)
glm_X1 <- glm(X_1 ~ X_0 + Z_0 + Z_1, family = binomial(), data = dat)
glm_X2 <- glm(X_2 ~ X_0 + X_1 + Z_0 + Z_1 + Z_2, family = binomial(), data = dat)
glm_X3 <- glm(X_3 ~ X_0 + X_1 + X_2 + Z_0 + Z_1 + Z_2 + Z_3, family = binomial(), data = dat)
glm_X4 <- glm(X_4 ~ X_0 + X_1 + X_2 + X_3 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4, family = binomial(), data = dat)
glm_X5 <- glm(X_5 ~ X_0 + X_1 + X_2 + X_3 + X_4 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4 + Z_5, family = binomial(), data = dat)
glm_X6 <- glm(X_6 ~ X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + Z_0 + Z_1 + Z_2 + Z_3 + Z_4 + Z_5 + Z_6, family = binomial(), data = dat)

ps_X0 <- predict(glm_X0, type = "response")
ps_X1 <- predict(glm_X1, type = "response")
ps_X2 <- predict(glm_X2, type = "response")
ps_X3 <- predict(glm_X3, type = "response")
ps_X4 <- predict(glm_X4, type = "response")
ps_X5 <- predict(glm_X5, type = "response")
ps_X6 <- predict(glm_X6, type = "response")

# Compute weights depending on actual values of X_t
w0 <- ifelse(dat$X_0 == 1, ps_X0, 1 - ps_X0)
w1 <- ifelse(dat$X_1 == 1, ps_X1, 1 - ps_X1)
w2 <- ifelse(dat$X_2 == 1, ps_X2, 1 - ps_X2)
w3 <- ifelse(dat$X_3 == 1, ps_X3, 1 - ps_X3)
w4 <- ifelse(dat$X_4 == 1, ps_X4, 1 - ps_X4)
w5 <- ifelse(dat$X_5 == 1, ps_X5, 1 - ps_X5)
w6 <- ifelse(dat$X_6 == 1, ps_X6, 1 - ps_X6)

dat$iptw_weights <- 1 / (w0 * w1 * w2 * w3 * w4 * w5 * w6)



msm_fit <- glm(Y ~ I(X_0 + X_1 + X_2 + X_3 + X_4 + X_5 + X_6) + C, 
               data = dat, weights = iptw_weights, family = Gamma(link = "log"))
sum_fit <- summary(msm_fit)

test_that("rmsm works as expected for longitudinal data and custom I(.) functions", {
  expect_lt(max(abs(sum_fit$coefficients[,1] - do_pars)/sum_fit$coefficients[,2]), 2.5)
})

