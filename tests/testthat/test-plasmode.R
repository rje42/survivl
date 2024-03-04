## first test for static covariates
set.seed(123)
n <- 1e4

df <- data.frame(Z1=rnorm(n), Z2=factor(sample(3, n, replace = TRUE)))

# msm_samp(T=3, formulas=forms, family=fams, pars=pars, )
#
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
