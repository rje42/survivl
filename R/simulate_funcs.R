##' Get log-density for rejection sampling
##'
##' @param family family for simulation
##' @param eta linear component
##' @param link link function
##' @param phi optional dispersion parameter
##'
##' @details Possible
##' families are Gaussian (=1), t (=2), Exponential (=3), beta (=4)
##' Bernoulli/categorical (=5) and log-normal (=6).
##'
##' Link functions are listed in `causl::links_list`.
##'
##' @return Returns a numeric vector of the same length as `eta`.
##'
# ##' @export
# glm_sim <- function(family, eta, link, phi) {
#
#   if (family %in% c(1,2,3,6) && missing(phi)) stop("Must specify a dispersion parameter for this family")
#   if (missing(link)) {
#     ## get the default link for this family
#     link <- causl::familyVals[causl::familyVals$val==family,2]
#   }
#
#   U <- runif(length(eta))
#
#   if (family == 1 || family == 6) {
#     if (link == "identity") {
#       out <- qnorm(U, mean = eta, sd=sqrt(phi))
#       lden <- dnorm(out, mean=eta, sd=sqrt(phi), log=TRUE)
#     }
#     else if (link == "log") {
#       out <- qnorm(U, mean = exp(eta), sd=sqrt(phi))
#       lden <- dnorm(out, mean=exp(eta), sd=sqrt(phi), log=TRUE)
#     }
#     else if (link == "inverse") {
#       out <- qnorm(U, mean = 1/eta, sd=sqrt(phi))
#       lden <- dnorm(out, mean=1/eta, sd=sqrt(phi), log=TRUE)
#     }
#     else stop("invalid link function for (log)-Gaussian distribution")
#
#     if (family == 6) {
#       out <- exp(out)
#     }
#   }
#   else if (family == 2) {
#     # out <- sqrt(phi)*qt(U, df=pars$par2) + eta
#
#     if (link == "identity") {
#       out <- sqrt(phi)*qt(U, df=pars$par2) + eta
#       lden <- dt((out - eta)/sqrt(phi), df=pars$par2, log=TRUE)-log(sqrt(phi))
#     }
#     else if (link == "log") {
#       out <- sqrt(phi)*qt(U, df=pars$par2) + exp(eta)
#       lden <- dt((out - exp(eta))/sqrt(phi), df=pars$par2, log=TRUE)-log(sqrt(phi))
#     }
#     else if (link == "inverse") {
#       out <- sqrt(phi)*qt(U, df=pars$par2) + 1/eta
#       lden <- dt((out - 1/eta)/sqrt(phi), df=pars$par2, log=TRUE)-log(sqrt(phi))
#     }
#     else stop("invalid link function for t-distribution")
#   }
#   else if (family == 3) {
#     # out <- qexp(U, rate = 1/(exp(eta)*sqrt(phi)))
#
#     if (link == "log") mu <- exp(eta)
#     else if (link == "identity") mu <- eta
#     else if (link == "inverse") mu <- 1/eta
#     else stop("invalid link function for Gamma distribution")
#
#     out <- qgamma(U, shape = 1/phi, rate = 1/(mu*phi))
#     lden <- dgamma(out, shape=1/phi, rate=1/(mu*phi), log=TRUE)
#
#   }
#   else if (family == 4) {
#     out <- qbeta(U, shape1 = 1, shape2 = 1)
#     lden <- dbeta(out, shape1 = 1, shape2 = 1, log=TRUE)
#   }
#   else if (family == 0 || family == 5) {
#     if (link == "probit") {
#       out <- 1*(eta + qnorm(U) > 0)
#       lden <- (out==1)*pnorm(eta, log.p = TRUE) + (out==0)*(pnorm(-eta, log.p = TRUE))
#     }
#     else if (link == "logit") {
#       out <- 1*(eta + qlogis(U) > 0)
#       lden <- (out==1)*eta - log1p(exp(eta))
#     }
#     else stop("invalid link function for binomial distribution")
#
#     # trunc <- pars$trunc
#     # trnc <- 1
#     #
#     # stop("Not finished family==5 yet")
#     # mat <- matrix(NA, length(U), length(trunc))
#     # for (j in seq_along(trunc[[trnc]])) {
#     #   mat[,j] <- 1*(U > trunc[[trnc]][j])
#     # }
#     # out <- rowSums(mat)
#   }
#   else stop("family must be between 0 and 6")
#
#   attr(out, "lden") <- lden
#
#   return(out)
# }

##' @describeIn glm_sim Get log-density for rejection sampling
##'
##' @param Y observed value
##' @param other_pars list containing parameters specific to this `family`
##'
##' @export
glm_ldens <- function(Y, family, eta, link, phi, other_pars) {

  if (family %in% c(1,2,3,6) && missing(phi)) stop("Must specify a dispersion parameter for this family")
  if (missing(link)) {
    ## get the default link for this family
    link <- causl::familyVals[causl::familyVals$val==family,2]
  }

  if (length(Y) != length(eta)) stop("Y and eta must have same length")

  if (family == 1 || family == 6) {
    if (family == 6) Y <- log(Y)

    if (link == "identity") {
      lden <- dnorm(Y, mean=eta, sd=sqrt(phi), log=TRUE)
    }
    else if (link == "log") {
      lden <- dnorm(Y, mean=exp(eta), sd=sqrt(phi), log=TRUE)
    }
    else if (link == "inverse") {
      lden <- dnorm(Y, mean=1/eta, sd=sqrt(phi), log=TRUE)
    }
    else stop("invalid link function for (log)-Gaussian distribution")
  }
  else if (family == 2) {

    if (link == "identity") {
      lden <- dt((Y - eta)/sqrt(phi), df=other_pars$par2, log=TRUE)-log(sqrt(phi))
    }
    else if (link == "log") {
      lden <- dt((Y - exp(eta))/sqrt(phi), df=other_pars$par2, log=TRUE)-log(sqrt(phi))
    }
    else if (link == "inverse") {
      lden <- dt((Y - 1/eta)/sqrt(phi), df=other_pars$par2, log=TRUE)-log(sqrt(phi))
    }
    else stop("invalid link function for t-distribution")
  }
  else if (family == 3) {
    if (link == "log") mu <- exp(eta)
    else if (link == "identity") mu <- eta
    else if (link == "inverse") mu <- 1/eta
    else stop("invalid link function for Gamma distribution")

    lden <- dgamma(Y, shape=1/phi, rate=1/(mu*phi), log=TRUE)

  }
  else if (family == 4) {
    lden <- dbeta(Y, shape1 = 1, shape2 = 1, log=TRUE)
  }
  else if (family == 0 || family == 5) {
    if (link == "probit") {
      lden <- (Y==1)*pnorm(eta, log.p = TRUE) + (Y==0)*(pnorm(-eta, log.p = TRUE))
    }
    else if (link == "logit") {
      lden <- (Y==1)*eta - log1p(exp(eta))
    }
    else stop("invalid link function for binomial distribution")

    # trunc <- other_pars$trunc
    # trnc <- 1
    #
    # stop("Not finished family==5 yet")
    # mat <- matrix(NA, length(U), length(trunc))
    # for (j in seq_along(trunc[[trnc]])) {
    #   mat[,j] <- 1*(U > trunc[[trnc]][j])
    # }
    # Y <- rowSums(mat)
  }
  else stop("family must be between 0 and 6")

  return(lden)
}
