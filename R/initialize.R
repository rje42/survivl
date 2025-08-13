## NEEDS TO ACCOUNT FOR LINK FUNCTIONS
## @param dat dataset in long format
## @param formulas formulas for variables to be fitted
## @param family vector of families of length 1 less than number of formulae
## @param link link functions
## @param full_form merged formulas
## @param kwd character string reserved for copula
## @param notInCop character strings for variables not included in the copula
##
## @details Computes the log-likelihood for a Cox MSM model
##
initialize_params_surv <- function(dat, formulas, family = rep(1, nv), pars, link, # init=FALSE,
                                   full_form, kwd, notInCop) {
  nv <- length(formulas) - 1

  if (length(family) != nv + 1) stop("Must have family parameter for every variable and copula")

  ## get terms labels
  trms <- terms(full_form$formula)
  labs <- if (attr(trms, "intercept") == 1) {
    c("(intercept)")
  } else {
    character(0)
  }
  labs <- c(labs, attr(trms, "term.labels"))

  wh <- full_form$wh

  ## define output matrix/vector for beta,phi
  LHS <- lhs(formulas)
  LHS <- LHS[-match(kwd, LHS)]
  if (length(LHS) != nv) stop("Formula for copula causing problems!")

  if (missing(notInCop)) {
    LHSc <- LHS
  } else {
    LHSc <- setdiff(LHS, notInCop)
  }
  nc <- length(LHSc)
  c2 <- utils::combn(nc, 2)
  colnms <- c(LHS, mapply(function(x, y) paste0(kwd, "_", LHSc[x], "_", LHSc[y]), c2[1, ], c2[2, ]))

  beta <- beta_m <- matrix(0, nrow = max(unlist(wh)), ncol = nv + choose(nc, 2), dimnames = list(labs, colnms))
  phi <- phi_m <- numeric(length(LHS))

  # FIGURE OUT HOW TO MAKE THIS WORK WITH NOT EVERYTHING IN THE COPULA

  # LHSs <- lhs(formulas)

  # if (init) {
  #   wts <- ipw_weights(dat, formulas[-length(formulas)])
  # }
  # else wts <- rep(1,nrow(dat))
  # dat <- cbind(dat, wts=wts)

  ## intialize parameters for each variable
  for (i in seq_along(phi)) {
    beta_m[wh[[i]], i] <- 1
    beta[wh[[i]], i] <- pars[[LHS[i]]]$beta

    if (family[i] <= 3) {
      # beta[1,i] <- mean(dat[[LHS[i]]])
      # phi[i] <- var(dat[[LHS[i]]])
      phi[i] <- pars[[LHS[i]]]$phi
      phi_m[i] <- 1
    }
    # else if (family[i] == 3) {
    #   mu_i <- mean(dat[[LHS[i]]])
    #   beta[1,i] <- exp(mu_i)
    #   phi[i] <- var(dat[[LHS[i]]])/mu_i^2
    #   phi_m[i] <- 1
    # }
    # else phi[i] <- NA
  }

  ## initialize copula parameters
  cp <- length(phi)
  beta_m[wh[[nv + 1]], nv + seq_len(choose(nc, 2))] <- 1
  beta[wh[[nv + 1]], nv + seq_len(choose(nc, 2))] <- c(pars[[kwd]]$beta)

  return(list(beta = beta, phi = phi, beta_m = beta_m, phi_m = phi_m))
}
