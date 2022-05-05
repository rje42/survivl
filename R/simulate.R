##' @import causl
NULL

##' Simulate from a survival model
##'
##' @param n number of samples
##' @param T number of time-points
##' @param formulas formulae for model
##' @param family families of distributions
##' @param par list of parameter values
##' @param link list of link functions
##' @param control list of options for the algorithm
##'
##' @details The formulas, family and (if specified) link arguments should be
##' lists with four entries; the first relates to any static baseline
##' covariates and the second to any (possibly) time-varying ones, the third
##' to any treatments, and the fourth to the outcome.  The four entries in
##' these lists should have the same length.
##'
##' The function gamma() that controls the acceleration defaults to just a
##' linear combination of the covariates, but it can be adjusted to an arbitrary
##' function of them.
##'
##' Link functions for the Gaussian, t and Gamma distributions can be the
##' identity, inverse or log functions.  Gaussian and t-distributions default to
##' the identity, and Gamma to the log link.  For the Bernoulli the logit and
##' probit links are available.
##'
##' The function currently doesn't support varying the parameters at each
##' time-point.
##'
##' @return A \code{data.frame} containing the variable on the left hand side
##' of each formula.
##'
##' @references Young, Hernan, Picciotto and Robins, 2008. Simulation from
##' Structural Survival Models under Complex Time-Varying Data Structures,
##' \emph{JSM Proceedings, Section on Statistics in Epidemiology: American Statistical Association}.
##'
##' @export
survSamp <- function (n, T, formulas, family, pars, link=NULL, control=list()) {

  con = list(max_wt = 1, warn = 1, cop="cop")
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))

  for (i in 1:3) if ("formula" %in% class(formulas[[i]])) formulas[[i]] <- list(formulas[[i]])

  LHS_C <- causl:::lhs(formulas[[1]])
  LHS_Z <- causl:::lhs(formulas[[2]])
  LHS_X <- causl:::lhs(formulas[[3]])
  LHS_Y <- causl:::lhs(formulas[[4]])
  nms_t <- c(LHS_Z, LHS_X, LHS_Y)
  nms <- c(LHS_C, outer(nms_t, seq_len(T)-1, paste, sep="_"))

  out <- data.frame(rep(list(rep(NA, n)), length(LHS_C) + T*length(nms_t)))
  names(out) <- nms

  RHS_Z <- causl:::rhs_vars(formulas[[1]])
  RHS_X <- causl:::rhs_vars(formulas[[2]])

  link <- causl:::linkSetUp(link, family = family, vars=list(LHS_C,LHS_Z,LHS_X,LHS_Y))

  if (family[[4]] == 3) out[[LHS_Y]] <- rgamma(n, shape=1, rate=pars[[LHS_Y]]$lambda)
  if (is.null(pars$gamma)) pars$gamma <- function(X, beta) X %*% beta
  else stop("Family should be 3")

  tmp <- standardize_formulae(formulas, static=LHS_C)
  var_nms_Z <- tmp$var_nms_Z
  var_nms_X <- tmp$var_nms_X
  var_nms_Y <- tmp$var_nms_Y

  ## simulate static covariates
  for (i in seq_along(LHS_C)) {
    ## now compute etas
    eta <- model.matrix(formulas[[1]][[i]][c(1,3)], data=out) %*% pars[[LHS_C[i]]]$beta
    out[[LHS_C[[i]]]] <- glm_sim(family[[1]][i], eta, link[[1]][i], pars[[LHS_C[[i]]]]$phi)
  }

  ## simulate time-varying covariates and survival
  for (t in seq_len(T)) {

    ## now simulate each variable according to formula and pars
    for (i in seq_along(LHS_Z)) {
      vZd <- var_nms_Z[[i]][!is.na(var_nms_Z[[i]]$lag), ]
      vZs <- var_nms_Z[[i]]$var[is.na(var_nms_Z[[i]]$lag)]
      out2 <- lag_data(out, t, vZd, static=vZs)

      ## now compute etas
      eta <- model.matrix(formulas[[2]][[i]][c(1,3)], data=out2) %*% pars[[LHS_Z[i]]]$beta
      var <- paste0(LHS_Z[[i]], "_", t-1)

      out[[var]] <- glm_sim(family[[2]][i], eta, link[[2]][i], pars[[LHS_Z[[i]]]]$phi)
    }
    for (i in seq_along(LHS_X)) {
      # out2 <- lag_data(out, t, var_nms_X[[i]], static=LHS_C)
      vXd <- var_nms_X[[i]][!is.na(var_nms_X[[i]]$lag), ]
      vXs <- var_nms_X[[i]]$var[is.na(var_nms_X[[i]]$lag)]
      out2 <- lag_data(out, t, vXd, static=vXs)


      ## now compute etas
      eta <- model.matrix(formulas[[3]][[i]][c(1,3)], data=out2) %*% pars[[LHS_X[i]]]$beta
      var <- paste0(LHS_X[[i]], "_", t-1)

      out[[var]] <- glm_sim(family[[3]][i], eta, link[[3]][i], pars[[LHS_X[[i]]]]$phi)
    }
  }

  # cum_haz <- rep(0, n)
  Y_res <- out[[LHS_Y]]

  out$T <- rep(0, n)
  gam_aft <- rep(NA, n)
  surv <- rep(TRUE, n)

  ## now rescale Y to match C, Z and Xs
  ## for now just assume simple time-homogeneous exponential hazard
  for (t in seq_len(T)) {
    na_vars <- grepl(paste0("_",t-1), names(out))
    out[na_vars][!surv,] <- NA

    vYd <- var_nms_Y[[i]][!is.na(var_nms_Y[[i]]$lag), ]
    vYs <- var_nms_Y[[i]]$var[is.na(var_nms_Y[[i]]$lag)]
    out2 <- lag_data(out, t, vYd, static=vYs)
    # out2 <- lag_data(out, t, var_nms_Y[[1]], static=LHS_C)

    ## now compute etas
    MM <- model.matrix(formulas[[4]][c(1,3)], data=out2)
    gam_aft[surv] <- pars$gamma(MM, pars[[LHS_Y[1]]]$beta)

    tmp <- out$T + surv*pmin(1, Y_res/exp(gam_aft))
    out$T[!is.na(tmp)] <- na.omit(tmp)
    tmp <- Y_res - exp(gam_aft)
    Y_res[!is.na(tmp)] <- na.omit(tmp)

    var <- paste0(LHS_Y[[1]], "_", t-1)
    out[[var]] <- ifelse(surv, 1*(Y_res < 0), NA)

    surv <- surv & Y_res > 0

    # cum_haz <- cum_haz + exp(gam_aft)
  }

  out <- cbind(id=seq_len(nrow(out)), out)
  out$status <- 1*surv
  out <- out[names(out) != "Y"]

  return(out)
}
