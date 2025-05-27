##' Simulate from Survival Model using Seamon-Keogh Bootstrap Method
##'
##' Obtain samples from a specified Cox Marginal Structural Model
##' using the frugal parameterization.
##'
##' @importFrom causl sim_variable
##' @importFrom causl `lhs<-`
##' @param MMs all the Model Matrices passed as arguments in a list.
##' @param pars list of parameter settings to use

##' @export
boot_sim <- function (out, proc_inputs, quantiles, kwd, MMs) {

  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  # dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X

  vars <- paste0(proc_inputs$vars_t, "_", proc_inputs$t)
  d <- lengths(formulas)

  ## simulate covariates and treatments
  for (j in 1:2) for (i in seq_len(d[j])) {
    vnm <- vars[i+(j-1)*length(formulas[[1]])]
    curr_link <- link[[j]][i]

    curr_form <- formulas[[j]][[i]]
    curr_fam <- family[[j]][[i]]

    MM <- MMs[[j]][[i]]
    eta <- MM %*% pars[[vnm]]$beta
    curr_phi <- pars[[vnm]]$phi
    tmp <- causl::glm_sim(fam=curr_fam, eta=eta, phi=curr_phi, link=curr_link,
                          other_pars=pars[[vnm]], quantiles = FALSE)
    out[[vnm]] <- tmp
    # quantiles[[vnm]] <- attr(tmp, "quantile")
  }

  # vnm <- lhs(formulas[[3]])
  # vnm_stm <- rmv_time(vnm)
  #
  # ## code to get Y quantiles conditional on different Zs
  # for (i in seq_along(formulas[[3]])) {
  #   ## simulate Y variable
  #
  #   forms <- list(formulas[[3]][[i]], formulas[[4]][[i]])
  #   fams <- list(family[[3]][[i]], family[[4]][[i]])
  #   prs <- list(c(pars[[vnm[i]]], list(x=proc_inputs$t)), pars[[kwd]][[vnm_stm[i]]])
  #   if (!is.null(prs[[1]]$lambda0)) prs[[1]]$phi <- 1 #prs[[1]]$phi <- prs[[1]]$lambda0
  #   lnk <- list(link[[3]][i], list()) # link[[4]][[i]])
  #
  #   out <- sim_variable(nrow(out), forms, fams, prs, lnk,
  #                              dat=out, quantiles=quantiles)
  #   quantiles <- attr(out, "quantiles")
  #   attr(out, "quantiles") <- NULL
  #
  # }

  return(list(dat=out, quantiles=quantiles))
}
