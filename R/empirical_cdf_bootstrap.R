##' Simulate from Survival Model using Seamon-Keogh Bootstrap Method
##'
##' Obtain samples from a specified Cox Marginal Structural Model
##' using the frugal parameterization.
##'
##' @importFrom causl sim_variable
##' @importFrom causl `lhs<-`
##' @param out Simulated data so far.
##' @param pars list of parameter settings to use

##' @export
empirical_cdf_bootstrap <- function(mod_inputs, out,
                                    riskH) {
  formulas <- mod_inputs$formulas
  pars <- mod_inputs$pars
  family <- mod_inputs$family
  link <- mod_inputs$link
  LHS_Z <- mod_inputs$LHSs$LHS_Z
  LHS_X <- mod_inputs$LHSs$LHS_X

  vars <- paste0(mod_inputs$vars_t, "_", mod_inputs$t)
  d <- lengths(formulas)
  ## simulate covariates and treatments
  X <- NULL
  for (j in 1:2) {
    for (i in seq_len(d[j])) {
      vnm <- vars[i + (j - 1) * length(formulas[[1]])]
      curr_link <- link[[j]][i]

      curr_form <- formulas[[j]][[i]]
      curr_fam <- family[[j]][[i]]
      MM <- model.matrix(delete.response(terms(curr_form)), data = out)

      eta <- MM %*% pars[[vnm]]$beta
      curr_phi <- pars[[vnm]]$phi
      tmp <- causl::glm_sim(
        fam = curr_fam, eta = eta, phi = curr_phi, link = curr_link,
        other_pars = pars[[vnm]], quantiles = FALSE
      )
      X <- cbind(X, tmp)
    }
  }
  H_j <- apply(X, 1, riskH)
  return(H_j)
}
