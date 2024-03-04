##' Log-likelihood for survival data
##'
##' Evaluate the log-likelihood for a survival dataset.
##'
##' @param dat data frame of class \code{survivl_dat}
##' @param control list of control parameters
##' @inheritParams cox_samp
##'
##' @details The only valid argument for \code{control} is \code{cop}, which
##' changes the keyword for the copula from its default (\code{"cop"}).
##'
##' @importFrom causl merge_formulas
##' @export
ll_surv <- function (dat, formulas, family, pars, link, control) {

  n <- nrow(dat)

  # get control parameters for optim or use defaults
  con <- list(cop="cop")
  matches = pmatch(names(control), names(con))
  con[matches[!is.na(matches)]] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ", paste(names(control[is.na(matches)]), sep = ", "))

  ## ensure copula keyword is not a variable name
  kwd <- con$cop
  if (kwd %in% names(dat)) stop(paste("Must not have a variable named '", kwd, "'", sep=""))

  ## tidy up the formulae
  forms <- tidy_formulas(formulas, kwd=kwd)
  fam_cop <- rje::last(family)
  link <- causl::link_setup(link, family = family[-length(family)])

  LHS <- lhs(forms[-length(forms)])
  full_form <- merge_formulas(forms)
  wh <- full_form$wh
  mm <- model.matrix(full_form$formula, data=dat)

  ## get initial parameters
  pars2 <- initialize_params_surv(dat, formulas=forms, family=family, pars=pars,
                                  link=link, full_form=full_form, kwd=kwd)
  # theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0],
  #               beta_start2$phi[beta_start2$phi_m > 0])
  #
  # theta_st <- c(beta_start2$beta[beta_start2$beta_m > 0], beta_start2$phi[beta_start2$phi_m > 0])

  ## get linear components
  eta <- mm %*% pars2$beta
  out <- numeric(n)

  for (i in seq_along(LHS)) {
    out <- out + glm_ldens(dat[[LHS[i]]], family=family[i], eta=eta[,i],
                           link=link[i], phi=pars2$phi[i])
  }

}
