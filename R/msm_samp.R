##' Simulate from Cox MSM Model
##'
##' Obtain samples from a specified Cox Marginal Structural Model
##' using the frugal parameterization.
##'
##' @param n number of samples
##' @param dat optional data frame for plasmode simulation
##' @param qtls associated quantiles from dat for plasmode simulation
##' @param T number of time points
##' @param formulas list of formulas to use
##' @param family list of families to use
##' @param pars list of parameter settings to use
##' @param link link functions for GLM-like models
##' @param control list of control parameters
##' @param method sampling method (defaults to `"inversion"`)
##'
##' @details Samples from a Marginal Structural Model specified
##' by a frugal parameterization; that is, one must specify the
##' marginal model (i.e. the dependence of the outcome on the
##' treatments), the distribution of the covariates, observed
##' confounders and treatments, and then a copula to join the distribution
##' of the outcome to that of the confounders.
##'
##' Among the left-hand sides of outcome variables, the variable 'Cen' has a
##' special meaning as censoring.  This keyword can be changed to something
##' else by using the argument \code{censor} in the \code{control} list.
##'
##' @return An object of class \code{survivl_dat} containing the simulated data.
##'
##' @importFrom survival Surv
##'
##' @export
msm_samp <- function(n, dat = NULL, qtls = NULL, T, formulas, family, pars, link = NULL,
                     method = "inversion", control = list()) {
  con <- list(verbose = FALSE, max_wt = 1, warn = 1, cop = "cop", censor = "Cen", start_at = 0)
  matches <- match(names(control), names(con))
  con[matches] <- control[!is.na(matches)]
  if (any(is.na(matches))) {
    warning(
      "Some names in control not matched: ",
      paste(names(control[is.na(matches)]),
        sep = ", "
      )
    )
  }
  verbose <- con$verbose
  kwd <- con$cop
  surv_model <- survivl_model(formulas, family, pars, link, T, dat, qtls, method,
    control = con, kwd = kwd
  )
  dat <- rmsm(n, surv_model)
  return(dat)
}

##' @export
rmsm <- function(n, surv_model, control = list()) {
  kwd <- surv_model$kwd
  dat <- surv_model$dat
  qtls <- surv_model$qtls
  ## infer sample size from dat if possible
  datNULL <- is.null(dat)
  if (datNULL && missing(n)) {
    stop("Must specify a sample size or supply a data frame")
  } else if (missing(n)) {
    n <- nrow(dat)
  } else if (!datNULL && n != nrow(dat)) {
    warning("Dimension of supplied data does not match n, using size of supplied data frame")
    n <- nrow(dat)
  }

  ## extract them again
  formulas <- surv_model$formulas
  pars <- surv_model$pars
  family <- surv_model$family
  link <- surv_model$link

  LHS_C <- surv_model$LHSs$LHS_C
  LHS_Z <- surv_model$LHSs$LHS_Z
  LHS_X <- surv_model$LHSs$LHS_X
  LHS_Y <- surv_model$LHSs$LHS_Y
  dC <- length(LHS_C)
  dZ <- length(LHS_Z)
  dX <- length(LHS_X)
  dY <- length(LHS_Y)
  surv_model$dC <- dC
  surv_model$dZ <- dZ
  surv_model$dX <- dX
  surv_model$dY <- dY
  var_nms_Z <- surv_model$std_form$var_nms_Z
  var_nms_X <- surv_model$std_form$var_nms_X
  var_nms_Y <- surv_model$std_form$var_nms_Y
  var_nms_cop <- surv_model$std_form$var_nms_cop
  order <- surv_model$ordering
  vars_t <- surv_model$vars_t
  vars <- surv_model$vars
  T <- surv_model$T
  # TODO: Fix censoring because the same as competing risk
  censoring <- surv_model$censoring
  # check for censoring
  out <- data.frame(rep(list(rep(NA, n)), length(LHS_C) + T * length(vars_t)), rep(0, n))
  names(out) <- vars
  out$T <- T
  # any columns that are in dat want to fill in new variable out
  # Copy matching columns from dat to out
  if (!datNULL) {
    out[, intersect(colnames(out), colnames(dat))] <- dat[, intersect(colnames(out), colnames(dat))]
  }


  ## simulate static covariates
  baseline_vars <- unlist(LHS_C)
  j <- 1
  for (i in seq_along(LHS_C)) {
    if (all(is.na(out[[LHS_C[[i]]]]))) {
      ## now compute etas
      eta <- model.matrix(update(formulas[[1]][[j]], NULL ~ .), data = out) %*% pars[[LHS_C[i]]]$beta
      tmp <- causl::glm_sim(
        family = family[[1]][[j]], eta = eta, phi = pars[[LHS_C[[i]]]]$phi,
        other_pars = pars[[LHS_C[[i]]]], link = link[[1]][j]
      )
      out[[LHS_C[[i]]]] <- tmp
      # qtls[[LHS_C[[i]]]] <- attr(tmp, "quantile")
      if (is.null(qtls)) {
        qtls <- setNames(data.frame(list(attr(tmp, "quantile"))), LHS_C[[i]])
        next
      }
      qtls <- cbind(
        qtls[seq_len(i - 1)],
        setNames(list(attr(tmp, "quantile")), LHS_C[[i]]),
        qtls[seq(i, ncol(qtls))]
      )
      j <- j + 1
    }
  }
  surv_model$qtls <- qtls


  method <- surv_model$method

  if (surv_model$survival_outcome) {
    method <- paste0(method, "-survivl")
  }

  if (method == "inversion-survivl") {
    out <- sim_inversion_survivl(out, surv_model)
  } else if (method == "inversion") {
    out <- sim_inversion_longitudinal(out, surv_model)
  } else if (method == "bootstrap-survivl") {
    out <- sim_bootstrap_survivl(out, surv_model)
  } else {
    stop("'method' should be \"inversion\" or \"bootstrap\" ")
  }


  out <- cbind(id = seq_len(nrow(out)), out)
  out <- out[, colSums(is.na(out)) < nrow(out)]
  # out$status <- 1*surv
  class(out) <- c("survivl_dat", "data.frame")

  return(out)
}

##' @describeIn msm_samp old name
coxSamp <- function(n, T, formulas, family, pars, link = NULL, control = list()) {
  .Deprecated(new = "msm_samp", msg = "This function is deprecated since version 0.4.0, please use msm_samp()")
  msm_samp(n = n, T = T, formulas = formulas, family = family, pars = pars, link = NULL, control = list())
}

##' @describeIn msm_samp old name
cox_samp <- function(n, T, formulas, family, pars, link = NULL, control = list()) {
  # deprecate_soft(when = "0.3.1", what = "cox_samp()", with = "msm_samp()")
  .Deprecated(new = "msm_samp", msg = "This function is deprecated since version 0.4.0, please use msm_samp()")
  msm_samp(n = n, T = T, formulas = formulas, family = family, pars = pars, link = NULL, control = list())
}
