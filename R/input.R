##' Process formulas, families and parameters
##'
## @param formulas list of lists of formulas
## @param pars list of lists of parameters
## @param family families for Z,X,Y and copula
## @param link list of link functions
## @param kwd keyword for copula
##' @inheritParams msm_samp
## @param control control variables
## @param method method to be used for sampling
## @param ordering logical: should an ordering of variables be computed?
##'
process_inputs <- function (formulas, pars, family, link, dat, T, method, control) {

  for (i in 1:5) if ("formula" %in% class(formulas[[i]])) formulas[[i]] <- list(formulas[[i]])

  ## check for censoring
  cens <- check_censoring(formulas[[4]], pars, cns_kwd = control$censor)
  censoring <- cens$censoring
  if (censoring) {
    formulas[[4]] <- cens$formulas
    pars <- cens$pars
  }

  LHSs <- lapply(formulas[1:4], causl::lhs)
  LHS_C <- LHSs[[1]]; LHS_Z <- LHSs[[2]]; LHS_X <- LHSs[[3]]; LHS_Y <- LHSs[[4]]
  # LHS_C <- causl::lhs(formulas[[1]])
  # LHS_Z <- causl::lhs(formulas[[2]])
  # LHS_X <- causl::lhs(formulas[[3]])
  # LHS_Y <- causl::lhs(formulas[[4]])

  dZ <- length(LHS_Z); dX <- length(LHS_X); dY <- length(LHS_Y)
  if (any(duplicated(unlist(LHSs)))) stop("Repeated variable names not allowed")


  ## check right number of parameters supplied
  forms <- lapply(formulas[1:4], function (x) lapply(x, terms))
  tms <- lapply(forms, function(x) lapply(x,  attr, "term.labels"))

  ## get response variables list
  RHS_vars <- rmv_lag(unlist(tms))

  if (!all(RHS_vars %in% unlist(LHSs))) {
    wh <- RHS_vars[!RHS_vars %in% unlist(LHSs)]
    wh <- unique.default(wh)
    stop(paste0("Variables ", paste(wh, collapse=", "), " appear on right hand side but are not simulated"))
  }

  ## introduce code from causl
  dims <- lengths(formulas)
  family <- causl::process_family(family=family, dims=dims, func_return=get_surv_family)

  ## now set up link functions
  link <- causl::link_setup(link, family = family[-(5)], vars=LHSs,
                            sources=c(links_list, surv_links_list))
  # link[[4]] <- "inverse"

  # LHSs <- list(LHS_C=LHS_C, LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y)
  dummy_dat <- causl::gen_dummy_dat(family=family, pars=pars, dat=dat, LHSs=LHSs, dims=dims)
  # causl::check_pars(formulas=formulas, family=family, pars=pars, dummy_dat=dummy_dat, LHSs=LHSs, kwd=control$cop, dims=dims)

  ## naive check for number of parameters
  for (j in seq_along(LHSs)) for (i in seq_along(formulas[[j]])) {
    npar <- length(tms[[j]][[i]]) + attr(forms[[j]][[i]], "intercept")
    if (length(pars[[LHSs[[j]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[j]][i], " does not match number of coefficients provided"))
  }

  ## useful summaries
  nms_t <- c(LHS_Z, LHS_X, LHS_Y)
  nms <- c(LHS_C, outer(nms_t, seq_len(T)-1, paste, sep="_"), "status")

  if (is.null(pars$gamma)) pars$gamma <- function(X, beta) X %*% beta

  ## get a variable ordering for simulation and 'standardized' formulas
  tmp <- var_order(formulas=formulas, nms_t=nms_t, LHSs=LHSs, dims=dims)
  ord <- tmp$ordering
  std_form <- tmp$std_form

  ## copula related things
  kwd <- control$cop

  if (method == 'inversion' || method == 'bootstrap') {

    tmp <- causl::pair_copula_setup(formulas=formulas[[5]], family=family[[5]], pars=pars[[kwd]],
                                     LHSs=LHSs, quans=character(0), ord=ord)
    formulas[[5]] <- tmp$formulas
    family[[5]] <- tmp$family
    pars[[kwd]] <- tmp$pars

  }


  ## check that outcomes are OK for time-to-event
  if (any(!is_surv_outcome(family[[4]]))) {
    whn <- which(!is_surv_outcome(family[[4]]))[1]
    stop(paste0("outcome '", LHSs[[4]][whn], "' must be of survival type (non-negative and continuous)"))
  }

  out <- list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_C=LHS_C, LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y),
              std_form=std_form, ordering=ord, vars=nms, vars_t=nms_t)
}


##' Obtain only time-varying inputs
##'
##' @param proc_inputs output of `process_inputs`
##'
modify_inputs <- function (proc_inputs) {
  proc_inputs$formulas <- proc_inputs$formulas[2:5]
  proc_inputs$family <- proc_inputs$family[2:5]
  proc_inputs$link <- proc_inputs$link[2:4]
  proc_inputs$dim <- lengths(proc_inputs$formulas[1:3])

  # proc_inputs$formulas[[4]] <- list(proc_inputs$formulas[[4]])

  return(proc_inputs)
}

modify_LHSs <- function (proc_inputs, t) {
  ## modify parameter names
  pnms <- names(proc_inputs$pars)
  ed <- pnms %in% proc_inputs$vars_t
  pnms[ed] <- paste0(pnms[ed], "_", t)
  names(proc_inputs$pars) <- pnms

  ## modify LHS expressions
  proc_inputs$LHSs$LHS_Z <- paste0(proc_inputs$LHSs$LHS_Z, "_", t)
  proc_inputs$LHSs$LHS_X <- paste0(proc_inputs$LHSs$LHS_X, "_", t)
  proc_inputs$LHSs$LHS_Y <- paste0(proc_inputs$LHSs$LHS_Y, "_", t)
  dZ <- proc_inputs$dim[1]
  dX <- proc_inputs$dim[2]

  ## list of variables
  proc_inputs$vars <- unlist(proc_inputs$LHSs[-1])

  return(proc_inputs)
}

##' Check for censoring among outcome formulas
##'
##' @param formulas formulas only for outcome variables
##' @param pars list of parameters
##' @param cns_kwd word used to denote censoring
##'
##' @return A list containing (possibly) edited `formulas`, `pars`, and a
##' logical `censoring` indicating the presence of censoring.
##'
check_censoring <- function (formulas, pars, cns_kwd) {

  LHS_Y <- causl::lhs(formulas, surv = TRUE)

  if (is.list(LHS_Y)) {
    censoring <- TRUE

    ## deal with censored variables
    lY <- lengths(LHS_Y)
    cenY <- which(lY == 2)
    LHS_Y <- c(cns_kwd, sapply(LHS_Y, function(x) x[[1]]))

    cens_form <- merge_formulas(formulas[cenY])
    formulas <- c(cens_form$formula, formulas)
    cenY <- cenY + 1
    causl::lhs(formulas) <- LHS_Y

    ## add in parameters
    if (is.null(pars[[cns_kwd]])) {
      if (length(cenY) == 1) {
        pars[[cns_kwd]] <- pars[[LHS_Y[cenY]]]
      }
      else if (length(cenY > 1)) {
        message("Multiple censoring specifications: attempting to construct merged parameter values")
        refs <- cens_form$reforms
        whF <- cens_form$wh
        trms <- terms(cens_form$formula)
        ntrms <- length(attr(trms, "order")) + attr(trms, "intercept")

        has_all <- sapply(whF, function (x) setequal(x,seq_len(ntrms)))
        if (any(has_all)) {
          wh_ha <- which(has_all)
          pars[[cns_kwd]] <- pars[[LHS_Y[cenY][wh_ha]]]
          message(paste0("Used parameters from outcome ", LHS_Y[cenY][wh_ha], " for censoring parameters"))
        }
        else stop("Censoring distribution not inferable from input")
      }
    }
  }
  else if (cns_kwd %in% LHS_Y) {
    censoring <- TRUE

    wh <- match(cns_kwd, LHS_Y)

    ## put censoring first if a competing risk
    LHS_Y <- LHS_Y[c(wh,seq_along(LHS_Y)[-wh])]
    formulas <- formulas[c(wh,seq_along(LHS_Y)[-wh])]
  }
  else censoring <- FALSE

  return(list(formulas=formulas, pars=pars, censoring=censoring))
}

var_order <- function (formulas, nms_t, LHSs, dims) {
  std_form <- standardize_formulae(formulas, static=LHSs[[1]])

  ## get ordering for terms
  A <- matrix(0,sum(dims[2:4]),sum(dims[2:4]))

  d_C <- dims[-1]

  for (j in 2:4) for (i in seq_along(LHSs[[j]])) {
    tab <- std_form[[j-1]][[i]]
    A[sum(d_C[seq_len(j-2)]) + i, ] <- 1*(nms_t %in% tab$var[tab$lag == 0])
  }

  ## get appropriate ordering
  ord <- causl:::topOrd(A)

  if (any(is.na(ord))) stop("Cyclic dependence in formulae provided")

  return(list(ordering=ord, std_form=std_form))
}

##' @importFrom causl gen_dummy_dat pair_copula_setup check_pars
NULL
