##' Process formulas, families and parameters
##'
##' @param formulas list of lists of formulas
##' @param pars list of lists of parameters
##' @param family families for Z,X,Y and copula
##' @param link list of link functions
## @param kwd keyword for copula
##' @param control control variables
## @param ordering logical: should an ordering of variables be computed?
##'
process_inputs <- function (formulas, pars, family, link, T, control) {

  for (i in 1:5) if ("formula" %in% class(formulas[[i]])) formulas[[i]] <- list(formulas[[i]])

  LHS_C <- causl:::lhs(formulas[[1]])
  LHS_Z <- causl:::lhs(formulas[[2]])
  LHS_X <- causl:::lhs(formulas[[3]])
  LHS_Y <- causl:::lhs(formulas[[4]])
  # LHS_cop <- causl:::lhs(formulas[[5]])
  dZ <- length(LHS_Z); dX <- length(LHS_X); dY <- length(LHS_Y)
  if (any(duplicated(c(LHS_C, LHS_Z, LHS_X, LHS_Y)))) stop("Repeated variable names not allowed")
  # if (any(duplicated(c(LHS_C, LHS_Z, LHS_X, LHS_Y, LHS_cop)))) stop("Repeated variable names not allowed")

  if (missing(family)) {
    ## assume everything is Gaussian
    stop("Families must be specified")
  }
  else if (is.list(family)) {
    if (!all(lengths(family[1:4]) == lengths(formulas[1:4]))) stop("mismatch in family and formulae specifications")
  }
  else if (length(family) == 5) family <- as.list(family)
  else stop("family should be a list or vector of length 5")

  ## extract families
  famC <- family[[1]]
  famZ <- family[[2]]
  famX <- family[[3]]
  famY <- family[[4]]
  famCop <- family[[5]]

  ## check right number of parameters supplied
  formsC <- lapply(formulas[[1]], terms)
  tmsC <- lapply(formsC, attr, "term.labels")
  formsZ <- lapply(formulas[[2]], terms)
  tmsZ <- lapply(formsZ, attr, "term.labels")
  formsX <- lapply(formulas[[3]], terms)
  tmsX <- lapply(formsX, attr, "term.labels")
  formsY <- lapply(formulas[[4]], terms)
  tmsY <- lapply(formsY, attr, "term.labels")

  for (i in seq_along(formulas[[1]])) {
    npar <- length(tmsC) + attr(formsC[[i]], "intercept")
    if (length(pars[[LHS_C[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_C[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[2]])) {
    npar <- length(tmsZ) + attr(formsZ[[i]], "intercept")
    if (length(pars[[LHS_Z[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[3]])) {
    npar <- length(tmsX) + attr(formsX[[i]], "intercept")
    if (length(pars[[LHS_X[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[4]])) {
    npar <- length(tmsY) + attr(formsY[[i]], "intercept")
    if (length(pars[[LHS_Y[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
  }

  ## check for censoring
  censoring <- (control$censor %in% LHS_Y)  # see if censoring variable is in list of Y variables
  if (censoring) {
    ## put censoring first if a competing risk
    LHS_Y <- c(control$censor, setdiff(LHS_Y, control$censor))
  }
  nms_t <- c(LHS_Z, LHS_X, LHS_Y)
  nms <- c(LHS_C, outer(nms_t, seq_len(T)-1, paste, sep="_"), "status")

  ## make sure family entries have a vector of integers
  if (any(sapply(family, is.null))) family[sapply(family, is.null)] <- list(integer(0))
  ## now set up link functions
  link <- causl:::link_setup(link, family = family[-(5)], vars=list(LHS_C,LHS_Z,LHS_X,LHS_Y))
  # link[[4]] <- "inverse"

  if (is.null(pars$gamma)) pars$gamma <- function(X, beta) X %*% beta

  std_form <- standardize_formulae(formulas, static=LHS_C)

  ## get ordering for terms
  A <- matrix(0,dZ+dX+dY,dZ+dX+dY)

  for (i in seq_along(LHS_Z)) {
    tab <- std_form$var_nms_Z[[i]]
    A[i,] <- 1*(nms_t %in% tab$var[tab$lag == 0])
  }
  for (i in seq_along(LHS_X)) {
    tab <- std_form$var_nms_X[[i]]
    A[i+dZ,] <- 1*(nms_t %in% tab$var[tab$lag == 0])
  }
  for (i in seq_along(LHS_Y)) {
    tab <- std_form$var_nms_Y[[i]]
    A[i+dZ+dX,] <- 1*(nms_t %in% tab$var[tab$lag == 0])
  }

  ordering <- causl:::topOrd(A)
  if (any(is.na(ordering))) stop("Cyclic dependence in formulae provided")

  out <- list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_C=LHS_C, LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y),
              std_form=std_form, ordering=ordering, var=nms, var_t=nms_t)
}

##' Modify inputs for simulation
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
  ed <- pnms %in% proc_inputs$var_t
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

