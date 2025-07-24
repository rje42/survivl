##' Define a `survivl_model` object
##'
##' This defines a `survivl_model` object, that can be used either for simulation
##' or inference.
##'
##' @param kwd word used for copula formula and parameters
##'
##' @details
##' The components `formulas` and `family` must both be specified, and have
##' matching lengths.  If `pars` is specified, then the model can be used for
##' simulation and inference, if not then only for inference.  `link` is optional,
##' and if not specified default links will be used.
##'
##' @export
survivl_model <- function (formulas, family, pars, link, T = T, dat=NULL, qtls = NULL, method="inversion",
                         kwd="cop", control=list()) {

  con = list(verbose=FALSE, max_wt=1, warn=1, cop="cop", censor="Cen", start_at=0,
             pm_cond = TRUE, pm_nlevs = 5, pm_cor_thresh = 0.25,quan_tol = 1e3*.Machine$double.eps)
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))
  out <- process_inputs(formulas=formulas, family=family, pars=pars,
                        link=link, dat=dat, control=con, method=method, T = T, 
                        qtls = qtls)
  class(out) <- "survivl_model"
  
  return(out)
}


##' Display output from `survivl_model`
##'
##' @param x an object of class `survivl_model`
##' @param ... additional arguments (not used)
##'
##' @exportS3Method base::print
print.survivl_model <- function (x, ...) {
  forms <- unlist(x$formulas)
  cat("Frugal survivl model with variables ", paste(unique(rmv_time(x$vars[1:(length(x$vars)-1)])), collapse=", "),"\n")
  
  invisible(x)
}

##' @export
modify <- function (x, ...) {
  UseMethod("modify")
}

##' @export
modify.default <- function (x, ...) {
  args <- list(...)
  for (i in seq_along(args)) x[[names(args)[i]]] <- args[[i]]
}

##' Modify `survivl_model` object
##'
##' Change one or more components of a `survivl_model` object.
##'
##' @inheritParams survivl_model
##' @param x an object of class `survivl_model`
##' @param over logical: should components be added/modified or entirely over-written?
##'
##' This function can be used to modify
##'
##' @exportS3Method
modify.survivl_model <- function (x, over=FALSE, formulas, family, pars, 
                                  T, link, dat, method,qtls,
                                kwd) {
  if (!is(x, "survivl_model")) stop("Must include an object of class 'survivl_model'")
  
  if (missing(formulas) && missing(family) && missing(pars) && missing(link) &&
      missing(dat) && missing(method) && missing(kwd)) return(x)
  if (missing(formulas)) formulas <- x$formulas
  if (missing(family)) family <- x$family
  if (missing(pars)) pars <- x$pars
  else {
    cpars <- x$pars
    cpars[names(pars)] <- pars
    pars <- cpars
  }
  if (missing(link)) link <- x$link
  if (missing(dat)) dat <- x$dat
  if (missing(method)) method <- x$method
  if (missing(kwd)) kwd <- x$kwd
  if(missing(T)) T <- x$T
  if(missing(qtls)) qtls <- x$qtls
  con = list(verbose=FALSE, max_wt=1, warn=1, cop=x$kwd, censor="Cen", 
             start_at=surv_model$start_at,
             pm_cond = TRUE, pm_nlevs = 5, pm_cor_thresh = 0.25,quan_tol = 1e3*.Machine$double.eps)
  
  out <- process_inputs(formulas=formulas, family=family, pars=pars,
                        link=link, dat=dat, method=method, T = T, control = con,
                        qtls = qtls)
  
  return(out)
}

# ##' @export
# modify_trt_prob.survivl_model <- function (cm, prob=0) {
#   if (length(cm$formulas[[2]]) > 1) stop("Must be a single treatment")
#   if (!(isTRUE(all.equal(cm$family[[2]], 5)) ||
#         isTRUE(all.equal(cm$family[[2]], 0)) ||
#         isTRUE(all.equal(cm$family[[2]], "binomial")) ||
#         isTRUE(all.equal(cm$family[[2]], binomial_survivl_fam())))) stop("Must be a single binary treatment")
#
#   cm0 <- cm
#   cm0$formulas[[2]][[1]] <- update.formula(cm0$formulas[[2]][[1]], ". ~ 1")
#
#   if (prob == 0) bd <- logit(1/(1e3*n))
#   else if (prob == 1) bd <- -logit(1/(1e3*n))
#   else stop("Probability should be 0 or 1")
#
#   ## set up new parameter values
#   new_pars0 <- list(list(beta = bd))
#   t_var <- lhs(cm$formulas[[2]])
#   names(new_pars0) <- t_var
#   cm0 <- modify.survivl_model(cm0, pars = new_pars0)
#
#   return(cm0)
# }
