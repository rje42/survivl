##' Families for survival distributions
##'
##' For each family, the function returns a list containing a density,
##' simulation function, quantile function, and cumulative distribution
##' function.  These have the generic names `ddist`, `rdist`, `qdist` and
##' `pdist` respectively.
##'
##' @param link link function to use
##'
##' @import stats
##' @name survivl_fams
NULL

##' @describeIn survivl_fams Weibull distribution
##' @export
weibull_surv_fam <- function (link) {
  if (missing(link)) link <- "log"

  dens <- function (x, scale, shape, med, log=FALSE, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(dweibull(x, shape=shape, scale=scale, log=log))
  }
  quan <- function (x, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(qweibull(x, shape=shape, scale=scale))
  }
  sim <- function (x, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(rweibull(x, shape=shape, scale=scale))
  }
  probs <- function (x, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(pweibull(x, shape=shape, scale=scale))
  }

  default <- function(theta) list(x=1, shape=1, scale=1, p=0.5)

  out <- list(name="weibull", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("shape", "scale"), default=default, link=link)
  class(out) <- c("surv_family", "causl_family")

  return(out)
}

##' Check if family can model a time-to-event outcome
##'
##' @param x vector or list of families
##'
is_surv_outcome <- function (x) {
  if (is.atomic(x)) {
    if(is.numeric(x)) return(x %in% c(3,6))
    else if (is.character(x)) return(x %in% c("Gamma", "weibull", "lognormal"))
  }
  else if (methods::is(x[[1]], "causl_family")) {
    nms <- sapply(x, function(y) y$name)
    return(x %in% c("Gamma", "weibull", "lognormal"))
  }
  else stop("Not a valid family representation")
}
