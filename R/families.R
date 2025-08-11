##' Numbers for parametric survival families
##'
##' Data frame containing
##' \itemize{
##' \item `val`: an integer
##' \item `family`: a vector giving the associated parametric family for that integer.
##' }
##' The integer `val` may be used in place of the name of the parametric family
##' when specifying the `family` object.  Note that the numbers are not currently
##' used, but this will be modified in a future version.
##'
##' @format `surv_family_vals` is a `data.frame` with 2 rows and 2 columns
##'
##' @export
surv_family_vals <- data.frame(val=1:3,
                               family=c("exp", "weibull", "lnorm"))

##' @describeIn surv_family_vals List of links for each family
##' @format `surv_links_list` is a list of length 2
##' @export
surv_links_list <- list(exp=c("inverse", "log", "identity"),
                        weibull=c("log"),
                        lnorm = c("identity", "log"))

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

##' Obtain families including survival outcomes
##'
##' @param val integer corresponding to distributional family
##'
##' @details
##' The functions `gaussian_causl_fam()` etc. represent the functions that are
##' returned by `get_family()`.  This function also returns families defined
##' by the `survivl` package.
##'
##'
##' @seealso [family_vals]
##'
get_surv_family <- function (val) {
  surv <- FALSE

  if (is.numeric(val)) {
    fm <- match(val, causl::family_vals$val)
    if (is.na(fm)) stop("Invalid family value")
  }
  else if (is.character(val)) {
    fm <- pmatch(val, causl::family_vals$family)
    if (is.na(fm)) {
      fm <- pmatch(val, surv_family_vals$family)
      surv <- TRUE
    }
    if (is.na(fm)) stop("Family not recognized")
  }

  if (surv) fmly <- surv_family_vals[fm,]$family
  else fmly <- causl::family_vals[fm,]$family
  if (is.numeric(val) && val > 5) stop("No function defined yet for this family")

  if (surv) out <- get(paste0(fmly, "_surv_fam"))
  else out <- get(paste0(fmly, "_causl_fam"))

  return(out)
}


##' @describeIn survivl_fams Weibull distribution
##' @export
weibull_surv_fam <- function (link) {
  if (missing(link)) link <- "log"

  dens <- function (x, mu, scale, shape, med, log=FALSE, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)*mu
      }
      else if (par == "R") {
        stop("Fill this in!")
      }
      else stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(dweibull(x, shape=shape, scale=scale, log=log))
  }
  quan <- function (p, x, mu, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)*mu
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    if (missing(x)) return(qweibull(p, shape=shape, scale=scale))
    else {
      qcweib <- function (p, x, shape, scale) {
        (x^shape - scale^shape*log(1-p))^(1/shape) - x
      }
      return(qcweib(p, x=x, shape=shape, scale=scale))
    }
  }
  # csim <- function (n, x, mu, shape, scale, par="PH") {
  #   if (!missing(scale)) {
  #     if (par == "AFT") {
  #       stop("Fill this in!")
  #     }
  #     else if (par == "PH") {
  #       scale <- scale^(-1/shape)*mu
  #     }
  #     else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
  #   }
  #   else stop("Must supply 'scale' parameter")
  #
  #   rcweib <- function (n, x, shape, scale) {
  #     u <- runif(n)
  #     (x^shape - scale^shape*log(u))^(1/shape)
  #   }
  #
  #   return(rcweib(n, x, shape, scale))
  # }
  sim <- function (n, x, mu, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)*mu
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    if (missing(x)) return(rweibull(n, shape=shape, scale=scale))
    else {
      rcweib <- function (n, x, shape, scale) {
        u <- runif(n)
        (x^shape - scale^shape*log(u))^(1/shape) - x
      }

      return(rcweib(n, x=x, shape=shape, scale=scale))
    }
  }
  probs <- function (x, mu, shape, scale, med, par="PH") {

    if (!missing(scale)) {
      if (par == "AFT") {
        stop("Fill this in!")
      }
      else if (par == "PH") {
        scale <- scale^(-1/shape)*mu
      }
      else if (par != "R") stop("'par' should be either \"PH\", \"AFT\", or \"R\"")
    }
    else if (!missing(med)) scale <- med/(log(2)^(1/shape))
    else stop("Must supply 'med' or 'scale' parameter")

    return(pweibull(x, shape=shape, scale=scale))
  }

  default <- function(theta) list(x=1, shape=1, scale=1, p=0.5)

  out <- list(name="weibull", ddist=dens, qdist=quan, rdist=sim, #rcdist=csim,
              pdist=probs, pars=c("mu", "shape", "scale"), default=default, link=link)
  class(out) <- c("surv_family", "causl_family")

  return(out)
}


##' @describeIn survivl_fams Exponential distribution
##' @export
exp_surv_fam <- function (link) {
  if (missing(link)) link <- "inverse"

  dens <- function (x, mu, log=FALSE) {
    return(dexp(x, rate = 1/mu, log=log))
  }
  quan <- function (p, x, mu) {
    return(qexp(p, rate = 1/mu))
  }
  sim <- function (n, x, mu) {
    return(rexp(n, rate = 1/mu))
  }
  probs <- function (x, mu) {
    return(pexp(x, rate = 1/mu))
  }

  default <- function(theta) list(x=1, mu=1)

  out <- list(name="exp", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu"), default=default, link=link)
  class(out) <- c("surv_family", "causl_family")

  return(out)
}

##' @describeIn survivl_fams Log-normal distribution
##' @export
lnorm_surv_fam <- function (link) {
  if (missing(link)) link <- "identity"
  
  dens <- function (x, mu, phi, log=FALSE) {
    return(dlnorm(x, meanlog = mu, sdlog = phi, log=log))
  }
  quan <- function (p, mu, phi) {
    return(qlnorm(p, meanlog =mu, sdlog = phi))
  }
  sim <- function (n, mu, phi) {
    return(rlnorm(n, meanlog = mu, sdlog =phi))
  }
  probs <- function (x, mu, phi) {
    return(plnorm(x, meanlog = mu, sdlog = phi))
  }
  
  default <- function(theta) list(x=1, mu=0, phi = 1)
  
  out <- list(name="lognormal", ddist=dens, qdist=quan, rdist=sim, pdist=probs,
              pars=c("mu", "phi"), default=default, link=link)
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
    else if (is.character(x)) return(x %in% c("exp", "Gamma", "weibull", "lognormal"))
  }
  else if (methods::is(x[[1]], "causl_family")) {
    nms <- sapply(x, function(y) y$name)
    return(nms %in% c("exp", "Gamma", "weibull", "lognormal"))
  }
  else stop("Not a valid family representation")
}
