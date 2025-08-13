##' Collect events
##'
##' @param dat data frame of observations
##' @param varY character vector of response variable stems
##' @param t time index
##' @param trunc logical indicating whether variables have already been truncated
##'
##' @export
collect_events <- function(dat, varY, t = 0L, trunc = FALSE) {
  ## obtain indices
  varYt <- paste0(varY, "_", t)
  indY <- match(varYt, names(dat))
  if (any(is.na(indY))) stop("Variable names not present")

  event <- dat$status

  if (!trunc) {
    ## obtain the details of the event
    T <- do.call(pmin, dat[indY])
    surv <- T >= 1

    event[!surv] <- max.col(-dat[!surv, indY, drop = FALSE])
    T[T >= 1] <- t + 1
    T[!surv] <- t + T[!surv]

    ## record the outcome in dat
    dat[, indY] <- 0L
    dat[cbind(which(!surv), indY[event])] <- 1L
  } else {
    surv <- do.call(pmax, dat[indY]) == 0
    T <- dat$T
    T[!surv] <- t + 0.5
    if (any(!surv)) event[!surv] <- max.col(-dat[!surv, indY, drop = FALSE])
  }

  dat$T <- T
  dat$status <- event

  return(list(dat = dat, surv = surv, event = event, T = t + T))
}
