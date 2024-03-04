##' Get lagged version of a dataset for a particular time point
##'
##' @param dat data set
##' @param t current time point
##' @param var_nms table output from `standardize_formulas`
##' @param static character vector of static variables
##'
lag_data <- function (dat, t, var_nms, static=character(0)) {
  nm <- paste(var_nms$var, t-1-var_nms$lag, sep="_")
  out2 <- do.call(data.frame, as.list(nm))
  names(out2) <- nm
  nmp <- nm[!grepl("-",nm)]
  nmn <- nm[grepl("-",nm)]
  out <- dat[nmp]
  out[nmn] <- 0

  nm2 <- paste(var_nms$var, var_nms$lag, sep="_l")
  names(out) <- nm2
  out <- cbind(dat[static], out)

  out
}
