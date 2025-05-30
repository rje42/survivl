tv_data <- function (dat, model=Y ~ X) {
  resp <- lhs(model)

  ## get number of timepoints
  nms <- names(dat)
  srch <- regexpr(pattern = paste0("^", resp, "_"), text = nms)
  kp <- nms[srch == 1]
  kp <- substr(kp, nchar(resp)+2, nchar(kp))
  T <- max(as.numeric(kp))
  strt <- min(as.numeric(kp))

  ## identify non-time dependent variables  MAKE MORE ROBUST!
  wh <- which(!grepl("_", nms))

  ##
  out <- survival::tmerge(dat[,wh], data2=dat, id=id, tstop=T)

  for (t in seq_len(T-strt-1)) {
    new_t <- strt + t
    out <- survival::tmerge(out, dat, id=id, Y=event(paste0(resp, "_", new_t)))
  }

  return(out)
}
