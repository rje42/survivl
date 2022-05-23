rmv_time_stamps <- function (nms) {
  wh2cut <- regexpr("_[0-9]+$", nms)
  stems <- substr(nms[wh2cut > 0], 1, wh2cut[wh2cut > 0]-1)
  tms <- as.numeric(substr(nms[wh2cut > 0], wh2cut[wh2cut > 0]+1, nchar(nms[wh2cut > 0])))

  mx_tms <- tapply(tms, factor(stems), FUN = max)
  if (any(mx_tms < 0)) stop("Error in format of names, negative number obtained")
  Tmax <- max(mx_tms)

  if (Tmax > 100) stop("too many time points")

  out <- c(nms[wh2cut < 0], names(mx_tms))
  attr(out, "times") <- c(rep(0,sum(wh2cut < 0)), mx_tms+1)
  out
}

##' @export
surv_to_long <- function(dat, lag=0) {#, formulas) {
  nms <- names(dat)
  nms <- nms[nms != "id" & nms != "status"]

  ## determine time-varying vs static covariates
  ts <- rmv_time_stamps(nms)
  Ts <- attr(ts, "times")
  Tmax <- max(Ts)
  static <- ts[Ts == 0]

  ## set up a blank data frame
  out <- c(list(id=rep(dat$id, each=Tmax),
           t=rep.int(seq_len(Tmax)-1, nrow(dat))),
           rep(list(rep.int(NA, Tmax*nrow(dat))), sum(Ts == 0)),
           rep(list(rep.int(NA, Tmax*nrow(dat))), (lag+1)*sum(Ts > 0)))
  # out <- vector(mode="list", length=2+length(ts))
  nms <- c("id", "t", ts[Ts==0], ts[Ts > 0])
  if (lag > 0) nms <- c(nms, outer(ts[Ts > 0], seq_len(lag), paste, sep="_l"))
  names(out) <- nms

  out <- data.frame(out)
  # out$t <- rep(seq_len(Tmax)-1, nrow(dat))

  # out <- cbind(out, Ts > 0)

  ## get locations of static and time-varying columns
  whStat <- c(FALSE, FALSE, rep(TRUE, sum(Ts==0)), rep(FALSE, ncol(out)-(sum(Ts==0)+2)))
  ## function to get (lagged) columns for time-varying covariates
  whTV <- function(l) {
    ret <- rep(FALSE, ncol(out))
    n_dyn <- sum(Ts > 0)
    ret[2+sum(Ts==0)+n_dyn*l+seq_len(n_dyn)] <- TRUE
    ret
  }

  ## now set values for each timepoint
  for (i in seq_len(Tmax)) {
    nms <- paste0(ts[Ts > 0], "_", i-1)
    out[(seq_len(nrow(dat))-1)*Tmax + i, whStat] <- dat[,ts[Ts == 0]]
    out[(seq_len(nrow(dat))-1)*Tmax + i, whTV(0)] <- dat[,nms]

    for (l in seq_len(min(i-1, lag))) {
      nms2 <- paste0(ts[Ts > 0], "_", i-l-1)
      out[(seq_len(nrow(dat))-1)*Tmax + i, whTV(l)] <- dat[,nms2]
    }
  }

  ## reconstruct factors among static covariates
  chk <- sapply(dat[ts[Ts==0]], is.factor)
  if (any(chk)) {
    for (i in which(chk)) out[[ts[Ts==0][i]]] <- factor(out[[ts[Ts==0][i]]], labels = levels(dat[[ts[Ts==0][i]]]))
  }

  status <- rep(0, nrow(out))
  status[Tmax*(seq_len(nrow(dat))-1) + ceiling(dat$T)] <- dat$status
  out$status <- status

  out <- out[!apply(out, 1, function(x) any(is.na(x))),]

  out$t_stop <- pmin(out$t+1, out$T)
  out2 <- out[setdiff(names(out), c("id","t","t_stop","T","status"))]
  out <- cbind(id=out$id, t=out$t, t_stop=out$t_stop, out2, T=out$T, status=out$status)
  # out <- mutate(out, t_stop=pmin(out$t+1, out$T), .after=t)
  out
}
