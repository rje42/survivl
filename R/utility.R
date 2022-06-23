##' Manipulate survival data
##'
##' @param dat survival dataset
##'
##' @details Functions to change the format of survival data as generated
##' by \code{coxSamp} or \code{coxSamp}.  The function \code{surv_to_long}
##' maps to one row for each time-point and individual, whereas \code{surv_to_wide}
##' has precisely one row for each individual.
##'
##' @name manipulate_survival
NULL

rmv_time_stamps <- function (nms) {
  wh2cut <- regexpr("_[0-9]+$", nms)
  stems <- substr(nms[wh2cut > 0], 1, wh2cut[wh2cut > 0]-1)
  tms <- as.numeric(substr(nms[wh2cut > 0], wh2cut[wh2cut > 0]+1, nchar(nms[wh2cut > 0])))

  mx_tms <- tapply(tms, factor(stems), FUN = max)
  mx_tms <- mx_tms[order(match(names(mx_tms), stems))]  # put in original order
  if (any(mx_tms < 0)) stop("Error in format of names, negative number obtained")
  Tmax <- max(mx_tms)

  if (Tmax > 100) stop("too many time points")

  out <- c(nms[wh2cut < 0], names(mx_tms))
  attr(out, "times") <- c(rep(0,sum(wh2cut < 0)), mx_tms+1)
  out
}

##' @describeIn manipulate_survival change to long format
##' @param lag number of earlier time-points to include
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

nms_from_zero <- function (x) {
  tab <- table(x)
  for (t in seq_along(tab)) {
    x[x == names(tab)[t]] <- paste(names(tab)[t], seq_len(tab[t])-1, sep="_")
  }
  x
}

##' @importFrom dplyr group_by ungroup summarise_if filter `%>%` full_join
##' @importFrom tidyr pivot_longer pivot_wider
##' @describeIn manipulate_survival change to wide format
##' @param tv_covs,fix_covs time-varying and fixed covariates
##' @export
surv_to_wide <- function(dat, tv_covs, fix_covs) {#, formulas) {
  nms <- names(dat)
  nms <- nms[nms != "t" & nms != "t_stop" & nms != "status" & nms != "T"]

  ## if not supplied, obtain lists of the time-varying and fixed covariates
  if (!missing(tv_covs) || !missing(fix_covs)) {
    if (missing(tv_covs)) tv_covs <- setdiff(nms, fix_covs)
    else if (missing(fix_covs)) fix_covs <- setdiff(nms, tv_covs)
    else if (length(intersect(tv_covs, fix_covs)) > 0) stop("Covariates must be either time-varying or static")
  }
  else {
    ## determine time-varying vs static covariates
    tv_covs_tmp <- dat[nms] %>%
      group_by(id) %>%
      summarise_if(.predicate=is.numeric, .funs=var, na.rm=TRUE)
    if (isTRUE(all.equal(names(tv_covs_tmp), "id"))) {
      tv_covs <- character(0)
    }
    else {
      tv_covs <- tv_covs_tmp %>%
        pivot_longer(-id) %>%
        ungroup %>% group_by(name) %>%
        summarise(var=var(value, na.rm=TRUE)) %>%
        filter(var > 0) %>%
        `$`("name")
    }
    ## check for time-varying factors
    tvf_covs_tmp <- dat[nms] %>%
      group_by(id) %>%
      summarise_if(.predicate=is.factor, .funs=function(x) var(as.numeric(x), na.rm=TRUE))
    if (isTRUE(all.equal(names(tvf_covs_tmp), "id"))) {
      tvf_covs <- character(0)
    }
    else {
      tvf_covs <- tvf_covs_tmp %>%
        pivot_longer(-id) %>%
        ungroup %>% group_by(name) %>%
        summarise(var=var(value, na.rm=TRUE)) %>%
        filter(var > 0) %>%
        `$`("name")
    }

    ## combine list of covariates
    tv_covs <- c(tv_covs, tvf_covs)

    fix_covs <- setdiff(nms, c(tv_covs, "id"))
  }

  ## number of timepoints
  tms <- max(dat$t)

  # spec <- tibble(.name = c("id", paste(tv_covs, rep(seq_len(tms+1)-1, each=length(tv_covs)), sep="_")),
  #                .value = c("id", rep(tv_covs, length(tms)+1)))
  # pivot_wider_spec(spec)

  # static <- vector(mode="list", length=length(fix_covs))
  ## this approach won't work for factors!
  static <- rep(list(median), length(fix_covs))
  names(static) <- fix_covs
  stat_vals <- dat %>%
    select(id, fix_covs) %>%
    group_by(id) %>%
    summarise(across(fix_covs, median, na.rm=TRUE))
  # %>%
  #   summarise(median)

  dat <- dat %>%
    pivot_wider(id_cols = id, names_from=t, values_from = tv_covs)

  dat <- dat[,c("id", paste(tv_covs, rep(seq_len(tms+1)-1, each=length(tv_covs)), sep="_"))]
  dat <- full_join(stat_vals, dat, by="id")

  dat
}


