##' Manipulate survival data
##'
##' @param dat survival dataset
##'
##' @details Functions to change the format of survival data as generated
##' by \code{coxSamp} or \code{cox_samp}.  The function \code{surv_to_long}
##' maps to one row for each time-point and individual, whereas \code{surv_to_wide}
##' has precisely one row for each individual.
##'
##' @name manipulate_survival
NULL

##' Drop left hand side from a formula
##'
##' @param formula object to remove LHS of
##'
##' @export
drop_LHS <- function (formula) {
  formula(delete.response(terms(formula)))
}

##' Add or remove time stamps to variable names
##'
##' @name time_stamps
NULL

##' @describeIn time_stamps remove time stamps
##' @param nms names to add or remove time stamps from
##'
rmv_time_stamps <- function (nms) {
  wh2cut <- regexpr("_[0-9]+$", nms)
  stems <- substr(nms[wh2cut > 0], 1, wh2cut[wh2cut > 0]-1)
  tms <- as.numeric(substr(nms[wh2cut > 0], wh2cut[wh2cut > 0]+1, nchar(nms[wh2cut > 0]))) # time stamps

  mx_tms <- tapply(tms, factor(stems), FUN = max)
  mx_tms <- mx_tms[order(match(names(mx_tms), stems))]  # put in original order
  if (any(mx_tms < 0)) stop("Error in format of names, negative number obtained")
  Tmax <- max(mx_tms)

  if (Tmax > 100) stop("too many time points")

  out <- c(nms[wh2cut < 0], names(mx_tms))
  attr(out, "times") <- c(rep(0,sum(wh2cut < 0)), mx_tms+1)
  out
}

##' Remove time or lag suffix or I(.) function
##'
##' @param x vector to trim
##'
rmv_time <- function (x) {
  rx <- regexpr("_[0-9]+$", x)
  x[rx > 0] <- substr(x[rx > 0], 1, rx[rx > 0]-1)
  return(x)
}

rmv_lag <- function (x) {
  rx <- regexpr("_l[0-9]+$", x)
  x[rx > 0] <- substr(x[rx > 0], 1, rx[rx > 0]-1)
  return(x)
}

##' Cleans tms object in `process_inputs`
##' 
##' @param tms a vector of tms (already flattened) that may have + ,_lk, or I(.) around them
##' 
clean_tms <- function(tms){
  #remove functionals in the I function
  tms <- unlist(strsplit(tms, "[\\s\\+\\-\\*\\/\\^\\(\\),]+"))
  tms <- trimws(tms) # gets rid of white space
  tms <- tms[tms != ""]  # remove empty strings
  
  tms <- tms[!tms %in% c("I", "exp", "sin", "cos", "sqrt", "log")]

  
  tms <- rmv_time(rmv_lag(tms))
  return(unique(tms))
  
}



##' @describeIn time_stamps replace lags with time stamps
##' @param t time point for adding
##' @param start_at first time point
##'
add_time_stamps <- function (nms, t, start_at=0) {
  wh2cut <- regexpr("_l[0-9]+$", nms)
  stems <- substr(nms[wh2cut > 0], 1, wh2cut[wh2cut > 0]-1)

  tms <- as.numeric(substr(nms[wh2cut > 0], wh2cut[wh2cut > 0]+2, nchar(nms[wh2cut > 0])))

  ## add in time points, setting to NA if prior to start_at
  stems <- paste0(stems, "_", t-tms)
  stems[which(t - tms < start_at)] <- NA_character_

  ## set variable names and return
  nms[wh2cut > 0] <- stems

  return(nms)
}


##' @describeIn manipulate_survival change to long format
##' @param lag number of earlier time-points to include
##' @export
surv_to_long <- function(dat, lag=0) {#, formulas) {
  nms <- names(dat)
  nms <- nms[nms != "id" & nms != "status"]
  n <- nrow(dat)

  ## determine time-varying vs static covariates
  ts <- rmv_time_stamps(nms)
  Ts <- attr(ts, "times")
  Tmax <- max(Ts)
  static <- ts[Ts == 0]

  ## set up a blank data frame
  out <- c(list(id=rep(dat$id, each=Tmax),
           t=rep.int(seq_len(Tmax)-1, n)),
           rep(list(rep.int(NA, Tmax*n)), sum(Ts == 0)),
           rep(list(rep.int(NA, Tmax*n)), (lag+1)*sum(Ts > 0)))
  # out <- vector(mode="list", length=2+length(ts))
  nms <- c("id", "t", ts[Ts==0], ts[Ts > 0])
  if (lag > 0) nms <- c(nms, outer(ts[Ts > 0], seq_len(lag), paste, sep="_l"))
  names(out) <- nms

  out <- data.frame(out)
  # out$t <- rep(seq_len(Tmax)-1, n)

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
    out[(seq_len(n)-1)*Tmax + i, whStat] <- dat[,ts[Ts == 0]]
    out[(seq_len(n)-1)*Tmax + i, whTV(0)] <- dat[,nms]

    for (l in seq_len(min(i-1, lag))) {
      nms2 <- paste0(ts[Ts > 0], "_", i-l-1)
      out[(seq_len(n)-1)*Tmax + i, whTV(l)] <- dat[,nms2]
    }
  }

  ## reconstruct factors among static covariates
  chk <- sapply(dat[ts[Ts==0]], is.factor)
  if (any(chk)) {
    for (i in which(chk)) out[[ts[Ts==0][i]]] <- factor(out[[ts[Ts==0][i]]], labels = levels(dat[[ts[Ts==0][i]]]))
  }

  ## record individual endpoints
  endpt <- rep(0, nrow(out))
  endpt[Tmax*(seq_len(n)-1) + ceiling(dat$T)] <- dat$status
  out$endpt <- endpt

  ## record overall status
  out$status <- rep(dat$status, each=Tmax)

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

##' @importFrom dplyr group_by ungroup summarise_if `%>%` full_join
##' @importFrom tidyr pivot_longer pivot_wider
##' @describeIn manipulate_survival change to wide format
##' @param tv_covs,fix_covs time-varying and fixed covariates
##' @export
surv_to_wide <- function(dat, tv_covs, fix_covs) {#, formulas) {
  nms <- names(dat)
  nms <- nms[nms != "t" & nms != "t_stop" & nms != "status" & nms != "T" & nms != "endpt"]

  ## if not supplied, obtain lists of the time-varying and fixed covariates
  if (!missing(tv_covs) || !missing(fix_covs)) {
    if (missing(tv_covs)) tv_covs <- setdiff(nms, fix_covs)
    else if (missing(fix_covs)) fix_covs <- setdiff(nms, tv_covs)
    else if (length(intersect(tv_covs, fix_covs)) > 0) stop("Covariates must be either time-varying or static")
  }
  else {
    ## determine time-varying vs static covariates
    tv_covs_tmp <- dat[nms] %>%
      dplyr::group_by(id) %>%
      dplyr::summarise_if(.predicate=is.numeric, .funs=var, na.rm=TRUE)
    if (isTRUE(all.equal(names(tv_covs_tmp), "id"))) {
      tv_covs <- character(0)
    }
    else {
      tv_covs <- tv_covs_tmp %>%
        tidyr::pivot_longer(-id) %>%
        dplyr::ungroup %>% group_by(name) %>%
        dplyr::summarise(var=var(value, na.rm=TRUE)) %>%
        dplyr::filter(var > 0) %>%
        `$`("name")
    }
    ## check for time-varying factors
    tvf_covs_tmp <- dat[nms] %>%
      dplyr::group_by(id) %>%
      dplyr::summarise_if(.predicate=is.factor, .funs=function(x) var(as.numeric(x), na.rm=TRUE))
    if (isTRUE(all.equal(names(tvf_covs_tmp), "id"))) {
      tvf_covs <- character(0)
    }
    else {
      tvf_covs <- tvf_covs_tmp %>%
        tidyr::pivot_longer(-id) %>%
        dplyr::ungroup %>% group_by(name) %>%
        dplyr::summarise(var=var(value, na.rm=TRUE)) %>%
        dplyr::filter(var > 0) %>%
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
    dplyr::select(id, fix_covs) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(dplyr::across(fix_covs, median, na.rm=TRUE))
  # %>%
  #   summarise(median)

  dat <- dat %>%
    tidyr::pivot_wider(id_cols = id, names_from=t, values_from = tv_covs)

  dat <- dat[,c("id", paste(tv_covs, rep(seq_len(tms+1)-1, each=length(tv_covs)), sep="_"))]
  dat <- dplyr::full_join(stat_vals, dat, by="id")

  dat
}

##' Extract patterns from a string
##'
##' @param pattern regular expression to match
##' @param x string for pattern extraction; should contain a single capture group
##' @param sub_patt string to use for replacement
##'
##' @details This function is supposed to extract all captured elements from a
##' string, and return them as a character vector.  It actually returns the
##' whole match
##'
regex_extr <- function (pattern, x, sub_patt="{") {

  subs <- gregexpr(pattern, text=x)
  ml <- lapply(subs, attr, "match.length")
  ed <- mapply(`+`, subs, ml, SIMPLIFY = FALSE)
  # ed <- subs + ml - 1
  out <- list()
  for (s in seq_along(x)) {
    out[[s]] <- character(length(subs[[s]]))
    for (i in seq_along(subs[[s]])) {
      out[[s]][i] <- substr(x[[s]], subs[[s]][i], ed[[s]][i]-1)
    }
  }
  out
}


# ## functions relating to testing that a name is valid
# is_valid_name <- function(string) {
#   grepl("^((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])$", string)
# }
#
# is_valid_unreserved <- function(string) {
#   make.names(string) == string
# }
#
# test_validity <- function(string) {
#   valid <- is_valid_name(string)
#   unreserved <- is_valid_unreserved(string)
#   reserved <- (valid & ! unreserved)
#   list("Valid"=valid,
#        "Unreserved"=unreserved,
#        "Reserved"=reserved)
# }

