standardize_formulae <- function (formulas, static=character(0)) {

  var_nms_Z <- var_nms_X <- var_nms_Y <- var_nms_cop <- list()

  var_nms_Z <- stand2(formulas[[2]], static=static)
  var_nms_X <- stand2(formulas[[3]], static=static)
  var_nms_Y <- stand2(formulas[[4]], static=static)
  if (length(formulas) > 4) {
    cop <- TRUE
    var_nms_cop <- stand2(formulas[[5]], static=static)
  }
  else cop <- FALSE

  out <- list(var_nms_Z=var_nms_Z, var_nms_X=var_nms_X, var_nms_Y=var_nms_Y)
  if (cop) out$var_nms_cop <- var_nms_cop

  return(out)
}

stand2 <- function (formulas, static=character(0)) {

  if ("formula" %in% class(formulas)) formulas <- list(formulas)
  var_nms <- vector("list", length(formulas))

  for (i in seq_along(formulas)) {
    tms <- terms(formulas[[i]])
    vars <- attr(tms, "term.labels")[attr(tms, "order") == 1]
    ends <- sapply(strsplit(vars, "_"), rje::last)
    lags <- rep(NA, length(vars))
    slg <- substr(ends, 1, 1) == "l"
    lags[slg] <- suppressWarnings(as.numeric(substr(ends[slg], 2, nchar(ends[slg]))))

    ## get the stem of each name
    var_stem <- vars
    var_stem[is.na(lags)] <- vars[is.na(lags)]
    var_lag <- vars[!is.na(lags)]
    var_stem[!is.na(lags)] <- substr(var_lag, 1, nchar(var_lag) - nchar(ends[!is.na(lags)])-1)

    ## now set variables without lags to be current
    lags[is.na(lags)] <- 0

    var_nms[[i]] <- data.frame(var=var_stem, lag=lags)

    ## now set lags for static variables to NA
    var_nms[[i]]$lag[var_nms[[i]]$var %in% static] <- NA
    # wh <- match(var_nms[[i]]$var, static, nomatch = 0L)
  }

  var_nms
}

## replace variables in a formula
replace_vars <- function(formula, replace) {

  tmp <- as.character(formula)
  R <- ifelse(length(tmp) == 2, 2, 3)

  for (i in seq_along(nrow(replace))) {
    old <- replace$old[[i]]
    new <- replace$new[[i]]

    ## need to be more careful about partial words here
    tmp[R] <- gsub(old, new, tmp[R])
  }
  if (R == 2) form <- paste(tmp[1], tmp[2])
  else form <- paste(tmp[2], tmp[1], tmp[3])

  as.formula(form)
  formula
}

# replace_var <- function(formula, old, new) {
#   tmp <- as.character(formula)
#   R <- ifelse(length(tmp) == 2, 2, 3)
#
#   tmp[R] <- gsub(old, new, tmp[R])
#
#   if (R == 2) form <- paste(tmp[1], tmp[2])
#   else form <- paste(tmp[2], tmp[1], tmp[3])
#
#   as.formula(form, env = NULL)
# }

curr_formulae <- function (formulas, t) {
  start_at <- 1
  mod_form <- function (form, t) {
    trms <- terms(form)
    trm_labs <- attr(trms, "term.labels")

    nos <- regex_extr("_l([0-9]+)", trm_labs)
    nos <- rapply(nos, function(x) substr(x, 3, nchar(x)))
    nos <- lapply(nos, as.numeric)
    drp <- sapply(nos, function(x) any(x > t-start_at))
    if (any(drp)) {
      if (all(drp)) return(update.formula(form, . ~ 1))
      trms <- drop.terms(trms, which(drp), keep.response=TRUE)
    }

    chr <- as.character(trms)[3]
    wh <- gregexpr("_l([0-9]+)", chr)[[1]]
    if (wh[1] < 0) {
      attributes(trms) <- NULL
      return(as.formula(trms))
    }
    ml <- attr(wh, "match.length")

    nos <- integer(length(wh))

    for (i in seq_along(wh)) {
      if (ml[i] < 3) stop("All matches should be at least three characters")
      nos[i] <- as.numeric(substr(chr, wh[i]+2, wh[i]+ml[i]-1))
    }
    for (i in seq_along(wh)) {
      chr <- sub("_l([0-9]+)", paste0("_", t-nos[i]), chr)
    }

    update.formula(form, paste0(". ~ ", chr))
  }

  rapply(formulas, mod_form, how = "replace", t=t)
}
