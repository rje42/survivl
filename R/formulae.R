##' @importFrom purrr pluck_depth list_flatten
standardize_formulae <- function (formulas, static=character(0)) {

  var_nms_Z <- var_nms_X <- var_nms_Y <- var_nms_cop <- list()

  var_nms_Z <- stand2(formulas[[2]], static=static)
  var_nms_X <- stand2(formulas[[3]], static=static)
  var_nms_Y <- stand2(formulas[[4]], static=static)
  if (length(formulas) > 4) {
    cop <- TRUE
    var_nms_cop <- stand2(unlist(formulas[[5]]), static=static)
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

##' Obtain reduced formulas and parameter vectors for earlier time-points
##'
##' @inheritParams process_inputs
##' @param ordering an ordering for simulation of variables
##' @param t the current time-point
##' @param done character vector of variables already simulated
##' @param vars_t stems for time-varying variables
##' @param kwd keyword for the copula
##'
##'
curr_inputs <- function (formulas, pars, t, ordering, done, vars_t, kwd) {
  start_at <- 0

  ## function to modify arguments
  mod_args <- function (form, beta, t, modLHS=FALSE) {
    trms <- terms(form)
    LHS <- lhs(form)
    trm_labs <- attr(trms, "term.labels")

    if (modLHS) {
      # LHSn <- paste0(regex_extr("_l([0-9]+)", LHS), "_", t)
      rgx <- regex_extr("_l([0-9]+)$", LHS)[[1]] # obtain lag
      if (nchar(rgx[[1]]) > 0) {
        LHS_lag <- as.numeric(substr(rgx, 3, nchar(rgx)))
        if (LHS_lag > t - start_at) {
          ## CODE TO DROP THIS FORMULA
          return(list(form=NA, beta=NA))
        }
        else {
          LHS_new <- sub("_l([0-9]+)$", paste0("_", t-LHS_lag), LHS)
        }
      }
      else LHS_new <- paste0(LHS, "_", t)
      form <- update.formula(form, paste0(LHS_new, " ~ ."))
    }
    # nos <- regex_extr("_l([0-9]+)$", trm_labs) # find lagged variables
    # form <- update.formula(form, as.formula(paste0(lhs(form), "_", t, " ~ .")))

    nos <- regex_extr("_l([0-9]+)$", trm_labs) # find lagged variables
    nos <- rapply(nos, function(x) substr(x, 3, nchar(x)))  # extract lags
    nos <- lapply(nos, as.numeric)
    drp <- sapply(nos, function(x) any(x > t-start_at))
    if (any(na.omit(drp))) {
      intc <- attr(trms, "intercept")
      if (!any(is.na(drp)) && all(drp)) {
        if (intc > 0) form <- update.formula(form, . ~ 1)
        else form <- update.formula(form, . ~ 0)
        beta <- beta[intc]
        # attr(form, "beta") <- pars[[LHS]]$beta
        return(list(form=form, beta=beta))
      }
      trms <- drop.terms(trms, which(drp), keep.response=TRUE)
      beta <- beta[c(intc, which(!drp | is.na(drp))+intc)]
    }

    chr <- as.character(trms)[3]
    wh <- gregexpr("_l([0-9]+)", chr)[[1]]
    if (wh[1] < 0) {
      # attributes(trms) <- list(class="formula")
      return(list(form=update.formula(form, paste0(". ~ ", chr)), beta=beta))
    }
    ml <- attr(wh, "match.length")
    if (any(ml < 3)) stop("All matches should be at least three characters")
    nos <- integer(length(wh))

    for (i in seq_along(wh)) {
      nos[i] <- as.numeric(substr(chr, wh[i]+2, wh[i]+ml[i]-1))
    }
    for (i in seq_along(wh)) {
      chr <- sub("_l([0-9]+)", paste0("_", t-nos[i]), chr)
    }

    list(form=update.formula(form, paste0(". ~ ", chr)), beta=beta)
  }

  # for (i in seq_along(formulas)) for (j in seq_along(formulas[[i]])) {
  #   LHS <- lhs(formulas[[i]][[j]])
  #   tmp <- mod_args(formulas[[i]][[j]], pars[[LHS]]$beta, t=t, modLHS = TRUE)
  #   formulas[[i]][[j]] <- tmp$form
  #   pars[[LHS]]$beta <- tmp$beta
  # }

  # mod_form <- function (form, beta, t, modLHS) mod_args(form, beta=beta, t, modLHS)[[1]]
  # mod_pars <- function (beta, form, t, modLHS) mod_args(form=form[1:3], beta=beta, t, modLHS)[[2]]
  #
  # betas <- purrr::transpose(pars)$beta
  # betas <- betas[names(betas) != kwd]
  #
  # forms <- rapply(formulas, mod_form, pars=pars, how = "replace", t=t, modLHS=TRUE)
  # pars2 <- purrr::flatten_list(rapply(betas, mod_pars, how = "list", t=t, modLHS=TRUE))
  # names(pars2) <- names(pars)[seq_along(pars2)]
  #
  # LHSs <- rapply(formulas, lhs, how="list")

  # ## deal with coefficient vectors
  # ## temporary solution for missing data, just assume 0 contribution
  # trms <- lapply(unlist(formulas), terms)
  dZ <- length(formulas[[1]])
  dX <- length(formulas[[2]])
  #
  # forms <- rapply(formulas, mod_form, pars=pars, how = "replace", t=t, modLHS=TRUE)
  # # forms[[4]] <- rapply(formulas[[4]], mod_form, pars=pars, how = "replace", t=t, modLHS=TRUE)
  #
  ## get list element and make sure cycles correctly
  for (i in seq_along(ordering)) {
    io <- ordering[i]
    wh <- 1 + 1*(io > dZ) + 1*(io > dZ+dX)
    i2 <- i - dZ*(io > dZ) - dX*(io > dZ+dX)
    norm <- dZ*(wh > 1) + dX*(wh > 2)

    frm <- formulas[[wh]][[i2]]
    LHS <- lhs(frm)
    prs <- pars[[LHS]]$beta
    tmp <- mod_args(frm, prs, t = t, modLHS = TRUE)
    formulas[[wh]][[i2]] <- tmp$form
    pars[[LHS]]$beta <- tmp$beta

    ## now look at copula parameters for this variable
    form_copY <- formulas[[4]][[LHS]]
    if (is.null(form_copY)) next
    pars_copY <- pars[[kwd]][[LHS]]
    LHSs <- lhs(form_copY)

    for (j in seq_along(form_copY)) {
      form_cops <- c()
      par_cops <- list()
      for(dt in 0:t){
        tmp <- mod_args(form_copY[[j]], pars_copY[[LHSs[j]]]$beta, t = dt, modLHS = TRUE)
        form_cops <- c(form_cops, tmp$form)
        par_cops <- c(par_cops, list(tmp$beta))
      }
      
      form_copY[[j]] <- form_cops
      pars_copY[[LHSs[j]]]$beta <- par_cops
    }
    formulas[[4]][[LHS]] <- form_copY
    pars[[kwd]][[LHS]] <- pars_copY
  }

  ## adjust parameter names
  nms <- names(pars)
  wh <- which(nms != kwd & nms %in% vars_t)
  nms[wh] <- paste0(nms[wh], "_", t)
  names(pars) <- nms

  # pars <- list_names(pars2, vars_t, t=t, kwd=kwd, start_at=start_at)

  # LHSs <- rapply(forms, lhs, how = "replace")
  # nms <- names(pars)
  # wh <- nms %in% vars_t # <- paste0(names(pars), "_", t)
  # nms[wh] <- paste0(nms[wh], "_", t)
  # names(pars) <- nms
  #
  # names(pars[[kwd]]) <- paste0(names(pars[[kwd]]), "_", t)
  # pars[[kwd]] <- lapply(pars[[kwd]], add_time_stamps)

  # LHSs <- rapply(forms, lhs, how = "replace")
  # forms2 <- mapply(function(x,y) `attr<-`(x, "beta", pars[[y]]$beta),
  #                 x=unlist(forms),
  #                 y=unlist(LHSs))
  # forms <- relist(forms2, skeleton = forms)
  # forms <- list_flatten(forms)
  # forms <- list_flatten(forms, skeleton = LHSs)
  # betas <- rapply(forms, function(x) `attr<-`(x, "beta", ))
  # pars[lhs] <- betas

  return(list(formulas=formulas, pars=pars))
}

# ##' Adds time-points to variable names given as names of a vector
# list_names <- function (x, vars_t, t, kwd, start_at=0) {
#   ## rename for univariate parameter vectors
#   wh_cop <- which(names(x) == kwd)
#   univ <- x[-wh_cop]
#   wh <- which(names(univ) %in% vars_t)
#   nms <- names(univ)
#   nms <- paste0(nms[wh], "_", t)
#   names(univ) <- nms
#
#   ## now consider copula parameters
#   cop_pars <- x[[kwd]]
#   names(cop_pars) <- paste0(names(cop_pars), "_", t)
#
#   ## apply time stamps to pair-copula parameters
#   cop_pars <- lapply(cop_pars, function(x) {
#     names(x) <- add_time_stamps(names(x), t=t, start_at=start_at)
#     return(x)
#   })
#
#   ## recreate pars object
#   out <- c(univ, list(cop_pars))
#   nms <- names(out)
#   nms[length(nms)] <- kwd
#   names(out) <- nms
#
#   return(out)
# }
