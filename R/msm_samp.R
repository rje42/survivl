##' Simulate from marginal structural (survival) model
##'
##' Obtain samples from a specified (survival-valued) marginal structural model,
##' using the frugal parameterization.
##'
##' @param n number of samples
##' @param dat optional data frame for plasmode simulation
##' @param T number of time points
##' @param formulas list of formulas to use
##' @param family list of families to use
##' @param pars list of parameter settings to use
##' @param link link functions for GLM-like models
##' @param control list of control parameters
##' @param method sampling method (defaults to `"inversion"`)
##'
##' @details Samples from a Marginal Structural Model specified
##' by a frugal parameterization; that is, one must specify the
##' marginal model (i.e. the dependence of the outcome on the
##' treatments), the distribution of the covariates, observed
##' confounders and treatments, and then a copula to join the distribution
##' of the outcome to that of the confounders.
##'
##' Among the left-hand sides of outcome variables, the variable `Cen` has a
##' special meaning as censoring.  This keyword can be changed to something
##' else by using the argument `censor` in the `control` list.
##'
##' Times for presimulated events are put at halfway between the two timepoints;
##' this is a hack.
##'
##' Setting the `control` argument `surv` to be false means that the data is
##' just longitudinal.
##'
##' @return An object of class `survivl_dat` containing the simulated data.
##'
##' @importFrom survival Surv
##'
##' @export
msm_samp <- function (n, dat=NULL, T, formulas, family, pars, link=NULL,
                      method="inversion", control=list()) {

  con = list(verbose=FALSE, max_wt=1, warn=1, cop="cop", censor="Cen", surv=TRUE, start_at=0)
  matches = match(names(control), names(con))
  con[matches] = control[!is.na(matches)]
  if (any(is.na(matches))) warning("Some names in control not matched: ",
                                   paste(names(control[is.na(matches)]),
                                         sep = ", "))
  verbose <- con$verbose
  kwd <- con$cop

  ## infer sample size from dat if possible
  datNULL <- is.null(dat)
  if (datNULL && missing(n)) stop("Must specify a sample size or supply a data frame")
  else if (missing(n)) n <- nrow(dat)
  else if (!datNULL && n != nrow(dat)) {
    warning("Dimension of supplied data does not match n, using size of supplied data frame")
    n <- nrow(dat)
  }

  ## process inputs
  proc_inputs <- process_inputs(n=n, formulas=formulas, pars=pars, family=family,
                                link=link, dat=dat, T=T, control=con, method=method)
  ## extract them again
  formulas <- proc_inputs$formulas; pars <- proc_inputs$pars; family <- proc_inputs$family; link <- proc_inputs$link

  LHS_C <- proc_inputs$LHSs$LHS_C; LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y
  dC <- length(LHS_C); dZ <- length(LHS_Z); dX <- length(LHS_X); dY <- length(LHS_Y)
  var_nms_Z <- proc_inputs$std_form$var_nms_Z
  var_nms_X <- proc_inputs$std_form$var_nms_X
  var_nms_Y <- proc_inputs$std_form$var_nms_Y
  var_nms_cop <- proc_inputs$std_form$var_nms_cop
  order <- proc_inputs$ordering
  vars_t <- proc_inputs$vars_t
  qtls <- proc_inputs$quantiles
  pres_var <- proc_inputs$pres_var

  # LHS_C <- causl::lhs(formulas[[1]])
  # LHS_Z <- causl::lhs(formulas[[2]])
  # LHS_X <- causl::lhs(formulas[[3]])
  # LHS_Y <- causl::lhs(formulas[[4]])
  # LHS_cop <- causl::lhs(formulas[[5]])
  # dZ <- length(LHS_Z); dX <- length(LHS_X); dY <- length(LHS_Y)
  # if (any(duplicated(c(LHS_C, LHS_Z, LHS_X, LHS_Y, LHS_cop)))) stop("Repeated variable names not allowed")

  ## check for censoring
  censoring <- (con$censor %in% LHS_Y)  # see if censoring variable is in list of Y variables
  if (censoring) {
    ## put censoring first if a competing risk
    LHS_Y <- c(con$censor, setdiff(LHS_Y, con$censor))
  }
  # nms_t <- c(LHS_Z, LHS_X, LHS_Y)
  # nms <- c(LHS_C, outer(vars_t, seq_len(T)-1, paste, sep="_"), "status")

  ## set up output
  out <- data.frame(matrix(NA, ncol=length(proc_inputs$vars), nrow=n))
  names(out) <- proc_inputs$vars
  out$status <- 0

  if (!datNULL) {
    to_back <- names(out)
    out[proc_inputs$dat_kp] <- dat[proc_inputs$dat_kp]
    out <- out[unique.default(c(proc_inputs$dat_kp, to_back))]
    # out <- cbind(dat, out)
  }
  if (anyDuplicated(na.omit(proc_inputs$vars))) stop("duplicated variable names")

  # if (datNULL) out <- data.frame(rep(list(rep(NA, n)), length(LHS_C) + T*length(vars_t)), rep(0,n))
  # else out <- cbind(dat, data.frame(rep(list(rep(NA, n)), length(LHS_C) + T*length(vars_t)), rep(0,n)))
  # names(out) <- vars
  out$T <- T - con$start_at

  ## simulate static covariates
  if (pres_var[1] < con$start_at && pres_var[2] < length(LHS_C)) {
    for (i in seq_along(LHS_C)[pres_var[2] + seq_len(length(LHS_C) - pres_var[2])]) {
      ## now compute etas
      eta <- model.matrix(update(formulas[[1]][[i]], NULL ~ .), data=out) %*% pars[[LHS_C[i]]]$beta
      tmp <- causl::glm_sim(family=family[[1]][[i]], eta=eta, phi=pars[[LHS_C[[i]]]]$phi,
                            other_pars=pars[[LHS_C[[i]]]], link=link[[1]][i])
      out[[LHS_C[[i]]]] <- tmp
      qtls[[LHS_C[[i]]]] <- attr(tmp, "quantiles")
    }
  }

  ## record who should still be simulated
  surv <- rep(TRUE, nrow(out))

  if (method == "inversion") {
    mod_inputs <- modify_inputs(proc_inputs)
    formulas <- mod_inputs$formulas
    done <- unlist(LHS_C)

    for (t in seq_len(T)-1) {

      if (t > pres_var[1] - con$start_at ||
          (t == pres_var[1] - con$start_at && pres_var[2] < length(vars_t))) {
        ## this time point not fully simulated already
        ## function to standardize formulae
        mod_inputs$t <- t
        cinp <- curr_inputs(formulas=formulas, pars=pars, ordering=order,
                            done=done, t=t, vars_t=vars_t, kwd=kwd, start_at=con$start_at)
        mod_inputs$formulas <- cinp$formulas
        if (t == pres_var[1] - con$start_at) mod_inputs$pres_var <- pres_var[2]
        else mod_inputs$pres_var <- 0
        # print(mod_inputs$formulas)
        mod_inputs$pars <- cinp$pars
        # tmp_pars <- rapply(mod_inputs$formulas, function(x) attr(x, "beta"), how="list")
        # while (pluck_depth(tmp_pars) > 2) {
        #   tmp_pars <- list_flatten(tmp_pars)
        # }

        mod2 <- modify_LHSs(mod_inputs, t=t)
        done <- c(done, paste0(vars_t, "_", t))

        ## use sim_inversion()
        tmp <- sim_block(out[surv,], mod_inputs, quantiles=qtls[surv,], kwd=kwd)
        out[surv,] <- tmp$dat; qtls[surv,] <- tmp$quantiles

        if (con$surv) {
          ## establish which survivors had events, and record
          tmp <- collect_events(dat=out[surv,], varY=LHS_Y, t=t)
        }
      }
      else if (con$surv) {
        ## fully simulated already, so just get results
        ## establish which survivors had events, and record
        ## (CHECK THIS IS OK FOR LONGITUDINAL CASE)
        tmp <- collect_events(dat=out[surv,], varY=LHS_Y, t=t, trunc=TRUE)
      }

      if (con$surv) {
        out[surv,] <- tmp$dat

        ## update list of survivors
        surv[surv] <- surv[surv] & tmp$surv

        ## if no-one has survived, then end the simulation
        if (!any(surv)) break
      }
    }
  }
  else if (method == "rejection") {
    ## simulate time-varying covariates and survival
    for (t in seq_len(T)) {
      OK <- rep(FALSE, nrow(out))
      OK[!surv] <- TRUE
      if (verbose) cat("T = ", t)

      while(!all(OK)) {

        lden <- matrix(0, sum(!OK), length(LHS_X))
        # max_lr <- rep(0, length(LHS_X))

        ## now simulate each treatment variable according to formula and pars
        for (i in seq_along(LHS_X)) {
          # out2 <- lag_data(out, t, var_nms_X[[i]], static=LHS_C)
          vXd <- var_nms_X[[i]][!is.na(var_nms_X[[i]]$lag), ]
          vXs <- var_nms_X[[i]]$var[is.na(var_nms_X[[i]]$lag)]
          out2 <- lag_data(out[!OK,,drop=FALSE], t, vXd, static=vXs)
          ## change this to replace with something closer to average value
          out2 <- data.frame(lapply(out2, function(x) {
            if (all(is.na(x))) return(rep(0,length(x)))
            else return(x)
          }))
          # out2[paste(LHS_Z, "0", sep="_l")] <- rep(colMeans(out2[paste(LHS_Z, "0", sep="_l")]), each=nrow(out2))

          ## now compute etas
          eta <- model.matrix(update(formulas[[3]][[i]], NULL ~ .), data=out2) %*% pars[[LHS_X[i]]]$beta
          var <- paste0(LHS_X[[i]], "_", t-1)

          tmp <- causl::glm_sim(family[[3]][i], eta, phi=pars[[LHS_X[[i]]]]$phi,
                                other_pars=pars[[LHS_X[[i]]]], link=link[[3]][i])
          lden[,i] <- attr(tmp, "quantile")
          out[[var]][!OK] <- tmp
          # max_lr[i] <- max(attr(tmp, "quantile"))
        }

        vcopd <- var_nms_cop[[i]][!is.na(var_nms_cop[[i]]$lag), ]
        vcops <- var_nms_cop[[i]]$var[is.na(var_nms_cop[[i]]$lag)]
        out2 <- lag_data(out[!OK,,drop=FALSE], t, vcopd, static=vcops)

        ## now compute model matrix
        copMM <- model.matrix(update(formulas[[5]][[1]], NULL ~ .), data=out2)
        resp <- paste0(c(LHS_Z, LHS_Y), "_", t-1)

        out[!OK,resp] <- causl::sim_copula(out[!OK,resp,drop=FALSE], family[[5]], par=pars$cop, par2=pars$cop$par2, model_matrix = copMM)

        for (i in seq_along(LHS_Z)) {
          vZd <- var_nms_Z[[i]][!is.na(var_nms_Z[[i]]$lag), ]
          vZs <- var_nms_Z[[i]]$var[is.na(var_nms_Z[[i]]$lag)]
          out2 <- lag_data(out[!OK,,drop=FALSE], t, vZd, static=vZs)

          ## now compute etas
          MM <- model.matrix(update(formulas[[2]][[i]], NULL ~ .), data=out2)
          var <- paste0(LHS_Z[[i]], "_", t-1)

          out[[var]][!OK] <- causl::rescale_var(out[[var]][!OK], X=MM,
                                                 family=family[[2]][i], pars=pars[[LHS_Z[i]]],
                                                 link=link[[2]][i])
        }
        for (i in seq_along(LHS_Y)) {
          vYd <- var_nms_Y[[i]][!is.na(var_nms_Y[[i]]$lag), ]
          vYs <- var_nms_Y[[i]]$var[is.na(var_nms_Y[[i]]$lag)]
          out2 <- lag_data(out[!OK,,drop=FALSE], t, vYd, static=vYs)

          ## now compute etas
          MM <- model.matrix(update(formulas[[4]][[i]], NULL ~ .), data=out2)
          var <- paste0(LHS_Y[[i]], "_", t-1)

          out[[var]][!OK] <- causl::rescale_var(out[[var]][!OK], X=MM,
                                                 family=family[[4]][i],
                                                 pars=c(pars[[LHS_Y[i]]], list(phi=1)),
                                                 link=link[[4]][i])
        }

        ## now compute new implied X densities
        lden2 <- matrix(0, sum(!OK), length(LHS_X))
        # max_lr <- rep(0, length(LHS_X))

        ## now simulate each treatment variable according to formula and pars
        for (i in seq_along(LHS_X)) {
          # out2 <- lag_data(out, t, var_nms_X[[i]], static=LHS_C)
          vXd <- var_nms_X[[i]][!is.na(var_nms_X[[i]]$lag), ]
          vXs <- var_nms_X[[i]]$var[is.na(var_nms_X[[i]]$lag)]
          out2 <- lag_data(out[!OK,,drop=FALSE], t, vXd, static=vXs)

          ## now compute etas
          eta <- model.matrix(update(formulas[[3]][[i]], NULL ~ .), data=out2) %*% pars[[LHS_X[i]]]$beta
          var <- paste0(LHS_X[[i]], "_", t-1)

          lden2[,i] <- glm_ldens(out[[var]][!OK], family[[3]][i], eta, link[[3]][i], pars[[LHS_X[[i]]]]$phi)
          # max_lr[i] <- max(attr(tmp, "lden"))
        }

        rat <- lden2 - lden
        rat <- rat - rep(apply(rat, 2, max), each=nrow(rat))

        ## now perform rejection sampling
        Us <- runif(sum(!OK))
        OK[!OK] <- (Us < exp(rat))
        if (verbose) rje::printCount(sum(OK))
      }

      ## update who is still alive and any censoring
      vars <- paste0(LHS_Y, "_", t-1)
      this_T <- apply(out[, vars, drop=FALSE], 1, min)
      surv[surv] <- (this_T[surv] > 1)
      new_fail <- !is.na(out[[vars[1]]]) & !surv
      out$T[new_fail] <- t - 1 + this_T[new_fail]
      out$status[new_fail] <- apply(out[new_fail,vars,drop=FALSE]==this_T[new_fail], 1, which.max) - censoring
      for (i in seq_along(vars)) out[[vars[i]]][!is.na(out[[vars[i]]])] <- 1*(out[[vars[i]]][!is.na(out[[vars[i]]])] < 1)

      # surv[surv] <- (out[[var]][surv] > 1)
      # out$T[!is.na(out[[var]]) & !surv] <- t - 1 + out[[var]][!is.na(out[[var]]) & !surv]
      # out[[var]][!is.na(out[[var]])] <- 1*(out[[var]][!is.na(out[[var]])] > 1)
    }
  }
  else stop("'method' should be \"inversion\" or \"rejection\"")

#   # cum_haz <- rep(0, n)
#   Y_res <- out[[LHS_Y]]

  # out$T <- rep(0, n)
  # gam_aft <- rep(NA, n)
  # surv <- rep(TRUE, n)
#
#   ## now rescale Y to match C, Z and Xs
#   ## for now just assume simple time-homogeneous exponential hazard
#   for (t in seq_len(T)) {
#     na_vars <- grepl(paste0("_",t-1), names(out))
#     out[na_vars][!surv,] <- NA
#
#     vYd <- var_nms_Y[[i]][!is.na(var_nms_Y[[i]]$lag), ]
#     vYs <- var_nms_Y[[i]]$var[is.na(var_nms_Y[[i]]$lag)]
#     out2 <- lag_data(out, t, vYd, static=vYs)
#     # out2 <- lag_data(out, t, var_nms_Y[[1]], static=LHS_C)
#
#     ## now compute etas
#     MM <- model.matrix(formulas[[4]][c(1,3)], data=out2)
#     gam_aft[surv] <- pars$gamma(MM, pars[[LHS_Y[1]]]$beta)
#
#     tmp <- out$T + surv*pmin(1, Y_res/exp(gam_aft))
#     out$T[!is.na(tmp)] <- na.omit(tmp)
#     tmp <- Y_res - exp(gam_aft)
#     Y_res[!is.na(tmp)] <- na.omit(tmp)
#
#     vars <- paste0(LHS_Y[[1]], "_", t-1)
#     out[[var]] <- ifelse(surv, 1*(Y_res < 0), NA)
#
#     surv <- surv & Y_res > 0
#
#     # cum_haz <- cum_haz + exp(gam_aft)
#   }

  out <- cbind(id=seq_len(nrow(out)), out)
  # out$status <- 1*surv
  class(out) <- c("survivl_dat", "data.frame")

  return(out)
}

##' @describeIn msm_samp old name
coxSamp <- function (n, T, formulas, family, pars, link=NULL, control=list()) {
  .Deprecated(new = "msm_samp", msg = "This function is deprecated and will be removed in version 0.4.0, please use msm_samp()")
  msm_samp(n=n, T=T, formulas=formulas, family=family, pars=pars, link=NULL, control=list())
}

##' @describeIn msm_samp old name
cox_samp <- function (n, T, formulas, family, pars, link=NULL, control=list()) {
  deprecate_soft(when = "0.3.1", what = "cox_samp()", with = "msm_samp()")
  # .Deprecated(new = "msm_samp", msg = "This function is deprecated, please use msm_samp()")
  msm_samp(n=n, T=T, formulas=formulas, family=family, pars=pars, link=NULL, control=list())
}

