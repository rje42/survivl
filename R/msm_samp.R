##' Simulate from Cox MSM Model
##'
##' Obtain samples from a specified Cox Marginal Structural Model
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
##' @return An object of class `survivl_dat` containing the simulated data.
##'
##' @importFrom survival Surv
##'
##' @export
msm_samp <- function (n, dat=NULL, T, formulas, family, pars, link=NULL,
                      method="inversion", control=list()) {

  con = list(verbose=FALSE, max_wt=1, warn=1, cop="cop", censor="Cen",
             start_at=0,
             risk_h = function(row) {
               treat <- row[1]; background <- row[2]; time_vary <- row[3];
               mean <- sum(treat) + sum(background) + sum(time_vary);
               mean
             },
             # eps = 0.05,
             bootsims = 5e3 - 1)

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
  proc_inputs <- process_inputs(formulas=formulas, pars=pars, family=family,
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
  vars_t <- proc_inputs$vars_t; vars <- proc_inputs$vars

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

  out <- data.frame(rep(list(rep(NA, n)), length(LHS_C) + T*length(vars_t)), rep(0,n))
  names(out) <- vars
  out$T <- T

  # RHS_Z <- causl::rhs_vars(formulas[[2]])
  # RHS_X <- causl::rhs_vars(formulas[[3]])

  # ## make sure family entries have a vector of integers
  # if (any(sapply(family, is.null))) family[sapply(family, is.null)] <- list(integer(0))
  # ## now set up link functions
  # link <- causl::link_setup(link, family = family[-(5)], vars=list(LHS_C,LHS_Z,LHS_X,LHS_Y))
  # # link[[4]] <- "inverse"
  #
  # if (is.null(pars$gamma)) pars$gamma <- function(X, beta) X %*% beta
  #
  # tmp <- standardize_formulae(formulas, static=LHS_C)
  # var_nms_Z <- tmp$var_nms_Z
  # var_nms_X <- tmp$var_nms_X
  # var_nms_Y <- tmp$var_nms_Y
  # var_nms_cop <- tmp$var_nms_cop

  qtls <- out[integer(0)]  # data frame of quantiles

  ## simulate static covariates

  for (i in seq_along(LHS_C)) {
    ## now compute etas
    eta <- model.matrix(update(formulas[[1]][[i]], NULL ~ .), data=out) %*% pars[[LHS_C[i]]]$beta
    tmp <- causl::glm_sim(family=family[[1]][i], eta=eta, phi=pars[[LHS_C[[i]]]]$phi,
                          other_pars=pars[[LHS_C[[i]]]], link=link[[1]][i])
    out[[LHS_C[[i]]]] <- tmp
    qtls[[LHS_C[[i]]]] <- attr(tmp, "quantiles")
  }

  surv <- rep(TRUE, nrow(out))

  if (method == "inversion") {
    mod_inputs <- modify_inputs(proc_inputs)
    formulas <- mod_inputs$formulas
    done <- unlist(LHS_C)
    for (t in seq_len(T)-1) {
      # this_time <- data.frame(rep(list(rep(NA,n)), dZ+dX+dY))
      # names(this_time) <- paste0(vars_t, "_", t)
      # out <- cbind(out, this_time)

      ## function to standardize formulae
      mod_inputs$t <- t
      cinp <- curr_inputs(formulas=formulas, pars=pars, ordering=order,
                          done=done, t=t, vars_t=vars_t, kwd=kwd)
      mod_inputs$formulas <- cinp$formulas
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

      ## determine if the individual had an event
      indYt <- dC + (t-con$start_at)*length(vars_t) + dZ + dX + seq_len(dY)  # indices of responses
      if (dY == 1) {
        surv_this <- out[surv, indYt] > 1
      }
      else {
        surv_this <- apply(out[surv, indYt] > 1, 1, all)
      }
      ## get time of event and which one
      out$T[surv][!surv_this] <- t + do.call(pmin, out[surv,][!surv_this, indYt, drop=FALSE])
      wh_fail <- max.col(-out[surv,][!surv_this, indYt, drop=FALSE])
      out$status[surv][!surv_this] <- wh_fail - (censoring)

      ## record 0 for intervals survived/censored, and i for failure due to ith competing risk
      out[surv, indYt] <- 0L
      out[cbind(which(surv)[!surv_this], indYt[wh_fail])] <- 1L

      ## update list of survivors
      surv[surv] <- surv[surv] & surv_this
      # out <- causl::sim_inversion(out, mod2)
      # qZ <- cbind(qZ, attr(out, "qZs"))

      ## if no-one has survived, then end the simulation
      if (!any(surv)) break
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
  else if (method == "bootstrap") {
    mod_inputs <- modify_inputs(proc_inputs)
    formulas <- mod_inputs$formulas
    done <- unlist(LHS_C)

    riskH <- con$risk_h
    # eps <- con$eps
    M <- con$bootsims
    rho <- unlist(proc_inputs$pars$cop)
    rho <- 2 * rje::expit(rho) - 1

    for (t in seq_len(T)-1) {
      ## function to standardize formulae
      mod_inputs$t <- t
      cinp <- curr_inputs(formulas=formulas, pars=pars, ordering=order,
                          done=done, t=t, vars_t=vars_t, kwd=kwd)
      mod_inputs$formulas <- cinp$formulas

      mod_inputs$pars <- cinp$pars

      mod2 <- modify_LHSs(mod_inputs, t=t)
      done <- c(done, paste0(vars_t, "_", t))

      ## simulate A_l, Z_l
      tmp <- sim_block(out[surv,], mod_inputs, quantiles=qtls[surv,], kwd=kwd)
      out[surv,] <- tmp$dat; qtls[surv,] <- tmp$quantiles

      # generate probabilites for Y_0
      lambda <- unlist(mod_inputs$pars[[paste0("Y_", t)]]$beta)
      sub_out <- out[surv,]
      n_surv <- nrow(sub_out)
      covariates <- cbind(sub_out[[paste0("X_", t)]], sub_out[['C']])
      g_t <- exp(-cbind(1,covariates) %*% lambda)
      # g_t[which(g_t <= 0)] <- eps
      # g_t[which(g_t >= 1)] <- 1 -eps

      # get scalar H
      covariates_0 <- cbind(sub_out[[paste0("X_", t)]], sub_out[['C']], sub_out[[paste0("Z_",t)]])
      H <- apply(covariates_0, 1, riskH)

      ## simulate M times

      ### First get the current model matrix:
      ## unpack mod_inputs
      vars <- paste0(mod_inputs$vars_t, "_", mod_inputs$t)
      d <- lengths(mod_inputs$formulas)
      MMs <- list(list(), list())

      ## simulate covariates and treatments
      for (j in 1:2) for (i in seq_len(d[j])) {
        vnm <- vars[i+(j-1)*length(formulas[[1]])]
        curr_link <-  mod_inputs$link[[j]][i]

        curr_form <- mod_inputs$formulas[[j]][[i]]
        curr_fam <-  mod_inputs$family[[j]][[i]]

        trm <- terms(curr_form)

        MM <- model.matrix(delete.response(trm), data=out[surv,])
        MMs[[j]][[i]] <- MM

        if (nrow(MM) != nrow(out[surv,])) {
          if (length(attr(trm, "factors")) == 0) {
            if (attr(trm, "intercept") == 1) MM <- matrix(1, nrow=nrow(out), ncol=1)
            else MM <- matrix(0, nrow=nrow(out), ncol=0)
          }
          else warning(paste0("Missing entries for ", vnm))
        }
      }

      ## simulate M times
      empirical_ys <- matrix(NA, nrow=n_surv, ncol=M)

      if (verbose) cat("Bootstrap iteration: ")
      for(j in seq_len(M)) {
        tmp_j <- boot_sim(out[surv,], mod_inputs, quantiles=qtls[surv,], kwd=kwd, MMs)
        dat <- tmp_j$dat

        covariates_j <- cbind(dat[[paste0("X_", t)]] ,dat[['C']], dat[[paste0("Z_",t)]])
        H_j <- apply(covariates_j, 1, riskH)
        empirical_ys[,j] <- H_j
        if (verbose) printCount(j, first=1, last=M)
      }

      Hs <- matrix(H, nrow=n_surv, ncol=M)

      # empirical_ys <- data.frame(empirical_ys)
      # colnames(empirical_ys) <- seq_len(M)

      # ecdfs <- apply(empirical_ys, 1, function(x) ecdf(x))

      # # apply ith ecdf to ith entry of H
      # us <- numeric(length(ecdfs))
      # for (i in seq_along(ecdfs)) {
      #   us[i] <- ecdfs[[i]](H[i])
      # }
      us <- (rowSums(empirical_ys < Hs) + 1)/(M + 1) - runif(n, max=1/(M+1))

      # # #options(warn=2, error=recover)
      # # then use empirical cdf to get the quantile
      # # truncate if above or below
      # us[which(us >= 1)] <- 1-eps
      # us[which(us <= 0)] <- eps

      # now with quantile and correlation coef, use copula, to get probability
      probs <- pnorm((qnorm(g_t) - rho * qnorm(us)) / (1 - rho^2))

      # now with probability, can sample next value
      y_k1 <- rbinom(length(probs), 1, probs)
      # y_k1 <- sapply(probs, function(x) rbinom(1, 1, x))
      tmp$dat[[paste0("Y_", t)]] <- y_k1
      out[surv,] <- tmp$dat

      ## determine if the individual had an event
      indYt <- dC + (t-con$start_at)*length(vars_t) + dZ + dX + seq_len(dY)  # indices of responses
      if (dY == 1) {
        surv_this <- out[surv, indYt] == 0
      }
      else {
        surv_this <- apply(out[surv, indYt] == 0, 1, all)
      }
      ## get time of event and which one
      out$T[surv][!surv_this] <- t + do.call(pmin, out[surv,][!surv_this, indYt, drop=FALSE])
      wh_fail <- max.col(-out[surv,][!surv_this, indYt, drop=FALSE])
      out$status[surv][!surv_this] <- wh_fail - (censoring)

      ## record 0 for intervals survived/censored, and i for failure due to ith competing risk
      out[surv, indYt] <- 0L
      out[cbind(which(surv)[!surv_this], indYt[wh_fail])] <- 1L

      ## update list of survivors
      surv[surv] <- surv[surv] & surv_this

      ## if no-one has survived, then end the simulation
      if (!any(surv)) break
    }
  }
  else stop("'method' should be \"inversion\", \"bootstrap\", or \"rejection\"")

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
