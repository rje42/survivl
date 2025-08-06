##' Simulate survivl outcome using rejection method 
##'
##' @param out data frame for output
##' @param surv_model the object surv_model
## @param control list of control parameters
##'
##'
##' @export
sim_rejection_survivl <- function(out, surv_model) {
  
  surv <- rep(TRUE, nrow(out))
  LHS_X <- surv_model$LHSs$LHS_X
  var_nms_X <- surv_model$std_form$var_nms_X
  
  formulas <- surv_model$formulas;pars <- surv_model$pars;
  family <- surv_model$family; link <- surv_model$link
  var_nms_cop <- surv_model$std_form$var_nms_cop
  ## simulate time-varying covariates and survival
  browser()
  for (t in seq_len(T)) {
    OK <- rep(FALSE, nrow(out))
    OK[!surv] <- TRUE
    
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
        
        tmp <- causl::glm_sim(family[[3]][[i]], eta, phi=pars[[LHS_X[[i]]]]$phi,
                              other_pars=pars[[LHS_X[[i]]]], link=link[[3]][i])
        lden[,i] <- attr(tmp, "quantile")
        out[[var]][!OK] <- tmp
        # max_lr[i] <- max(attr(tmp, "quantile"))
      }
      
      vcopd <- var_nms_cop[[i]][!is.na(var_nms_cop[[i]]$lag), ]
      vcops <- var_nms_cop[[i]]$var[is.na(var_nms_cop[[i]]$lag)]
      out2 <- lag_data(out[!OK,,drop=FALSE], t, vcopd, static=vcops)
      browser()
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
  return(out)
}
