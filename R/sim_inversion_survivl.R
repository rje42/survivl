##' Simulate survivl outcome using inversion method 
##'
##' @param out data frame for output
##' @param surv_model the object surv_model
## @param control list of control parameters
##'
##'
##' @export
sim_inversion_survivl <- function(out, surv_model) {
  mod_inputs <- modify_inputs(surv_model)
  formulas <- mod_inputs$formulas

  T <- surv_model$T
  pars <- surv_model$pars; order <- surv_model$ordering; 
  LHS_C <- surv_model$LHSs$LHS_C; done <- unlist(LHS_C)
  vars_t <- surv_model$vars_t; kwd <- surv_model$kwd
  qtls <- surv_model$qtls; surv <- rep(TRUE, nrow(out))
  dC <- surv_model$dC; dZ <- surv_model$dZ; dX <- surv_model$dX; dY <- surv_model$dY
  censoring <- surv_model$censoring

  for (t in seq_len(T)-1) {
    rje::printCount(t)
    ## function to standardize formulae
    mod_inputs$t <- t
    
    cinp <- curr_inputs(formulas=formulas, pars=pars, ordering=order,
                        done=done, t=t, vars_t=vars_t, kwd=kwd)
    
    
    mod_inputs$formulas <- cinp$formulas

    mod_inputs$pars <- cinp$pars

    done <- c(done, paste0(vars_t, "_", t))
    if(t > 0){
      if(t > 1){
        qtls <- dplyr::select(qtls, -contains("prev"))
      }
      colnames(qtls)[(dC + 1):ncol(qtls)] <- sapply(colnames(qtls)[(dC+1):ncol(qtls)], function(z) paste0(z,"_prev"))
    }
    tmp <- sim_block(out[surv,], mod_inputs, quantiles=qtls[surv,, drop = FALSE], kwd=kwd)
    
    out[surv, ] <- tmp$dat
    if(is.null(qtls)){
      qtls <- tmp$quantiles
    }else{
      for(name in names(tmp$quantiles)){
        qtls[surv, name] <- tmp$quantiles[[name]]
      }
    }
    
    
    
    
    
    ## determine if the individual had an event
    indYt <- dC + (t-surv_model$start_at)*length(vars_t) + dZ + dX + seq_len(dY)  # indices of responses
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
    
    if (mean(surv) < 0.05) {
      warning(paste0("Time: ", t, " Many samples lost. 
                       May result in unstable calculations."))
    }
    if(sum(surv) < 20){
      stop("Time: ", t, " Too many samples lost. Stopping simulation.")
    }
    
    ## if no-one has survived, then end the simulation
    if (!any(surv)) break
  }
  return(out)
}
