##' Simulate longitudinal non survival outcome using inversion method
##'
##' @param out data frame for output
##' @param surv_model the object surv_model
## @param control list of control parameters
##'
##'
##' @export
sim_inversion_longitudinal <- function(out, surv_model) {
  mod_inputs <- modify_inputs(surv_model)
  formulas <- mod_inputs$formulas

  pars <- surv_model$pars; order <- surv_model$ordering; 
  LHS_C <- surv_model$LHSs$LHS_C; done <- unlist(LHS_C)
  LHS_Y <- surv_model$LHSs$LHS_Y
  vars_t <- surv_model$vars_t; kwd <- surv_model$kwd
  qtls <- surv_model$qtls; surv <- rep(TRUE, nrow(out))
  dZ <- surv_model$dZ
  family <- surv_model$family; link <- surv_model$link
  T <- surv_model$T
  for (t in seq_len(T)-1) {
    mod_inputs$t <- t
    cinp <- curr_inputs(formulas=formulas, pars=pars, ordering=order,
                        done=done, t=t, vars_t=vars_t, kwd=kwd)
    mod_inputs$formulas <- cinp$formulas
    mod_inputs$pars <- cinp$pars
    
    done <- c(done, paste0(vars_t, "_", t))


    tmp <- sim_block(out[surv,], mod_inputs, quantiles=qtls[surv,, drop = FALSE], kwd=kwd)
    out[surv, ] <- tmp$dat
    if(is.null(qtls)){
      qtls <- tmp$quantiles
    }else{
      for(name in names(tmp$quantiles)){
        qtls[surv, name] <- tmp$quantiles[[name]]
      }
    }
  }
  # start with unif(0,1)
  vt <- runif(nrow(out), 0, 1)
  
  # make a (k+1) x dZ upper matrix (list) of copula formulas, pars, and families
  cop_fams <- array(vector("list", (T)*dZ), dim = c(dZ, T))
  cop_forms <- array(vector("list", (T)*dZ), dim = c(dZ, T))
  cop_pars <- array(vector("list", (T)*dZ), dim = c(dZ, T))
  idx <- 1
  for (j in seq_len(T)) {
    for (l in seq_len(dZ)) {
      cop_fams[[l, j]] <- family[[5]][[1]][[l]]
      cop_forms[[l, j]] <- mod_inputs$formulas[[4]][[1]][[l]][[idx]]
      cop_pars[[l, j]] <- list(beta = mod_inputs$pars$cop[[1]][[l]]$beta[[idx]], 
                               df =  mod_inputs$pars$cop[[1]][[l]]$df)
      
    }
    idx <- idx + 1
    
  }
  
  time_vars <- surv_model$LHSs$LHS_Z 
  for(j in rev(seq_len(T)-1)){
    for(i in rev(seq_len(dZ)[-1])){
      L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], 
                                               collapse = ""), "_",j)
      qs <- cbind(
        qtls[[L_col]],
        vt
      )
      
      idxs <- c(i, j+1)
      form <- do.call("[[", c(list(cop_forms), as.list(idxs)))
      X <- model.matrix(form, data = out)
      
      vt <- compute_copula_quantiles(qs, X = X, cop_fams, 
                                     cop_pars, idxs, inv = TRUE)
    }
    L_col <- paste0(time_vars[1], "_", j)
    qs <- cbind(qtls[[L_col]], vt)
    
    idxs <- c(1, j+1)
    form <- do.call("[[", c(list(cop_forms), as.list(idxs)))
    X <- model.matrix(form, data = out)
    
    vt <- compute_copula_quantiles(qs, X = X, cop_fams, 
                                   cop_pars, idxs, inv = TRUE)
    
  }
  
  qY <- vt
  X_do <- model.matrix(delete.response(terms(formulas[[3]][[1]])), data = out)
  Y <- rescale_var(qY, X=X_do, family=family[[4]][[1]], pars=pars[[LHS_Y]], link=link[[4]])
  out[[LHS_Y]] <- Y
  return(out)
}
