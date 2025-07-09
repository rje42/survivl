##' @importFrom causl sim_variable
##' @importFrom causl `lhs<-`
##' @importFrom rje printCount expit
sim_block <- function (out, proc_inputs, quantiles, kwd) {
  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y
  # famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]

  vars <- paste0(proc_inputs$vars_t, "_", proc_inputs$t)
  time_vars <- proc_inputs$vars_t[1:dZ]
  outcome_vars <- proc_inputs$vars_t[(dZ + dX + 1):(dZ+dX+dY)] # TODO: fix this
  survival_vars <- proc_inputs$survival_outcome
  d <- lengths(formulas)
  k = proc_inputs$t
  ## simulate covariates and treatments
  for (j in 1:2) for (i in seq_len(d[j])) {
    vnm_q <- vnm <- vars[i+(j-1)*length(formulas[[1]])]
    if(!all(is.na(out[[vnm]]))){
      next;
    }
    insert_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k)

    if(j ==1 & i > 1){
        vnm_q <- paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k)
    }
    curr_link <- link[[j]][i]

    curr_form <- formulas[[j]][[i]]
    curr_fam <- family[[j]][[i]]

    trm <- terms(curr_form)
    # curr_form2 <- delete.response(terms(curr_form))
    MM <- model.matrix(delete.response(trm), data=out)
    if (nrow(MM) != nrow(out)) {
      if (length(attr(trm, "factors")) == 0) {
        if (attr(trm, "intercept") == 1) MM <- matrix(1, nrow=nrow(out), ncol=1)
        else MM <- matrix(0, nrow=nrow(out), ncol=0)
      }
      else warning(paste0("Missing entries for ", vnm))
    }
    eta <- MM %*% pars[[vnm]]$beta
    curr_phi <- pars[[vnm]]$phi
    tmp <- causl::glm_sim(fam=curr_fam, eta=eta, phi=curr_phi, link=curr_link,
                          other_pars=pars[[vnm]])
    out[[vnm]] <- tmp
    if(is.null(quantiles)){
      quantiles <- data.frame(attr(tmp, "quantile"))
      names(quantiles) = vnm_q
    }else{
      quantiles[[vnm_q]] <- attr(tmp, "quantile")
    }
  }


  vnm <- lhs(formulas[[3]])
  vnm_stm <- rmv_time(vnm)
  # vnm_t <- paste0(vnm, "_", proc_inputs$t)
  first <- TRUE
  prev <- paste0(time_vars, collapse = "")
  for(outcome_vnm in outcome_vars) {
      quantiles[[paste0(outcome_vnm, "|", prev, "_", k)]] <- runif(nrow(out))

      prev <- paste0(prev, outcome_vnm)
    
  }
  
  if(k > 0 & length(time_vars) > 1){
    for(i in length(time_vars):2){
      insert_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k)
      quantiles[[insert_col]] <- runif(nrow(out))
    }
    quantiles[[paste0(time_vars[1], "_", k)]] <- runif(nrow(out))
  }


  ## code to get Y quantiles conditional on different Zs

  for (i in seq_along(formulas[[3]])) {
    ## simulate Y variable
    # qY <- runif(n)
    # print(wh)
    forms <- list(formulas[[3]][[i]], formulas[[4]][[i]])
    fams <- list(family[[3]][[i]], family[[4]][[i]])
    prs <- list(c(pars[[vnm[i]]], list(x=proc_inputs$t)), pars[[kwd]][[vnm_stm[i]]])
    if (!is.null(prs[[1]]$lambda0)) prs[[1]]$phi <- 1 #prs[[1]]$phi <- prs[[1]]$lambda0
    lnk <- list(link[[3]][i], list()) # link[[4]][[i]])
    cop_pars <- pars[grepl("^cop", names(pars))]
    cop_pars[["doPar"]] <- pars[[vnm[i]]]
    if(i == 1){ #TODO: fix logic here
      type_event <- "primary"
    }else{
      type_event <- "competing"
    }
    out <- survivl::sim_variable(nrow(out), forms, fams, cop_pars, lnk,
                                 dat = out, quantiles=quantiles, type_event)
    browser()
    collect_events(out, "Y")
    # out <- causl::sim_variable(nrow(out), forms, fams, prs, lnk,
    #                            dat=out, quantiles=quantiles)
    quantiles <- attr(out, "quantiles")
    attr(out, "quantiles") <- NULL

  }


  return(list(dat=out, quantiles=quantiles))
  }

