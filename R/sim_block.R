##' @importFrom causl sim_variable
##' @importFrom causl `lhs<-`

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

  d <- lengths(formulas)

  ## simulate covariates and treatments
  for (j in 1:2) for (i in seq_len(d[j])) {

    vnm <- vars[i+(j-1)*length(formulas[[1]])]
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
    quantiles[[vnm]] <- attr(tmp, "quantile")
    ## just for recreating XI experiment can fix more later
    if(proc_inputs$t == 0 & vnm == "L_0"){

      out[["L_0"]] = rgamma(dim(out)[1], shape = 25 + 0.8*out$B, scale = 2.5)
      quantiles[["L_0"]] = pgamma(out[["L_0"]], shape = 25 + 0.8*out$B, scale = 2.5)
    }

  }


  vnm <- lhs(formulas[[3]])
  vnm_stm <- rmv_time(vnm)
  # vnm_t <- paste0(vnm, "_", proc_inputs$t)


  ## code to get Y quantiles conditional on different Zs
  k = proc_inputs$t
  quantiles[[paste0("Y_L_", k)]]<- runif(nrow(out))
  for (i in seq_along(formulas[[3]])) {
    ## simulate Y variable
    # qY <- runif(n)
    # print(wh)

    forms <- list(formulas[[3]][[i]], formulas[[4]][[i]])
    fams <- list(family[[3]][[i]], family[[4]][[i]])
    prs <- list(c(pars[[vnm[i]]], list(x=proc_inputs$t)), pars[[kwd]][[vnm_stm[i]]])
    if (!is.null(prs[[1]]$lambda0)) prs[[1]]$phi <- 1 #prs[[1]]$phi <- prs[[1]]$lambda0
    lnk <- list(link[[3]][i], list()) # link[[4]][[i]])


    out <- survivl::sim_variable(nrow(out), forms, fams, prs, lnk,
                                 dat = out, quantiles=quantiles)
    # out <- causl::sim_variable(nrow(out), forms, fams, prs, lnk,
    #                            dat=out, quantiles=quantiles)
    quantiles <- attr(out, "quantiles")
    attr(out, "quantiles") <- NULL

  }

  return(list(dat=out, quantiles=quantiles))
}
