##' @importFrom causl sim_variable

sim_block <- function (out, proc_inputs, quantiles, kwd) {

  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y
  # famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]

  vars <- paste0(proc_inputs$var_t, "_", proc_inputs$t)

  d <- lengths(formulas)

  ## simulate covariates and treatments
  for (j in 1:2) for (i in seq_len(d[j])) {
    vnm <- vars[i+(j-1)*length(formulas[[1]])]
    curr_link <- link[[j]][i]

    curr_form <- formulas[[j]][[i]]
    curr_fam <- family[[j]][i]

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
                          par2=pars[[vnm]]$par2)
    out[[vnm]] <- tmp
    quantiles[[vnm]] <- attr(tmp, "quantile")
  }

  vnm <- lhs(formulas[[3]])
  vnm_stm <- remove_time(vnm)
  # vnm_t <- paste0(vnm, "_", proc_inputs$t)

  ## code to get Y quantiles conditional on different Zs
  for (i in seq_along(formulas[[3]])) {
    ## simulate Y variable
    # qY <- runif(n)
    # print(wh)

    forms <- list(formulas[[3]][[i]], formulas[[4]][[i]])
    fams <- list(family[[3]][i], family[[4]][[i]])
    prs <- list(pars[[vnm[i]]], pars[[kwd]][[vnm_stm[i]]])
    if (!is.null(prs[[1]]$lambda0)) prs[[1]]$phi <- 1 #prs[[1]]$phi <- prs[[1]]$lambda0
    lnk <- list(link[[3]][i], list()) # link[[4]][[i]])

    out <- causl::sim_variable(nrow(out), forms, fams, prs, lnk,
                               dat=out, quantiles=quantiles)
    quantiles <- attr(out, "quantiles")
    attr(out, "quantiles") <- NULL
      # out[[vars[dZ+dX+i]]] <- causl:::sim_Y(n, formulas=formulas[[4]][[i]],
      #                                       family=family[[4]][[i]],
      #                                       pars=pars[[kwd]]$beta[[i]],
      #                                       formY = formulas[[3]][[i]],
      #                                       famY=family[[3]][i], parsY=pars[[LHS_Y[i]]],
      #                                       linkY=link[[3]][i], qZ=quantiles, vars=vars,
      #                                       dat=out)
  }
      # for (j in seq_len(dZ)) {
      #   curr_qZ <- qZs[[vars[j]]]
      #   X <- model.matrix(formulas[[4]][[wh]][[j]], data=out)
      #   curr_fam <- family[[4]][wh,j]
      #   curr_par <- pars[[kwd]]$beta[[wh]][[j]]
      #   # eta <- X %*% curr_par
      #   qY <- rescaleCop(cbind(curr_qZ,qY), X=X, pars=curr_par, family=curr_fam) #, link=link[[4]][i,j])
      # }
      # ##
      # X <- model.matrix(formulas[[3]][[wh]], data=out)
      # qY <- rescaleVar(qY, X=X, family=famY[wh], pars=pars[[LHS_Y[wh]]],
      #                  link=link[[3]][wh])
      #
      # out[[vars[order[i]]]] <- qY

  return(list(dat=out, quantiles=quantiles))
}
