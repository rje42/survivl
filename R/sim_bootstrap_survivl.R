##' Simulate survivl outcome using inversion method
##'
##' @param out data frame for output
##' @param surv_model the object surv_model
## @param control list of control parameters
##'
##'
##' @export
sim_bootstrap_survivl <- function(out, surv_model) {
  mod_inputs <- modify_inputs(surv_model)
  formulas <- mod_inputs$formulas
  LHS_C <- surv_model$LHSs$LHS_C
  done <- unlist(LHS_C)


  rho <- unlist(surv_model$pars$cop)
  rho <- 2 * rje::expit(rho) - 1
  formulas <- mod_inputs$formulas
  pars <- surv_model$pars
  order <- surv_model$ordering
  kwd <- surv_model$kwd
  vars_t <- surv_model$vars_t
  link <- surv_model$link
  qtls <- surv_model$qtls
  surv <- rep(TRUE, nrow(out))
  riskH <- surv_model$risk$risk_h
  risk_form <- surv_model$risk$risk_form
  # # eps <- con$eps
  M <- surv_model$bootsims
  direct_risk <- mod_inputs$family[[3]][[1]]$name == "binomial"
  T <- surv_model$T
  vnm <- surv_model$LHSs$LHS_Y
  dC <- surv_model$dC
  for (t in seq_len(T) - 1) {
    ## function to standardize formulae
    mod_inputs$t <- t
    cinp <- curr_inputs(
      formulas = formulas, pars = pars, ordering = order,
      done = done, t = t, vars_t = vars_t, kwd = kwd
    )
    mod_inputs$formulas <- cinp$formulas

    mod_inputs$pars <- cinp$pars

    done <- c(done, paste0(vars_t, "_", t))

    ## simulate A_l, Z_l
    tmp <- sim_block(out[surv, ], mod_inputs,
      quantiles = qtls[surv, , drop = FALSE], kwd = kwd,
      sim_Y = FALSE
    )
    out[surv, ] <- tmp$dat

    # generate probabilites for Y_0
    lambda <- unlist(mod_inputs$pars[[paste0("Y_", t)]]$beta)

    trm <- terms(terms(mod_inputs$formulas[[3]][[1]]))
    MM <- model.matrix(delete.response(trm),
      data = out[surv, ]
    )
    if (nrow(MM) != nrow(out[surv, ])) {
      if (length(attr(trm, "factors")) == 0) {
        if (attr(trm, "intercept") == 1) {
          MM <- matrix(1, nrow = nrow(out), ncol = 1)
        } else {
          MM <- matrix(0, nrow = nrow(out), ncol = 0)
        }
      } else {
        warning(paste0("Missing entries for ", vnm))
      }
    }
    eta <- MM %*% pars[[vnm]]$beta
    curr_phi <- pars[[vnm]]$phi
    curr_link <- mod_inputs$link[[3]]

    if (direct_risk) {
      g_t <- link_apply(eta, curr_link, "binomial")
    } else {
      Y_fam <- mod_inputs$family[[3]][[1]]
      # TODO: continuos time for seamon and keogh
      if (Y_fam$name == "Gamma") {
        lambda <- causl::glm_sim(Y_fam,
          eta = eta, phi = curr_phi,
          other_pars = pars[[vnm]]
        )
        g_t <- 1 - exp(-lambda)
      } else if (Y_fam$name == "weibull") {
        lambda <- causl::glm_sim(Y_fam,
          eta = eta, phi = curr_phi,
          other_pars = pars[[vnm]]
        )
        # g_t <- (1/lambda^shape) ((t+1)^shape - t^shape)
        shape <- pars[[vnm]]$shape
        # sometimes really small so pin to eps
        exp_lambda_t <- pmax(exp(-(lambda^(-shape) * ((t + 1)^shape - t^shape))), 0.005)
        g_t <- 1 - exp_lambda_t
      } else if (Y_fam$name == "lognormal") {
        stop("Numerical integration required for lognormal hazard. \n
             Not supported for Seaman + Keogh method.")
      } else {
        (stop("Distribution must be one a survival family."))
      }
    }

    # create risk function model matrix interventions, covariates, baseline covariates
    risk_form_t <- delete.response(terms(mod_args(risk_form,
      beta = rep(1, length(rhs_vars(risk_form)[[1]])),
      t = t, modLHS = FALSE
    )$form))
    MM <- model.matrix(risk_form_t, data = out[surv, ])
    MM <- MM[, 2:ncol(MM), drop = FALSE] # get rid of intercept only want ranking
    # get scalar H
    H <- apply(MM, 1, riskH)

    ## simulate M times

    ### First get the current model matrix:
    ## unpack mod_inputs
    vars <- paste0(mod_inputs$vars_t, "_", mod_inputs$t)
    d <- lengths(mod_inputs$formulas)
    ## simulate M times

    empirical_ys <- matrix(NA, nrow = sum(surv), ncol = M)
    # if (surv_model$verbose) cat("Bootstrap iteration: ")


    for (j in seq_len(M)) {
      H_j <- empirical_cdf_bootstrap(
        mod_inputs,
        out[surv, ], riskH
      )
      empirical_ys[, j] <- H_j
      printCount(j, first = 1, last = M)
    }

    Hs <- matrix(H, nrow = sum(surv), ncol = M)
    us <- (rowSums(empirical_ys < Hs) + 1) / (M + 1) - runif(sum(surv), max = 1 / (M + 1))


    # now with quantile and correlation coef, use copula, to get probability
    zs <- qnorm(us)
    zY <- rnorm(sum(surv), rho * zs, 1 - rho^2)
    probs_sim <- pnorm(zY)
    y_k1 <- as.numeric(probs_sim < g_t)


    tmp$dat[[paste0("Y_", t)]] <- y_k1
    surv_this <- y_k1 == 0
    if (!direct_risk) {
      cont_t <- t + log(1 - probs_sim[!surv_this]) / log(1 - g_t[!surv_this])
    } else {
      cont_t <- t + 1
    }

    out[surv, ] <- tmp$dat

    # surv[surv] <- surv[surv] & surv_this

    # ## determine if the individual had an event
    # indYt <- dC + (t-con$start_at)*length(vars_t) + dZ + dX + seq_len(dY)  # indices of responses
    # if (dY == 1) {
    #   surv_this <- out[surv, indYt] == 0
    # }
    # else {
    #   surv_this <- apply(out[surv, indYt] == 0, 1, all)
    # }
    # ## get time of event and which one
    out$T[surv][!surv_this] <- cont_t

    ## update list of survivors
    surv[surv] <- surv[surv] & surv_this

    ## if no-one has survived, then end the simulation
    if (!any(surv)) break
  }
  return(out)
}
