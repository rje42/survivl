##' Simulate a single variable using the inversion method. We note that
##' in Xi et al paper notation in the paper never has Y_0, treatment at time
##' t = 0 results in outcome at t = 1, Y_1. For organizational reasons, we
##' assume treatment at time t = 0 results in outcome at time t = 0 and so forth.
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##' @param type_event character denoting if primary or competing
##' @param num_time_varying number of time varying covariates
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, 
                          dat, quantiles, type_event, num_time_varying) {
  ## get variables

  vnm_t <- lhs(formulas[[1]])
  if (length(vnm_t) != 1L) stop("Unable to extract variable name")

  LHS_cops <- sapply(formulas[[2]], \(x) lhs(x))
  k <- as.numeric(substr(vnm_t, nchar(vnm_t), nchar(vnm_t)))

  time_vars <- unique(sapply(as.character(unlist(formulas[[2]])), function(x) sub("_.*", "", x)))
  p <- num_time_varying

  vnm <- rmv_time(vnm_t)
  # make a (k+1) x p upper matrix (list) of copula formulas, pars, and families
  cop_fams <- array(vector("list", (k+1)*p), dim = c(p, k+1))
  cop_forms <- array(vector("list", (k+1)*p), dim = c(p, k+1))
  cop_pars <- array(vector("list", (k+1)*p), dim = c(p, k+1))
  idx <- 1

    for (j in seq_len(k+1)) {
      for (l in seq_len(p)) {
        cop_fams[[l, j]] <- family[[2]][l]
        cop_forms[[l, j]] <- formulas[[2]][[l]][[idx]]
        cop_pars[[l, j]] <- list(beta = pars$cop[[1]][[l]]$beta[[idx]], 
                                    df =  pars$cop[[1]][[l]]$df)
        
      }
      idx <- idx + 1
      
  }

  if(type_event == "primary"){
    # Pass all needed variables to survival_loop
    quantiles <- survival_loop(k=k, p=p, vnm=vnm, time_vars=time_vars, 
                             quantiles=quantiles, family=cop_fams, 
                             pars=cop_pars, formulas = cop_forms, data = dat)
    
  }else if(type_event == "competing"){
    quantiles <- competing_loop(k=k, p=p, vnm=vnm, time_vars=time_vars, 
                                quantiles=quantiles, family=family[[2]], 
                                pars=pars$cop[[1]], X = X)
  }else{
    stop("Must be one of primary or competing")
  }

  qY <- quantiles[[vnm]]
  ## now rescale to correct margin
  X_do <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)

  Y <- rescale_var(qY, X=X_do, family=family[[1]], pars=pars[["doPar"]], link=link[[1]])
  
  eta <- X_do %*% pars[["doPar"]]$beta
  phi <- pars[["doPar"]]$phi
  lambda <- (link_apply(eta, link[[1]], family = family[[1]]$name)) # not working with ordinal or categorical
  for(i in 1:p){
    pdist_pars <- pars[["doPar"]]; pdist_pars$x <- 1; 
    pdist_pars$beta <- NULL; pdist_pars$mu <- lambda;
    quantiles[[paste0("q", k, "_", 0, "_", i)]] = do.call(family[[1]]$pdist, 
                                                          pdist_pars)
  }
  
  
  if(k > 0){
    for(j in 0:(k-1)){
      for(i in 1:p){
        if(i > 1){
          L_col <- paste0(time_vars[[i]], "|", time_vars[[1:(i-1)]], "_", (j))
        }
        else{
          L_col <- paste0(time_vars[[1]], "_", (j))
        }
        
        Y_col <- paste0("q", k, "_", (j), "_", i) # q_{k,j,i} where i is the multiple time-varying
        insert_col <- paste0("q", k, "_", (j+1), "_", i)
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        idxs <- c(i, (j+1)) 
        form <- do.call("[[", c(list(cop_forms), as.list(idxs)))
        X <- model.matrix(form, data = dat)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, cop_fams, 
                                                            cop_pars, idxs, inv = FALSE)
        }
        
        
      }
  }
 
    



  # Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm_t]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}



survival_loop <- function(k, p, vnm, time_vars, quantiles, 
                          family, pars, formulas, data){
  # Get the Quantiles of Y and L

  if (k >0){
    #get the distributions of Ls in the survivors
    for (j in rev(0:(k-1))){
        for(i in rev(seq_len(p)[-1])){
          L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_",j, "_prev")
          Y_col = paste0(paste0("q", (k-1), "_", j, "_", i, "_prev"))
          
          qs <- cbind(
            quantiles[[L_col]],
            quantiles[[Y_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          idxs <- c(i, (j+1))
          form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
          numerator <- quantiles[[L_col]] - compute_copula_quantiles(qs, X, family,
                                                                     pars, idxs, inv = F, p = TRUE)
          

          quantiles[[insert_col]] <-numerator / (1 - quantiles[[Y_col]]) 
        }
        L_col = paste0(time_vars[1], "_", j, "_prev")
        Y_col = paste0("q", (k-1), "_", j, "_", "1_prev")
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        idxs <- c(1, (j+1))
        form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
        numerator <- quantiles[[L_col]] - compute_copula_quantiles(qs, X, family,
                                                                   pars, idxs, inv = F, p = TRUE)
        
        
        quantiles[[insert_col]] <-numerator / (1 - quantiles[[Y_col]]) 
      }
    }
  
  # now get the quantiles of the Ys
  for (j in rev(0:k)){ # in paper use k+1, we're not so go from 1, k
    if(j > 0){
      Y_col = paste0(vnm, "|", paste0(time_vars[1:p], collapse = ""), "_", j)
      for(i in rev(seq_len(p)[-1])){
        L_col = paste0(time_vars[p], "|", paste0(time_vars[1:(p-1)], collapse = ""), "_", j)
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        if(j == 0){
          insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                              "_", j)
        }else{
          insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                              "_", j, "_", time_vars[i],"_", j-1)
        }
        idxs <- c(i, (j+1))
        form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                            pars, idxs, inv = T)
        
        Y_col = insert_col
    }
    
    L_col = paste0(time_vars[1], "_", j)
    # Y_col is still the last column we calculated
    qs = cbind(
      quantiles[[L_col]],
      quantiles[[Y_col]]
    )
    insert_col = insert_col = paste0(vnm, "|", paste0(time_vars, collapse = ""), "_", j-1)
    idxs <- c(1, (j+1))
    form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
    quantiles[[insert_col]] = compute_copula_quantiles(qs, X, family, pars, idxs, inv = T)
    }else{
      
    for(i in rev(seq_len(p)[-1])){
      L_col = paste0(time_vars[p], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0")
      Y_col = paste0(vnm, "|", paste0(time_vars[1:i], collapse = ""), "_0")
      qs <- cbind(
        quantiles[[L_col]],
        quantiles[[Y_col]]
      )
      insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                          "_", j)
      idxs <- c(i, (j+1))
      form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
      quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                          pars, idxs, inv = T)
    }
      L_col <- paste0(time_vars[1], "_0")
      Y_col <- paste0(vnm, "|", time_vars[1], "_", j)
      insert_col <- vnm
      qs <- cbind(
        quantiles[[L_col]],
        quantiles[[Y_col]]
      )
      idxs <- c(1, (j+1))
      form <- do.call("[[", c(list(formulas), as.list(idxs))); X <- model.matrix(form, data = data)
      quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                          pars, idxs, inv = T)
      

    }
    
    
  }
  return(quantiles)
  }



competing_loop <- function(k, p, vnm, time_vars, quantiles, family, pars, X){
  # Get the Quantiles of Y and L
  if (k >0){
    #get the distributions of Ls in the survivors
    for (j in (1:k)){
      if (j<k){
        
        for(i in rev(seq_len(p)[-1])){
 
          Y_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                         "_", k-j, "_", time_vars[i],"_", k-j-1, "_prev")
          L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k-j, "_prev")
          qs <- cbind(
            quantiles[[Y_col]],
            quantiles[[L_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                              pars, i, inv = F)
        }
        Y_col = paste0(vnm, "|", paste0(time_vars, collapse = ""), "_", k-j-1, "_prev")
        L_col = paste0(time_vars[1], "_",k-j,"_prev")
        qs <- cbind(
          quantiles[[Y_col]],
          quantiles[[L_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family, pars, 1, F)
        
      } else{
        
        for(i in rev(seq_len(p)[-1])){
          L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0", "_prev")
          Y_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0", "_prev")
          qs <- cbind(
            quantiles[[Y_col]],
            quantiles[[L_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                              pars, i, inv = F)
        }
        L_col = paste0(time_vars[1], "_0_prev")
        Y_col = "Y_prev"
        qs <- cbind(
          quantiles[[Y_col]],
          quantiles[[L_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family, pars, 1, F)
        
      }
    }}
  
  
  for (j in (0:k)){
    
    if (j<k){
      Y_col = paste0(vnm, "|", paste0(time_vars[1:p], collapse = ""), "_", k-j)
      for(i in rev(seq_len(p)[-1])){
        L_col = paste0(time_vars[p], "|", paste0(time_vars[1:(p-1)], collapse = ""), "_", k-j)
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                            "_", k-j, "_", time_vars[i],"_", k-j-1)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                            pars, i, inv = T)
        
        Y_col = insert_col
      }
      L_col = paste0(time_vars[1], "_", k-j)
      # Y_col is still the last column we calculated
      qs = cbind(
        quantiles[[L_col]],
        quantiles[[Y_col]]
      )
      insert_col = paste0(vnm, "|", paste0(time_vars, collapse = ""), "_", k-j-1)
      quantiles[[insert_col]] = compute_copula_quantiles(qs, X, family, pars, 1, inv = T)
      
      
    }else{
      
      for(i in rev(seq_len(p)[-1])){
        L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0")
        Y_col = paste0(vnm, "|", paste0(time_vars[1:i], collapse = ""), "_0")
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k-j)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, X, family,
                                                            pars, i, inv = T)
      }
      L_col <- paste0(time_vars[1], "_0")
      Y_col <- paste0(vnm, "|", time_vars[1], "_", k-j)
      insert_col <- vnm
      ## rescale quantiles for pair-copula
      qs <- cbind(
        quantiles[[L_col]],
        quantiles[[Y_col]]
      )

      qY <- compute_copula_quantiles(qs, X, family, pars, 1, inv = T)
      quantiles[[insert_col]] <- qY
      
    }
  }
  return(quantiles)
}

##' A wrapper function to compute h or inverse h functions. Also can apply a
##' copula CDF if p = TRUE
##' @importFrom copula cCopula pCopula
##' @importFrom copula normalCopula tCopula
##' @importFrom copula claytonCopula gumbelCopula frankCopula joeCopula fgmCopula
##' @param qs a matrix of two quantiles qs = (q1, q2) as we are using bivariate copulas.
##' @param family A list of numeric copula family 1,2,3,4,5 see `copula_vals` for more.
##' @param pars Parameter for the copula.
##' @param i Use i as an index for which copula we are working with if there are multiple
##' @param inv Boolean if inverse or regular h function
##' @param p Boolean if we want to compute a cdf instead of h function required for survival outcomes.
##' @export
compute_copula_quantiles <- function(qs, X, family, pars, idxs, inv, p = FALSE) {
  

  # pin from 1 to 1-eps for eps = 0.0005
  fam <- family[[idxs[1], idxs[2]]]
  copPars <- pars[[idxs[1], idxs[2]]]
  beta <- copPars$beta
  df <- copPars$df
  if(fam == 11){stop("Not supported for fgm")}

  qY <- rescale_cop(qs, X, beta, family = fam, df = df, inv = inv, cdf = p)

  return(qY)
}


# pCopNorm <- function(U, param){
#   VGAM::pbinormcop(U[,1], U[,2], rho = param)
# }

