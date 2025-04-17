##' Simulate a single variable using the inversion method
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##' @param survival character denoting if primary or competing
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, dat, quantiles, type_event) {
  ## get variables

  vnm_t <- lhs(formulas[[1]])
  if (length(vnm_t) != 1L) stop("Unable to extract variable name")

  LHS_cop <- lhs(formulas[[2]])
  k <- as.numeric(substr(vnm_t, nchar(vnm_t), nchar(vnm_t)))

  time_vars <- sapply(as.character(unlist(formulas[[2]])), function(x) sub("_.*", "", x))
  p <- length(time_vars)
  vnm <- rmv_time(vnm_t)
  X <- model.matrix(formulas[[1]], data=dat)

  if(type_event == "primary"){
    # Pass all needed variables to survival_loop
    quantiles <- survival_loop(k=k, p=p, vnm=vnm, time_vars=time_vars, 
                             quantiles=quantiles, family=family, pars=pars$cop[[1]])
    
  }else if(type_event == "competing"){
    quantiles <- competing_loop(k=k, p=p, vnm=vnm, time_vars=time_vars, 
                                quantiles=quantiles, family=family, pars=pars$cop[[1]])
  }else{
    stop("Must be one of primary or competing")
  }
  
  qY <- quantiles[[vnm]]
  ## now rescale to correct margin
  X <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)
  Y <- rescale_var(qY, X=X, family=family[[1]], pars=pars[["doPar"]], link=link[[1]])
  

  eta <- X %*% pars[["doPar"]]$beta
  lambda <- 1/ (link_apply(eta, link[[1]], family[[1]]))
  quantiles[[paste0("q1_0")]] <- pexp(1, lambda)


  for(j in 0:k){
    for(i in 1:p){
      if(i > 1){
        L_col <- paste0(time_vars[[i]], "|", time_vars[[1:(i-1)]], "_", j)
      }
      else{
        L_col <- paste0(time_vars[[1]], "_", j)
      }
      Y_col <- paste0("q", i, "_", j)
      qs <- cbind(
        quantiles[[L_col]],
        quantiles[[Y_col]]
      )
      insert_col <- paste0("q", i+1, "_", j)
      if(i+1 > p){
        insert_col <- paste0("q1", "_", j+1)
      }
      quantiles[[insert_col]] <- compute_copula_quantiles(qs, family, pars$cop[[1]], i, inv = FALSE)
    }
  }



  # Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm_t]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}



survival_loop <- function(k, p, vnm, time_vars, quantiles, family, pars){
  # Get the Quantiles of Y and L
  if (k >0){
    #get the distributions of Ls in the survivors
    for (j in (1:k)){
        for(i in p:(min(p, 2))){
          if(p == 1){
            break;
          }else{
            L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_",k-j, "_prev")
            Y_col = paste0(paste0("q", i, "_", k-j, "_prev"))
          }
          
          qs <- cbind(
            quantiles[[L_col]],
            quantiles[[Y_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          # quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
          #                                                     pars, i, inv = F)
        
          numerator <- quantiles[[L_col]] - compute_copula_quantiles(qs, family,
                                                                     pars, i, inv = F, p = TRUE)
          

          quantiles[[insert_col]] <-numerator / (1 - quantiles[[Y_col]]) 
        }
        L_col = paste0(time_vars[1], "_", k-j, "_prev")
        Y_col = paste0("q1_", k-j, "_prev")
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        numerator <- quantiles[[L_col]] - compute_copula_quantiles(qs, family,
                                                                   pars, 1, inv = F, p = TRUE)
        
        
        quantiles[[insert_col]] <-numerator / (1 - quantiles[[Y_col]]) 
      }
    }
  
  
  for (j in (0:k)){
    
    if (j<k){
      Y_col = paste0(vnm, "|", paste0(time_vars[1:p], collapse = ""), "_", k-j)
      for(i in p:2){
        if(p == 1){
          break;
        }else{
          L_col = paste0(time_vars[p], "|", paste0(time_vars[1:(p-1)], collapse = ""), "_", k-j)
        }
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                            "_", k-j, "_", time_vars[i],"_", k-j-1)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
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
      quantiles[[insert_col]] = compute_copula_quantiles(qs, family, pars, 1, inv = T)
      
      
    }else{
      
      for(i in p:(min(p, 2))){
        if(p == 1){
          break;
        }else{
          L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0")
          Y_col = paste0(vnm, "|", paste0(time_vars[1:i], collapse = ""), "_0")
        }
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k-j)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
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
      qY <- compute_copula_quantiles(qs, family, pars, 1, inv = T)
      quantiles[[insert_col]] <- qY
      
    }
  }
  return(quantiles)
}


competing_loop <- function(k, p, vnm, time_vars, quantiles, family, pars){
  # Get the Quantiles of Y and L
  if (k >0){
    #get the distributions of Ls in the survivors
    for (j in (1:k)){
      if (j<k){
        
        for(i in p:(min(p,2))){
          if(p == 1){
            break;
          }else{
            Y_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                           "_", k-j, "_", time_vars[i],"_", k-j-1, "_prev")
            L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k-j, "_prev")
          }
          qs <- cbind(
            quantiles[[Y_col]],
            quantiles[[L_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
                                                              pars, i, inv = F)
        }
        Y_col = paste0(vnm, "|", paste0(time_vars, collapse = ""), "_", k-j-1, "_prev")
        L_col = paste0(time_vars[1], "_",k-j,"_prev")
        qs <- cbind(
          quantiles[[Y_col]],
          quantiles[[L_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family, pars, 1, F)
        
      } else{
        
        for(i in p:(min(p, 2))){
          if(p == 1){
            break;
          }else{
            L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0", "_prev")
            Y_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0", "_prev")
          }
          
          qs <- cbind(
            quantiles[[Y_col]],
            quantiles[[L_col]]
          )
          insert_col = gsub("_prev", "", L_col)
          quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
                                                              pars, i, inv = F)
        }
        L_col = paste0(time_vars[1], "_0_prev")
        Y_col = "Y_prev"
        qs <- cbind(
          quantiles[[Y_col]],
          quantiles[[L_col]]
        )
        insert_col = gsub("_prev", "", L_col)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family, pars, 1, F)
        
      }
    }}
  
  
  for (j in (0:k)){
    
    if (j<k){
      Y_col = paste0(vnm, "|", paste0(time_vars[1:p], collapse = ""), "_", k-j)
      for(i in p:2){
        if(p == 1){
          break;
        }else{
          L_col = paste0(time_vars[p], "|", paste0(time_vars[1:(p-1)], collapse = ""), "_", k-j)
        }
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""),
                            "_", k-j, "_", time_vars[i],"_", k-j-1)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
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
      quantiles[[insert_col]] = compute_copula_quantiles(qs, family, pars, 1, inv = T)
      
      
    }else{
      
      for(i in p:(min(p, 2))){
        if(p == 1){
          break;
        }else{
          L_col = paste0(time_vars[i], "|", paste0(time_vars[1:(i-1)], collapse = ""), "_0")
          Y_col = paste0(vnm, "|", paste0(time_vars[1:i], collapse = ""), "_0")
        }
        qs <- cbind(
          quantiles[[L_col]],
          quantiles[[Y_col]]
        )
        insert_col = paste0(vnm, "|", paste0(time_vars[1:(i-1)], collapse = ""), "_", k-j)
        quantiles[[insert_col]] <- compute_copula_quantiles(qs, family,
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
      qY <- compute_copula_quantiles(qs, family, pars, 1, inv = T)
      quantiles[[insert_col]] <- qY
      
    }
  }
  return(quantiles)
}


##' @importFrom copula cCopula
##' @importFrom copula pCopula
compute_copula_quantiles <- function(qs, family, pars, i, inv, p = FALSE) {
  library(copula)

  copula_functions_p <- list(
    #function(U, param, par2, inv) pnorm(qnorm(U[, 2]) * sqrt(1 - param^2) + param * qnorm(U[, 1])),
    function(U, param, par2) pCopula(U, copula = tCopula(param = param, dim = 2, dispstr = "un")),#pCopNorm(U, param),
    function(U, param, par2) pCopula(U, copula = tCopula(param = param, df = par2, dim = 2, dispstr = "un")),
    function(U, param, par2) pCopula(U, copula = claytonCopula(param = param, dim = 2)),
    function(U, param, par2) pCopula(U, copula = gumbelCopula(param = param, dim = 2)),
    function(U, param, par2) pCopula(U, copula = frankCopula(param = param, dim = 2)),
    function(U, param, par2) pCopula(U, copula = joeCopula(param = param, dim = 2))
  )
  copula_functions <- list(    
    #function(U, param, par2, inv) pnorm(qnorm(U[, 2]) * sqrt(1 - param^2) + param * qnorm(U[, 1])),
    function(U, param, par2, inv) cCopula(U, copula = normalCopula(param = param, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = tCopula(param = param, df = par2, dim = 2, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = claytonCopula(param = param, dim = 2), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = gumbelCopula(param = param, dim = 2), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = frankCopula(param = param, dim = 2), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = joeCopula(param = param, dim = 2), inverse = inv)
  )
  # pin from 1 to 1-eps for eps = 0.0005
  fam <- family[[2]][[i]]
  copPars <- unlist(pars[[i]])
  beta <- copPars[1]
  par2 <- copPars[2]
  if(fam != 4){
    beta <- 2 * expit(beta) - 1
  }
  
  if (fam >= 1 && fam <= 6) {
    if(p){

      qY <- copula_functions_p[[fam]](qs, beta, par2 = par2)
    }else{
      qY <- copula_functions[[fam]](qs, beta, par2 = par2, inv)[,2]
      
    }
  } else {
    stop("family must be between 0 and 5")
  }
  return(qY)
}


pCopNorm <- function(U, param){
  VGAM::pbinormcop(U[,1], U[,2], rho = param)
}

