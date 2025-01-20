##' Simulate a single variable using the inversion method
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, dat, quantiles) {
  ## get variables
  vnm <- lhs(formulas[[1]])
  if (length(vnm) != 1L) stop("Unable to extract variable name")

  LHS_cop <- lhs(formulas[[2]])
  k <- as.numeric(substr(vnm, nchar(vnm), nchar(vnm)))


  for (i in rev(seq_along(formulas[[2]]))) {
    X <- model.matrix(formulas[[2]][[i]], data=dat)
    # eta <- X %*% pars[[2]][[i]]$beta

    # Get the Quantiles of Y and L
    if (k >0){
      #get the distributions of Ls in the survivors
      for (j in (1:k)){
        if (j<k){

          Y_L_col <- paste0("Y_L_", k-j-1,"_prev")
          L_col <- paste0("L_", k-j, "_prev")
          L_col_insert <- paste0("L_",k-j)
          qs <- as.matrix(cbind(quantiles[[Y_L_col]], quantiles[[L_col]]))
          quantiles[[L_col_insert]] <- compute_copula_quantiles(qs, family, pars, i, FALSE)
        } else{
          L_col <- paste0("L_", k-j,"_prev")
          L_col_insert <- paste0("L_",k-j)
          Y_col <- paste0("Y","_prev")

          qs <- as.matrix(cbind(quantiles[[Y_col]],quantiles[[L_col]]))
          quantiles[[L_col_insert]] <- compute_copula_quantiles(qs, family, pars, i, FALSE)

        }
      }
    }
    for (j in (0:k)){
      if (j<k){
        Y_L1_column = paste0("Y_L_",k-j-1)
        Y_L_column = paste0("Y_L_", k-j)
        L_column = paste0("L_", k -j)
        qs <-  as.matrix(cbind(
          quantiles[[L_column]],
          quantiles[[Y_L_column]]
        ))

        quantiles[[Y_L1_column]] = compute_copula_quantiles(qs, family, pars, i, TRUE)

      }else{
        Y_column = "Y"
        ## rescale quantiles for pair-copula
        qs <- as.matrix(
          cbind(quantiles[["L_0"]],
                quantiles[["Y_L_0"]])
        )

        qY <- compute_copula_quantiles(qs, family, pars, i, inv = TRUE)


        quantiles[[Y_column]] <- qY
      }
    }
  }

  ## now rescale to correct margin
  X <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)

  Y <- rescale_var(qY, X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  # Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}


compute_copula_quantiles <- function(qs, family, pars, i, inv) {
  library(copula)
  #change for now don't know how to add df

  copula_functions <- list(
    function(U, param, par2, inv) pnorm(qnorm(U[, 2]) * sqrt(1 - param^2) + param * qnorm(U[, 1])),
    function(U, param, par2, inv) cCopula(U, copula = tCopula(param = param, df = par2, dim = 2, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = claytonCopula(param = param, dim = 2, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = gumbelCopula(param = param, dim = 2, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = frankCopula(param = param, dim = 2, dispstr = "un"), inverse = inv),
    function(U, param, par2, inv) cCopula(U, copula = joeCopula(param = param, dim = 2, dispstr = "un"), inverse = inv)
  )
  fam <- family[[2]][[i]]
  beta <- pars[[2]][[i]]$beta
  beta <- 2 * expit(beta) - 1
  if (fam >= 1 && fam <= 6) {
    qY <- copula_functions[[fam]](qs, beta, par2 = 5, inv)
  } else {
    stop("family must be between 0 and 5")
  }
  return(qY[,2])
}
