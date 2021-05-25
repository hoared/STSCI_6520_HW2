#' Vector 2-Norm
#'
#' @param x A vector 
#'
#' @return The norm of x
#' 
#'
#'
vec_norm <- function(x){
  # Convenience function for computing vectorn norm
  return(sqrt(sum(x^2)))
}

#' Elastic Net
#' 
#' @description Fit a linear model with an elastic net penalty for a sequence of values for the regularization coefficient lambda.
#'
#' @param X Matrix of inputs
#' @param Y Vector of outputs
#' @param lambdas A sequence of lambda values for the elastic net algorithm
#' @param alpha The elastic net mixing parameter. Alpha should be between 0 and 1
#'
#' @return A matrix of the fitted coefficients. The coefficients in column i correspond to lambda_i
#' @export
#'
#' 
elnet_coord <- function(X, Y, lambdas = seq(0.01, 3, length.out = 25), alpha = 0.5){
  if(alpha<0){
    stop("alpha should be between 0 and 1")
  }
  if(alpha>1){
    stop("alpha should be between 0 and 1")
  }
  
  p <- ncol(X)
  n <- nrow(X)
  
  # Standardize X
  x<-X
  for(i in 1:p){
    x[,i] <- x[,i]/vec_norm(x[,i])
  }
  
  nlambda = length(lambdas)
  beta_hat = matrix(0, nrow = p, ncol = nlambda) #Initialize beta
  
  
  for(lidx in 1:nlambda){
    lambda = lambdas[lidx] # Set Lambda
    
    # Elastic Net Algorithm
    niter = 100
    
    for(iter in 1:niter){
      
      for(j in 1:p){
        # Compute the cross term
        cross <- 0
        for(i in 1:n){
          for(k in 1:p){
            cross <- cross +  x[i,j]*x[i,k]*beta_hat[k, lidx]*(k!=j)
          }
        }
        cond = (2/n)*(x[,j] %*% Y)-2/n*cross
        if(cond >= lambda*alpha){
          beta_hat[j, lidx] <- (cond - lambda*alpha)/(2/n + 2*lambda*(1-alpha))
        } else if(cond <= -lambda*alpha){
          beta_hat[j, lidx] <- (cond + lambda*alpha)/(2/n + 2*lambda*(1-alpha))
        } else{
          beta_hat[j, lidx] <- 0
        }
      }
    }
  }
  
  return(list(betas = beta_hat))
}