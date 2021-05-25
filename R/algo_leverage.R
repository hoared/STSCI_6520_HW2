#' Weighted Leverating
#' 
#' @description A method for computing the coefficients in a linear regression using subsampling by selecting samples either uniformly or with probabilities weighted by the leverage scores.
#'
#' @param X A matrix of input values
#' @param Y A vector of outputs
#' @param r The subsample size
#' @param method Either 'uniform' or 'leverage' methods. The uniform method
#'
#' @return Returns the coefficients computed using subsampling by selecting samples either uniformly or with probabilities weighted by the leverage scores.
#' @export
#'
algo_leverage <- function(X, Y, r, method = "uniform"){
  X <- as.matrix(X)
  n <- nrow(X)
  
  if(method=="uniform"){
    idx <- sample(1:n, size = r, replace = FALSE)
    beta_unif <- stats::lm(Y[idx]~X[idx,]+0)$coefficients
    
    return(beta_unif)
  } else if(method=="leverage"){
    H = X %*% solve(t(X)%*%X) %*% t(X)
    scores = diag(H)/sum(diag(H)) # These are the leverage scores
    
    idx <- sample(1:n, size = r, replace = TRUE, prob = scores) # Sample
    probs <- scores[idx] # Probabilities
    beta_blev <- stats::lm(Y[idx]~X[idx,]+0, weights = probs)$coefficients
    
    return(beta_blev)
  } else{
    stop("Invalid Method. Please use 'uniform' or 'leverage'")
  }
  return(0)
}
