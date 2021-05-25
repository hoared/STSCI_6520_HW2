#' Ordinary Least Squares
#'
#' @description Computes the Ordinary Least Squares solution using either the Gauss-Seidel or the Jacobi iterative methods. 
#'
#' @param A A matrix of coefficients
#' @param b A vector
#' @param method Either "Gauss" for the Gauss-Seidel method or "Jacobi" for the Jacobi Method. Implements the parallel-Jacobi method if numcores is more than 1 and the sequential Jacobi method if numcores is 1.
#' @param numcores An integer for the numebr of cores to use for computing the Jacobi method
#' @param niter The number of iterations to perform.
#'
#' @return A solution x to the equation Ax = b
#' @export
#'
#'
solve_ols <- function(A, b, method = "Gauss", numcores = 1, niter = 10000){
  # Some conditions
  if(numcores < 1){
    stop("Please input a positive value for numcores")
  }
  
  x <- rep(0, length(b))
  n <- length(b)
  
  if(method=="Gauss"){
    for(i in 1:niter){
      x[1] <- 1/A[1,1]*(b[1] + x[2])
      for(j in 2:(n-1)){
        # Update x[j]
        x[j] <- 1/A[j,j]*(b[j] + x[j-1] + x[j+1])
      }
      # Last step
      x[n] <- 1/A[n,n]*(b[n]+x[n-1])
    }
    return(x)
    
  } else if(method=="Jacobi" && numcores==1){
    # If the number of cores is 1 perform sequential jacobi
    for(i in 1:niter){
      temp = x
      #First Step
      x[1] <- 1/A[1,1]*(b[1] + temp[2])
      
      for(j in 2:(n-1)){
        x[j] <- 1/A[j,j]*(b[j] + temp[j-1]+temp[j+1])
      }
      
      #Last Step
      x[n] <- 1/A[n,n]*(b[n]+temp[n-1])
    }
    
    return(x)
    
  } else if(method=="Jacobi" && numcores >1){
    # If the number of cores is greater than 1, perform parallel jacobi
    for(i in 1:niter){
      #Do the computation
      
      x <- unlist(parallel::mclapply(1:n, function(j){
        if(j==1){
          tmp <- 1/A[1,1]*(b[1] + x[2])
        } else if(j==100){
          tmp <- 1/A[100,100]*(b[100]+x[99])
        } else{
          tmp <- 1/A[j,j]*(b[j] + x[j-1]+x[j+1])
        }
        tmp
      }, mc.cores = numcores))
      
    }
    
    return(x)
    
  } else{
    stop("Please give either Gauss or Jacobi for method")
  }
  
  
  
}