
##################################################################################################
########## Simulating multivarite singular Gaussian random numbers using Eigen value-decompostion  ##########
##################################################################################################
#### Note that the cholesky deocompostion is faster than the eigen decompostion, so when the covariance matric or
## precision matrix is non-singular, then always use choleskey-factorization 

#' Truncated Eigen Square Root (or inverse square root) 
#'
#' Designed for Covaraice or Precision matricies - throws error if 
#' there is any large negative Eigenvalues. 
#'
#' @param Sigma p x p matrix to compute truncated matrix square root
#' @param inv (boolean) whether to compute inverse of square root
#' @return p x k matrix (where k is rank of Sigma) 
trunc_eigen_sqrt <- function(Sigma, inv){
  es <- eigen(Sigma)
  
  # Small negative eigenvalues can occur just due to numerical error and should
  # be set back to zero. 
  es$values[(es$values < 1e-12) & (es$values > -1e-12)] <- 0
  
  # If any eigen values are large negative throw error (something wrong)
  if (any(es$values < -1e-12)) stop("Non-trivial negative eigenvalues present")
  
  # calculate square root and reveal rank (k)
  k <- sum(es$values > 0)
  if (!inv){
    L <- es$vectors %*% diag(sqrt(es$values))  
  } else if (inv) {
    L <- es$vectors %*% diag(1/sqrt(es$values))  
  }
  return(L[,1:k, drop=F])
}

#' Covariance parameterization (Eigen decomposition)
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma covariance matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormC_eigen <- function(n, mu, Sigma){
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv=FALSE)
  k <- ncol(L)
  Z <- matrix(rnorm(k*n), k, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}

#' Precision parameterization (Eigen decomposition)
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma precision matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormP_eigen <- function(n, mu, Sigma){
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv=TRUE) # only difference vs. cov parameterization, this is done in terms of covariance matrix because eigen values of a matrix and its inverse is the same
  k <- ncol(L)
  Z <- matrix(rnorm(k*n), k, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}
# n <- 10000
# mu <- 1:3
# Omega <- MASS::ginv(Sigma) # use pseudoinverse 
# x1 <- rMVNormC_eigen(n, mu, Sigma)
# x2 <- rMVNormP_eigen(n, mu, Omega)
# 
# # Create function that tests for equality with high tolerance due to 
# # random number generation
# weak_equal <- function(x, y) all.equal(x, y, tolerance=.05)
# 
# # check row means match up with mu and agree with eachother
# weak_equal(rowMeans(x1), mu) & weak_equal(rowMeans(x2), mu)
# 
# # check empirical covariances agree
# weak_equal(var(t(x1)), var(t(x2)))

