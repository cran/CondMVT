CondMVT=function (mean, sigma, df, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  #library(mvtnorm)
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,dependent.ind])))
  #This give the marginal mean and variance
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  phi=solve(D)
  CDinv <- C %*% solve(D)
  c_df=df+length(X.given)
  KK=t(X.given - mean[given.ind])%*% solve(D)%*%(X.given - mean[given.ind])
  R=as.numeric((df+KK)/c_df)
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  CIJK=(B - CDinv %*% t(C))
  cVar <- R*CIJK
  return(list(condMean = cMu,CIJK=CIJK ,condVar = cVar,cond_df=c_df,phi=solve(D),R=R))
}

