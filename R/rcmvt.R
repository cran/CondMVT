rcmvt=function (n, mean, sigma, df,dependent.ind, given.ind, X.given,
                check.sigma = TRUE,type = c("Kshirsagar", "shifted") ,method = c("eigen", "svd", "chol"))
{
  #library(mvtnorm)
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(mvtnorm::rmvt(n, delta= mean[dependent.ind], sigma = as.matrix(sigma[dependent.ind,
                                                                                dependent.ind]),df=df))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  ret <- CondMVT(X.given = X.given, mean = mean, sigma = sigma, df=df,
                 dependent.ind = dependent.ind, given.ind = given.ind,
                 check.sigma = check.sigma)
  
  mvtnorm::rmvt(n, delta = ret$condMean, sigma = ret$condVar, df=ret$cond_df,type = type)
}
