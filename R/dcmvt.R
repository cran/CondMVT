dcmvt=function (x, mean, sigma,df, dependent.ind, given.ind, X.given,
                check.sigma = TRUE, log = FALSE)
{
  #library(mvtnorm)
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (length(x) != length(dependent.ind))
    stop("lengths of `x' and `dependent.ind' must be same")
  if (missing(given.ind) & missing(X.given))
    return(mvtnorm::dmvt(x, delta = mean[dependent.ind], sigma = as.matrix(sigma[dependent.ind,dependent.ind]),df=df, log = log))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  ret <- CondMVT(X.given = X.given, mean = mean, sigma = sigma,df=df,
                 dependent.ind = dependent.ind, given.ind = given.ind,
                 check.sigma = check.sigma)
  mvtnorm::dmvt(x, delta = ret$condMean, sigma = ret$condVar,df=ret$cond_df, log = log)
}

