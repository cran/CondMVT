\name{MISS}
\alias{MISS}
\title{Creating Missing Values at Random in Multivariate Datasets

}
\description{
This function randomly creates missing values in a multivariate dataset. The resultant missing data mechanism is missing at random (MAR). The percentage of missingness has to be specified. This percentage is computed as a proportion of the sample size. In addition, the function allows for more than one missing value in any given case. It is set such that in a \eqn{p}-variate dataset, for any \eqn{i^{th}} case, the maximum allowable number of missing values is \eqn{p-1}. This helps avoid a situation where a case has no observed value.
}
\usage{
MISS (TT, Percent)
}
\arguments{
  \item{TT}{n×p complete dataset.}
  \item{Percent}{the proportion of missing values, which must be specified.}
  }
\value{Data \code{Y} of size \code{n×p}  with missing values (\code{NA}) created at random. The missing values are logical in nature. 
}
\examples{
# 3-dimensional multivariate t distribution
n <- 10
p=3
df=3
mu=c(1:3)
A <- matrix(rt(p^2,df), p, p)
A <- tcrossprod(A,A) #A %*% t(A)

Y7 <-mvtnorm::rmvt(n, delta=mu, sigma=A, df=df)
Y7
TT=Y7 #Complete Dataset

#Introduce MAR Data
Y8= MISS(TT,20) #The newly created incomplete dataset.
Y8  
}
