\name{dcmvt}
\alias{dcmvt}
\title{Conditional Multivariate t Density and Random Deviates
}
\description{This function provides the density function for the conditional multivariate t distribution, [Y given X], where Z = (X,Y) is the fully-joint multivariate t distribution with location vector (or mode) equal to mean and covariance matrix sigma.}
\usage{
dcmvt(x, mean, sigma,df, dependent.ind, given.ind, X.given, check.sigma=TRUE, log = FALSE)
}
\arguments{
  \item{x}{vector or matrix of quantiles of Y. If x is a matrix, each row is taken to be a quantile.}
  \item{mean}{location vector, which must be specified.}
  \item{sigma}{a symmetric, positive-definte matrix of dimension n x n, which must be specified.}
  \item{df}{degrees of freedom, which must be specified.}
  \item{dependent.ind}{a vector of integers denoting the indices of dependent variable Y.}
  \item{given.ind}{a vector of integers denoting the indices of conditoning variable X.}
 \item{X.given}{a vector of reals denoting the conditioning value of X. When both given.ind and X.given are missing, the distribution of Y becomes Z[dependent.ind]}
 \item{check.sigma}{logical; if TRUE, the scatter matrix is checked for appropriateness (symmetry, positive-definiteness). This could be set to FALSE if the user knows it is appropriate.}
 \item{log}{logical; if TRUE, densities d are given as log(d).}
}
\value{
\code{numeric}
}
\references{Genz, A. and Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195. Springer-Verlag, Heidelberg.

S. Kotz and S. Nadarajah (2004), Multivariate t Distributions and Their Applications. Cambridge University Press. Cambridge.}
\examples{
# 10-dimensional multivariate t distribution
n <- 10
df=3
A <- matrix(rt(n^2,df), n, n)
A <- tcrossprod(A,A) #A %*% t(A)

# density of Z[c(2,5)] given Z[c(1,4,7,9)]=c(1,1,0,-1)
dcmvt(x=c(1.2,-1), mean=rep(1,n), sigma=A,dependent.ind=c(2,5),df=df, given.ind=c(1,4,7,9),
X.given=c(1,1,0,-1))

dcmvt(x=-1, mean=rep(1,n), sigma=A,df=df, dep=3, given=c(1,4,7,9,10), X=c(1,1,0,0,-1))
dcmvt(x=c(1.2,-1), mean=rep(1,n), sigma=A,df=df, dep=c(2,5))

# gives an error since `x' and `dep' are incompatibe
#dcmvt(x=-1, mean=rep(1,n), sigma=A,df=df, dep=c(2,3),
#      given=c(1,4,7,9,10), X=c(1,1,0,0,-1))

rcmvt(n=10, mean=rep(1,n), sigma=A,df=df, dep=c(2,5),
         given=c(1,4,7,9,10), X=c(1,1,0,0,-1),type="shifted",
         method="eigen")

rcmvt(n=10, mean=rep(1,n), sigma=A,df=df, dep=3,
         given=c(1,4,7,9,10), X=c(1,1,0,0,-1),type="Kshirsagar",
         method="chol")
}

