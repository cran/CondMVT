\name{CondMVT}
\alias{CondMVT}
\title{Conditional Location Vector, Scatter Matrix, and Degrees of Freedom of Multivariate t Distribution}
\usage{
CondMVT(mean, sigma, df, dependent.ind, given.ind, X.given,
  check.sigma = TRUE)
}
\arguments{
\item{mean}{location vector, which must be specified.}

\item{sigma}{a symmetric, positive-definte matrix of dimension n x n, which must be specified.}

\item{df}{degrees of freedom, which must be specified.}

\item{dependent.ind}{a vector of integers denoting the indices of dependent variable Y.}

\item{given.ind}{a vector of integers denoting the indices of conditoning variable X.}

\item{X.given}{a vector of reals denoting the conditioning value of X. When both given.ind and X.given are missing, the distribution of Y becomes Z[dependent.ind]}

\item{check.sigma}{logical; if TRUE, the scatter matrix is checked for appropriateness (symmetry, positive-definiteness). This could be set to FALSE if the user knows it is appropriate.}
}
\value{
Returns the conditional location vector (\code{condMean}), conditional scatter matrix (\code{condVar}), and the conditional degrees of freedom (\code{cond_df}) for the multvariate t distribution. 
}
\description{
These functions provide the conditional location vector, scatter matrix, and degrees of freedom of [Y given X], where Z = (X,Y)
is the fully-joint multivariate t distribution with location vector equal to mean, scatter matrix sigma, and degrees of freedom df. For more details on the computation of the parameters and their respective formulae, see Roth (2013).
}
\references{Roth, M. (2013). On the multivariate t-distribution, Tech Rep.
}
\examples{
 # 10-dimensional multivariate normal distribution
 n <- 10
 df=3
 A <- matrix(rt(n^2,df), n, n)
 A <- tcrossprod(A,A) #A %*% t(A)
 CondMVT(mean=rep(1,n), sigma=A, df=df, dependent=c(2,3,5), given=c(1,4,7,9),X.given=c(1,1,0,-1))
 CondMVT(mean=rep(1,n), sigma=A, df=df, dep=3, given=c(1,4,7,9), X=c(1,1,0,-1))

}
