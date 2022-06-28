\name{pcmvt}
\alias{pcmvt}
\title{Conditional Multivariate t Distribution
}
\description{Computes the distribution function of the conditional multivariate t, [Y given X], where Z = (X,Y) is the fully-joint multivariate t distribution with mean equal to location vector, df equal to degrees of freedom and scatter matrix sigma. Computations are based on algorithms by Genz and Bretz.
}
\usage{
pcmvt(lower = -Inf, upper = Inf, mean, sigma, df, dependent.ind, given.ind, X.given,
check.sigma = TRUE, algorithm = GenzBretz(), ...)
}
\arguments{
  \item{lower}{
the vector of lower limits of length n.
}
  \item{upper}{
the vector of upper limits of length n.
}
  \item{mean}{
the mean vector of length n.
}
  \item{sigma}{
a symmetric, positive-definte matrix, of dimension n x n, which must be specified.
}
  \item{df}{
degrees of freedom, which must be specified.
}
  \item{dependent.ind}{
a vector of integers denoting the indices of the dependent variable Y.
}
  \item{given.ind}{
a vector of integers denoting the indices of the conditioning variable X.

}
  \item{X.given}{
a vector of reals denoting the conditioning value of X. When both given.ind and X.given are missing, the distribution of Y becomes Z[dependent.ind]
}
  \item{check.sigma}{
logical; if TRUE, the variance-covariance matrix is checked for appropriateness (symmetry, positive-definiteness). This could be set to FALSE if the user knows it is appropriate.
}
  \item{algorithm}{
an object of class GenzBretz, Miwa or TVPACK specifying both the algorithm to be used as well as the associated hyper parameters.
}
  \item{\dots}{
additional parameters (currently given to GenzBretz for backward compatibility issues).
}
}
\details{
This program involves the computation of multivariate t probabilities with arbitrary correlation matrices.
}
\value{
The evaluated distribution function is returned with attributes

 \item{error}{estimated absolute error and}
  \item{msg }{Normal Completion}
}

\references{
Genz, A. and Bretz, F. (1999), Numerical computation of multivariate t-probabilities with application to power calculation of multiple contrasts. Journal of Statistical Computation and Simulation, 63, 361–378.

Genz, A. and Bretz, F. (2002), Methods for the computation of multivariate t-probabilities. Journal of Computational and Graphical Statistics, 11, 950–971.

Genz, A. (2004), Numerical computation of rectangular bivariate and trivariate normal and t-probabilities, Statistics and Computing, 14, 251–260.

Genz, A. and Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195. Springer-Verlag, Heidelberg.
}
\seealso{\code{\link{dcmvt}()},\code{\link{rcmvt}()},\code{\link{pmvt}()},\code{\link{GenzBretz}()}}

\examples{
n <- 10
df=3
A <- matrix(rt(n^2,df), n, n)
A <- tcrossprod(A,A) #A %*% t(A)

pcmvt(lower=-Inf, upper=1, mean=rep(1,n), sigma=A, df=df,	dependent.ind=3, 
      given.ind=c(1,4,7,9,10), X.given=c(1,1,0,0,-1))


pcmvt(lower=-Inf, upper=c(1,2), mean=rep(1,n),
      sigma=A,df=df, dep=c(2,5), given=c(1,4,7,9,10),
      X=c(1,1,0,0,-1))

pcmvt(lower=-Inf, upper=c(1,2), mean=rep(1,n), sigma=A,df=df,
      dep=c(2,5))
}
