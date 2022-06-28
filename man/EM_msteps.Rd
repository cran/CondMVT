\name{EM_msteps}
\alias{EM_msteps}
\title{Data Imputation Using EM (Multiple Iterations; Degrees of Freedom Known)
}
\description{
The sub-package contains subroutines for imputation of missing values as well as parameter estimation (for the location vector and the scatter matrix) in multivariate t distribution using the Expectation Maximization (EM) algorithm when the degrees of freedom are known. EM algorithm iteratively imputes the missing values and computes the estimates for the multivariate t parameters in two steps (E-step and M-step) as explained in Kinyanjui et al. (2021). For a single iteration, the function \code{EM_onestep} is run. For multiple iterations, the function \code{EM_msteps} is run. The function \code{LIKE} (specifying the likelihood) facilitates the setting of tolerance level for convergence of the algorithm (that is \eqn{L(\theta^{t+1})-L(\theta^{t})\leq{\delta}}, where \eqn{\delta} is a set tolerance level and t denotes the number of iterations).
}
\usage{
EM_msteps(Y,mu,Sigma,df,K,error)
}
\arguments{
    \item{Y}{the multivariate t dataset}
  \item{mu}{the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
  \item{K}{the number of iterations, which must be specified.}
   \item{error}{tolerance level for convergence of the EM algorithm.}
}
 \value{
 Completed dataset (with imputed values), updated location vector, and scatter matrix. All outputs are numeric. 
 }
\references{
Kinyanjui, P.K., Tamba, C.L., & Okenye, J.O. (2021). Missing Data Imputation in a t -Distribution with Known Degrees of Freedom Via Expectation Maximization Algorithm and Its Stochastic Variants. International Journal of Applied Mathematics and Statistics.
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

#Initializing Values
mu_stat=c(0.5,1,2)
Sigma_stat=matrix(c(0.33,0.31,0.3,0.31,0.335,0.295,0.3,0.295,0.32),3,3)

#Imputing Missing Values and Updating Parameter Estimates
#Single Iteration (EM)
EM1=EM_onestep(Y=Y8,mu=mu_stat,Sigma=Sigma_stat,df=df)

#Multiple Iterations (EM)
EM=EM_msteps(Y=Y8,mu=mu_stat,Sigma=Sigma_stat,df=3,K=1000,error=0.00001)

#Results for Newly Completed Dataset (EM)
EM$IMP   #Newly completed Dataset (with imputed values)
EM$mu	   #updated location vector
EM$Sigma #updated scatter matrix
EM$K1	#number of iterations the algorithm takes to converge

}
