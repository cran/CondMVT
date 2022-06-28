\name{EM_onestep}
\alias{EM_onestep}
\title{Data Imputation Using EM (Single Iteration; Degrees of Freedom Known)
}
\description{
The sub-package contains subroutines for imputation of missing values as well as parameter estimation (for the location vector and the scatter matrix) in multivariate t distribution using the Expectation Maximization (EM) algorithm when the degrees of freedom are known. EM algorithm imputes the missing values and computes the estimates for the multivariate t parameters in two steps (E-step and M-step) as explained in Kinyanjui et al. (2021). For a single iteration, the function EM_onestep is run.Arbitrary starting values are supplied to initiate the algorithm.
}
\usage{
EM_onestep(Y,mu,Sigma,df)
}
\arguments{
   \item{Y}{the multivariate t dataset}
  \item{mu}{the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
 algorithm.
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

#Results for Newly Completed Dataset (EM)
EM1$Y2    #Newly completed Dataset (with imputed values)
EM1$mu	  #updated location vector
EM1$Sigma #updated scatter matrix

}
