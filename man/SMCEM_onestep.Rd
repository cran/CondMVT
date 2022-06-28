\name{SMCEM_onestep}
\alias{SMCEM_onestep}
\title{Data Imputation Using SEM and MCEM (Single Iteration, Degrees of Freedom Known)
}
\description{
This sub-package contains the subroutines for iterative imputation of missing values as well as parameter estimation (for the location vector and the scatter matrix) in multivariate t distribution using Stochastic EM (SEM) and Monte Carlo EM (MCEM). In this case, the degrees of freedom for the distribution are known or fixed a priori. SEM is implemented when the analyst specifies a single draw in the E-step. In case we have multiple draws in the E-step, the algorithm changes to MCEM. In both algorithms, the function SMCEM_onestep is run when we are only interested in the imputed values and the parameter updates in a single iteration.
}
\usage{SMCEM_onestep(Y,mu,Sigma,df,nob)
}
\arguments{
 \item{Y}{the multivariate t dataset}
  \item{mu}{	the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
  \item{nob}{number of draws in the E-step}
}
\value{
 Completed dataset, updated location vector, and scatter matrix when employing the SEM and MCEM algorithms. All outputs are numeric. 
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

#Initializing Values
mu_stat=c(0.5,1,2)
Sigma_stat=matrix(c(0.33,0.31,0.3,0.31,0.335,0.295,0.3,0.295,0.32),3,3)

#Imputing Missing Values and Updating Parameter Estimates
#Single Iteration (SEM)
SEM1=SMCEM_onestep(Y=Y8,mu= mu_stat,Sigma=A,df=df,nob=1)

#Single Iteration (MCEM)
MCEM1=SMCEM_onestep(Y=Y8,mu= mu_stat,Sigma=A,df=df,nob=100)

#Results for Newly Completed Dataset (SEM)
SEM1$Y2    #Newly completed Dataset (with imputed values)
SEM1$mu	   #updated location vector
SEM1$Sigma #updated scatter matrix

#Results for Newly Completed Dataset (MCEM)
MCEM1$Y2    #Newly completed Dataset (with imputed values)
MCEM1$mu	   #updated location vector
MCEM1$Sigma #updated scatter matrix
}
