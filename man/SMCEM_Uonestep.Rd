\name{SMCEM_Uonestep}
\alias{SMCEM_Uonestep}
\title{Data Imputation Using SEM and MCEM (Single Iteration; Degrees of Freedom Unknown)
}
\description{
This sub-package provides subroutines for implementation of SEM and MCEM techniques in imputing missing values as well as estimating multivariate t parameters when the degrees of freedom are unknown. It has 4 functions namely fun1, dfun1, Bisec, and SMCEM_Uonestep. The functions fun1 and dfun1 in the sub-package constitute the equation for the degrees of freedom and its derivative respectively. The Bisec function contains the bisection method subroutines to facilitate the iterative estimation of the degrees of freedom using fun1 and dfun1. The function SMCEM_Uonestep constitute the SEM and MCEM algorithms for single-iteration data imputation and parameter estimation for multivariate t data with unknown degrees of freedom. The functions represent SEM when the number of draws in the E-step (denoted by nob) is 1 and MCEM when we have more than one draw in the E-step.Details of how SEM and MCEM impute missing values and estimate parameters in multivariate t context (unknown degrees of freedom) are explained by Kinyanjui et al. (2020).
}
\usage{
SMCEM_Uonestep(Y,mu,Sigma,df,nob,e)
}
\arguments{
 \item{Y}{the multivariate t dataset}
  \item{mu}{	the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
  \item{nob}{number of draws in the E-step}
  \item{e}{tolerance level for convergence of the bisection method for estimation of df. }
}
\value{
 Completed dataset, updated location vector,scatter matrix, and degrees of freedom when employing the SEM and MCEM algorithms. All outputs are numeric. 
 }
\references{
Kinyanjui, P. K., Tamba, C. L., Orawo, L. A. O., & Okenye, J. O. (2020). Missing data imputation in multivariate t distribution with unknown degrees of freedom using expectation maximization algorithm and its stochastic variants. Model Assisted Statistics and Applications, 15(3), 263-272.
}
\examples{
# 3-dimensional multivariate t distribution
n <- 25
p=3
df=3
mu=c(10,20,30)
A=matrix(c(14,10,12,10,13,9,12,9,18), 3,3)
Y7 <-mvtnorm::rmvt(n, delta=mu, sigma=A, df=df)
Y7
TT=Y7 #Complete Dataset

#Introduce MAR Data
Y8= MISS(TT,20) #The newly created incomplete dataset.
Y8

#Initializing Values
mu_stat=c(0.5,1,2)
Sigma_stat=matrix(c(0.33,0.31,0.3,0.31,0.335,0.295,0.3,0.295,0.32),3,3)
df_stat=6

#Imputing Missing Values and Updating Parameter Estimates

#Single Iteration (SEM)
SEMU1=SMCEM_Uonestep(Y=Y8,mu=mu,Sigma=Sigma_stat,df= df_stat,nob=1,e=0.00001)

#Single Iteration (MCEM)
MCEMU1=SMCEM_Uonestep(Y=Y8,mu=mu,Sigma=Sigma_stat,df= df_stat,nob=1000,e=0.00001)

#Results for Newly Completed Dataset (SEM)
SEMU1$Y2     #Newly completed Dataset (with imputed values)
SEMU1$mu	   #updated location vector
SEMU1$Sigma  #updated scatter matrix

#Results for Newly Completed Dataset (MCEM)
MCEMU1$Y2     #Newly completed Dataset (with imputed values)
MCEMU1$mu	    #updated location vector
MCEMU1$Sigma  #updated scatter matrix
MCEMU1$df     #updated degrees of freedom
}
