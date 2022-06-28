\name{SMCEM_msteps}
\alias{SMCEM_msteps}
\title{Data Imputation Using SEM and MCEM (Multiple Iterations, Degrees of Freedom Known)

}
\description{
This sub-package contains the subroutines for iterative imputation of missing values as well as parameter estimation (for the location vector and the scatter matrix) in multivariate t distribution using Stochastic EM (SEM) and Monte Carlo EM (MCEM). In this case, the degrees of freedom for the distribution are known or fixed a priori. SEM is implemented when the analyst specifies a single draw in the E-step. In case we have multiple draws in the E-step, the algorithm changes to MCEM. In both algorithms, the function SMCEM_onestep is run when we are only interested in the imputed values and the parameter updates in a single iteration. The function SMCEM_msteps is run when we are interested in multiple iterations (this is usually the case). Essentially, the first iterations (for instance, 10 percent of all iterations) is usually burnt-in in order to ward off the effects of initial values. Details of how SEM and MCEM operate can be found in among others Kinyanjui et al. (2021), Nielsen (2000), Levine and Casella (2001) Jank (2005) and Karimi et al. (2019).
}
\usage{SMCEM_msteps(Y,mu,Sigma,df, nob,K)
}
\arguments{
  \item{Y}{the multivariate t dataset}
  \item{mu}{	the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
  \item{nob}{number of draws in the E-step}
  \item{K}{the number of iterations, which must be specified.}
}
\value{
 Completed dataset, updated location vector, and scatter matrix when employing the SEM and MCEM algorithms. All outputs are numeric. 
 }
\references{
Karimi, B., Lavielle, M., and Moulines, É. (2019). On the Convergence Properties of the Mini-Batch EM and MCEM Algorithms.

Kinyanjui, P.K., Tamba, C.L., & Okenye, J.O. (2021). Missing Data Imputation in a t -Distribution with Known Degrees of Freedom Via Expectation Maximization Algorithm and Its Stochastic Variants. International Journal of Applied Mathematics and Statistics.

Levine, R. A. and Casella, G. (2001). Implementations of the Monte Carlo EM algorithm. Journal of Computational and Graphical Statistics, 10(3), 422-439.

Nielsen, S.F. (2000). The stochastic EM algorithm: estimation and asymptotic results. Bernoulli, 6(3), 457-489.
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
#Single Iteration (SEM)
SEM1=SMCEM_onestep(Y=Y8,mu= mu_stat,Sigma=Sigma_stat,df=df,nob=1)

#Single Iteration (MCEM)
MCEM1=SMCEM_onestep(Y=Y8,mu= mu_stat,Sigma=Sigma_stat,df=df,nob=100)

#Multiple Iterations (SEM)
SEM=SMCEM_msteps(Y=Y8,mu= mu_stat,Sigma= Sigma_stat,df=df,nob=1,K=500)

#Results for Newly Completed Dataset (Burning in first 50 iterations in SEM)
T_mu=rep(0,3)
T_Sigma=matrix(rep(0,3*3),nrow=3)
T_Data=matrix(rep(0,3*10), nrow =10)
for (l in 51:500){
  T_mu = T_mu + SEM$muchain[l,]
  T_Sigma = T_Sigma + SEM$SigmaChain[,,l]
  T_Data= T_Data+ SEM$YChain[,,l]
}
#updated location vector
round((T_mu/450),4) 
#updated scatter matrix  
round((T_Sigma/450),4)
#complete dataset as an average of (K-50) complete datasets for the various iterations.  
T_Data1=  T_Data/450	
T_Data1

#Multiple Iterations (MCEM)
MCEM=SMCEM_msteps(Y=Y8,mu=mu_stat,Sigma=Sigma_stat,df=df,nob=100,
K=500)

#Results for Newly Completed Dataset (Burning in first 50 iterations in MCEM)
T_mu=rep(0,3)
T_Sigma=matrix(rep(0,3*3),nrow=3)
T_Data=matrix(rep(0,3*10), nrow =10)
for (l in 51:500){
  T_mu = T_mu + MCEM$muchain[l,]
  T_Sigma = T_Sigma + MCEM$SigmaChain[,,l]
  T_Data= T_Data+ MCEM$YChain[,,l]
}
#updated location vector
round((T_mu/450),4) 
#updated scatter matrix
round((T_Sigma/450),4) 
#complete dataset as an average of (K-50) complete datasets for the various iterations.
T_Data1= T_Data/450 
T_Data1                    

}
