\name{EM_Umsteps}
\alias{EM_Umsteps}
\title{Data Imputation Using EM (Multiple Iterations, Degrees of Freedom Unknown)
}
\description{
This sub-package constitutes the subroutines for EM algorithm (for multiple iterations). It has 2 functions namely \code{LIKE} and \code{EM_Umsteps}.  The function \code{EM_Umsteps} carries out missing data imputation as well as parameter estimation in multivariate t distribution in multiple iterations; assuming that the degrees of freedom are unknown. In addition to updating the location vector and the scatter matrix, therefore, the function also finds an estimate for the degrees of freedom. The bisection method is employed in the algorithm to iteratively update the degrees of freedom.The function \code{LIKE} (specifying the likelihood) facilitates the setting of tolerance level for convergence of the EM algorithm (that is \eqn{L(\theta^{t+1})-L(\theta^{t})\leq{\delta}}, where \eqn{\delta} is a set tolerance level and t denotes the number of iterations).Details of how EM works in light of unknown degrees of freedom can be found in Kinyanjui et al. (2020) and Liu and Rubin (1995).
}
\usage{EM_Umsteps(Y,mu,Sigma,df,K,e,error)
}
\arguments{
    \item{Y}{the multivariate t dataset}
  \item{mu}{	the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{Scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
   \item{e}{tolerance level for convergence of the bisection method for estimation of df. }
  \item{error}{tolerance level for convergence of the EM algorithm. }
  \item{K}{the number of iterations, which must be specified.}

}
\value{
 Completed dataset, updated location vector,scatter matrix, and degrees of freedom. All outputs are numeric. 
 }
\references{
Kinyanjui, P. K., Tamba, C. L., Orawo, L. A. O., & Okenye, J. O. (2020). Missing data imputation in multivariate t distribution with unknown degrees of freedom using expectation maximization algorithm and its stochastic variants. Model Assisted Statistics and Applications, 15(3), 263-272.

Liu, C. and Rubin, D. B. (1995). ML estimation of the t distribution using EM and its extensions, ECM and ECME. Statistica Sinica, 19-39.
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
#Single Iteration (EM)
EMU1=EM_Uonestep (Y=Y8,mu=mu,Sigma= Sigma_stat,df= df_stat,e=0.00001)

#Multiple Iterations (EM)
EMU=EM_Umsteps(Y=Y8,mu=mu_stat,Sigma=Sigma_stat,df=df_stat,K=1000,e=0.00001,error=0.00001)

#Results for Newly Completed Dataset (EM)
EMU$IMP #Newly completed Dataset (with imputed values)
EMU$mu   #updated location vector
EMU$Sigma #updated scatter matrix
EMU$df	    #Updated degrees of freedom.
EMU$K1	#number of iterations the algorithm takes to converge

}
