\name{EM_Uonestep}
\alias{EM_Uonestep}
\title{Data Imputation Using EM (Single Iteration, Degrees of Freedom Unknown)
}
\description{
This sub-package constitutes the subroutines for EM algorithm (for a single iteration). It has 4 functions namely \code{fun1}, \code{dfun1}, \code{Bisec}, and \code{EM_Uonestep}. The function EM_Uonestep carries out missing data imputation as well as parameter estimation in multivariate t distribution in one iteration; assuming that the degrees of freedom are unknown. In addition to updating the location vector and the scatter matrix, therefore, the function also finds an estimate for the degrees of freedom. The bisection method is employed in the algorithm to iteratively update the degrees of freedom. In this respect, function \code{fun1} specifies the degrees of freedom equation to be solved. \code{dfun1} is its derivative. The two functions (\code{fun1} and \code{dfun1}) are then solved numerically using the bisection method as specified in the function \code{Bisec}.Details of how EM works in light of unknown degrees of freedom can be found in Kinyanjui et al. (2020) and Liu and Rubin (1995).
}
\usage{
EM_Uonestep(Y,mu,Sigma,df,e)
}
\arguments{
 \item{Y}{the multivariate t dataset}
  \item{mu}{	the location vector, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{Sigma}{Scatter matrix, which must be specified. In cases where it is unknown, starting values are provided.}
  \item{df}{degrees of freedom, which must be specified.}
   \item{e}{tolerance level for convergence of the bisection method for estimation of df. }
}
\value{
 Completed dataset, updated location vector,scatter matrix, and degrees of freedom. All outputs are numeric. 
 }
\references{
Kinyanjui, P. K., Tamba, C. L., Orawo, L. A. O., & Okenye, J. O. (2020). Missing data imputation in multivariate t distribution with unknown degrees of freedom using expectation maximization algorithm and its stochastic variants. Model Assisted Statistics and Applications, 15(3), 263-272.

Liu, C., & Rubin, D. B. (1995). ML estimation of the t distribution using EM and its extensions, ECM and ECME. Statistica Sinica, 19-39.
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

#Results for Newly Completed Dataset (EM)
EMU1$Y2    #Newly completed Dataset (with imputed values)
EMU1$mu	   #updated location vector
EMU1$Sigma #updated scatter matrix
EMU1$df
}
