EM_onestep=function(Y,mu,Sigma,df){
  n=nrow(Y)
  p=ncol(Y)
  Weight=rep(NA,n)
  bn=seq(1:p)
  R=ifelse(is.na(Y),0,1)
  Y1=matrix(rep(NA,n*p),n,p)
  FGH=which(R==0,arr.ind = T)
  W=sort(unique(FGH[,1]))# incomplete observations
  PHI=matrix(rep(0,p*p),p)
  PHI_2=array(rep(NA,p*p*n),dim=c(p,p,n))
  for (i in 1:n){
    if (is.element(i, W))
    {
      bb=bn[R[i,]==0]
      Y1[i,-bb]=Y[i,-bb]
      K1=CondMVT(mean=mu, sigma=Sigma,df=df, dependent.ind =bb, given.ind =bn[-bb],
                 X.given=Y[i,-bb],check.sigma = F)

      Y1[i,bb]=K1$condMean

      #M=mahalanobis(Y1[i,-bb],mu[-bb],K1$phi)
      Weight[i]=1/K1$R #(df+pp)/(df+M)

      MPHI=matrix(rep(0,p*p),p)
      MPHI[bb,bb]=K1$CIJK

    }
    else{Y1[i,]=Y[i,]
    M=mahalanobis(Y1[i,],mu,Sigma)
    Weight[i]=(df+p)/(df+M)
    MPHI=matrix(rep(0,p*p),p)

    }
    PHI_2[,,i]=MPHI
  }

  mu1=rep(0,p)
  Sigma1=matrix(rep(0,p*p),p,p)
  PHI1=matrix(rep(0,p*p),p,p)
  for (i in 1:n){
    mu1=mu1+Weight[i]*Y1[i,]
    PHI1=PHI1+PHI_2[,,i]
  }
  mu2=mu1/sum(Weight)
  for (i in 1:n){
    mu1=mu1+Weight[i]*Y1[i,]
    DFG=((Y1[i,]-mu2)%*%t(Y1[i,]-mu2))
    Sigma1=Sigma1+Weight[i]*DFG

  }



  PHI_3=(1/(n))*PHI1
  Sigma=(1/(n))*Sigma1+PHI_3


  return(list(Y2=Y1,mu=mu2,Sigma=Sigma,df=df,PHI1=PHI1,P2=PHI_2))

}
