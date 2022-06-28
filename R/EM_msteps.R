EM_msteps=function (Y,mu,Sigma,df,K,error){
  LIKE=function(Y,delta1,Sigma1,df1){
    L=prod(dmvt(Y,delta = delta1,sigma= Sigma1,df=df1))
    return(L)
  }
  n=nrow(Y)
  p=ncol(Y)
  K1=1

  muchain=matrix(rep(NA,p*K),nrow=p)

  SigmaChain=array(rep(NA,p*p*K),dim = c(p,p,K))
  YChain=array(rep(0,n*p*K),dim = c(n,p,K))
  Y1=matrix(rep(0,n*p),nrow=n)
  muchain[,K1]=mu
  YChain[,,K1]=Y1
  SigmaChain[,,K1]=Sigma
  ERR1=LIKE(Y1,mu,Sigma,df)

  Cond=TRUE
  while(Cond==TRUE)
  { ERR2=ERR1
  K1=K1+1
  GH= EM_onestep (Y,mu,Sigma,df)
  muchain[,K1]=GH$mu

  YChain[,,K1]=GH$Y2
  SigmaChain[,,K1]=GH$Sigma
  mu=GH$mu
  Sigma=GH$Sigma

  Y23=abs(GH$Y2-Y1)
  Y1=GH$Y2
  ERR1=LIKE(Y1,mu,Sigma,df)

  Cond=(abs(ERR1-ERR2)>error & K1<K)

  }

  IMP=Y1

  Y2=GH$Y2
  return(list(Y2=Y2,IMP=IMP,mu=mu,Sigma=Sigma,df=df,K1=K1,YChain=YChain[,,1:K1],muchain=muchain[,1:K1],SigmaChain=SigmaChain[,,1:K1]))
}
