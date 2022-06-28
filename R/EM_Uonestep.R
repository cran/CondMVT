EM_Uonestep=function(Y,mu,Sigma,df,e){
  fun1=function(df,Weight,p){
    SS1=0
    SS2=0
    n=length(Weight)
    for (i in 1:n){
      SS1=SS1+log(Weight[i])-Weight[i]
      SS2=SS2+digamma((df+p[i])/2)-log((df+p[i])/2)
    }

    ty=-digamma(df/2)+log(df/2)+(1/n)*SS1+1+(1/n)*SS2
    return(ty)
  }


  dfun1=function(n,p,df){
    -0.5*trigamma(df/2)+(1/df)+(1/n)*sum(.5*trigamma((df+p)/2)-(1/(df+p)))
  }

  Bisec=function(a,b,fun1,e){
    soln=c()
    k=2
    d=(a+b)/2
    soln[1:2]=c(a,b)
    while (abs(a-b)>e) {
      if((fun1(a)*fun1(d))<0){
        a=a;b=d }else{a=d;b=b}
      d=(a+b)/2
      k=k+1
      soln=c(soln,d)
    }
    return(list(root=d,soln=soln,k=k))
  }
  n=nrow(Y)
  p=ncol(Y)
  pp=rep(NA,n)
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
      #pp=length(R[i,]==0)
      pp[i]=sum(R[i,]==1)
      #for (j in 1:p){
      Y1[i,-bb]=Y[i,-bb]
      K1=CondMVT(mean=mu, sigma=Sigma,df=df, dependent.ind =bb, given.ind =bn[-bb],
                 X.given=Y[i,-bb],check.sigma = F)

      Y1[i,bb]=K1$condMean

      #M=mahalanobis(Y1[i,-bb],mu[-bb],K1$phi)
      Weight[i]<-1/K1$R #(df+pp)/(df+M)

      MPHI=matrix(rep(0,p*p),p)
      MPHI[bb,bb]=K1$CIJK

    }
    else{Y1[i,]=Y[i,]
    pp[i]=p
    #print(pp[i])
    M=mahalanobis(Y1[i,],mu,Sigma)
    Weight[i]=(df+pp[i])/(df+M)
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

  #df1=df
  ww=Weight
  pp1=pp

  ff=function(df1){

    fun1(df1,ww,pp1)

  }

  df_1=Bisec(0.001,200,ff,.00001)


  df=df_1$root
  return(list(Y2=Y1,mu=mu2,Sigma=Sigma,df=df,PHI1=PHI1,P2=PHI_2))

}
