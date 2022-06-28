MISS=function(TT, Percent){
  n=dim(TT)[1]
  p=dim(TT)[2]
  X=sample(1:n,size = round(Percent*n/100,0),
           replace = FALSE)
  TT1=TT #matrix(rep(NA,n*p),nrow=n)
  for (w in X){
    S=sample(1:(p-1),1,replace = T)
    S1=sample(1:p,size = S,replace = FALSE)
    TT1[w,]=replace(TT[w,],S1,NA)
  }
  return(TT1)
}
