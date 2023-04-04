MDSIC <- function(distances, X_bmds_out) {
  n=nrow(distances)
  m=as.integer(n*(n-1)/2)
  maxp=length(X_bmds_out)

  ssr=LR=mdsic=llike=penalty= c(1:maxp)*0
  s=r=matrix(0,maxp,maxp);

  for(p in 1:maxp)  {
    xst=as.matrix(X_bmds_out[[p]])  # xst=as.matrix(xst.sv[, 1:p,p])
    ssr[p] = sum((distances-distRcpp(xst))^2)/2
    for(j in 1:p)  s[j,p]= t(xst[,j])%*%xst[,j]
  }

  mdsic[1]= (m-2)* log(ssr[1])

  if (maxp > 1) {

    for(p in 1:(maxp-1)){
      for(j in 1:p) r[j,p+1]= s[j,p+1]/s[j,p]
    }

    for(p in 1:(maxp-1)){
      sr=0;
      for(j in 1:p)
        sr=sr+log(r[j,p+1]*(n+1)/( n+r[j,p+1]))
      llike[p] = (m-2)*(log( ssr[p+1])-log( ssr[p]))
      penalty[p]=(n+1)*sr  + (n+1)*log(n+1)
      LR[p] = llike[p]+penalty[p]
    }

    for(p in 2 :maxp) mdsic[p] = mdsic[1] + sum( LR[1:(p-1)])
  }

  mdsic
}
