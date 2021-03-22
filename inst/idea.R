#install.packages("StratifiedSampling")
#install.packages("transport")
#install.packages("proxy")

#library(transport)
#library(sampling)
#library(StratifiedSampling)

rm(list=ls())

N=100
p=2
n1=30
n2=50
X=array(rnorm(N*p,5,1),c(N,p))

pik1=sampling::inclusionprobabilities(runif(N),n1)
pik2=sampling::inclusionprobabilities(runif(N),n2)
s1=sampling::UPpivotal(pik1)
s2=sampling::UPpivotal(pik2)

EPS=0.000001

X1=X[s1>1-EPS,]
X2=X[s2>1-EPS,]
dk1=1/pik1[s1>1-EPS]
dk2=1/pik2[s2>1-EPS]
dk1=dk1/sum(dk1)*N
dk2=dk2/sum(dk2)*N

w1=dk1
w2=dk2

statistical.linkage.random<-function(X1,X2,w1,w2,method = "Euclidean",EPS=1e-7)
{
  # distance
  D <- proxy::dist(X1,X2,method = method)
  
  # optimal transport
  res <- transport::transport(w1,w2,D)
  #XR=cbind(as.matrix(res[,1:3]) , X1[res[,1],] , X2[res[,2],] )
  
  # new p
  q <- rep(0,nrow(res))
  for(i in 1:n1){
    q[res[,1]==i] <- res[,3][res[,1]==i]/sum(res[,3][res[,1]==i])
  } 
  
  # 
  Z <- cbind( res[,3]*X2[res[,2],] )
  
  #intializing for sampling
  strata <- res[,1][q<1-EPS]
  XXX <- Z[q<1-EPS,]
  qqq <- q[q<1-EPS]
  
  # select stratified sample
  s <- StratifiedSampling::stratifiedcube(X=XXX,strata=strata,pik=qqq)
  
  # finalization
  ss <- rep(1,nrow(res))
  ss[q<1-EPS] <- s
  return(res[ss==1,1:3])
}


statistical.linkage.expected<-function(X1,X2,w1,w2,EPS=0.000001)
{
  D=proxy::dist(X1,X2,method = "Euclidean")
  res <- transport::transport(w1,w2,D)
  XR=cbind(X1[res[,1],] , X2[res[,2],] )
  q=rep(0,nrow(res))
  for(i in 1:n1) q[res[,1]==i]=res[,3][res[,1]==i]/sum(res[,3][res[,1]==i])
  Z=cbind( res[,3]*X2[res[,2],] )
  strata=res[,1][q<1-EPS]
  XXX=Z[q<1-EPS,]
  qqq=q[q<1-EPS]
  t(q*disjunctive(res[,1]))%*%XR
}



statistical.linkage.full<-function(X1,X2,w1,w2,EPS=0.000001)
{
  D=proxy::dist(X1,X2,method = "Euclidean")
  res <- transport::transport(w1,w2,D)
  cbind(as.matrix(res[,1:3]) , X1[res[,1],] , X2[res[,2],] )
}




T=statistical.linkage.random(X1,X2,dk1,dk2,EPS=0.000001)
MAT1=cbind(X1,X2[T[,2],])

MAT2=statistical.linkage.full(X1,X2,dk1,dk2,EPS=0.000001)

MAT3=statistical.linkage.expected(X1,X2,dk1,dk2,EPS=0.000001)


colSums(X)
colSums(dk1*X1)
colSums(dk2*X2)
colSums(dk1*MAT1)
colSums(MAT2[,3]*MAT2)
colSums(dk1*MAT3)









