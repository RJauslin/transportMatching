context("Test for linkage function")


test_that("linakge --- All identificators appear.",{
 
  rm(list = ls())
  N=1000
  p=5
  X=array(rnorm(N*p),c(N,p))
  EPS= 1e-9

  n1=100
  n2=100
  
  s1=srswor(n1,N)
  s2=srswor(n2,N)

  id1=(1:N)[s1==1]
  id2=(1:N)[s2==1]

  d1=rep(N/n1,n1)
  d2=rep(N/n2,n2)

  X1=X[s1==1,]
  X2=X[s2==1,]

  re=harmonize(X1,d1,id1,X2,d2,id2)
  w1=re$w1
  w2=re$w2
  object = linkage(X1,id1,X2,id2,w1,w2)
 
  expect_equal( nlevels(factor(object$conc[,1])),nlevels(factor(id1)))
  expect_equal( nlevels(factor(object$conc[,2])),nlevels(factor(id2)))
})


test_that("linkage --- Respect totals.",{
  
  rm(list = ls())
  N=1000
  p=5
  X=array(rnorm(N*p),c(N,p))
  EPS= 1e-9
  
  n1=100
  n2=100
  
  s1=srswor(n1,N)
  s2=srswor(n2,N)
  
  id1=(1:N)[s1==1]
  id2=(1:N)[s2==1]
  
  d1=rep(N/n1,n1)
  d2=rep(N/n2,n2)
  
  X1=X[s1==1,]
  X2=X[s2==1,]
  
  re=harmonize(X1,d1,id1,X2,d2,id2)
  w1=re$w1
  w2=re$w2
  object = linkage(X1,id1,X2,id2,w1,w2)
  
  test1 <- as.numeric(round(colSums(object$conc$weight*object$conc)[4:ncol(object$conc)],5))
  test2 <- as.numeric(round(colSums(cbind(w1*X1,w2*X2)),5))
  
  expect_equal(any(abs(test1 - test2) > 1e-4),FALSE) 
  
})


test_that("linkage --- No intersect.",{
  
  rm(list = ls())
  N=1000
  p=5
  X=array(rnorm(N*p),c(N,p))
  EPS= 1e-9
  
  n1=100
  n2=100
  
  s1=srswor(n1,N)
  s2=srswor(n2,N)
  
  i <- intersect(which(s1 == 1),which(s2 == 1))
  
  s1[i] <- 0
  s2[i] <- 0
  
  n1 <- sum(s1)
  n2 <- sum(s2)
  
  
  id1=(1:N)[s1==1]
  id2=(1:N)[s2==1]
  
  d1=rep(N/n1,n1)
  d2=rep(N/n2,n2)
  
  X1=X[s1==1,]
  X2=X[s2==1,]
  
  re=harmonize(X1,d1,id1,X2,d2,id2)
  w1=re$w1
  w2=re$w2
  object = linkage(X1,id1,X2,id2,w1,w2)
  
  test1 <- as.numeric(round(colSums(object$conc$weight*object$conc)[4:ncol(object$conc)],5))
  test2 <- as.numeric(round(colSums(cbind(w1*X1,w2*X2)),5))
  
  expect_equal(any(abs(test1 - test2) > 1e-4),FALSE) 
  
})



test_that("linkage --- Character for identificators.",{
  
  rm(list = ls())
  N=1000
  p=5
  X=array(rnorm(N*p),c(N,p))
  EPS= 1e-9
  
  n1=100
  n2=100
  
  s1=srswor(n1,N)
  s2=srswor(n2,N)
  
  # i <- intersect(which(s1 == 1),which(s2 == 1))
  # 
  # s1[i] <- 0
  # s2[i] <- 0
  # 
  # n1 <- sum(s1)
  # n2 <- sum(s2)
  
  
  id1 <- paste0("A",(1:N)[s1==1])
  id2 <- paste0("B",(1:N)[s2==1])
  
  d1=rep(N/n1,n1)
  d2=rep(N/n2,n2)
  
  X1=X[s1==1,]
  X2=X[s2==1,]
  
  re=harmonize(X1,d1,id1,X2,d2,id2)
  w1=re$w1
  w2=re$w2
  object = linkage(X1,id1,X2,id2,w1,w2)
  
  
  expect_equal( nlevels(factor(object$conc[,1])),nlevels(factor(id1)))
  expect_equal( nlevels(factor(object$conc[,2])),nlevels(factor(id2)))
  
  test1 <- as.numeric(round(colSums(as.numeric(object$conc$weight)*apply(object$conc[,4:ncol(object$conc)],MARGIN = 2,as.numeric)),5))
  test2 <- as.numeric(round(colSums(cbind(w1*X1,w2*X2)),5))
  
  expect_equal(any(abs(test1 - test2) > 1e-4),FALSE) 
  
})





