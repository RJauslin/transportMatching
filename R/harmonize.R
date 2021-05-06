
#' Harmonisation by calibration
#'
#' @param X1 dataset 1 
#' @param d1 weights 1
#' @param id1 identifiers 1 
#' @param X2 dataset 2
#' @param d2 weigths 2
#' @param id2 identifiers 2
#' @param totals 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' N <- 10000
#' X <- data.frame(x1 = rnorm(N,0,1), x2 = rnorm(N,0,1))
#' 
#' n1=1000
#' n2=3000
#' s1=srswor(n1,N)
#' s2=srswor(n2,N)
#' id1=(1:N)[s1==1]
#' id2=(1:N)[s2==1]
#' d1=rep(N/n1,n1)
#' d2=rep(N/n2,n2)
#' X1 = X[s1==1,]
#' X2 = X[s2==1,]
#' 
#' W <- harmonize(X1,d1,id1,X2,d2,id2)
#' 
#' colSums(W$w1*X1)
#' colSums(W$w2*X2)
#' 
#' 
#' ############## if knowing the true totals
#' 
#' totals <- c(N,colSums(X))
#' W <- harmonize(X1,d1,id1,X2,d2,id2,totals)
#' 
#' colSums(W$w1*X1)
#' colSums(W$w2*X2)
#' colSums(X)
#' 
harmonize <- function(X1,d1,id1,X2,d2,id2,totals)
{
  
  # number of units in each sample
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  # add constant vector to ensure same sum
  XX1=cbind(rep(1,n1),X1)
  XX2=cbind(rep(1,n2),X2)
  
  # we can specify the desired total (for example if we know the totals of the population)
  if(missing(totals)){
    n12=length(intersect(id1,id2))
    a=(n1-n12)/(n1+n2-2*n12)  
    totals=a*colSums(d1*XX1)+(1-a)*colSums(d2*XX2)
  }

  # calibration with sampling package
  w1=d1*calibR(as.matrix(XX1),d1,totals,q = rep(1,length(d1)))
  w2=d2*calibR(as.matrix(XX2),d2,totals,q = rep(1,length(d2)))
  
  # w1=d1*calibRaking(as.matrix(XX1),d1,totals,q = rep(1,length(d1)))
  # w2=d2*calibRaking(as.matrix(XX2),d2,totals,q = rep(1,length(d2)))

  return(list(w1=w1,w2=w2))
}

