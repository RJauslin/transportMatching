#' Title
#'
#' @param object output from the function \code{linkage}
#'
#' @return
#' @export
#'
#' @examples
#' rm(list = ls())
#' N=10000
#' p=5
#' X=array(rnorm(N*p),c(N,p))
#' EPS= 1e-9
#' 
#' n1=1000
#' n2=1000
#' s1=srswor(n1,N)
#' s2=srswor(n2,N)
#' 
#' 
#' id1=(1:N)[s1==1]
#' id2=(1:N)[s2==1]
#' 
#' d1=rep(N/n1,n1)
#' d2=rep(N/n2,n2)
#' 
#' X1=X[s1==1,]
#' X2=X[s2==1,]
#' 
#' re=harmonize(X1,d1,id1,X2,d2,id2)
#' w1=re$w1
#' w2=re$w2
#' object = linkage(X1,id1,X2,id2,w1,w2)
#' out <- statMatch_expected(object,Z)
#' 
#' round(colSums(object$conc$weight*object$conc[,4:ncol(object$conc)]),3)
#' round(colSums(cbind(w1*X1,w2*X2)),3)
#' round(colSums(w1*out$conc[,2:ncol(out$conc)]),3)
#' 
statMatch_expected <- function(object,Z){
  
  q <- rep(0,nrow(object$conc))
  n1 <- nlevels(factor(object$conc$id1))
  l1 <- levels(factor(object$conc$id1))
  
  q_l <- split(object$conc$weight,f =object$conc$id1)
  q <- as.numeric(do.call(c,lapply(q_l, function(x){x/sum(x)})))
  
  
  out <- object
  out$conc <- cbind(id1 = levels(factor(object$conc$id1)), t(q*disjunctive(object$conc$id1))%*%as.matrix((object$conc[,4:ncol(object$conc)])))
  out$conc <- apply(out$conc,2,as.numeric)
  out$weight <- do.call(c,lapply(q_l,sum))
  out$Z <- t(q*disjunctive(object$conc$id1))%*%disjunctive(Z[as.character(object$conc$id2),])
  return(out)
  
}