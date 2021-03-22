#' Title
#'
#' @param X1 
#' @param id1 
#' @param X2 
#' @param id2 
#' @param w1 
#' @param w2 
#' @param EPS 
#' @param CST 
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
#' nlevels(factor(object$conc[,1]))
#' nlevels(factor(id1))
#' nlevels(factor(object$conc[,2]))
#' nlevels(factor(id2))
#' 
#' round(colSums(object$conc$weight*object$conc[,4:ncol(object$conc)]),3)
#' round(colSums(cbind(w1*X1,w2*X2)),3)
#' 
#' 
linkage <- function(X1,id1,X2,id2,w1,w2,method = "Euclidean",EPS= 1e-9){
  
  
  # to increase precision failed from calib function in harmonize : TODO check that
  # w1 <- (w1/sum(w1))*N
  # w2 <- (w2/sum(w2))*N
  
  # distance
  D <- proxy::dist(X1,X2,method = "Euclidean")
  
  
  # to find units that are in intersection between the two sample
  inter <- base::intersect(id1,id2)
  inter1 <- do.call(c,lapply(inter,function(x){which(x == id1)}))
  inter2 <- do.call(c,lapply(inter,function(x){which(x == id2)}))
  
  # select the minimum weights between the two sample
  wtmp <- pmin(w1[inter1],w2[inter2])
  
  # change the weight of the intersected value to be the weight - min(w), some weights are put equal to 0.
  ww1 <- w1
  ww2 <- w2 
  ww1[inter1] <- ww1[inter1] - wtmp
  ww2[inter2] <- ww2[inter2] - wtmp
  
  # weights that are put equal to 0, we keep only weights that are greater than 0
  w01 <- ww1>EPS
  w02 <- ww2>EPS
  
  # keep track of the right indices (some weights are put to 0 and then the indices have changes)
  id1_tmp <- id1[w01]
  id2_tmp <- id2[w02]
  
  # keep weights greter than 0
  www1 <- ww1[w01]
  www2 <- ww2[w02]
  
  # adapt distance to optimal transport
  DD <-D[w01,w02]
  
  # optimal transport
  res <- transport::transport(www1,www2,DD)
  
  # 
  U <- cbind(id1 = id1_tmp[res[,1]],id2 = id2_tmp[res[,2]],weight = res[,3] , X1[w01,][res[,1],] , X2[w02,][res[,2],])
  
  # for common unit 
  D1 <- data.frame(w1,id1,X1 = X1)
  D2 <- data.frame(w2,id2,X2 = X2)
  V <- merge(D1,D2,by.x="id1",by.y="id2")
  V <- cbind(id1=V$id1,id2=V$id1,weight=pmin(V$w1,V$w2),V[,-which(names(V)%in% c("w1","w2","id1"))])
  
  colnames(U)[4:ncol(U)] <- colnames(V)[4:ncol(V)]
  
  
  conc <- rbind(V,U)
  conc <- conc[order(conc[,1]),]
  
  
  out <- list(res = res,
             conc = conc)
  class(out) <- "match"
  
  return(out)
  
  # E <- cbind(id1=id1[res[,1]],id2=id2[res[,2]],res[,3] , X1[ww1>EPS,][res[,1],] , X2[ww2>EPS,][res[,2],])
  #
 
  #
  # rbind(final,as.matrix(E))
}





