
library(parallel)
library(ppcor)
rm(list = ls())

numCores <- detectCores()
numCores
N <- 30000
mu <- c(4,9)
Sigma <- matrix(c(1,0.3,0.3,1),ncol = 2)
Xm <- mvrnorm(N,mu,Sigma)
colnames(Xm) <- c("x1","x2")
Xm <- as.data.frame(Xm)

Y <- data.frame(y1 = Xm$x1^2  + rnorm(N,0,1),y2 = exp(Xm$x1) + rnorm(N,0,0.1))
Z <- data.frame(z1 = (1/30)*Xm$x2^3 + rnorm(N,0,1),z2 = sqrt(abs(Xm$x2)) + rnorm(N,0,0.1))

# Y <- data.frame(y1 = 3*X$x1 + 4*X$x2 + rnorm(N,0,1),y2 = X$x1 + rnorm(N,0,0.1))
# Z <- data.frame(z1 = 8*X$x1 - 3*X$x2 + rnorm(N,0,1),z2 = abs(X$x2) + rnorm(N,0,0.1))


round(pcor(cbind(Xm,Y,Z))$estimate,4)

n1=1000
n2=3000



sim <- function(n,Xm,Y,Z,n1,n2){
  N <- nrow(Xm)
  s1=srswor(n1,N)
  s2=srswor(n2,N)
  
  id1=(1:N)[s1==1]
  id2=(1:N)[s2==1]
  
  d1=rep(N/n1,n1)
  d2=rep(N/n2,n2)
  
  X1 = Xm[s1==1,]
  X2 = Xm[s2==1,]
  Y1 <- data.frame(Y[s1 == 1,])
  Z2 <- data.frame(Z[s2 == 1,])
  
  re <- harmonize(X1,d1,id1,X2,d2,id2)
  
  w1 = re$w1
  w2 = re$w2
  
  # optimal transport full
  object = linkage(X1,id1,X2,id2,w1,w2)
  Y_opt <- Y[as.character(object$conc$id1),]
  Z_opt <- Z[as.character(object$conc$id2),]
  Y_opt_c <- apply(Y_opt,MARGIN = 2, function(x){x - mean(x)})
  Z_opt_c <- apply(Z_opt,MARGIN = 2, function(x){x - mean(x)})
  
  
  full <- t(Y_opt_c)%*%Z_opt_c/(nrow(object$conc))
  
  # optimal transport balanced
  out_r <- statMatch_random(object)
  Y_opt_r <- Y[as.character(out_r$conc$id1),]
  Z_opt_r <- Z[as.character(out_r$conc$id2),]
  Y_opt_r_c <- apply(Y_opt_r,MARGIN = 2, function(x){x - mean(x)})
  Z_opt_r_c <- apply(Z_opt_r,MARGIN = 2, function(x){x - mean(x)})
  
  balanced <- t(Y_opt_r_c)%*%Z_opt_r_c/(nrow(out_r$conc))
  
  # Renssen
  f1 <- lm(as.matrix(Y1) ~ X1[,1] + X1[,2] -1)
  Y2_pred <- as.matrix(X2)%*%f1$coefficients
  f2 <- lm(as.matrix(Z2) ~ X2[,1] + X2[,2] -1)
  Z2_pred <- as.matrix(X2)%*%f2$coefficients
  
  Y1_pred <- as.matrix(X1)%*%f1$coefficients
  Z1_pred <- as.matrix(X1)%*%f2$coefficients
  
  
  Y2_pred_c <- apply(Y2_pred,MARGIN = 2, function(x){x - mean(x)})
  Z2_pred_c <- apply(Z2_pred,MARGIN = 2, function(x){x - mean(x)})
  Y1_pred_c <- apply(Y1_pred,MARGIN = 2, function(x){x - mean(x)})
  Z1_pred_c <- apply(Z1_pred,MARGIN = 2, function(x){x - mean(x)})
  
  ren_1 <- t(Y1_pred_c)%*%Z1_pred_c/nrow(Y1_pred_c)
  ren_2 <- t(Y2_pred_c)%*%Z2_pred_c/nrow(Y2_pred_c)
  
  return(list( true = cov(Y,Z),
               full = full,
               balanced = balanced,
               ren_1 = ren_1,
               ren_2 = ren_2
               ))
  
}



cl <- makeCluster(detectCores())
clusterEvalQ(cl,{
  library(devtools)
  devtools::load_all("C:/Users/jauslinr/switchdrive/matching_optimal_transport/transportMatching")
})

SIM = 100
l_sim <- parLapply(cl = cl,
          X = 1:SIM,
          fun = sim,
          Xm = X,
          Y = Y,
          Z = Z,
          n1 = n1,
          n2 = n2)

stopCluster(cl)

true <- full <- balanced <- ren_1 <- ren_2 <- matrix(rep(0,4),ncol = 2)
for(i in 1:length(l_sim)){
  true <- true + l_sim[[i]]$true
}
true <- true/SIM

mse_full <- matrix(rep(0,4),ncol = 2)
full <- matrix(rep(0,4),ncol = 2)
mse_ren_1 <- matrix(rep(0,4),ncol = 2)
ren_1 <- matrix(rep(0,4),ncol = 2)
for(i in 1:length(l_sim)){
  mse_full <- mse_full + (l_sim[[i]]$full-true)^2
  full <- full + l_sim[[i]]$full
  mse_ren_1 <- mse_ren_1 + (l_sim[[i]]$ren_1-true)^2
  ren_1 <- ren_1 + l_sim[[i]]$ren_1
  # full <- (full-true)^2 + l_sim[[i]]$full
  # balanced <- balanced + l_sim[[i]]$balanced
  # ren_1 <- ren_1 + l_sim[[i]]$ren_1
  # ren_2 <- ren_2 + l_sim[[i]]$ren_2
}


mse_full <- mse_full/SIM 
full <- full/SIM
b_full <- full-true
b_full^2
mse_full

mse_ren_1 <- mse_ren_1/SIM
ren_1 <- ren_1/SIM
b_ren_1 <- ren_1-true
b_ren_1^2
mse_ren_1

true/100 - balanced/100
true/100 - ren_1/100
true/100 - ren_2/100

mat <- do.call(rbind,l_sim)

# for(i in 1:ncol(mat)){

apply(mat,MARGIN = 2,FUN = function(x){
  
  tmp <- x
  tmp2 <- x[[1]]
  print(tmp2)
  for(j in 2:length(tmp2)){
    tmp2 <- tmp2 + x[[j]]
  }
  # print(tmp2)
  return(tmp2)
},simplify = FALSE)

i = 1
  tmp <- mat[,i][[1]]
  for(j in 2:nrow(mat)){
    tmp <- tmp + mat[,i][[j]]
    
  }
  tmp/10
# }

          