
library(parallel)
library(ppcor)
rm(list = ls())

numCores <- detectCores()
numCores
N <- 30000
mu <- c(4,9,20,30,100,200)


tmp <- matrix(rnorm(36),nrow=6)
Sigma <- tmp%*%t(tmp)

# tmp <- c(1,0.3,0.4,0.5,0.6,0.7)
# Sigma <- tmp%*%t(tmp)
# diag(Sigma) <- rep(1,6)

dat <- mvrnorm(N,mu,Sigma)
Xm <- dat[,1:2]
colnames(Xm) <- c("x1","x2")
Xm <- as.data.frame(Xm)

Y <- dat[,3:4]
colnames(Y) <- c("y1","y2")
Y <- as.data.frame(Y)
Z <- dat[,5:6]
colnames(Z) <- c("z1","z2")
Z <- as.data.frame(Z)
# Y <- data.frame(y1 = Xm$x1^2  + rnorm(N,0,1),y2 = exp(Xm$x1) + rnorm(N,0,0.1))
# Z <- data.frame(z1 = (1/30)*Xm$x2^3 + rnorm(N,0,1),z2 = sqrt(abs(Xm$x2)) + rnorm(N,0,0.1))

# Y <- data.frame(y1 = 3*X$x1 + 4*X$x2 + rnorm(N,0,1),y2 = X$x1 + rnorm(N,0,0.1))
# Z <- data.frame(z1 = 8*X$x1 - 3*X$x2 + rnorm(N,0,1),z2 = abs(X$x2) + rnorm(N,0,0.1))


# round(pcor(cbind(Xm,Y,Z))$estimate,4)

n1=1000
n2=3000


# function to call for parLapply

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


# cluster creation

cl <- makeCluster(detectCores())
clusterEvalQ(cl,{
  library(devtools)
  devtools::load_all("C:/Users/jauslinr/switchdrive/matching_optimal_transport/transportMatching")
})


# call parLapply
SIM = 100
l_sim <- parLapply(cl = cl,
                   X = 1:SIM,
                   fun = sim,
                   Xm = Xm,
                   Y = Y,
                   Z = Z,
                   n1 = n1,
                   n2 = n2)

# close cluster
stopCluster(cl)


true <- lapply(l_sim, FUN = function(x){x[[1]]})[[1]]
l_full <- lapply(l_sim, FUN = function(x){x[[2]]})
l_balanced <- lapply(l_sim, FUN = function(x){x[[3]]})
l_ren_1 <- lapply(l_sim, FUN = function(x){x[[4]]})
l_ren_2 <- lapply(l_sim, FUN = function(x){x[[5]]})

# biais

full <- balanced <- ren_1 <- ren_2 <- matrix(rep(0,4),ncol = 2)

for(i in 1:length(l_sim)){
  full <- full + l_full[[i]]/SIM
  balanced <- balanced + l_balanced[[i]]/SIM
  ren_1 <- ren_1 + l_ren_1[[i]]/SIM
  ren_2 <- ren_2 + l_ren_2[[i]]/SIM
}


b_full <- full - true
b_balanced <- balanced - true
b_ren_1 <- ren_1 - true
b_ren_2 <- ren_2 - true

# variance
var_full <- var_balanced <- var_ren_1 <- var_ren_2 <- matrix(rep(0,4),ncol = 2)

for(i in 1:length(l_sim)){
  var_full <- var_full + (l_full[[i]]- full)^2/SIM
  var_balanced <- var_balanced + (l_balanced[[i]]-balanced)^2/SIM
  var_ren_1 <- var_ren_1 + (l_ren_1[[i]]-ren_1)^2/SIM
  var_ren_2 <- var_ren_2 + (l_ren_2[[i]]-ren_2)^2/SIM
}


# mse
mse_full <- mse_balanced <- mse_ren_1 <- mse_ren_2 <- matrix(rep(0,4),ncol = 2)

for(i in 1:length(l_sim)){
  mse_full <- mse_full + (l_full[[i]]- true)^2/SIM
  mse_balanced <- mse_balanced + (l_balanced[[i]]-true)^2/SIM
  mse_ren_1 <- mse_ren_1 + (l_ren_1[[i]]-true)^2/SIM
  mse_ren_2 <- mse_ren_2 + (l_ren_2[[i]]-true)^2/SIM
}

b_full^2 + var_full
mse_full

tabNOCIA <- rbind(b_full^2,var_full,mse_full)
tabNOCIA <- cbind(tabNOCIA, rbind(b_balanced^2,var_balanced,mse_balanced))
tabNOCIA <- cbind(tabNOCIA, rbind(b_ren_1^2,var_ren_1,mse_ren_1))
# tab <- cbind(tab, rbind(b_ren_2^2,var_ren_2,mse_ren_2))


saveRDS(tabNOCIA,file = "C:/Users/jauslinr/switchdrive/matching_optimal_transport/paper/tabNOCIA.rds")

