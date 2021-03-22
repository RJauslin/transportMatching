rm(list = ls())
set.seed(1)

library(survey)
library(StatMatch)
library(MASS)
# setup

N <- 10000
X <- data.frame(x1 = rnorm(N,0,1), x2 = rnorm(N,0,1))

mu <- c(4,9)
Sigma <- matrix(c(1,0.3,0.3,1),ncol = 2)
X <- mvrnorm(N,mu,Sigma)
colnames(X) <- c("x1","x2")
X <- as.data.frame(X)

Y <- data.frame(y1 = X$x1^2 + X$x2^2 + rnorm(N,0,1),y2 = exp(X$x1) + rnorm(N,0,0.1))
Z <- data.frame(z1 = log(abs(X$x1)) + rnorm(N,0,1),z2 = sqrt(abs(X$x2)) + rnorm(N,0,0.1))

# Y <- data.frame(y1 = 3*X$x1 + 4*X$x2 + rnorm(N,0,1),y2 = X$x1 + rnorm(N,0,0.1))
# Z <- data.frame(z1 = 8*X$x1 - 3*X$x2 + rnorm(N,0,1),z2 = abs(X$x2) + rnorm(N,0,0.1))



n1=1000
n2=3000

s1=srswor(n1,N)
s2=srswor(n2,N)

id1=(1:N)[s1==1]
id2=(1:N)[s2==1]

d1=rep(N/n1,n1)
d2=rep(N/n2,n2)

X1 = X[s1==1,]
X2 = X[s2==1,]
Y1 <- data.frame(Y[s1 == 1,])
Z2 <- data.frame(Z[s2 == 1,])




re <- harmonize(X1,d1,id1,X2,d2,id2)

w1=re$w1
w2=re$w2


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






XX.w1 = t(as.matrix(X1))%*% (as.matrix(X1) * w1)
XX.w2 = t(as.matrix(X2))%*% (as.matrix(X2) * w2)

gamma.p <- n1/(n1 + n2)
XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
YZ.CIA <- t(f1$coefficients) %*% XX.pool %*% f2$coefficients
YZ.CIA/N
t(Y2_pred_c)%*%Z2_pred_c/nrow(Y2_pred_c)
t(Y1_pred_c)%*%Z1_pred_c/nrow(Y1_pred_c)



# # f1c <- apply(f1$fitted.values,MARGIN = 2, function(x){x - mean(x)})
# f2c <- apply(f2$fitted.values,MARGIN = 2, function(x){x - mean(x)})
# f1c%*%f2c


Yc <- apply(Y,MARGIN = 2, function(x){x - mean(x)})
Zc <- apply(Z,MARGIN = 2, function(x){x - mean(x)})

t(Yc)%*%Zc/(N-1)
cov(Y,Z,method = "pearson")





##############


# cst <- sum(w1)
# w1 <- w1/sum(w1)*cst
# w2 <- w2/sum(w2)*cst
object = linkage(X1,id1,X2,id2,w1,w2)
# object <- statMatch_random(object)

Y_opt <- Y[as.character(object$conc$id1),]
Z_opt <- Z[as.character(object$conc$id2),]
Y_opt_c <- apply(Y_opt,MARGIN = 2, function(x){x - mean(x)})
Z_opt_c <- apply(Z_opt,MARGIN = 2, function(x){x - mean(x)})


t(Y_opt_c)%*%Z_opt_c/(nrow(object$conc))
cov(Y,Z,method = "pearson")
t(Y2_pred_c)%*%Z2_pred_c/nrow(Y2_pred_c)
t(Yc)%*%Zc/(N-1)

