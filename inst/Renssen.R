rm(list = ls())
# set.seed(6)

library(survey)
library(StatMatch)

# setup

N <- 10000
X <- as.data.frame(cbind(sex = rbinom(N,3,0.5), age = rbinom(N,6,prob = 0.5)))
tapply(rep(1,nrow(X)),list(X$sex,X$age),FUN = sum)


# Y <- X[,2]^2

# Y <- rbinom(N,4,0.4)+ 1
# Z <- rbinom(N,3,0.5)+ 10

Y <- ceiling(log(X[,2]+1) + rbinom(N,2,0.3))
# Y <- 2*X[,1]^2 + X[,2] + ceiling(rnorm(N,0,0.1))
Z <- 0.5*X[,1]^1/2 + ceiling(rnorm(N,0,0.01))


# YZ <- addmargins(table(Y,Z))

YZ <- tapply(rep(1,nrow(X)),list(Y,Z),FUN = sum)
YZ[is.na(YZ)] <- 0
YZ <- addmargins(YZ)
YZ


for(i in 1:ncol(X)){
  X[,i] <-  factor(X[,i])
}
Y <- data.frame(y = factor(Y))
Z <- data.frame(Z = factor(Z))




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



Y1 <- data.frame(y = factor(Y[s1 == 1,]))
rownames(Y1) <- id1
Z2 <- data.frame(z = factor(Z[s2 == 1,]))
rownames(Z2) <- id2

X1_dis <- model.matrix(~age:sex-1,data = X1)
X2_dis <- model.matrix(~age:sex-1,data = X2)
Y1_dis <- model.matrix(~y-1,data = Y1)
Z2_dis <- model.matrix(~z-1,data = Z2)



X1_ren <- cbind(X1,Y1)
X2_ren <- cbind(X2,Z2)


# creates svydesign objects
svy.X1 <- svydesign(~1, weights=~d1, data=X1_ren)
svy.X2 <- svydesign(~1, weights=~d2, data=X2_ren)


out.hz <- harmonize.x(svy.A=svy.X1,
                      svy.B=svy.X2,
                      form.x=~age:sex-1,
                      # x.tot =  colSums(model.matrix(~age:sex-1,data = X)),
                      cal.method="raking")





#------------------------------------- optimal transport


re <- harmonize(X1_dis,d1,id1,X2_dis,d2,id2,method = "raking")

# population totals
# totals <- colSums(model.matrix(~age:sex-1,data = X))
# totals <- c(N,totals)
# re <- harmonize(X1_dis,d1,id1,X2_dis,d2,id2,method = "raking",totals = as.numeric(totals))


w1=re$w1
w2=re$w2

w1 <- out.hz$weights.A
w2 <- out.hz$weights.B

# cst <- sum(w1)
# w1 <- w1/sum(w1)*cst
# w2 <- w2/sum(w2)*cst
object = linkage(X1_dis,id1,X2_dis,id2,w1,w2)
out <- statMatch_random(object) # statmatch random
out_test <- statMatch_random(object,Z2) # statmatch random

colSums(w2*Z2_dis)
colSums((out$conc$weight/out$q)*disjunctive(Z2[as.character(out$conc$id2),]))
colSums((out_test$conc$weight/out_test$q)*disjunctive(Z2[as.character(out_test$conc$id2),]))


sum((colSums(w2*Z2_dis) - colSums((out$conc$weight/out$q)*disjunctive(Z2[as.character(out$conc$id2),])))^2)
sum((colSums(w2*Z2_dis) - colSums((out_test$conc$weight/out_test$q)*disjunctive(Z2[as.character(out_test$conc$id2),])))^2)

# out <- out_test


out_exp <- statMatch_expected(object,Z2) # statmatch random
round(colSums(w1*out_exp$conc[,2:ncol(out_exp$conc)]),3)

#---------------------------------  check
diag(as.matrix(t(object$conc[,4:13]))%*% (as.matrix(object$conc[,14:ncol(object$conc)]) * object$conc[,3]))
out_p <- tapply(object$conc[,3],list(X1[as.character(object$conc[,1]),]$sex,X1[as.character(object$conc[,1]),]$age),sum)
out_p
tapply(w1,list(X1$sex, X1$age),sum)
tapply(w2,list(X2$sex, X2$age),sum)
colSums(w1*X1_dis) # totals from X1
colSums(w2*X2_dis) # totals from X2
tapply(rep(1,N),list(X$sex, X$age),sum) # from population


#---------------------------------- Data integration from statmatch random
# do.call(c,lapply(out$conc$id2,FUN = function(x){which(x == id2)})) # OLD


Y1_random <- cbind(X1[as.character(out$conc$id1),],y = Y1[as.character(out$conc$id1),])
Z2_random <- cbind(X2[as.character(out$conc$id2),],z = Z2[as.character(out$conc$id2),])
YZ_random <- tapply(out$conc$weight/out$q,list(Y1_random$y,Z2_random$z),sum)
YZ_random[is.na(YZ_random)] <- 0
YZ_random <- addmargins(YZ_random)


#---------------------------------- Data integration from optimal transport

Y1_optimal <- cbind(X1[as.character(object$conc$id1),],y = Y1[as.character(object$conc$id1),])
Z2_optimal <- cbind(X2[as.character(object$conc$id2),],z = Z2[as.character(object$conc$id2),])
YZ_optimal <- tapply(object$conc[,3],list(Y1_optimal$y,Z2_optimal$z),sum)
YZ_optimal[is.na(YZ_optimal)] <- 0
YZ_optimal <- addmargins(YZ_optimal)


#---------------------------------- Data integration from prediction


YZ_exp <-  addmargins( t(disjunctive(as.matrix(Y1)))%*%( out_exp$Z * out_exp$weight))
YZ_exp


#---------------------------------- Renssen via StatMatch



tapply(out.hz$weights.A,list(X1$sex, X1$age),sum)
tapply(out.hz$weights.B,list(X2$sex, X2$age),sum)


outStat <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                        svy.C=NULL, y.lab="y", z.lab="z",
                        form.x=~age:sex-1,
                        micro = TRUE)

YZ_ren <- addmargins(outStat$yz.CIA)



### HOW THE YZ_ren is calculated

svy.A=out.hz$cal.A
svy.B=out.hz$cal.B
svy.C=NULL
y.lab="y"
z.lab="z"
form.x=~age:sex-1



f1 <- lm(disjunctive(as.matrix(Y1))~X1_dis-1)
lm(Y1_dis ~ X1_dis-1)
f1
f2 <- lm(disjunctive(as.matrix(Z2))~X2_dis-1)
f2

XX.w1 = t(X1_dis)%*% (X1_dis * w1)
XX.w2 = t(X2_dis)%*% (X2_dis * w2)

gamma.p <- n1/(n1 + n2)
XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
YZ.CIA <- t(f1$coefficients) %*% XX.pool %*% f2$coefficients
YZ.CIA
YZ_ren


# 
# 
# outStat
# outStat$yz.CIA
# YZ2 <- addmargins(outStat$yz.CIA)




#---------------------------------- Comparison results


# YZ_random
YZ_optimal
round(YZ_ren,10)
# YZ_exp
YZ

sum((YZ_random[-nrow(YZ_random),-ncol(YZ_random)] - YZ[-nrow(YZ),-ncol(YZ)])^2)
sum(abs(YZ_optimal[-nrow(YZ_optimal),-ncol(YZ_optimal)] - YZ[-nrow(YZ),-ncol(YZ)])^2)
sum(abs(YZ_ren[-nrow(YZ_ren),-ncol(YZ_ren)] - YZ[-nrow(YZ),-ncol(YZ)])^2)

