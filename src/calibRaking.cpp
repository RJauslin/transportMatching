#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title 
//'
//' @description
//' 
//'
//' @param X 
//' @param k 
//' @param tore 
//' @param toreBound 
//'
//' @details
//'
//'
//'
//' @return 
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//'
//' @seealso
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calibRaking(arma::mat Xs,
                    arma::vec d,
                    arma::vec total,
                    arma::vec q,
                    int max_iter = 500,
                    double tol = 1e-9){
  
  // intiializing
  int n = Xs.n_rows;
  int p = Xs.n_cols;
  
  arma::vec lambda(p,arma::fill::zeros);
  arma::vec w1 = d % exp(Xs * lambda % q);
  
  // tol and iteration counter
  double crit = 1;
  int iter = 0;
  
  // convergence loop
  do{
    arma::vec phi = Xs.t() * w1 - total;
    arma::mat phiprim = Xs.t() * arma::diagmat(w1) * Xs;
    lambda = lambda - arma::pinv(phiprim) * phi;
    w1 = d % exp(Xs * lambda % q);
    arma::vec tr = Xs.t() * w1;
    crit = max(abs(tr - total)/total);
    iter = iter + 1;
  } while (crit > tol & iter < max_iter);
  
  // return(w1/d);
  Rcpp::NumericVector out = Rcpp::wrap(w1/d);
  out.attr("dim") = R_NilValue;
  return(out);
}



/*** R

rm(list = ls())
library(MASS)
N <- 10000
X <- data.frame(x1 = rnorm(N,0,1), x2 = rnorm(N,0,1))

mu <- c(4,9)
Sigma <- matrix(c(1,0.3,0.3,1),ncol = 2)
X <- mvrnorm(N,mu,Sigma)
colnames(X) <- c("x1","x2")
X <- as.data.frame(X)

Y <- data.frame(y1 = 3*X$x1^2 + 4*X$x2^2 + rnorm(N,0,1),y2 = exp(X$x1) + rnorm(N,0,0.1))
Z <- data.frame(z1 = 8*X$x1^2 - 3*X$x2^2 + rnorm(N,0,1),z2 = sqrt(abs(X$x2)) + rnorm(N,0,0.1))

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


# number of units in each sample
n1 <- nrow(X1)
n2 <- nrow(X2)

# add constant vector to ensure same sum
XX1=cbind(rep(1,n1),X1)
XX2=cbind(rep(1,n2),X2)

# we can specify the desired total (for example if we know the totals of the population)
n12=length(intersect(id1,id2))
a=(n1-n12)/(n1+n2-2*n12)  
totals=a*colSums(d1*XX1)+(1-a)*colSums(d2*XX2)

q <- rep(1,length(d1))
# calibration with sampling package
# w1=d1*calib(XX1,d1,totals,method=method)

d1*calibRaking(as.matrix(XX1),d1,totals,q)


*/




