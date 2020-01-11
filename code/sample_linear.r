####################################################################################
## Estimate parameters using ADMM for high-dimensional linear regression 
## influence beta_2 and beta_3
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta_N)+rho v^T (beta-theta)+0.5rho ||beta-theta||_2^2
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v: the Langrange multiplier
######################################################
library(SIS)
library(ncvreg)

# calculate the weight given a sub-matrix
# return mean(sub_X[,j]^2)-cov(sub_X[,j],sub_X[,-j])^TCov(sub_X[,-j],sub_X[,-j])^{-1}cov(sub_X[,j],sub_X[,-j])
# j is the index indicating which variable needed for inference
omega <- function(sub_X, j){
  # extract the dimension
  sub_n <- dim(sub_X)[1]
  sub_p <- dim(sub_X)[2]
  
  if (sub_p==1){
    reg <- NULL
    sig2 <- mean(sub_X[,j]^2)
  }
  else{
    rho <- crossprod(sub_X[,j],sub_X[,-j])
    reg <- solve(crossprod(sub_X[,-j],sub_X[,-j]), t(rho))
    sig2 <- mean(sub_X[,j]^2)- rho%*%reg/sub_n
  }
  
  return(list(reg=reg, sig2=sig2))
}

# ADMM step: update beta
# beta=(rho I+ X^T X/n)^{-1} (X^T Y/n+ rho(theta-v))
CI <- function(X, Y, sn=20){
  n <- length(Y)
  p <- dim(X)[2]
  # the design matrix including intercept
  #  X0 <- cbind(1, X)
  
  # initialize
  weight <- matrix(0, n, 4)
  pres <- matrix(0, n, 4)
  
  # initial estimator for beta
  re <- SIS(X, Y)
  support0 <- rep(FALSE, p)
  support0[re$ix] <- TRUE
  beta0 <- rep(0, p)
  beta0[re$ix] <- re$coef.est[-1]
  
  for (i in (sn+1):n){
    # sure screening
    re <- SIS(X[1:(i-1),], Y[1:(i-1)])
    
    # weight for the third coefficient
    support <- rep(FALSE, p)
    support[re$ix] <- TRUE
    support.copy <- support
    support.copy[3] <- TRUE
    X_sub <- as.matrix(X[, support.copy])
    re1 <- omega(X_sub, j=(1+sum(support.copy[1:2])))
    
    # refit for the third coefficient
    support.copy[3] <- FALSE
    pres[i-sn, 1] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re1)){
      weight[i-sn, 1] <- X[i,3]/sqrt(re1$sig2)
    }
    else{
      support.copy[3] <- FALSE
      weight[i-sn, 1] <- (X[i,3] - X[i,support.copy] %*% re1$reg)/sqrt(re1$sig2)
    }
    
    # weight for the fourth coefficient
    support.copy <- support
    support.copy[4] <- TRUE
    X_sub <- as.matrix(X[, support.copy])
    re2 <- omega(X_sub, j=(1+sum(support.copy[1:3])))
    
    # refit for the fourth coefficient
    support.copy[4] <- FALSE
    pres[i-sn, 2] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re2)){
      weight[i-sn, 2] <- X[i,4]/sqrt(re2$sig2)
    }
    else{
      support.copy[4] <- FALSE
      weight[i-sn, 2] <- (X[i,4] - X[i,support.copy] %*% re2$reg)/sqrt(re2$sig2)
    }
    
    # weight for the fourth coefficient
    support.copy <- support
    support.copy[5] <- TRUE
    X_sub <- as.matrix(X[, support.copy])
    re3 <- omega(X_sub, j=(1+sum(support.copy[1:4])))
    
    # refit for the fourth coefficient
    support.copy[5] <- FALSE
    pres[i-sn, 3] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re3)){
      weight[i-sn, 3] <- X[i,5]/sqrt(re3$sig2)
    }
    else{
      support.copy[5] <- FALSE
      weight[i-sn, 3] <- (X[i,5] - X[i,support.copy] %*% re3$reg)/sqrt(re3$sig2)
    }
    
    # weight for the fifth coefficient
    support.copy <- support
    support.copy[6] <- TRUE
    X_sub <- as.matrix(X[, support.copy])
    re4 <- omega(X_sub, j=(1+sum(support.copy[1:5])))
    
    # refit for the fifth coefficient
    support.copy[6] <- FALSE
    pres[i-sn, 4] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re4)){
      weight[i-sn, 4] <- X[i,6]/sqrt(re4$sig2)
    }
    else{
      support.copy[6] <- FALSE
      weight[i-sn, 4] <- (X[i,6] - X[i,support.copy] %*% re4$reg)/sqrt(re4$sig2)
    }
  }
  
  # for the third few observations
  re <- SIS(X[(sn+1):n,], Y[(sn+1):n])
  support <- rep(FALSE, p)
  support[re$ix] <- TRUE
  
  # estimate the parameter for the third element
  support.copy <- support
  support.copy[3] <- TRUE
  X_sub <- as.matrix(X[, support.copy])
  re1 <- omega(X_sub, j=(1+sum(support.copy[1:2])))
  
  # refit for the third coefficient
  for (i in 1:sn){
    support.copy[3] <- FALSE
    pres[i+n-sn, 1] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re1)){
      weight[i+n-sn, 1] <- X[i,3]/sqrt(re1$sig2)
    }
    else{
      support.copy[3] <- FALSE
      weight[i+n-sn, 1] <- (X[i,3] - X[i,support.copy] %*% re1$reg)/sqrt(re1$sig2)
    }
  }
  
  # estimate the parameter for the fourth element
  support.copy <- support
  support.copy[4] <- TRUE
  X_sub <- as.matrix(X[, support.copy])
  re2 <- omega(X_sub, j=(1+sum(support.copy[1:3])))
  
  # refit for the fourth coefficient
  for (i in 1:sn){
    support.copy[4] <- FALSE
    pres[i+n-sn, 2] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re2)){
      weight[i+n-sn, 2] <- X[i,4]/sqrt(re2$sig2)
    }
    else{
      support.copy[4] <- FALSE
      weight[i+n-sn, 2] <- (X[i,4] - X[i,support.copy] %*% re2$reg)/sqrt(re2$sig2)
    }
  }
  
  # estimate the parameter for the fifth element
  support.copy <- support
  support.copy[5] <- TRUE
  X_sub <- as.matrix(X[, support.copy])
  re3 <- omega(X_sub, j=(1+sum(support.copy[1:4])))
  
  # refit for the fifth coefficient
  for (i in 1:sn){
    support.copy[5] <- FALSE
    pres[i+n-sn, 3] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re3)){
      weight[i+n-sn, 3] <- X[i,5]/sqrt(re3$sig2)
    }
    else{
      support.copy[5] <- FALSE
      weight[i+n-sn, 3] <- (X[i,5] - X[i,support.copy] %*% re3$reg)/sqrt(re3$sig2)
    }
  }
  
  # estimate the parameter for the sixth element
  support.copy <- support
  support.copy[6] <- TRUE
  X_sub <- as.matrix(X[, support.copy])
  re4 <- omega(X_sub, j=(1+sum(support.copy[1:5])))
  
  # refit for the sixth coefficient
  for (i in 1:sn){
    support.copy[6] <- FALSE
    pres[i+n-sn, 4] <- Y[i] - X[i,support.copy] %*% beta0[support.copy]
    
    if (is.null(re4)){
      weight[i+n-sn, 4] <- X[i,6]/sqrt(re4$sig2)
    }
    else{
      support.copy[6] <- FALSE
      weight[i+n-sn, 4] <- (X[i,6] - X[i,support.copy] %*% re4$reg)/sqrt(re4$sig2)
    }
  }
  
  ## aggregate the results
  beta <- rep(0, 4)
  se <- rep(0, 8)
  
  beta[1] <- mean((weight[,1]*pres[,1]))/mean((weight[,1]*X[c((sn+1):n,1:sn),3]))
  beta[2] <- mean((weight[,2]*pres[,2]))/mean((weight[,2]*X[c((sn+1):n,1:sn),4]))
  beta[3] <- mean((weight[,3]*pres[,3]))/mean((weight[,3]*X[c((sn+1):n,1:sn),5]))
  beta[4] <- mean((weight[,4]*pres[,4]))/mean((weight[,4]*X[c((sn+1):n,1:sn),6]))
  
  ## refitted cross-validation to estimate sigma
  indice1 <- sample(1:n, n/2)
  indice2 <- setdiff(1:n, indice1)
  
  re <- SIS(X[indice1,], Y[indice1])
  sig21 <- sigma(lm(Y[indice2]~X[indice2,re$ix]-1))
  
  re <- SIS(X[indice2,], Y[indice2])
  sig22 <- sigma(lm(Y[indice1]~X[indice1,re$ix]-1))
  
  sig <- sqrt((sig21^2+sig22^2)/2)
  
  ## standard error
  se[1] <- (sqrt(n)/sum(weight[,1]*X[c((sn+1):n,1:sn),3]))*sig
  se[2] <- (sqrt(n)/sum(weight[,2]*X[c((sn+1):n,1:sn),4]))*sig
  se[3] <- (sqrt(n)/sum(weight[,3]*X[c((sn+1):n,1:sn),5]))*sig
  se[4] <- (sqrt(n)/sum(weight[,4]*X[c((sn+1):n,1:sn),6]))*sig
  
  se[5] <- (sqrt(sum(weight[,1]^2))/sum(weight[,1]*X[c((sn+1):n,1:sn),3]))*sig
  se[6] <- (sqrt(sum(weight[,2]^2))/sum(weight[,2]*X[c((sn+1):n,1:sn),4]))*sig
  se[7] <- (sqrt(sum(weight[,3]^2))/sum(weight[,3]*X[c((sn+1):n,1:sn),5]))*sig
  se[8] <- (sqrt(sum(weight[,4]^2))/sum(weight[,4]*X[c((sn+1):n,1:sn),6]))*sig
  
  return(list(beta=beta, se=se))
}

##############################################################################
## Generating simulation models
## beta_1=1.5, beta_2=-1.5, beta_j=0 for j>=3
## rho: Cov(X_j, X_i)=rho^{|i-j|}
## n: sample size
## p: dimension 
## sig: Var(Y|X)
##############################################################################
Gen_Model <- function(n=100, p=500, sig=1, rho=0){
  if (rho==0)
    X <- matrix(rnorm(n*p), n, p)
  else{
    X <- matrix(0, n, p)
    for (i in 1:n){
      X[i, ] <- arima.sim(model=list(ar=rho), sd=sqrt(1-rho^2), n=p)
    }
  }
  
  Y <- X[,1]+X[,2]+X[,3]+X[,4]+X[,5]+sig*rnorm(n)
  
  return(list(X=X, Y=Y))
}

################################################################################
## Simulation code for constructing confidence intervals for beta_3
## L is the simulation replication 
################################################################################
  # cp initialization
beta.online <- matrix(0, 500, 4)
se.online <- matrix(0, 500, 8)

for (l in 1:25){
  print(l)
  
  set.seed(12345*l)
  Model <- Gen_Model(n=200, p=1000, rho=0.5)
  X <- Model$X
  Y <- Model$Y
  n <- length(Y)
  result1 <- CI(X, Y, sn=floor(2*n/log(n)))
  
  beta.online[l,] <- result1$beta; se.online[l,] <- result1$se

 save.image(file="simu_online_nonsparse_rho05_new_seed1.RData")
} 
