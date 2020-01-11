####################################################################################
## Estimate parameters using ADMM for high-dimensional logistic regression model
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
# hess: the diagonal of the negative Hessian matrix
omega <- function(sub_X, j, hess){
  # extract the dimension
  sub_n <- dim(sub_X)[1]
  sub_p <- dim(sub_X)[2]
  
  if (sub_p==1){
    reg <- NULL
    sig2 <- mean(sub_X[,j]^2*hess)
  }
  else{
    rho <- crossprod(sub_X[,j], hess*sub_X[,-j])
    reg <- solve(crossprod(sub_X[,-j], hess*sub_X[,-j]), t(rho))
    sig2 <- mean(sub_X[,j]^2*hess)-rho%*%reg/sub_n
  }
  
  return(list(reg=reg, sig2=sig2))
}

# select model using sure screening
MD <- function(X, Y, sn=20){
  n <- length(Y)
  p <- dim(X)[2]
  
  # variable selection based on the first few observations
  md <- list()
  for (i in ((sn+1):n)){
    re <- try(SIS(X[1:(i-1),], Y[1:(i-1)], family="binomial"), silent = T)
    if (!inherits(re, "try-error")){
      md[[i-sn]] <- re$ix
    }
    else{
      print("use SCAD only")
      re <- cv.ncvreg(X[1:(i-1),], Y[1:(i-1)], family="binomial", penalty="SCAD")
      md[[i-sn]] <- as.vector(which(re$fit$beta[-1,re$min]!=0))
    }
  }
  # variable selection based on the last couple observations
  re <- try(SIS(X[(sn+1):n,], Y[(sn+1):n], family="binomial"), silent = T)
  if (!inherits(re, "try-error")){
    md[[n-sn+1]] <- re$ix
  }
  else{
    print("use SCAD only")
    re <- cv.ncvreg(X[(sn+1):n,], Y[(sn+1):n], family="binomial", penalty="SCAD")
    md[[n-sn+1]] <- as.vector(which(re$fit$beta[-1,re$min]!=0))
  }
  
  return (md)
}

CI <- function(X, Y, sn=20, j0, md, K=5){
  n <- length(Y)
  p <- dim(X)[2]
  
  # initial estimator for beta
  re <- try(SIS(X, Y, family="binomial"), silent = T)
  if (!inherits(re, "try-error")){
    support0 <- rep(FALSE, p)
    support0[re$ix] <- TRUE
    beta0 <- rep(0, p)
    beta0[re$ix] <- re$coef.est[-1]
  }
  else{
    print("use SCAD only")
    re <- cv.ncvreg(X, Y, family="binomial", penalty="SCAD")
    beta0 <- re$fit$beta[-1,re$min]
  }
  
  # weight, pres and hess0
  weight <- rep(0, n)
  pres <- rep(0, n)
  hess0 <- rep(0, n)
  
  se <- rep(0, 2)
  
  for (k in 1:K){
    
    # calculate hess
#    pi <- plogis(X%*%beta0) 
#    hess <- as.vector(pi*(1-pi))
    
    for (i in (sn+1):n){
      # calculate the weight
      support <- rep(FALSE, p)
      support[md[[i-sn]]] <- TRUE
      support[j0] <- TRUE
      X_sub <- as.matrix(X[, support])
	    hess <- as.vector(plogis(X[, support] %*% beta0[support]))
      re <- omega(X_sub, j=(1+support[j0-1]), hess)
      
      if (is.null(re$reg)){
        weight[i-sn] <- X[i,j0]/sqrt(re$sig2)
      }
      else{
        support[j0] <- FALSE
        weight[i-sn] <- (X[i,j0] - X[i,support] %*% re$reg)/sqrt(re$sig2)
      }
    }
  
    # weight for the first few observations
    support <- rep(FALSE, p)
    support[md[[n-sn+1]]] <- TRUE
    support[j0] <- TRUE
    X_sub <- as.matrix(X[, support])
	  hess <- as.vector(plogis(X[, support] %*% beta0[support]))
    re <- omega(X_sub, j=(1+support[j0-1]), hess)
  
    for (i in 1:sn){
      if (is.null(re$reg)){
        weight[i+n-sn] <- X[i,j0]/sqrt(re$sig2)
      }
      else{
        support[j0] <- FALSE
        weight[i+n-sn] <- (X[i,j0] - X[i,support] %*% re$reg)/sqrt(re$sig2)
      }
    }
  
    for (i in (sn+1):n){
      # calculate pres and hess0
      support <- rep(FALSE, p)
      support[md[[i-sn]]] <- TRUE
      support[j0] <- TRUE
      
      pres[i-sn] <- Y[i] - plogis(X[i,support] %*%beta0[support])
      hess0[i-sn] <- plogis(X[i,support] %*%beta0[support])*
        (1-plogis(X[i,support]%*%beta0[support]))
    }
    
    # for the first few observations
    support <- rep(FALSE, p)
    support[md[[n-sn+1]]] <- TRUE
    support[j0] <- TRUE
    
    for (i in 1:sn){
      # refit for the second coefficient
      pres[i+n-sn] <- Y[i] - plogis(X[i,support] %*% beta0[support])
      hess0[i+n-sn] <- plogis(X[i,support] %*%beta0[support])*
        (1-plogis(X[i,support]%*%beta0[support]))
    }
    
    # for the second coefficient
    beta0[j0]<-beta0[j0]+mean((weight*pres))/mean((weight*X[c((sn+1):n,1:sn),j0]*hess0))
    se[1] <- (sqrt(n)/sum(weight*X[c((sn+1):n,1:sn),j0]*hess0))
    se[2] <- (sqrt(n)/sum(weight*X[c((sn+1):n,1:sn),j0]*hess))
  }
  
  return(list(beta=beta0[j0], se=se))
}

##############################################################################
## Generating simulation models
## beta_1=1.5, beta_2=-1.5, beta_j=0 for j>=3
## rho: Cov(X_j, X_i)=rho^{|i-j|}
## n: sample size
## p: dimension 
## sig: Var(Y|X)
##############################################################################
Gen_Model <- function(n=400, p=1000, beta=1.5, rho=0){
  if (rho==0)
    X <- matrix(rnorm(n*p), n, p)
  else{
    X <- matrix(0, n, p)
    for (i in 1:n){
      X[i, ] <- arima.sim(model=list(ar=rho), sd=sqrt(1-rho^2), n=p)
    }
  }
  
  prob <- plogis(beta*X[,1]+beta*X[,2])
  Y <- rbinom(n, 1, prob)
  
  return(list(X=X, Y=Y))
}

################################################################################
## Simulation code for constructing confidence intervals for beta_3
## L is the simulation replication 
################################################################################
# cp initialization
beta.online <- matrix(0, 500, 2)
se.online <- matrix(0, 500, 4)

for (l in 1:25){
  print(l)
  
  set.seed(12345*l)
  Model <- Gen_Model(n=400, p=1000, beta=1)
  X <- Model$X
  Y <- Model$Y
  
  md <- MD(X, Y, sn=floor(3*sqrt(400)))
  
  result <- CI(X, Y, sn=floor(3*sqrt(400)), j0=2, md=md, K=5)
  
  beta.online[l,1] <- result$beta
  se.online[l,c(1,3)] <- result$se
  
  result <- CI(X, Y, sn=floor(3*sqrt(400)), j0=3, md=md, K=5)
  
  beta.online[l,2] <- result$beta
  se.online[l,c(2,4)] <- result$se
  
 save.image(file="v5_online_5steps_n400_p1000_beta1_seed1.RData")
} 
