## Simulation code for high-dimensional linear hypothesis testing project 
## Chengchun Shi, Jun. 19, 2017
####################################################################################
####################################################################################
## Estimate parameters using ADMM for high-dimensional linear regression 
## with the constraints: beta_N=0
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta)+rho v^T (beta-theta)+0.5rho ||beta-theta||_2^2
## subject to beta_N=0
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v: the Langrange multiplier
######################################################

# ADMM step, update beta
# beta=(rho I+ X^T X/n )^{-1} (X^T Y/n + rho(theta-v))
update_beta <- function(X, Y, theta, v, rho=1){
  n <- length(Y)
  p <- dim(X)[2]
  beta <- crossprod(X, Y)/n+rho*(theta-v)
  
  if (p<n){
    beta <- solve( crossprod(X,X)/n+rho*diag(p), beta )
  }
  else{
    # using woodbury formula when p is large
    Xbeta <- solve( X%*%t(X)+rho*n*diag(n), X%*%beta )
    Xbeta <- crossprod(X, Xbeta)
    beta <- (beta-Xbeta)/rho
  }
  
  return (beta)
}

# ADMM step: update theta
# SCAD thresholding with rho=1
update_theta <- function(beta, v, lambda, rho=1, a=3.7){
  theta <- v+beta
  # when absolute value smaller than lambda
  theta[abs(theta)<=lambda] <- 0
  # when larger than lambda but smaller than 2lambda
  index <- abs(theta)>lambda & abs(theta)<=2*lambda
  theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
  # when larger than 2lambda but smaller than a lambda
  index <- abs(theta)>2*lambda & abs(theta)<=a*lambda
  theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2) 
  
  return (theta)
}

# calculate the primal and dual residuals
res <- function(beta_new, theta_new, beta_old, theta_old, rho=1){
  # the primal residual
  r <- beta_new-theta_new
  # the dual residual 
  s <- rho*(theta_new-theta_old)
  
  return (list(r=r, s=s))
}

# ADMM algorithm
SCAD_ADMM_re <- function(X, Y, N, beta0, lambda, rho=1, a=3.7, err=0.5e-3){
  
  # initial value of parameters
  beta_old <- beta0[!N]
  theta_old <- beta0[!N]
  v <- 0
  
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err || norm(as.matrix(eps$s), "F")>=err)&&iter<=300)
  {
    # ADMM update
    # update v
    v <- v+(beta_old-theta_old)
    
    # update beta
    beta_new <- update_beta(X=X[,!N], Y=Y, theta=theta_old, v=v, rho=rho)
    
    # update theta
    theta_new <- update_theta(beta=beta_new, v=v, lambda=lambda, rho=rho, a=a)
    
    # calculate the residuals
    eps <- res(beta_new, theta_new, beta_old, theta_old, rho=rho)
    
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter+1
  }
  
  thetaf <- rep(0, dim(X)[2])
  thetaf[!N] <- theta_new
  
  return(thetaf)
}

##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM_re <- function(X, Y, N, beta0, K=5, rho=1, a=3.7, err=0.5e-3, tune="cv"){
  
  # potential tuning parameters
  lambdalist <- exp(seq(-2.5, 1.0, 0.07))
  lambdan <- length(lambdalist)
  
  # data splitting
  n <- length(Y)
  p <- dim(X)[2]
  
  if (tune=="cv"){
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    # calculate MSE for each folds
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      X0 <- X[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      
      # training dataset
      X1 <- X[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM_re(X=X1, Y=Y1, N=N, beta0=beta0, lambda=lambdalist[j], 
                              rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum((Y0-X0%*%beta0)^2)
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
    beta <- SCAD_ADMM_re(X=X, Y=Y, N=N, beta0=beta0, lambda=lambda, rho=rho, a=a, 
                         err=err)
  }
  else{
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta0 <- SCAD_ADMM_re(X=X, Y=Y, N=N, beta0=beta0, lambda=lambdalist[j], 
                            rho=rho, a=a, err=err)
      BIC[j] <- sum((Y-X%*%beta0)^2)+2*sum(beta0!=0)*max(log(n),log(log(n))*log(p))
    }
    
    lambda <- lambdalist[which.min(BIC)]
    
    beta <- SCAD_ADMM_re(X=X, Y=Y, N=N, beta0=beta0, lambda=lambda, rho=rho, a=a, 
                         err=err)
  }
  
  return(beta)
}

####################################################################################
## Estimate parameters using ADMM for high-dimensional linear regression 
## without the constraints: beta_N=0
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

# ADMM step: update beta
# beta=(rho I+ X^T X/n)^{-1} (X^T Y/n+ rho(theta-v))
update_beta_unre <- function(X, Y, theta, v, rho=1){
  n <- length(Y)
  p <- dim(X)[2]
  beta <- crossprod(X, Y)/n+rho*(theta-v)
  
  if (p<n){
    beta <- solve( crossprod(X,X)/n+rho*diag(p), beta )
  }
  else{
    # using woodbury formula when p is large
    Xbeta <- solve( X%*%t(X)+rho*n*diag(n), X%*%beta )
    Xbeta <- crossprod(X, Xbeta)
    beta <- (beta-Xbeta)/rho
  }
  
  return (beta)
}

# ADMM step: update theta
# SCAD thresholding
update_theta_unre <- function(beta, N, v, lambda, rho=1, a=3.7){
  theta <- v+beta
  # when absolute value smaller than lambda
  theta[(!N) & abs(theta)<=lambda] <- 0
  # when larger than lambda but smaller than 2lambda
  index <- (!N) & abs(theta)>lambda & abs(theta)<=2*lambda
  theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
  # when larger than 2lambda but smaller than a lambda
  index <- (!N) & abs(theta)>2*lambda & abs(theta)<=a*lambda
  theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2) 
  
  return (theta)
}

# calculate the primal and dual residuals
res_unre <- function(beta_new, theta_new, beta_old, theta_old, rho=1){
  # the primal residual
  r <- beta_new-theta_new
  # the dual residual 
  s <- rho*(theta_new-theta_old)
  
  return (list(r=r, s=s))
}

# ADMM algorithm
SCAD_ADMM_unre <- function(X, Y, N, beta0, lambda, rho=1, a=3.7, err=0.5e-3){
  
  # initial value of parameters
  beta_old <- beta0
  theta_old <- beta0
  v <- 0
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err || norm(as.matrix(eps$s), "F")>=err)&&iter<=300)
  {
    # ADMM update
    # update v
    v <- v+(beta_old-theta_old)
    
    # update beta
    beta_new <- update_beta_unre(X=X, Y=Y, theta=theta_old, v=v, rho=rho)
    
    # update theta
    theta_new <- update_theta_unre(beta=beta_new, N=N, v=v, lambda=lambda, rho=rho, a=a)
    
    # calculate the residuals
    eps <- res_unre(beta_new, theta_new, beta_old, theta_old, rho=rho)
    
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter+1
  }
  
  return(theta_new)
}

##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM_unre <- function(X, Y, N, beta0, K=5, rho=1, a=3.7, err=0.5e-3, tune="cv"){
  
  # potential tuning parameters
  lambdalist <- exp(seq(-2.5, 1.0, 0.07))
  lambdan <- length(lambdalist)
  
  # data splitting
  n <- length(Y)
  p <- dim(X)[2]
  
  if (tune=="cv"){
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    # calculate MSE for each folds
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      X0 <- X[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      
      # training dataset
      X1 <- X[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM_unre(X=X1, Y=Y1, N=N, beta0=beta0, lambda=lambdalist[j], 
                                rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum((Y0-X0%*%beta0)^2)
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
    beta <- SCAD_ADMM_unre(X=X, Y=Y, N=N, beta0=beta0, lambda=lambda, rho=rho, a=a, err=err)
  }
  else{
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta0 <- SCAD_ADMM_unre(X=X, Y=Y, N=N, beta0=beta0, lambda=lambdalist[j], 
                              rho=rho, a=a, err=err)
      BIC[j] <- sum((Y-X%*%beta0)^2)+2*sum(beta0!=0)*max(log(n),log(log(n))*log(p))
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(BIC)]
    beta <- SCAD_ADMM_unre(X=X, Y=Y, N=N, beta0=beta0, lambda=lambda, rho=rho, a=a, err=err)
  }
  
  return(beta)
}

##############################################################################
## Generating simulation models
## beta_1=1,beta_2=-1,beta_3=beta_4=h, beta_j=0 for j>=5
## rho: Cov(X_j, X_i)=rho^{|i-j|}
## n: sample size
## p: dimension 
## sig: Var(Y|X)
##############################################################################
Gen_Model <- function(n=100, p=10, h=0, sig=1, beta=2, rho=0){
  if (rho==0)
    X <- matrix(rnorm(n*p), n, p)
  else{
    X <- matrix(0, n, p)
    for (i in 1:n){
      X[i, ] <- arima.sim(model=list(ar=rho), n=p)
    }
  }
  
  Y <- beta*X[,1]-beta*X[,2]+h*(X[,3]+X[,4])+sig*rnorm(n)
  
  return(list(X=X, Y=Y))
}

################################################################################
## Simulation code for testing the null hypothesis: beta_N=0
## constructing the partial penalized Wald, score and likelihood ratio statistics
## L is the simulation replication 
################################################################################
simu.ex2 <- function(L=200, n=100, p=50, h=0, beta=2, sig=1, rho=0){
  
  # p-value initialization
  pv <- rep(0, 6)
  sig.all <- rep(0, L)
  
  # location of zero components
  N <- c(FALSE, FALSE, TRUE, TRUE, rep(FALSE, p-4))
  
  Tall <- matrix(0, L, 6)
  beta.al <- matrix(0, L, 2)
  
  for (l in 1:L){
    # each iteration
    Model <- Gen_Model(n=n, p=p, h=h, beta=beta, sig=sig, rho=rho)
    
    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC")
    indice.unre <- beta.unre!=0
    
    # estimate the constrained estimator
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC")
    indice.re <- beta.re!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)*n/(n-sum(beta.unre!=0))
    
    # construct the likelihood ratio statistic
    TL <- sum((Model$X%*%beta.re-Model$Y)^2)-sum((Model$X%*%beta.unre-Model$Y)^2)
    
    # construct the Wald statistic
    B_0 <- solve(crossprod(Model$X[,indice.unre|N], Model$X[,indice.unre|N]))
    d_0 <- sum(which(indice.unre|N)<=3)
    TW <- crossprod(beta.unre[N], solve(B_0[d_0:(d_0+1),d_0:(d_0+1)], beta.unre[N])) 
    
    # construct the score statistic
    eps <- Model$Y-Model$X%*%beta.re
    Xeps <- crossprod(Model$X[,indice.re|N], eps)
    TS <- crossprod(Xeps, solve(crossprod(Model$X[,indice.re|N], Model$X[,indice.re|N]), Xeps))
    
    if (TL>=sig2*qchisq(0.95, df=2))
      pv[1] <- pv[1]+1/L
    if (TW>=sig2*qchisq(0.95, df=2))
      pv[2] <- pv[2]+1/L
    if (TS>=sig2*qchisq(0.95, df=2))
      pv[3] <- pv[3]+1/L
    if (TL>=qchisq(0.95, df=2))
      pv[4] <- pv[4]+1/L
    if (TW>=qchisq(0.95, df=2))
      pv[5] <- pv[5]+1/L
    if (TS>=qchisq(0.95, df=2))
      pv[6] <- pv[6]+1/L
    sig.all[l] <- sig2
    
    Tall[l, ] <- c(TL/sig2, TW/sig2, TS/sig2, TL, TW, TS)
    beta.al[l, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
  }
  
  return(list(pv=pv, sig2=sig.all, TS=Tall, beta=beta.al))
}

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
## Estimate parameters using ADMM for high-dimensional linear regression 
## without the linear constraints: C beta=b
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta)+rho v^T(beta-theta)+0.5rho ||beta-theta||_2^2
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v: the Langrange multiplier
######################################################
# ADMM algorithm
SCAD_ADMM <- function(X, Y, beta0, lambda, rho=1, a=3.7, err=0.5e-3){
  
  # initial value of parameters
  beta_old <- beta0
  theta_old <- beta0
  v <- 0
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err || norm(as.matrix(eps$s), "F")>=err)&&iter<=300)
  {
    # ADMM update
    # update v
    v <- v+(beta_old-theta_old)
    
    # update beta
    beta_new <- update_beta(X=X, Y=Y, theta=theta_old, v=v, rho=rho)
    
    # update theta
    theta_new <- update_theta(beta=beta_new, v=v, lambda=lambda, rho=rho, a=a)
    
    # calculate the residuals
    eps <- res(beta_new, theta_new, beta_old, theta_old, rho=rho)
    
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter +1
  }
  
  return(theta_new)
}

##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM <- function(X, Y, beta0, K=5, rho=1, a=3.7, err=0.5e-3, tune="cv"){
  
  # potential tuning parameters
  lambdalist <- exp(seq(-2.5, 0.5, 0.05))
  lambdan <- length(lambdalist)
  n <- length(Y)
  p <- dim(X)[2]
  
  if (tune=="cv"){
    # data splitting
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    # calculate MSE for each folds
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      X0 <- X[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      
      # training dataset
      X1 <- X[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM(X=X1, Y=Y1, beta0=beta0, lambda=lambdalist[j], 
                           rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum((Y0-X0%*%beta0)^2)
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
  }
  
  if (tune=="BIC"){
    # calculate BIC for each tuning parameter
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta0 <- SCAD_ADMM(X=X, Y=Y, beta0=beta0, lambda=lambdalist[j], 
                         rho=rho, a=a, err=err)
      BIC[j] <- BIC[j]+sum((Y-X%*%beta0)^2)+2*sum(beta0!=0)*max(log(n),log(log(n))*log(p))
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(BIC)]
  }
  
  beta <- SCAD_ADMM(X=X, Y=Y, beta0=beta0, lambda=lambda, rho=rho, a=a, err=err)
  
  return(list(beta=beta, lambda=lambda))
}

####################################################################################
## Estimate parameters using ADMM for high-dimensional linear regression 
## with the linear constraints: C beta=b
####################################################################################
## Objective function:
## f(Y,X,beta)+p_a(lambda,theta_N^c)+rho v_1^T (beta-theta)+rho v_2^T (Cbeta-b)+0.5rho ||beta-theta||_2^2
## subject to C beta=b
## f: the loglikelihood
## p: The SCAD penalty function 
## a: parameter in the SCAD penalty, default=3.7
## Y: the response vector
## X: the design matrix
## lambda: the regularization parameter
## rho: regularization parameter in ADMM, we set rho=1
## v_1,v_2: the Langrange multiplier
## N={1,2}
#####################################################################################
# ADMM step: update beta
# beta=(rho I+ X^T X/n+ rho C^T C)^{-1} (X^T Y/n+ rho(theta-v1)+ rho C^T (b-v2))
update_beta_con <- function(X, Y, theta, C, b, v1, v2, rho=1){
  n <- length(Y)
  p <- dim(X)[2]
  d <- dim(C)[1]
  beta <- crossprod(X, Y)/n+rho*(theta-v1)+rho*as.vector(t(C)%*%(b-v2))
  
  if (p<n){
    beta <- solve( crossprod(X,X)/n+rho*diag(p)+rho*crossprod(C,C), beta )
  }
  else{
    # using woodbury formula when p is large
    Xbeta <- solve( rbind(X,C)%*%t(rbind(X/n, rho*C))+rho*diag(n+d), rbind(X,C)%*%beta )
    Xbeta <- crossprod(rbind(X/n, rho*C), Xbeta)
    beta <- (beta-Xbeta)/rho
  }
  
  return (beta)
}

# ADMM step: update theta
# SCAD thresholding
update_theta_con <- function(beta, v1, lambda, rho=1, a=3.7){
  theta <- v1+beta
  N <- c(TRUE, TRUE, rep(FALSE, length(theta)-2))
  # when absolute value smaller than lambda
  theta[(!N) & abs(theta)<=lambda] <- 0
  # when larger than lambda but smaller than 2lambda
  index <- (!N) & abs(theta)>lambda & abs(theta)<=2*lambda
  theta[index] <- sign(theta[index])*(abs(theta[index])-lambda)
  # when larger than 2lambda but smaller than a lambda
  index <- (!N) & abs(theta)>2*lambda & abs(theta)<=a*lambda
  theta[index] <- ((a-1)*theta[index]-sign(theta[index])*a*lambda)/(a-2) 
  
  return (theta)
}

# calculate the primal and dual residuals
res_con <- function(beta_new, theta_new, beta_old, theta_old, C, b, rho=1){
  # the primal residual
  r <- c(beta_new-theta_new, C%*%beta_new-b)
  # the dual residual 
  s <- rho*(theta_new-theta_old)
  
  return (list(r=r, s=s))
}

# ADMM algorithm
SCAD_ADMM_con <- function(X, Y, C, b, beta0, lambda, rho=1, a=3.7, err=0.5e-3){
  
  # initial value of parameters
  beta_old <- beta0
  theta_old <- beta0
  v1 <- 0
  v2 <- 0
  
  iter <- 0
  
  eps <- list(r=1, s=1)
  while ((norm(as.matrix(eps$r), "F")>=err || norm(as.matrix(eps$s), "F")>=err)&&iter<=300)
  {
    # ADMM update
    # update v
    v1 <- v1+(beta_old-theta_old)
    v2 <- v2+(C%*%beta_old-b)
    
    # update beta
    beta_new <- update_beta_con(X=X, Y=Y, theta=theta_old, C=C, b=b, v1=v1, v2=v2, rho=rho)
    
    # update theta
    theta_new <- update_theta_con(beta=beta_new, v1=v1, lambda=lambda, rho=rho, a=a)
    
    # calculate the residuals
    eps <- res_con(beta_new, theta_new, beta_old, theta_old, C=C, b=b, rho=rho)
    
    # beta_old and theta_old
    beta_old <- beta_new
    theta_old <- theta_new
    
    iter <- iter+1
  }
  
  return(theta_new)
}

##############################################################################
## K-folded cross validation, K=5 by default
##############################################################################
cv.SCAD_ADMM_con <- function(X, Y, C, b, beta0, K=5, rho=1, a=3.7, err=0.5e-3, tune="cv"){
  
  # potential tuning parameters
  lambdalist <- exp(seq(-2.5, 1.0, 0.07))
  lambdan <- length(lambdalist)
  
  # data splitting
  n <- length(Y)
  p <- dim(X)[2]
  
  # calculate MSE for each folds
  if (tune=="cv"){
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    
    MSE <- rep(0, lambdan)
    for (k in 1:K){
      # testing dataset
      X0 <- X[folds[[k]], ]
      Y0 <- Y[folds[[k]]]
      
      # training dataset
      X1 <- X[-folds[[k]], ]
      Y1 <- Y[-folds[[k]]]
      
      # est the training dataset
      for (j in 1:lambdan){
        beta0 <- SCAD_ADMM_con(X=X1, Y=Y1, C=C, b=b, beta0=beta0, lambda=lambdalist[j], 
                               rho=rho, a=a, err=err)
        MSE[j] <- MSE[j]+sum((Y0-X0%*%beta0)^2)
      }
    }
    
    # take the minimum lambda
    lambda <- lambdalist[which.min(MSE)]
    beta <- SCAD_ADMM_con(X=X, Y=Y, C=C, b=b, beta0=beta0, lambda=lambda, rho=rho, a=a, 
                          err=err)
  }
  else{
    BIC <- rep(0, lambdan)
    for (j in 1:lambdan){
      beta0 <- SCAD_ADMM_con(X=X, Y=Y, C=C, b=b, beta0=beta0, lambda=lambdalist[j], 
                             rho=rho, a=a, err=err)
      BIC[j] <- BIC[j]+sum((Y-X%*%beta0)^2)+2*sum(beta0!=0)*max(log(n),log(log(n))*log(p))
    }
    lambda <- lambdalist[which.min(BIC)]
    beta <- SCAD_ADMM_con(X=X, Y=Y, C=C, b=b, beta0=beta0, lambda=lambda, rho=rho, a=a, 
                          err=err)
  }
  
  return(beta)
}

##############################################################################
## Generating simulation models
## beta_1=1,beta_2=-1+h, beta_j=0 for j>=3
## rho: Cov(X_j, X_i)=rho^{|i-j|}
## n: sample size
## p: dimension 
## sig: Var(Y|X)
##############################################################################
Gen_Model <- function(n=100, p=10, h=0, sig=1, beta=2.0, rho=0){
  if (rho==0)
    X <- matrix(rnorm(n*p), n, p)
  else{
    X <- matrix(0, n, p)
    for (i in 1:n){
      X[i, ] <- arima.sim(model=list(ar=rho), n=p)
    }
  }
  
  Y <- beta*X[,1]-(beta+h)*X[,2]+sig*rnorm(n)
  
  return(list(X=X, Y=Y))
}

################################################################################
## Simulation code for testing the null hypothesis: C beta = b
## constructing the partial penalized Wald, score and likelihood ratio statistics
## L is the simulation replication 
################################################################################
simu.ex2 <- function(L=200, n=100, p=50, h=0, beta=1, sig=1, rho=0){
  
  # p-value initialization
  pv <- rep(0, 6)
  sig.all <- rep(0, L)
  Tall <- matrix(0, L, 6)
  beta.al <- matrix(0, L, 2)
  
  C <- t(as.matrix(c(1,1, rep(0,p-2))))
  
  for (l in 1:L){
    # each iteration
    Model <- Gen_Model(n=n, p=p, h=h, beta=beta, sig=sig, rho=rho)
    
    # estimate the uncontrained estimator
    beta.uncon <- cv.SCAD_ADMM(X=Model$X, Y=Model$Y, beta0=rep(0, dim(Model$X)[2]), err=1e-4,
                               tune="BIC")
    lambda.uncon <- beta.uncon$lambda
    beta.uncon <- beta.uncon$beta
    indice.uncon <- beta.uncon!=0
    
    # estimate of the variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.uncon)^2)
    
    # estimate the constrained estimator
    beta.con <- cv.SCAD_ADMM_con(X=Model$X, Y=Model$Y, C=t(as.matrix(c(1,1, rep(0,dim(Model$X)[2]-2)))), 
                                 b=0, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC")
    indice.con <- beta.con!=0
    
    # estimate of the variance
    sig2 <- mean((Model$Y-Model$X%*%beta.uncon)^2)*n/(n-sum(beta.uncon!=0))
    
    # construct the likelihood ratio statistic
    TL <- sum((Model$X%*%beta.con-Model$Y)^2)-sum((Model$X%*%beta.uncon-Model$Y)^2)
    
    # construct the Wald statistic
    TW <- (C%*%beta.uncon)^2 / (C[,indice.uncon]%*% solve(crossprod(Model$X[,indice.uncon], 
                                                                    Model$X[,indice.uncon]), as.matrix(C[,indice.uncon])))
    
    # construct the score statistic
    eps <- Model$Y-Model$X%*%beta.con
    Xeps <- crossprod(Model$X[,indice.con], eps)
    TS <- crossprod(Xeps, solve(crossprod(Model$X[,indice.con], Model$X[,indice.con]), Xeps))
    
    if (TL>=sig2*qchisq(0.95, df=1))
      pv[1] <- pv[1]+1/L
    if (TW>=sig2*qchisq(0.95, df=1))
      pv[2] <- pv[2]+1/L
    if (TS>=sig2*qchisq(0.95, df=1))
      pv[3] <- pv[3]+1/L
    if (TL>=qchisq(0.95, df=1))
      pv[4] <- pv[4]+1/L
    if (TW>=qchisq(0.95, df=1))
      pv[5] <- pv[5]+1/L
    if (TS>=qchisq(0.95, df=1))
      pv[6] <- pv[6]+1/L
    
    sig.all[l] <- sig2
    Tall[l, ] <- c(TL/sig2, TW/sig2, TS/sig2, TL, TW, TS)
    beta.al[l, ] <- c(sum(beta.con!=0), sum(beta.uncon!=0))
  }
  
  return(list(pv=pv, sig2=sig.all, TS=Tall, beta=beta.al))
}