# simulation setting CE
# using subagging estimation
# n=500

dyn.load("subagging.so")

# generate data
Gen_Data <- function(n=300){
  # generate the covariates
  X <- cbind(rbinom(n, 1, 0.5), runif(n, min=-2, max=2))
  # generate treatment
  piA <- 0.5+0.1*X[,1]
  A <- rbinom(n, 1, piA)
  # generate response
  Y <- X[,2]^2+A*X[,1]*X[,2]^2+0.5*rnorm(n)
  
  return(list(X=X, A=A, Y=Y))
}

# Estimate the response using B-spline
Est.bs <- function(Xtest, X, Y, K=5){
  n <- length(Y)
  n0 <- c(n, length(Xtest))
  if (n<=5){
    Ytest <- rep(0, length(Xtest))
  }
  else if (n<=15){
    X0 <- cbind(X, X^2, X^3)
    beta0 <- as.matrix(lm(Y~X0)$coef)
    Ytest <- cbind(1, Xtest, Xtest^2, Xtest^3)%*%beta0
  }
  else{
    Ytest <- rep(0, max(n0[1], n0[2]))
    # 5-folded cross-validation
    ntrain <- rep(0, K+1)
    ntest <- rep(0, K)
    folds <- split(sample(n, n, replace=FALSE), as.factor(1:K))
    ITrain <- NULL
    ITest <- NULL
    for (i in 1:K){
      ITrain <- c(ITrain, setdiff(1:n, folds[[i]]))
      ITest <- c(ITest, folds[[i]])
      ntrain[i] <- n - length(folds[[i]])
      ntest[i] <- length(folds[[i]])
    }
    ntrain[K+1] <- max(ntrain[1:K])
    ITrain <- ITrain-1
    ITest <- ITest-1
    
    result <- .C("CV_BS", as.double(Ytest), as.double(Xtest), as.double(X), as.double(Y), as.integer(ITrain), 
                 as.integer(ITest), as.integer(ntrain), as.integer(ntest), as.integer(n0), as.integer(5))
    Ytest <- result[[1]][1:n0[2]]
  }
  
  return (Ytest)
}

# Estimate contrast function, propensity score and conditional mean functions
Est.tau <- function(Xtest, X, A, Y){
  
  m <- dim(Xtest)[1]
  pi.est <- rep(0, m)
  h0.est <- rep(0, m)
  h1.est <- rep(0, m)
  tau.est <- rep(0, m)
  
  pi.est[Xtest[,1]==0] <- Est.bs(Xtest[Xtest[,1]==0,2], X[X[,1]==0,2], A[X[,1]==0])
  pi.est[Xtest[,1]==1] <- Est.bs(Xtest[Xtest[,1]==1,2], X[X[,1]==1,2], A[X[,1]==1])
  pi.est[pi.est<=0.05] <- 0.05
  pi.est[pi.est>=0.95] <- 0.95

  h0.est[Xtest[,1]==0] <- Est.bs(Xtest[Xtest[,1]==0,2], X[(X[,1]==0)&(A==0),2], Y[(X[,1]==0)&(A==0)])
  h0.est[Xtest[,1]==1] <- Est.bs(Xtest[Xtest[,1]==1,2], X[(X[,1]==1)&(A==0),2], Y[(X[,1]==1)&(A==0)])
  
  h1.est[Xtest[,1]==0] <- Est.bs(Xtest[Xtest[,1]==0,2], X[(X[,1]==0)&(A==1),2], Y[(X[,1]==0)&(A==1)])
  h1.est[Xtest[,1]==1] <- Est.bs(Xtest[Xtest[,1]==1,2], X[(X[,1]==1)&(A==1),2], Y[(X[,1]==1)&(A==1)])
  
  tau.est <- h1.est-h0.est
  d.est <- (tau.est>0)+0
  
  return(list(pi.est=pi.est, h0.est=h0.est, h1.est=h1.est, d.est=d.est))
}

# Value function estimation
V.est <- function(X, A, Y, est){
  
  # all estimators
  pi.est <- est$pi.est
  h0.est <- est$h0.est
  h1.est <- est$h1.est
  d.est <- est$d.est
  
  # AIPWE
  AIPWE <- (d.est*A/pi.est+(1-d.est)*(1-A)/(1-pi.est))*Y
  AIPWE <- AIPWE-(d.est*A/pi.est+(1-d.est)*(1-A)/(1-pi.est)-1)*(h1.est*d.est+h0.est*(1-d.est))
}

### simulation
L <- 1000
v0.all <- rep(0, L)
sd.all <- rep(0, L)
B <- 4000
n <- 500
sn <- floor(4*n/log(n))
num <- rep(0, 4)


	for (l in 1:10){
	  print(l)
	  # generate model
	  set.seed(1234567*l)
	  Md <- Gen_Data(n=n)
	  X <- Md$X
	  A <- Md$A
	  Y <- Md$Y
	  # subsampling
	  iter <- 0
	  v0 <- rep(0, n)
	  num0 <- rep(0, n)
	  while (iter<B){

		Ind0 <- sample(n, sn, replace=FALSE)
		num[1] <- sum((A[Ind0]==0)&(X[Ind0,1]==0))
		num[2] <- sum((A[Ind0]==0)&(X[Ind0,1]==1))
		num[3] <- sum((A[Ind0]==1)&(X[Ind0,1]==0))
		num[4] <- sum((A[Ind0]==1)&(X[Ind0,1]==1))
		if (any(num<=10)){
		  next
		}
		else{
		  # estimate the treatment regime
		  d.est <- Est.tau(X, X[Ind0,], A[Ind0], Y[Ind0])$d.est
		  
		  # randomly split the samples
		  Ind0c <- setdiff(1:n, Ind0)
		  Ind1 <- sample(Ind0c, floor((n-sn)/2), replace=FALSE)
		  Ind2 <- setdiff(Ind0c, Ind1)
		  num0[Ind1] <- num0[Ind1]+1
		  num0[Ind2] <- num0[Ind2]+1
		  
		  # estimate the propensity score and conditional mean functions
		  Est1 <- Est.tau(X[Ind2,], X[union(Ind0, Ind1),], A[union(Ind0, Ind1)], Y[union(Ind0, Ind1)])
		  Est1$d.est <- d.est[Ind2]
		  Est2 <- Est.tau(X[Ind1,], X[union(Ind0, Ind2),], A[union(Ind0, Ind2)], Y[union(Ind0, Ind2)])
		  Est2$d.est <- d.est[Ind1]
		  
		  # AIPWE
		  v0[Ind1] <- v0[Ind1]+V.est(X[Ind1,], A[Ind1], Y[Ind1], Est2)
		  v0[Ind2] <- v0[Ind2]+V.est(X[Ind2,], A[Ind2], Y[Ind2], Est1)
		  
		  # the next iteration
		  iter <- iter+1
		}
	  }
	  v0.all[l] <- (sum(v0))/(sum(num0))
	  sd.all[l] <- sqrt(var((v0/num0)))/sqrt(n)
	  
		save.image(file="sb_CE_I_n500_sn4_seed1.RData") 
	}
