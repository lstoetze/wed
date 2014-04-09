#####################
# wed.R
#####################

# This provides functions to estimate the nonseperable preference model
# as outlined in the article for the ESS

# There are two seperate Liklihood functions
# one for the standard model and one for the nonsep-model

#  A function estM opitimizes the two Liklihood functions
# returning ml estimates of the parameters; se and so on.



#######################
# lIKLIHOOD FUNCTIONS

lik.mnl.normal <- function(theta,y,S,P,X){
  
  require(dummies)
  
  N <- length(y)
  # Vote Matrix
  V <- dummy(y) == 1
  # Number of parties
  J <- ncol(V)
  #
  num.cov <- ncol(X)
  
  # Weigth Matrix
  L <- matrix(0, ncol=2,nrow=2)
  L[1,1] <- theta[1]
  L[2,2] <- theta[2]
  L[2,1] <- 0
  
  # Transformation
  A <- L %*% t(L)
  
  
  # Betas
  beta <- matrix(0,ncol=num.cov+1,nrow=J)
  beta[2:J,] <- theta[3:(2+(num.cov+1)*(J-1))]
  
  
  
  # Utility
  U <- matrix(NA,ncol=J,nrow=N)
  U[,1] <- 0
  
  for (j in 2:J){
    U[,j] <- beta[j,1] - sqrt(  A[1,1]* (P[(j-1),1] - S[,1])^2     + 2*A[1,2]*(P[(j-1),1] - S[,1])*(P[(j-1),2] - S[,2])   + A[2,2] * (P[(j-1),2] - S[,2])^2) +  X %*% beta[j,-1] 
  }
  
  # Probabilities
  P <- matrix(NA,ncol=J,nrow=N)
  sU <- apply(exp(U),1,sum)
  for(j in 1:J){
    P[,j] <- exp(U[,j])/sU
  }
  
  # liklihood
  loglik <- sum(log(P[V]))
  
  return(loglik)
  
}

lik.mnl.nonsep <- function(theta,y,S,P,X){
  
  require(dummies)
  
  N <- length(y)
  # Vote Matrix
  V <- dummy(y) == 1
  # Number of parties
  J <- ncol(V)
  #
  num.cov <- ncol(X)
  
  # Weigth Matrix
  L <- matrix(0, ncol=2,nrow=2)
  L[1,1] <- theta[1]
  L[2,2] <- theta[2]
  L[2,1] <- theta[3]
  
  # Transformation
  A <- L %*% t(L)
  
  
  # Betas
  beta <- matrix(0,ncol=num.cov+1,nrow=J)
  beta[2:J,] <- theta[4:(3+(num.cov+1)*(J-1))]
  
  
  
  # Utility
  U <- matrix(NA,ncol=J,nrow=N)
  U[,1] <- 0
  
  for (j in 2:J){
    U[,j] <- beta[j,1] - sqrt(  A[1,1]* (P[(j-1),1] - S[,1])^2     + 2*A[1,2]*(P[(j-1),1] - S[,1])*(P[(j-1),2] - S[,2])   + A[2,2] * (P[(j-1),2] - S[,2])^2) +  X %*% beta[j,-1] 
  }
  
  # Probabilities
  P <- matrix(NA,ncol=J,nrow=N)
  sU <- apply(exp(U),1,sum)
  for(j in 1:J){
    P[,j] <- exp(U[,j])/sU
  }
  
  # liklihood
  loglik <- sum(log(P[V]))
  
  return(loglik)
  
}




##################
# Likelihood Maximization routine

wed <- function(list_data,  sep=FALSE ,rep = 2,se.A=FALSE) {
  
  if(sep==TRUE){
  lik.func <- lik.mnl.nonsep
  } else {
  lik.func <-   lik.mnl.normal
  }
  # Iterate over maximization with different starting values

  # Gte number of paramters
  k <- ncol(list_data$S) # number of dimensions
  nX <- ncol(list_data$X)+1 # Number of Covariates
  J <- length(unique(list_data$Rat))
  off.diag <- sum(lower.tri(diag(rep(0,k))))
    
  num.par <- ifelse(sep==TRUE
                     ,off.diag+k+nX*(J-1)
                     ,k+nX*(J-1))
  
  cat("Maximize Likelihood \n\n")
  
  res1 <- optim(rep(0,num.par), lik.func
                ,y=list_data$Rat,S=list_data$S,P=list_data$P
                ,X=as.matrix(list_data$X)
                ,control=list(fnscale=-1,trace=10)
                ,method="BFGS"
                ,hessian=TRUE)
  
  
  
  cat(paste("\n\n Check Results using", rep-1, "different starting values\n\n"))
  
  t <- 1
  while(t < rep){
    cat(paste("########\n",t, "ROUND\n\n"))
    res <- optim(rnorm( num.par ), lik.func
                 ,y=list_data$Rat,S=list_data$S,P=list_data$P
                 ,X=as.matrix(list_data$X)
                 ,control=list(fnscale=-1,trace=10)
                 ,method="BFGS"
                 ,hessian=TRUE)
    
    
    if (res$value > res1$value) {
      cat("\n RESULT: Changed value \n\n") ; res1 <- res
    } else {
      cat("\n RESULT: No change \n\n")
    }
    t <- t +1
  }
  
  
  
  # Get Estimates
  par <- res1$par
  
  # Get Se
  se <- sqrt(diag(solve(-res1$hessian)))
  
  # For A
  L <- matrix(0,2,2)
  L[1,1] <- par[1]
  L[2,2] <- par[2]
  L[2,1] <- ifelse(sep==TRUE,par[3],0)
  
  A <- L %*% t(L)
  
  # Se for A
  if(se.A==TRUE){
    require("MASS")
    V <- solve(-res1$hessian)
    sp <- ifelse(nonsep==TRUE,3,2)
    S <- mvrnorm(1000,par,V)
    S <- S[,1:sp]
    
    se.S <- function(S,nonsep){
      T <- matrix(0,2,2)
      T[1,1] <- S[1]
      T[2,2] <- S[2]
      T[2,1] <- ifelse(nonsep==TRUE,S[3],0)
      res <- T %*% t(T)
      return(c(res))
    }
    As <- t(apply(S,1,se.S,nonsep))
    nonsep <- As[,2]/sqrt(As[,1]*As[,4])
    se.A.res <- apply(As,2,quantile,c(0.025,0.5,0.975))
    d.nonsep <- quantile(nonsep,c(0.025,0.5,0.975))
    
  } else{
    
    se.A.res  <- NA
    d.nonsep <- NA
  }
  
  
  
  return(list("est"=par
              , "SE"=se
              ,"A"=A
              ,"se.A" = se.A.res
              ,"d.nonsep" = d.nonsep
              , "lik"=res1$value))
}



