library(MASS)

###-------- Functions for sampling from the generalized Lehmann model with Gumbel's bivariate exponential distribution as our baseline --------
S1 <- function(t1,a1,alpha,Z){
  return(exp(-a1*t1*exp(alpha*Z)))
}


S2 <- function(t2,a2,beta,Z){
  return(exp(-a2*t2*exp(beta*Z)))
}


exp_A <- function(t1,t2,a1,a2,alpha,beta, gamma,Z){
  return(exp( -a1*a2*t1*t2*exp(gamma*Z)))
}


S <- function(t1,t2,a1,a2,alpha,beta, gamma,Z){ 
  return(S1(t1,a1,alpha,Z)*S2(t2,a2,beta,Z)*(exp_A(t1,t2,a1,a2,alpha,beta, gamma,Z)))
}


F1_inv <- function(u,a1,alpha,Z){
  return(-log(1-u)*exp(-alpha*Z)/a1)
}

F21 <- function(t1,t2,a1,a2,alpha,beta, gamma,Z){
  S21 <- S2(t2,a2,beta,Z)*exp_A(t1,t2,a1,a2,alpha,beta, gamma,Z)*(1+a2*t2*exp((gamma-alpha)*Z))
  return(1-S21)
}

f <- function(x,v,t1i,a1,a2,alpha,beta, gamma,Zi){
  return(F21(t2=x,t1=t1i,a1,a2,alpha,beta, gamma,Zi)-v)
}

data_gen_NOD_no_cens <- function(a1=1, a2=1, alpha, beta, gamma, Z, setseed=FALSE, seed=123){
  Failure_times <- data.frame(T1=NA,T2=NA,covar=NA)
  n <- length(Z)
  if (isTRUE(setseed)){
    set.seed(seed)
  }
  u <- runif(n)
  w <- runif(n)
  for (j in 1:n){
    Zi <- Z[j]
    t1i <- F1_inv(u[j],a1,alpha,Zi)
    root <- uniroot(f, c(0,10^2),v=w[j],t1=t1i,a1=a1,a2=a2,alpha=alpha,beta=beta, gamma=gamma,Zi=Zi)
    t2i <- root$root
    Failure_times[j,] <- c(t1i,t2i,Zi)
  }
  return(Failure_times)
}

###-------- Functions for sampling from the generalized Lehmann model with Frank's copula as our baseline --------
exp_A_Frank <- function(t1,t2,a1,a2, gamma,Z, Frank_param){
  exp_marg_1 <- exp(-a1*t1) #u=S1(t1,a1,alpha=0,Z=0)
  exp_marg_2 <- exp(-a2*t2) #v=S2(t2,a2,alpha=0,Z=0)
  log_term <- log(1+(exp(-Frank_param*exp_marg_1)-1)*(exp(-Frank_param*exp_marg_2)-1)/(exp(-Frank_param)-1))
  exp_A <- (-((Frank_param*exp_marg_1*exp_marg_2)^{-1})*log_term)^exp(gamma*Z)
  return(exp_A)
}

S_Frank <-  function(t1,t2,a1,a2, alpha,beta, gamma,Z, Frank_param){
  return(S1(t1,a1,alpha,Z)*S2(t2,a2,beta,Z)*exp_A_Frank(t1,t2,a1,a2, gamma,Z, Frank_param))
}

Frank_cop <- function(a1,a2,t1,t2,Frank_param){
  u <- S1(t1,a1,alpha=0,Z=0)
  v <- S2(t2,a2,beta=0,Z=0)
  return(-log(1+(exp(-Frank_param*u)-1)*(exp(-Frank_param*v)-1)/(exp(-Frank_param)-1))/Frank_param)
}

Frank_deriv <- function(a1,a2,t1,t2,Frank_param){
  u <- S1(t1,a1,alpha=0,Z=0)
  v <- S2(t2,a2,beta=0,Z=0)
  return((exp(-Frank_param*u)*(exp(-Frank_param*v)-1))/((exp(-Frank_param)-1)+(exp(-Frank_param*u)-1)*(exp(-Frank_param*v)-1)))
}


F21_Frank <- function(t1,t2,a1,a2,alpha,beta, gamma,Z,Frank_param){
  u <- S1(t1,a1,alpha=0,Z=0)
  frac_Cu_C <- Frank_deriv(a1,a2,t1,t2,Frank_param)/Frank_cop(a1,a2,t1,t2,Frank_param)
  S21 <- S2(t2,a2,beta,Z)*exp_A_Frank(t1,t2,a1,a2, gamma,Z,Frank_param)*(exp(alpha*Z)-exp(gamma*Z)+exp(gamma*Z)*frac_Cu_C*u)/exp(alpha*Z)
  return(1-S21)
}

f_Frank <- function(x,v,t1i,a1,a2,alpha,beta, gamma,Zi,Frank_param){
  return(F21_Frank(t2=x,t1=t1i,a1,a2,alpha,beta, gamma,Zi,Frank_param)-v)
}

data_gen_Frank_no_cens <- function(a1=1, a2=1, alpha, beta, gamma, Z, Frank_param, setseed=FALSE, seed=123){
  Failure_times <- data.frame(T1=NA,T2=NA,covar=NA)
  n <- length(Z)
  if (isTRUE(setseed)){
    set.seed(seed)
  }
  u <- runif(n)
  w <- runif(n)
  for (j in 1:n){
    Zi <- Z[j]
    t1i <- F1_inv(u[j],a1,alpha,Zi)
    root <- uniroot(f_Frank, c(0,20),v=w[j],t1=t1i,a1=a1,a2=a2,alpha=alpha,beta=beta, gamma=gamma,Zi=Zi,Frank_param=Frank_param)
    t2i <- root$root
    Failure_times[j,] <- c(t1i,t2i,Zi)
  }
  return(Failure_times)
}



###--- Functions for sampling from the generalized Lehmann model with a Clayton copula as our baseline (any alpha)--------
S1_PODa <- function(lambda1,t1,alpha,Z){
  return(exp(-lambda1*t1*exp(alpha*Z)))
}


S2_PODa <- function(lambda2,t2,beta,Z){
  return(exp(-lambda2*t2*exp(beta*Z)))
}


exp_A_PODa <- function(lambda1,lambda2,t1,t2, gamma,Z,clay_param){
  u <- S1_PODa(lambda1,t1,alpha=0,Z=0)
  v <- S2_PODa(lambda2,t2,beta=0,Z=0)
  return((v^(clay_param-1)+u^(clay_param-1)-(u*v)^(clay_param-1))^(exp(gamma*Z)/(1-clay_param)))
}


S_PODa <- function(lambda1,lambda2,t1,t2,alpha,beta, gamma,Z,clay_param){ 
  return(S1_PODa(lambda1,t1,alpha,Z)*S2_PODa(lambda2,t2,beta,Z)*(exp_A_PODa(lambda1,lambda2,t1,t2,gamma,Z,clay_param)))
}


F1_inv_PODa <- function(u,lambda1,alpha,Z){
  return(-(log(1-u)/lambda1)*exp(-alpha*Z))
}

F21_PODa <- function(lambda1,lambda2,t1,t2,alpha,beta, gamma,Z,clay_param){
  S21_POD <- S2_PODa(lambda2,t2,beta,Z)*exp_A_PODa(lambda1,lambda2,t1,t2, gamma,Z,clay_param)*(1-exp((gamma-alpha)*Z)*exp(-lambda1*t1*(clay_param-1))*(1-exp(-lambda2*t2*(clay_param-1)))*(exp(-lambda1*t1*(clay_param-1))+exp(-lambda2*t2*(clay_param-1))-exp(-(lambda1*t1*(clay_param-1)+lambda2*t2*(clay_param-1))))^(-1))
  return(1-S21_POD)
}

f_PODa <- function(x,v,t1i,lambda1,lambda2,alpha,beta, gamma,Zi,clay_param){
  return(F21_PODa(lambda1,lambda2,t2=x,t1=t1i,alpha,beta, gamma,Zi,clay_param)-v)
}



data_gen_PODa_no_cens <- function(lambda1,lambda2,alpha, beta, gamma, Z,clay_param, setseed=FALSE, seed=123){
  Failure_times <- data.frame(T1=NA,T2=NA,covar=NA)
  n <- length(Z)
  if (isTRUE(setseed)){
    set.seed(seed)
  }
  u <- runif(n)
  T1 <- F1_inv_PODa(u,lambda1,alpha,Z)
  w <- runif(n)
  for (j in 1:n){
    Zi <- Z[j]
    t1i <- T1[j]
    root <- uniroot(f_PODa, c(0,10^3),v=w[j],t1=t1i,lambda1=lambda1,lambda2=lambda2,alpha=alpha,beta=beta, gamma=gamma,Zi=Zi,clay_param=clay_param)
    t2i <- root$root
    Failure_times[j,] <- c(t1i,t2i,Zi)
  }
  return(Failure_times)
}



###-------Bivariate lognormal failure times-----

gen_LN_data <- function(Sigma ,gamma1=0.2,gamma2=1, Z,seed){
  n <- length(Z) #Z is a numeric vector of size n (a univariate covariate)
  mu <- data.frame(mu1=gamma1*Z,mu2=gamma2*Z)
  set.seed(seed)
  Y2 <- t(apply(X=mu,MARGIN = 1,FUN = mvrnorm,n=1,Sigma=Sigma))
  Failure_times <- exp(Y2)
  #add univariate censoring
  set.seed(1000+seed)
  C <- rexp(n, rate=0.3)
  obs <- data.frame(obs1=pmin(Failure_times[,1],C),delta1=as.numeric(Failure_times[,1]<=C),obs2=pmin(Failure_times[,2],C), delta2=as.numeric(Failure_times[,2]<=C), C=C, Z=Z)
  obs$id <- 1:n
  return(obs)
}
