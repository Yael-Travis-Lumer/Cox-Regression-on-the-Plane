library(pseudo)
library(mhazard)
library(Rfast)

###----- Dabrowska for simple models (1D theta) ---------
PO_func_1D <- function(obs,t0){
  n <- nrow(obs)
  old_name <- colnames(obs)
  m <- ncol(obs)
  l <- nrow(t0)
  obs_temp <- obs
  obs_temp[,(m+1):(m+l)] <- NA
  S_dabrowska <- npSurv2(Y1 = obs$obs1, Delta1 = obs$delta1,Y2 = obs$obs2, Delta2 = obs$delta2,newT1 = t0[,1],newT2 = t0[,2])
  S_hat <- diag(S_dabrowska$Fhat_est)
  S_hat_mat  <- do.call("rbind", replicate(n,S_hat,simplify = FALSE))
  #compute leave-one-out estimate - using for loop or apply with function
  for (i in 1:n){
    obs_i <- obs[-i,]
    S_dabrowska_i <- npSurv2(Y1 = obs_i$obs1, Delta1 = obs_i$delta1,Y2 = obs_i$obs2, Delta2 = obs_i$delta2,newT1 = t0[,1],newT2 = t0[,2])
    obs_temp[i,(m+1):(m+l)] <- as.list(diag(S_dabrowska_i$Fhat_est))
  }
  #compute POs
  obs[,(m+1):(m+l)] <- n*S_hat_mat-(n-1)*obs_temp[,(m+1):(m+l)]
  
  colnames(obs) <- c(old_name,paste0("PO_t",1:l,""))
  obs$id <- 1:n
  return(obs)
}




###---- Dabrowska for generalized Lehmann model (3D theta) ----
PO_func_3D_v2 <- function(obs,times){
  n <- nrow(obs)
  obs$id <- 1:n
  old_name <- colnames(obs)
  k <- nrow(times)
  t0 <- matrix(0,3*k,2)
  obs_rep <- do.call("rbind", replicate(k, obs, simplify = FALSE))
  obs_rep <- obs_rep[sort(obs_rep$id),]
  obs_rep$times <- rep(1:k,n)
  for (j in 1:k){
    t1 <- times[j,1]
    t2 <- times[j,2]
    t0[(3*(j-1)+1):(3*j),1] <- c(t1,0,t1)
    t0[(3*(j-1)+1):(3*j),2] <- c(0,t2,t2)
  }
  l <- nrow(t0)
  S_dabrowska <- npSurv2(Y1 = obs$obs1, Delta1 = obs$delta1,Y2 = obs$obs2, Delta2 = obs$delta2,newT1 = t0[,1],newT2 = t0[,2])
  S_hat <- diag(S_dabrowska$Fhat_est)
  theta <- matrix(0,k,3)
  for (j in 1:k){
    theta[j,] <- c(S_hat[3*(j-1)+1],S_hat[3*(j-1)+2],S_hat[3*(j-1)+3])
  }
  theta_mat  <- do.call("rbind", replicate(n,theta,simplify = FALSE))
  
  thetas <-data.frame(matrix(0,n*k,3))
  for (i in 1:n){
    obs_i <- obs[-i,]
    S_dabrowska_i <- npSurv2(Y1 = obs_i$obs1, Delta1 = obs_i$delta1,Y2 = obs_i$obs2, Delta2 = obs_i$delta2,newT1 = t0[,1],newT2 = t0[,2])
    S_hat_i<- diag(S_dabrowska_i$Fhat_est)
    theta_i <- matrix(0,k,3)
    for (j in 1:k){
      theta_i[j,] <- c(S_hat_i[3*(j-1)+1],S_hat_i[3*(j-1)+2],S_hat_i[3*(j-1)+3])
    }
    thetas[(k*(i-1)+1):(k*i),] <- theta_i
  }
  POs <- n*theta_mat-(n-1)*thetas
  obs_rep <- cbind(obs_rep,POs)
  colnames(obs_rep) <- c(old_name,"times","PO1","PO2","PO3")
  return(obs_rep)
}
