calc_var <- function(fit1,fit2,obs,obs_long,covar_names,POD=FALSE,exchange=FALSE){ 
  n <- length(unique(obs$id))
  l <- cov_U1_U2(obs_long, obs, covar_names, fit1,fit2,POD,exchange)
  Sigma1 <- l$S1/n
  Sigma2 <- l$S2/n
  Sigma12 <- l$S12/n
  if (isTRUE(POD)){
    B2_inv <- n*solve(l$B2)
  } else {
    B2_inv <- n*fit2$vbeta.naiv
  }
  B1 <- B1_func2( obs_marg_unique=obs_long,
                  obs=obs,
                  fit1=fit1,
                  fit2=fit2,
                  form_marg = fit1$formula,
                  POD = POD)
  Sig21_invA1_B1 <- t(Sigma12) %*%   fit1$vbeta.naiv %*% t(B1) 
  middle1 <- Sigma2 +B1 %*% fit1$vbeta.naiv %*% Sigma1 %*% fit1$vbeta.naiv %*% t(B1)- Sig21_invA1_B1- t(Sig21_invA1_B1)
  middle2 <- Sigma2-t(Sigma12)%*%solve(Sigma1,Sigma12)
  V1 <- B2_inv %*% middle1 %*% B2_inv
  V2 <- B2_inv %*% middle2 %*% B2_inv
  return(list(V1,V2))
}

cov_U1_U2 <- function(obs_long, obs, covar_names, fit1,fit2, POD=FALSE, exchange=FALSE){
  n <- length(unique(obs$id))
  p <- length(covar_names)
  k1 <- length(unique(obs_long$time_margin))
  k2 <- length(unique(obs$times))
  if(isTRUE(exchange)){
    k1 <- length(unique(obs_long$new_time_names))
    k2 <- length(unique(obs$new_time_names))
  }
  
  clusters_1 <- split(obs_long, obs_long$id)
  clusters_2 <- split(obs, obs$id)
  
  xi1_hat <- fit1$beta
  xi2_hat <- fit2$beta
  
  B2 <- matrix(data = 0, nrow = p+k2, ncol = p+k2)
  if (isTRUE(exchange)){
    S12 <- matrix(data = 0, nrow = p+k1, ncol = p+k2)
    S11 <-  matrix(data = 0, nrow = p+k1, ncol = p+k1)
  } else {
  S12 <- matrix(data = 0, nrow = 2*p+k1, ncol = p+k2)
  S11 <-  matrix(data = 0, nrow = 2*p+k1, ncol = 2*p+k1)
  }
  S22 <- matrix(data = 0, nrow = p+k2, ncol = p+k2)
  
  nf <- gaussian()
  mu_1 <- rep(0,k1)
  mu_2 <- rep(0,k2)
  var_1 <- nf$variance(mu_1)
  var_2 <- nf$variance(mu_2)
  
  form1 <- fit1$formula
  form2 <- fit2$formula
  
  for (i in 1:n) {
    clus1 <- clusters_1[[i]]
    clus2 <- clusters_2[[i]]
    
    X_1i <- model.matrix(form1, data = clus1)
    X_2i <- model.matrix(form2, data = clus2)
    
    eta_1i <- as.vector(X_1i %*% xi1_hat)
    eta_2i <- as.vector(X_2i %*% xi2_hat)
    
    mu_1i     <- 1 - exp(-exp(eta_1i))              # mean of Y = 1 - PO
    dmu1_deta1 <- exp(eta_1i) * exp(-exp(eta_1i))     # dμ/dη
    
    resid_1i  <- (1 - clus1$value) - mu_1i                 # must match model response
    
    if (isTRUE(POD)){
      eta_2i_clip <- pmin(eta_2i, 5)
      mu_2i       <- pmax(exp(exp(eta_2i_clip)), .Machine$double.eps)
      dmu2_deta2  <- pmax(exp(eta_2i_clip) * exp(exp(eta_2i_clip)), .Machine$double.eps)
      
      resid_2i <- clus2$PO3_new - mu_2i
    } else {
      mu_2i <- 1-exp(-exp(eta_2i))
      dmu2_deta2 <- exp(eta_2i) * exp(-exp(eta_2i)) 
      resid_2i <- (1-clus2$PO3_new) - mu_2i
    }
    
    
    D_1i <- X_1i * dmu1_deta1  # each column scaled by dmu/deta
    D_2i <- X_2i * dmu2_deta2  # each column scaled by dmu/deta
    
    V_i_inv <- diag(1 / var_1, nrow = length(mu_1i))
    W_i_inv <- diag(1 / var_2, nrow = length(mu_2i))
    
    U_1i <- t(D_1i) %*% V_i_inv %*% resid_1i
    U_2i <- t(D_2i) %*% W_i_inv %*% resid_2i
    B2_i <- t(D_2i) %*% W_i_inv %*% D_2i
    #Sigma_12 <- Sigma_12 + U_1i %*% t(U_2i)
    S11 <- S11 + U_1i %*% t(U_1i)
    S22 <- S22 + U_2i %*% t(U_2i)
    S12 <- S12 + U_1i %*% t(U_2i)
    B2 <- B2 +B2_i
  }
  return(list(S1=S11,S2=S22,S12=S12, B2=B2))
}


B1_func2 <- function(obs_marg_unique, obs, fit1, fit2, form_marg, POD = FALSE) {
  # Handles both POD=FALSE (cloglog, geese) and POD=TRUE (double-log, geem)
  #
  # obs_marg_unique : obs_PO_v2_marginals_unique
  # obs             : obs_PO_v2 (must contain PO3_new)
  # fit1            : fit_v2_marg (step-1 fit)
  # fit2            : fit_2_step (step-2 fit: geese or geem)
  # form_marg       : step-1 formula used to build marginal model matrix
  # POD             : logical; TRUE -> POD (double log / geem), FALSE -> standard (cloglog / geese)
  #
  # Returns: B1 matrix (q2 x q1) with sign convention B1 = - dU2/dxi1 (so returns -B1)
  
  # Split clusters
  cl_marg <- split(obs_marg_unique, obs_marg_unique$id)
  cl_2    <- split(obs, obs$id)
  
  xi1 <- fit1$beta
  xi2 <- fit2$beta
  
  # Dimensions
  q1 <- length(xi1)  # step-1 params
  q2 <- length(xi2)  # step-2 params
  
  # Check marginal model matrix ordering matches fit1
  X_marg_cols <- colnames(model.matrix(form_marg, obs_marg_unique))
  if (!is.null(names(xi1))) {
    stopifnot(all(X_marg_cols == names(xi1)))
  }
  
  # Initialize B1 (q2 x q1)
  B1 <- matrix(0, nrow = q2, ncol = q1)
  rownames(B1) <- names(xi2)
  colnames(B1) <- names(xi1)
  
  form2 <- fit2$formula
  
  for (id in names(cl_2)) {
    clus2 <- cl_2[[id]]
    clusm <- cl_marg[[id]]
    if (is.null(clusm)) next
    
    # use only rows used in step-2 (non-marginal)
    clus2_nm <- dplyr::filter(clus2, t1 != 0, t2 != 0)
    if (nrow(clus2_nm) == 0) next
    
    # Step-2 model matrix and linear predictor
    X2  <- model.matrix(form2, data = clus2_nm)
    eta2 <- as.vector(X2 %*% xi2)
    
    # Choose mean and derivative depending on POD
    if (isTRUE(POD)) {
      # double-log link: mu = exp(exp(eta)), clip eta to avoid overflow
      eta2_clip <- pmin(eta2, 5)
      mu2  <- pmax(exp(exp(eta2_clip)), .Machine$double.eps)
      dmu2 <- pmax(exp(eta2_clip) * exp(exp(eta2_clip)), .Machine$double.eps)  # dmu/deta
      # response is PO3_new under POD branch
      y2 <- clus2_nm$PO3_new
      resid2 <- y2 - mu2
    } else {
      # cloglog link: mu = 1 - exp(-exp(eta))
      mu2  <- 1 - exp(-exp(eta2))
      dmu2 <- exp(eta2) * exp(-exp(eta2))
      # response in non-POD case is 1 - PO3_new
      y2 <- 1 - clus2_nm$PO3_new
      resid2 <- y2 - mu2
    }
    
    # D2 and W_inv (Gaussian working variance => identity)
    D2 <- X2 * dmu2
    W_inv <- diag(1, nrow = length(mu2))
    
    # Prepare dy/dxi1 (rows = observations in clus2_nm, cols = q1)
    dy_dxi1 <- matrix(0, nrow = nrow(clus2_nm), ncol = q1)
    
    # Marginal model matrix for this subject
    Xm <- model.matrix(form_marg, data = clusm)
    # prepare lookup: time_margin -> row index in clusm
    if (!("time_margin" %in% names(clusm))) {
      # If create_marg_df used a different naming, attempt to reconstruct
      clusm$time_margin <- paste0(clusm$name, "_", ifelse(clusm$name == "PO1", clusm$t1, clusm$t2))
    }
    lookup <- setNames(seq_len(nrow(clusm)), clusm$time_margin)
    
    for (j in seq_len(nrow(clus2_nm))) {
      t1j <- clus2_nm$t1[j]
      t2j <- clus2_nm$t2[j]
      
      r1 <- lookup[[paste0("PO1_", t1j)]]
      r2 <- lookup[[paste0("PO2_", t2j)]]
      
      if (is.null(r1) || is.null(r2)) next
      
      x1 <- Xm[r1, , drop = FALSE]
      x2 <- Xm[r2, , drop = FALSE]
      
      # linear predictors for marginals
      eta1 <- as.vector(x1 %*% xi1)
      eta2m <- as.vector(x2 %*% xi1)
      
      # derivative of log(pred) depending on marginal link (marginal used cloglog in your pipeline)
      # marginals are fit with cloglog: pred = exp(-exp(eta)) => log(pred) = -exp(eta)
      # d log(pred) / d xi1 = -exp(eta) * x
      dlogp1 <- -exp(eta1)  * x1
      dlogp2 <- -exp(eta2m) * x2
      
      # PO3_new = PO3 / (p1 * p2)
      # d PO3_new = - PO3_new * (d log p1 + d log p2)
      # For non-POD: y = 1 - PO3_new => dy = + PO3_new * (d log p1 + d log p2)
      # For POD:   y = PO3_new       => dy = - PO3_new * (d log p1 + d log p2)
      if (isTRUE(POD)) {
        # POD: y = PO3_new => dy = - PO3_new * (dlogp1 + dlogp2)
        dy_dxi1[j, ] <- - as.numeric(clus2_nm$PO3_new[j]) * as.numeric(dlogp1 + dlogp2)
      } else {
        # non-POD: y = 1 - PO3_new => dy = + PO3_new * (dlogp1 + dlogp2)
        dy_dxi1[j, ] <-  as.numeric(clus2_nm$PO3_new[j]) * as.numeric(dlogp1 + dlogp2)
      }
    } # end j loop
    
    # Contribution to B1: sum_i t(D2_i) W^{-1} (dy/dxi1)_i
    term <- t(D2) %*% W_inv %*% dy_dxi1
    B1 <- B1 + term
  } # end id loop
  
  # Return with same sign convention as before: B1 = - dU2/dxi1
  -B1
}
