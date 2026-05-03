rm(list = ls())

source("run_sim.R")
source("data_generating.R")
source("three_PO_functions.R")
source("GEE_functions.R")
source("geem source code.R")
source("fit_predict_gjrm.R")
source("variance and coverage functions.R")

#' An example for a sim_fn function
#' @param setting_param_list parameters that changes at each run
#' @param other_param_list static parameters for all the runs
#' @return a tibble with the results
sim_fn <- function(setting_param_list,fixed_param_list){
  seed <-  setting_param_list$seed
  n <-  setting_param_list$n
 
  a1 <- fixed_param_list$a1
  a2 <- fixed_param_list$a2
  alpha <- fixed_param_list$alpha
  beta <- fixed_param_list$beta
  gamma <- fixed_param_list$gamma
  Frank_param <- fixed_param_list$Frank_param
  covar <- fixed_param_list$covar
  Z <- covar[[which(names(covar) == n)]]
  t0 <- fixed_param_list$t0
  cop_names <- fixed_param_list$cop_names
  f.list <- fixed_param_list$f.list
  
  ###--- Sample simulated data----
  obs_uncensored <- data_gen_Frank_no_cens(a1=a1,a2=a2,alpha=alpha, beta=beta, gamma=gamma, Z=Z, Frank_param=Frank_param, setseed=TRUE, seed=seed)
  C <- rexp(n, rate=0.3)
  obs <- data.frame(obs1=pmin(obs_uncensored$T1,C),delta1=as.numeric(obs_uncensored$T1<=C),obs2=pmin(obs_uncensored$T2,C), delta2=as.numeric(obs_uncensored$T2<=C), C=C, Z=Z)
  obs$id <- 1:n  
  
  ###--- Calculate the POs------
  obs_PO_v2 <- PO_func_3D_v2(obs=obs,times=t0)
  t0$times <- 1:nrow(t0)
  obs_PO_v2 <- obs_PO_v2 %>%
    left_join(t0, by = "times")  
  
  
  ###--- Simple Lehmann model---
  # Fit simple model
  k <- nrow(t0)
  if (k==1) {
    form_simple <- formula(1-PO3~Z)
  } else {
    form_simple <- formula(1-PO3~as.factor(times)+Z)
  }
  fit_simple <- gee_d_variate(pseudo_obs= obs_PO_v2, k=k,form=form_simple,d=1)
  betas_df_simple <- fit_simple$beta
  mat_simple <- model.matrix(form_simple, obs_PO_v2)
  pred_simple_i <- exp(-exp(mat_simple%*%fit_simple$beta))
  pred_simple_df <- data.frame(matrix(pred_simple_i,n,k,byrow = TRUE))
  colnames(pred_simple_df) <- paste0("S_",1:nrow(t0),"")
  pred_simple <- pred_simple_df
  var_betas_sim <- data.frame(fit_simple$vbeta)
  
  ###--- Generalized Lehmann model-------
  # First step
  marg_list <- create_marg_df(obs_PO_v2) 
  obs_PO_v2_marginals_unique <- marg_list$df_unique
  k1 <- marg_list$k1
  k2 <- marg_list$k2
  
  form_marg <- formula(1-value~factor(time_margin)+factor(name):Z-1)
  
  fit_v2_marg <- gee_d_variate_new(pseudo_obs= obs_PO_v2_marginals_unique, k_vec = c(k1, k2),form=form_marg,d=2)
  var_marg <- data.frame(fit_v2_marg$vbeta) 
  
  # Marginal predictions
  mat_v2_marg <- model.matrix(form_marg,obs_PO_v2_marginals_unique)
  obs_PO_v2_marginals_unique$pred_v2_marg <-  exp(-exp(mat_v2_marg%*%fit_v2_marg$beta))
  
  obs_PO_v2 <- create_df_for_2step(obs_PO_v2, obs_PO_v2_marginals_unique)
  obs_PO_v2_no_marginals <- obs_PO_v2 %>%  filter(t1 != 0, t2 != 0)
  k <- length(unique(obs_PO_v2_no_marginals$times))
  
  if (mean(obs_PO_v2$PO3_new)>1){
    pod=TRUE
  } else {
    pod=FALSE
  }
  
  if (isTRUE(pod)){
    if (k==1) {
      form_2step <- formula(PO3_new~Z)
    } else {
      form_2step <- formula(PO3_new~as.factor(times)+Z)
    }
  } else {
    if (k==1) {
      form_2step <- formula(1-PO3_new~Z)
    } else {
      form_2step <- formula(1-PO3_new~as.factor(times)+Z)
    }
  }
  
  fit_2_step <- tryCatch({
    if (isTRUE(pod)){
      geem2(form_2step, id=id, data=obs_PO_v2_no_marginals, family = FunList, corstr = "independence")
    } else{
      gee_d_variate(pseudo_obs= obs_PO_v2_no_marginals, k=k,form=form_2step, d=1)
    }
  }, error = function(e) {
    message(sprintf("Error at iteration %d: %s", i, e$message))
    return(NULL)
  })
  if (!is.null(fit_2_step)) {
    betas_df_gen_v2 <- c(fit_v2_marg$beta, fit_2_step$beta)
    names(betas_df_gen_v2) <- c(names(fit_v2_marg$beta),c(paste0("gamma0_",1:k,""),"gamma1"))
    mat_v2 <- model.matrix(form_2step, obs_PO_v2_no_marginals)
    gamma_hat <- fit_2_step$beta
    if (isTRUE(pod)){
      pred_v2_vec <- exp(exp(mat_v2%*%gamma_hat))
      var_betas_gen_gee <- as.data.frame(as.matrix(fit_2_step$var))
    } else {
      pred_v2_vec <- exp(-exp(mat_v2%*%gamma_hat))
      var_betas_gen_gee <- data.frame(fit_2_step$vbeta)
    }
   
    obs_PO_v2_no_marginals <- obs_PO_v2 %>%
      filter(t1 != 0, t2 != 0) %>%
      mutate(pred = pred_v2_vec) %>%
      dplyr::select(id, t1, t2, pred)  # keep only columns to join
    
    # Join back only the pred column
    obs_PO_v2 <- obs_PO_v2 %>%
      left_join(obs_PO_v2_no_marginals, by = c("id", "t1","t2"))
    
    obs_PO_v2 <- obs_PO_v2 %>%
      mutate(
        pred = case_when(
          !is.na(pred) ~ pred*pred_PO1*pred_PO2,        # product of three terms
          t2 == 0       ~ pred_PO1,   # if t2==0, use pred_PO1
          t1 == 0       ~ pred_PO2,   # if t1==0, use pred_PO2
          TRUE          ~ NA_real_     # fallback, if needed
        )
      )
    
    pred_v2_df <- data.frame(matrix(obs_PO_v2$pred,n, nrow(t0),byrow = TRUE))
    colnames(pred_v2_df) <- paste0("S_",1:nrow(t0),"")

    
    var_beta <- calc_var(fit1=fit_v2_marg,fit2=fit_2_step,obs=obs_PO_v2,obs_long=obs_PO_v2_marginals_unique,covar_names="Z",POD = pod)
    var_betas_schur <- as.data.frame(as.matrix(var_beta[[2]] / n))
    var_betas_gen <- as.data.frame(as.matrix(var_beta[[1]] / n))
  } else {
    betas_df_gen_v2 <- c(fit_v2_marg$beta, rep(NA,1+k))
    names(betas_df_gen_v2) <- c(names(fit_v2_marg$beta),c(paste0("gamma0_",1:k,""),"gamma1"))
    pred_v2_df <- NA
    var_betas_gen<- NA
    var_betas_gen_gee <- NA
  }
  
  ###--- Copulas---
  fit_best <- fit_and_predict_gjrm(data=obs,formula_list = f.list, grid=t0,copulas = cop_names)
  mod <- fit_best$best_model
  var <- diag(mod$Vb)
  var_cop <- var[grep("^z|^\\(Intercept\\)", names(var))]
  coef_cop <- mod$coefficients
  betas_df_copulas <- unname(coef_cop[names(coef_cop) == "Z"])
  names(betas_df_copulas) <- c("beta1","beta2","beta3")
  # Predicted survival probabilities
  pred <- fit_best$predicted_survival
  pred_wide <- pred %>%
    mutate(grid_label = paste0("t1_", t1, "_t2_", t2)) %>%
    dplyr::select(subject, grid_label, S12) %>%
    pivot_wider(
      names_from = grid_label,
      values_from = S12
    ) %>%
    dplyr::select(-subject)
  
  colnames(pred_wide) <- paste0("S_",1:nrow(t0),"")
  pred_wide <- as.data.frame(lapply(pred_wide, as.numeric))
  pred_copulas <- pred_wide

  return(tibble::tibble(
    n = n,
    
    betas_simple = list(betas_df_simple),
    var_betas_simple = list(var_betas_sim),
    pred_simple = list(pred_simple),
    
    betas_gen = list(betas_df_gen_v2),
    var_marg = list(var_marf),
    var_betas_gen = list(var_betas_gen),
    var_betas_schur = list(var_betas_schur),
    var_betas_gen_gee = list(var_betas_gen_gee),
    pred_gen = list(pred_v2_df),
    
    betas_cop = list(betas_df_copulas),
    var_betas_cop = list(var_cop),
    pred_cop = list(pred_copulas)
  ))
}

n=c(400,800)
covars <- list()
for (i in 1:length(n)){
  sam_size <- n[i]
  set.seed(11111+i)
  covars[[i]] <- runif(sam_size,0,1) 
}
names( covars) <- n

a <- c(0.5,0.7)
b <- c(0.6,0.7,0.8)
t0 <- expand.grid(a,b)
colnames(t0) <- c("t1","t2")
t0_merged <- paste0(t0[,1],",",t0[,2], sep = "")

eq1 <- obs1 ~ Z+s(log(obs1), bs = "mpi")
eq2 <- obs2 ~ Z+s(log(obs2), bs = "mpi")
eq3 <- ~ Z
f.list <- list(eq1, eq2, eq3)

setting_param_list <- list(
  seed = 1:100,         
  n = n 
)

fixed_param_list <- list(a1=1,a2=1, alpha=1,beta=0.7,gamma=0.3, covar=covars,Frank_param=-5,t0=t0,cop_names="F",
                         f.list=f.list)


run_sim(sim_fn, setting_param_list, fixed_param_list ,nworkers = 60)