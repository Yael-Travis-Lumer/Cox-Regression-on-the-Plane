###---- GEE with d-variate cloglog link-------
gee_d_variate <- function(pseudo_obs, k,form,d){
  link <- rep("cloglog",d)
  link_rep <- rep(link, k)
  # Running the model fit
  fit <- geese(formula = form,
               data = pseudo_obs,
               id = id,
               mean.link = link_rep,
               variance = rep("gaussian", k * d),
               sca.link = rep("identity", k * d),
               scale.value = 1,
               scale.fix = 1,
               link.same = T,
               gm = 1)
  return(fit)
} 



###---- Create data frame from marginal analysis------
create_marg_df <- function(obs_PO_v2){
  
  # Drop PO3 and PO3_new safely (if present)
  obs_PO_v2_marginals <- obs_PO_v2 %>%
    dplyr::select(-dplyr::any_of(c("PO3", "PO3_new")))
  
  # Long format
  obs_PO_v2_marginals_long <- obs_PO_v2_marginals %>%
    pivot_longer(cols = starts_with("PO"))
  
  # Add the time margin variable
  obs_PO_v2_marginals_long <- obs_PO_v2_marginals_long %>%
    mutate(time_margin = case_when(
      name == "PO1" ~ paste0(name, "_", t1),
      name == "PO2" ~ paste0(name, "_", t2),
      TRUE ~ NA_character_
    ))
  
  # Count zero-margins BEFORE removal
  n_PO1_zero <- sum(obs_PO_v2_marginals_long$name == "PO1" & obs_PO_v2_marginals_long$t1 == 0)
  n_PO2_zero <- sum(obs_PO_v2_marginals_long$name == "PO2" & obs_PO_v2_marginals_long$t2 == 0)
  
  # Print warnings if needed
  if (n_PO1_zero > 0) {
    warning(paste0("Removed ", n_PO1_zero, " PO1 marginal(s) with t1 = 0."))
  }
  if (n_PO2_zero > 0) {
    warning(paste0("Removed ", n_PO2_zero, " PO2 marginal(s) with t2 = 0."))
  }
  
  # Remove rows where the margin corresponds to a zero time
  obs_PO_v2_marginals_long <- obs_PO_v2_marginals_long %>%
    filter(!(name == "PO1" & t1 == 0)) %>%
    filter(!(name == "PO2" & t2 == 0))
  
  # Keep unique margins
  obs_PO_v2_marginals_unique <- obs_PO_v2_marginals_long %>%
    distinct(id, time_margin, .keep_all = TRUE)
  
  # Count distinct marginal times for PO1
  k1 <- obs_PO_v2_marginals_unique %>%
    filter(name == "PO1") %>%
    distinct(t1) %>%
    nrow()
  
  # Count distinct marginal times for PO2
  k2 <- obs_PO_v2_marginals_unique %>%
    filter(name == "PO2") %>%
    distinct(t2) %>%
    nrow()
  
  return(list(df_unique = obs_PO_v2_marginals_unique,
              k1 = k1, k2 = k2))
}


###---- Create data frame from the second step estimation procedure ------
create_df_for_2step <- function(obs_PO_v2, obs_PO_v2_marginals_unique){
  
  # Split into two data frames: PO1 and PO2 predictions
  pred_PO1 <- obs_PO_v2_marginals_unique %>%
    filter(name == "PO1") %>%
    dplyr::select(id, t1, pred_v2_marg) %>%
    rename(pred_PO1 = pred_v2_marg)
  
  pred_PO2 <- obs_PO_v2_marginals_unique %>%
    filter(name == "PO2") %>%
    dplyr::select(id, t2, pred_v2_marg) %>%
    rename(pred_PO2 = pred_v2_marg)
  
  # Join PO1 predictions by id, t1
  obs_PO_v2 <- obs_PO_v2 %>%
    left_join(pred_PO1, by = c("id", "t1"))
  
  # Join PO2 predictions by id, t2
  obs_PO_v2 <- obs_PO_v2 %>%
    left_join(pred_PO2, by = c("id", "t2"))
  
  
  # ----- Assign pred = 1 when the margin is zero -----
  
  # Number of zero t1 replaced
  n_PO1_zero_replaced <- sum(obs_PO_v2$t1 == 0)
  # Number of zero t2 replaced
  n_PO2_zero_replaced <- sum(obs_PO_v2$t2 == 0)
  
  # Set pred_PO1 = 1 for t1 = 0
  obs_PO_v2 <- obs_PO_v2 %>%
    mutate(pred_PO1 = ifelse(t1 == 0, 1, pred_PO1))
  
  # Set pred_PO2 = 1 for t2 = 0
  obs_PO_v2 <- obs_PO_v2 %>%
    mutate(pred_PO2 = ifelse(t2 == 0, 1, pred_PO2))
  
  # Warning messages
  if (n_PO1_zero_replaced > 0) {
    warning(paste0("pred_PO1 was set to 1 for ", n_PO1_zero_replaced, " observation(s) with t1 = 0."))
  }
  if (n_PO2_zero_replaced > 0) {
    warning(paste0("pred_PO2 was set to 1 for ", n_PO2_zero_replaced, " observation(s) with t2 = 0."))
  }
  
  
  # ----- Compute PO3_new -----
  obs_PO_v2$PO3_new <- obs_PO_v2$PO3 / (obs_PO_v2$pred_PO1 * obs_PO_v2$pred_PO2)
  
  
  return(obs_PO_v2)
}


###---Two cloglog links and an additional log link---------
gee_trivariate_log_link <-  function(pseudo_obs, k,form){
  link <- c("cloglog", "cloglog", "log")
  link <- rep(link,k)
  fit <- geese(formula = form, 
               data = pseudo_obs,
               id = id,
               mean.link = link,
               variance = rep("gaussian", k * 3),
               sca.link = rep("identity", k * 3),
               scale.value = 1,
               scale.fix = 1,
               link.same = F,
               gm = 1)
  return(fit)
} 

###-------- loglog link ------------
linkfun = function(mu) log(log(mu))  # g(x)=loglog(x)
variance <- function(mu) { rep.int(1, length(mu)) }
linkinv = function(eta) {pmax(exp(exp(eta)), .Machine$double.eps)}  # Inverse: x = exp(exp(y))
mu.eta = function(eta) {# Derivative of inverse
  eta <- pmin(eta, 5)
  pmax(exp(eta) * exp(exp(eta)), .Machine$double.eps)
}  
FunList <- list(LinkFun=linkfun, VarFun=variance, InvLink=linkinv, InvLinkDeriv=mu.eta)
