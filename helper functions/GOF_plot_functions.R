# Goodness of fit
gof_general <- function(obs_w_PO_pred,PO_name,pred_name, cont_covar_name, categ_name_color, categ_name_shape,
                            times_of_interest=NULL,x_label,y_label, ylim_low=NULL,ylim_up=NULL, xlim_low=NULL,xlim_up=NULL, shape_size=2){
  if (!is.null( times_of_interest)){
    obs_w_PO_pred <- subset(obs_w_PO_pred,
                          times %in% times_of_interest)
  }
  obs_w_PO_pred$pearson_resid <- unlist((obs_w_PO_pred[PO_name]-obs_w_PO_pred[pred_name])/sqrt(obs_w_PO_pred[pred_name]*(1-obs_w_PO_pred[pred_name])))
  clean_data <- obs_w_PO_pred[is.finite(obs_w_PO_pred$pearson_resid), ]
 if (is.null(ylim_low) & is.null(ylim_up)){
   ylim_low <- min(clean_data$pearson_resid)
   ylim_up <- max(clean_data$pearson_resid)
 }
  if (is.null(xlim_low) & is.null(xlim_up)){
    xlim_low <- min(clean_data[[cont_covar_name]])
    xlim_up <- max(clean_data[[cont_covar_name]])
  }
   g <- ggplot(obs_w_PO_pred, aes(x = as.numeric(.data[[cont_covar_name]]), 
                               y = as.numeric(pearson_resid),
                               color = .data[[categ_name_color]],
                               shape = .data[[categ_name_shape]])) +
  facet_wrap(~times_name) +
  geom_point(size = shape_size, alpha = 0.7) +
  ylim(ylim_low, ylim_up) +
  xlim(xlim_low, xlim_up) +
  #geom_smooth(method = "loess", se = FALSE, colour = "black") +  # Removed group = times
  geom_smooth(method = "loess", se = FALSE, aes(group = times), colour="black")+
     labs(title = "Pearson residuals based on POs",
       x = "Age (in months)",
       y = y_label) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

return(g)
}

# Conditional probabilities
get_index_matrix <- function(k, n) {
  outer(1:k, 1:n, function(j, i) k * i - (k - j))
}

get_conditional_prob <- function(obs_orig,obs_w_pred, pred_name,times_of_interest=1:6) {
  k <- length(times_of_interest)
  n <- nrow(obs_orig)
  ind_mat <- get_index_matrix(k,n)
  obs_w_pred <- subset(obs_w_pred, times %in% times_of_interest)
  
  obs_orig$T1_gt_24_given_T2_leq_12 <-  pmin(pmax(0,
                                        (obs_w_pred[[pred_name]][ind_mat[2,]]-obs_w_pred[[pred_name]][ind_mat[3,]])/(1-obs_w_pred[[pred_name]][ind_mat[1,]])),1)
  obs_orig$T1_gt_24_given_T2_gt_12 <-pmin(pmax(0,obs_w_pred[[pred_name]][ind_mat[3,]]/obs_w_pred[[pred_name]][ind_mat[1,]]),1)
  obs_orig$T2_gt_24_given_T1_leq_12 <- pmin(pmax(0,(obs_w_pred[[pred_name]][ind_mat[5,]]-obs_w_pred[[pred_name]][ind_mat[6,]])/(1-obs_w_pred[[pred_name]][ind_mat[4,]])),1)
  obs_orig$T2_gt_24_given_T1_gt_12 <-pmin(pmax(0,obs_w_pred[[pred_name]][ind_mat[6,]]/obs_w_pred[[pred_name]][ind_mat[4,]]),1)
  return(obs_orig)
}


plot_conditional <- function(df, title,y_name,xlim_low =0, xlim_up = 125) {
  ggplot(df, aes(x = age, y = .data[[y_name]],
                 color = income_level, shape = tumor_stage)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = title,
         x = "Age (in months)", y = "Survival probability") +
    ylim(0, 1) + xlim(xlim_low, xlim_up) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

