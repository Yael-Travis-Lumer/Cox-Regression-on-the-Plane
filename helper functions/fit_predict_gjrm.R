library(copula)
library(GJRM)


###--- Copula predictions per patient ------
predict_biv_surv_vectorized <- function(fit, Z_row, grid, copula_name = "N") {

  n_grid <- nrow(grid)

  # Marginal survival for obs1
  newdata1 <- Z_row[rep(1, n_grid), , drop = FALSE]
  newdata1$obs1 <- grid[,1]
  eta1 <- as.numeric(predict(fit, newdata = newdata1, eq = 1, type = "link"))
  s1 <- exp(-exp(eta1))

  # Marginal survival for obs2
  newdata2 <- Z_row[rep(1, n_grid), , drop = FALSE]
  newdata2$obs2 <- grid[,2]
  eta2 <- as.numeric(predict(fit, newdata = newdata2, eq = 2, type = "link"))
  s2 <- exp(-exp(eta2))

  # Extract copula parameter (theta) — vector of length n_grid (same Z replicated)
  eta <- predict(fit, newdata = newdata1, eq = 3, type = "link")
  # Transform eta to the association partameter theta using the corresponding inverse link function used in eq3
  if (copula_name == "N") {
    theta <- tanh(eta)
    par <- pmax(pmin(as.numeric(theta), 0.999), -0.999)
    cop_list <- lapply(par, function(r) normalCopula(param = r, dim = 2))
  } else if (copula_name == "F") {
    par <- as.numeric(eta)
    cop_list <- lapply(par, function(p) frankCopula(param = p, dim = 2))
  } else if (copula_name == "C0") {
    theta <- exp(eta)
    par <- pmax(as.numeric(theta), 0.001)
    cop_list <- lapply(par, function(p) claytonCopula(param = p, dim = 2))
  } else if (copula_name == "G0") {
    theta <- 1+exp(eta)
    par <- pmax(as.numeric(theta), 1.001)
    cop_list <- lapply(par, function(p) gumbelCopula(param = p, dim = 2))
  } else {
    stop("Unsupported copula: ", copula_name)
  }

  # Compute bivariate survival probabilities
  s12 <- numeric(n_grid)
  for (i in seq_len(n_grid)) {
    u1 <- 1 - s1[i]
    u2 <- 1 - s2[i]
    cval <- pCopula(c(u1, u2), copula = cop_list[[i]])
    s12[i] <- s1[i] + s2[i] - 1 + cval  # Inclusion-exclusion
  }

  # Return grid with S1, S2, and S12
  return(cbind(grid, S1 = s1, S2 = s2, S12 = s12, par=par))
}

####----- Copula predictions for all patients ------
predict_biv_surv_all <- function(fit, newdata, grid, copula_name = "N") {
  n <- nrow(newdata)

  all_results <- lapply(1:n, function(i) {
    subject_data <- newdata[i, , drop = FALSE]
    out <- data.frame(predict_biv_surv_vectorized(fit, subject_data, grid, copula_name))
    out$subject <- i
    return(out)
  })

  res <- do.call(rbind, all_results)
}


####---- Function that fits several copula models and selects the best one ----------
fit_and_predict_gjrm <- function(data, formula_list, margins = c("-cloglog", "-cloglog"),
                                 copulas = c("N", "C0", "G0","F"),
                                 grid=t0, individuals = NULL, verbose = TRUE) {
  
  fits <- list()
  bic_vals <- numeric(length(copulas))
  
  # === Fit all candidate copula models ===
  for (i in seq_along(copulas)) {
    cop <- copulas[i]
    if (verbose) cat("Fitting copula:", cop, "\n")
    fit <- tryCatch({
      gjrm(
        formula = formula_list,
        copula = cop,
        data = data,
        margins = margins,
        cens1 = data$delta1,
        cens2 = data$delta2,
        model = "B"
      )
    }, error = function(e) {
      message("Failed to fit copula: ", cop)
      return(NULL)
    })
    fits[[cop]] <- fit
    bic_vals[i] <- if (!is.null(fit)) BIC(fit) else Inf
  }
  
  copula_comparison <- data.frame(Copula = copulas, BIC = bic_vals)
  copula_comparison <- copula_comparison[order(copula_comparison$BIC), ]
  best_copula <- copula_comparison$Copula[1]
  best_fit <- fits[[best_copula]]
  if (verbose) cat("Best copula selected:", best_copula, "\n")
  
  # === Identify individuals to predict for ===
  if (is.null(individuals)) {
    individuals <- seq_len(nrow(data))
  }
  if (!is.null(fit)){
    grid_result <- predict_biv_surv_all(fit=best_fit, newdata=data, grid=grid, copula_name = best_copula)
  } else {
    grid_reult <- NULL
  }

  return(list(
    copula_comparison = copula_comparison,
    best_model = best_fit,
    best_copula = best_copula,
    predicted_survival = grid_result
  ))
}

