run_sim <- function(sim_fn, setting_param_list, fixed_param_list, nworkers = 40) {
  # Load required packages
  library(furrr)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(tibble)
  
  # 1. Set up parallel workers
  future::plan("multisession", workers = min(nworkers, future::availableCores() - 1))
  
  # 2. Safely wrap sim_fn to capture errors
  sim_fn_safe <- purrr::safely(sim_fn)
  
  # 3. Create grid of parameter combinations
  param_grid <- do.call(tidyr::expand_grid, setting_param_list) %>%
    mutate(scenario = row_number())
  
  # 4. Convert each row into a named list of parameters
  setting_list <- purrr::transpose(param_grid %>% select(-scenario))
  
  # 5. Run simulations in parallel
  sim_results <- furrr::future_map(setting_list, ~ sim_fn_safe(.x, fixed_param_list), .progress = TRUE)
  
  # 6. Store results in a structured format
  results_final <- param_grid %>%
    mutate(
      setting = setting_list,
      res = map(sim_results, "result"),
      err = map(sim_results, "error")
    )
  
  # 7. Save results with timestamp
  timestamp <- str_replace_all(format(Sys.time(), "%Y%m%d_%H%M%S"), "[:\\s-]", "")
  filename <- paste0("res_", timestamp, ".RDS")
  saveRDS(results_final, file = filename, version = 2)
  
  # 8. Reset the parallel plan
  future::plan("sequential", .cleanup = TRUE)
  
  # Return the final results tibble
  return(results_final)
}
