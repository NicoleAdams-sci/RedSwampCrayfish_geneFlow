library(abc)
library(tidyverse)

# Set working directory
setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2/rsc.modSel")

# Define files and parameters
files_info <- data.frame(
  file = c("rsc_modsel_100ksubset_MNLOG-0.005.rda",
           "rsc_modsel_100ksubset_MNLOG-0.01.rda", 
           "rsc_modsel_100ksubset_MNLOG-0.025.rda",
           "rsc_modsel_subset_NNET-0.005.rda",
           "rsc_modsel_subset_NNET-0.01.rda",
           "rsc_modsel_subset_NNET-0.025.rda"),
  method = c("mnlogistic", "mnlogistic", "mnlogistic",
             "neuralnet", "neuralnet", "neuralnet"),
  tolerance = c(0.005, 0.01, 0.025, 0.005, 0.01, 0.025)
)

# Function to load and extract results
## need new.env bc all loaded .rda files have the same named objects
# Function to safely load and extract results
extract_results_safe <- function(file, method, tolerance) {
  temp_env <- new.env()
  load(file, envir = temp_env)
  
  if(method == "mnlogistic") {
    obj <- temp_env$MNLOG
  } else {
    obj <- temp_env$NNET
  }
  
  obj_summary <- summary(obj)
  
  # Calculate correct accepted reps based on tolerance
  accepted_reps <- tolerance * 1000000
  
  # Extract results based on method
  if(method == "mnlogistic") {
    pred_probs <- obj_summary$mnlogistic$Prob
  } else if(method == "neuralnet") {
    pred_probs <- obj_summary$neuralnet$Prob
  }
  
  # Also get rejection results
  rejection_probs <- obj_summary$rejection$Prob
  
  pred_list <- as.list(pred_probs)
  rejection_list <- as.list(rejection_probs)
  
  # Create method results
  method_results <- data.frame(
    Algorithm = method,
    Tolerance = tolerance,
    AcceptedReps = accepted_reps
  )
  
  for(model_name in names(pred_list)) {
    method_results[[model_name]] <- pred_list[[model_name]]
  }
  
  # Create rejection results  
  rejection_results <- data.frame(
    Algorithm = "rejection",
    Tolerance = tolerance,
    AcceptedReps = accepted_reps  # Same calculation
  )
  
  for(model_name in names(rejection_list)) {
    rejection_results[[model_name]] <- rejection_list[[model_name]]
  }
  
  return(rbind(method_results, rejection_results))
}

# Extract all results
all_results <- list()
for(i in 1:nrow(files_info)) {
  all_results[[i]] <- extract_results_safe(
    files_info$file[i], 
    files_info$method[i], 
    files_info$tolerance[i]
  )
}

# Combine all results
combined_results <- bind_rows(all_results)

# Remove duplicate rejection results (they're the same across methods for each tolerance)
combined_results_clean <- combined_results %>%
  distinct(Algorithm, Tolerance, .keep_all = TRUE)

# Round numeric columns to 4 significant figures
combined_results_rounded <- combined_results_clean %>%
  mutate(across(where(is.numeric) & !c(AcceptedReps, Tolerance), ~ signif(.x, 4)))

# Save results
write.csv(combined_results_rounded, "table2A_model_selection_results_table.csv", row.names = FALSE, quote = FALSE)

# View the results in a clean format
print(combined_results_rounded)