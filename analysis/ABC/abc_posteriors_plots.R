#Script to make combined posterior distribution plots for the best RSC ABC models
# Figure S8

library(abc)

setwd("~/Documents/crayfish_lab/RedSwampCrayfish_geneFlow/")

parms <- read.csv("output/ABC_output/rsc_100k_params4PEplots_new.csv", header = TRUE)


########### Plot only neural net (no loclinear) ###########

#### Group 2 models ####
load("output/ABC_output/rsc_PE_subset_new_Grp2_step_indGrp1-5_neuralnet-0.05.rda", Grp2.step.indGrp1_5.nnet.05 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp2_bridge_indGrp1-5_neuralnet-0.05.rda", Grp2.bridge.indGrp1_5.nnet.05 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp2_step_indGrp1-5_neuralnet-0.1.rda", Grp2.step.indGrp1_5.nnet.1 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp2_bridge_indGrp1-5_neuralnet-0.1.rda", Grp2.bridge.indGrp1_5.nnet.1 <- new.env())


#Now make the plots...
parmorder = c(8, 5, 6, 3, 4, 7, 9, 12, 11, 13);

parmnames = Grp2.step.indGrp1_5.nnet.05$PE$names$parameter.names;
longname = c(
  "Ne - LA", # all are contemporary Ne
  "Ne - UNK",
  "Ne - Group 3_SW",
  "Ne - Group 3_NE",
  "Ne - Group 2_Hotel", 
  "Ne - Group 2_WestGC",
  "Ne - Group 4",
  "Ne - Group 1",
  "Ne - Group 5",
  "Migration Rate within native range",
  "Migration Rate within Group 3",
  "Migration Rate within Group 2",
  "Bottleneck severity",
  "Time of Divergence - UNK", 
  "Time of Divergence - Group 3_SW",
  "Time of Divergence - Group 3_NE",
  "Time of Divergence - Group2_Hotel",
  "Time of Divergence - Group 4",
  "Time of Divergence - Group2_WestGC",
  "Time of Divergence - Group 1",
  "Time of Bottleneck - Group 5");

bws = c(123, 123, 23, 23, 23, 50, 23, 23, 50, 0.0006563, 0.0001229, 0.0001229, 0.003051, 66, 0.7877, 0.5085, 0.779, 0.5509, 0.4662, 0.8136, 0.5509);


upy=c();
for(j in 1:length(parmnames)) {
  upy[j] = 1.1*max(c(max(density(parms[,j], bw = bws[j])$y), max(density(Grp2.step.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j])$y), max(density(Grp2.bridge.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j])$y)))
}

today <- format(Sys.time(), format = "%Y-%m-%d")
pdf.nam <- paste0("figures/figureS8A_PE_Grp2_plots_nNet_2tols_", today, ".pdf")
pdf(pdf.nam)

par(mfrow=c(4,3), mar = c(4,4,2,2));
for(j in parmorder[c(1:10)]) {
  
  plot(density(parms[,j], bw = bws[j]), ylim = c(0,upy[j]), lty = "dotted", lwd = 0.75, main = longname[j], xlab = "Parameter Estimate", ylab = "Posterior Density");
  lines(density(Grp2.step.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j]), col= "blue3", lty = "solid", lwd = 1.5);
  lines(density(Grp2.bridge.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j]), col = "skyblue", lty = "solid", lwd = 0.7);
  lines(density(Grp2.step.indGrp1_5.nnet.1$PE$adj.values[,j], bw = bws[j]), col= "gray40", lty = "solid", lwd = 1.5);
  lines(density(Grp2.bridge.indGrp1_5.nnet.1$PE$adj.values[,j], bw = bws[j]), col = "gray", lty = "solid", lwd = 0.7);
  lines(density(parms[,j], bw = bws[j]), col = "black", lty = "dotted", lwd = 0.75)
}

plot.new()
legend("center", legend = c("Prior","Grp2.step NNet_0.05", "Grp2.bridge NNet_0.05",
                            "Grp2.step NNet_0.1", "Grp2.bridge NNet_0.1"),
       lty = c("dotted","solid", "solid", "solid", "solid"),
       col=c("black", "blue3", "skyblue", "gray40", "gray"),
       lwd = c(0.75, 1.5, 0.7,  1.5, 0.7), cex = 0.85, y.intersp=0.9)

dev.off()


#### Group 3 models ####
load("output/ABC_output/rsc_PE_subset_new_Grp3_step_indGrp1-5_neuralnet-0.05.rda", Grp3.step.indGrp1_5.nnet.05 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp3_bridge_indGrp1-5_neuralnet-0.05.rda", Grp3.bridge.indGrp1_5.nnet.05 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp3_step_indGrp1-5_neuralnet-0.1.rda", Grp3.step.indGrp1_5.nnet.1 <- new.env())
load("output/ABC_output/rsc_PE_subset_new_Grp3_bridge_indGrp1-5_neuralnet-0.1.rda", Grp3.bridge.indGrp1_5.nnet.1 <- new.env())

#identical(Grp2.step.indGrp1_5.nnet.05$PE$names$parameter.names, Grp3.step.indGrp1_5.nnet.05$PE$names$parameter.names)

fc.upy=c();
for(j in 1:length(parmnames)) {
  fc.upy[j] = 1.1*max(c(max(density(parms[,j], bw = bws[j])$y), max(density(Grp3.step.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j])$y), max(density(Grp3.bridge.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j])$y)))
}

pdf.nam2 <- paste0("figures/figureS8B_PE_Grp3_plots_nNet_2tols_", today, ".pdf")
pdf(pdf.nam2)
par(mfrow=c(4,3), mar = c(4,4,2,2));

for(j in parmorder[c(1:10)]) {
  
  plot(density(parms[,j], bw = bws[j]), ylim = c(0,fc.upy[j]), lty = "dotted", lwd = 0.75, main = longname[j], xlab = "Parameter Estimate", ylab = "Posterior Density");
  lines(density(Grp3.step.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j]), col= "forestgreen", lty = "solid", lwd = 1.5);
  lines(density(Grp3.bridge.indGrp1_5.nnet.05$PE$adj.values[,j], bw = bws[j]), col = "palegreen2", lty = "solid", lwd = 0.7);
  lines(density(Grp3.step.indGrp1_5.nnet.1$PE$adj.values[,j], bw = bws[j]), col= "gray40", lty = "solid", lwd = 1.5);
  lines(density(Grp3.bridge.indGrp1_5.nnet.1$PE$adj.values[,j], bw = bws[j]), col = "gray", lty = "solid", lwd = 0.7);
  lines(density(parms[,j], bw = bws[j]), col = "black", lty = "dotted", lwd = 0.75)
}

plot.new()
legend("center", legend = c("Prior","Grp3.step NNet_0.05", "Grp3.bridge NNet_0.05",
                            "Grp3.step NNet_0.1", "Grp3.bridge NNet_0.1"),
       lty = c("dotted","solid", "solid", "solid", "solid"),
       col=c("black", "forestgreen", "palegreen2", "gray40", "gray"),
       lwd = c(0.75, 1.5, 0.7, 1.5, 0.7), cex = 0.85, y.intersp=0.9)

dev.off()



##### To get point estimates of parameters -- Table S8
#summary(Grp2.bridge.indGrp1_5.nnet.05$PE)

# Function to extract summary statistics from ABC results
extract_PE_summary <- function(abc_result, model_name) {
  summ <- summary(abc_result$PE)
  
  # Extract weighted median (row 3) and 95% CI (rows 2 and 6)
  median_vals <- summ[3, ]
  lower_ci <- summ[2, ]  # 2.5%
  upper_ci <- summ[6, ]  # 97.5%
  
  # Get parameter names for conditional formatting
  param_names <- names(median_vals)
  
  # Create formatted strings with different precision based on parameter type
  formatted_vals <- character(length(median_vals))
  
  for(i in 1:length(median_vals)) {
    # Use 3 decimal places for migration rates and bottleneck severity
    if(grepl("mig|bottleneck", param_names[i], ignore.case = TRUE)) {
      formatted_vals[i] <- paste0(
        sprintf("%.3f", median_vals[i]), 
        " (", 
        sprintf("%.3f", lower_ci[i]), 
        ", ", 
        sprintf("%.3f", upper_ci[i]), 
        ")"
      )
    } else {
      # Use 2 decimal places for population sizes and other parameters
      formatted_vals[i] <- paste0(
        sprintf("%.2f", median_vals[i]), 
        " (", 
        sprintf("%.2f", lower_ci[i]), 
        ", ", 
        sprintf("%.2f", upper_ci[i]), 
        ")"
      )
    }
  }
  
  # Create data frame
  result_df <- data.frame(
    Parameter = param_names,
    Value = formatted_vals,
    stringsAsFactors = FALSE
  )
  
  # Rename the value column to the model name
  names(result_df)[2] <- model_name
  
  return(result_df)
}

# Extract summaries for each model
grp2_bridge <- extract_PE_summary(Grp2.bridge.indGrp1_5.nnet.05, "Grp2_bridge_indGrp1_5")
grp2_step <- extract_PE_summary(Grp2.step.indGrp1_5.nnet.05, "Grp2_step_indGrp1_5")
grp3_bridge <- extract_PE_summary(Grp3.bridge.indGrp1_5.nnet.05, "Grp3_bridge_indGrp1_5")
grp3_step <- extract_PE_summary(Grp3.step.indGrp1_5.nnet.05, "Grp3_step_indGrp1_5")

# Combine all results
combined_table <- merge(grp2_bridge, grp2_step, by = "Parameter", all = TRUE)
combined_table <- merge(combined_table, grp3_bridge, by = "Parameter", all = TRUE)
combined_table <- merge(combined_table, grp3_step, by = "Parameter", all = TRUE)

# Create parameter name mapping for prettier display
param_mapping <- c(
  "N.Grp1" = "Group 1 Ne",
  "N.Grp2A" = "Group 2_Hotel Ne", 
  "N.Grp2B" = "Group 2_WestGC Ne",
  "N.Grp3B" = "Group 3_NE Ne",
  "N.Grp3A" = "Group 3_SW Ne",
  "N.Grp4" = "Group 4 Ne",
  "N.Grp5" = "Group 5 Ne",
  "migGrp2" = "Group 2 migration",
  "migGrp3" = "Group 3 migration", 
  "bottleneckSeverity" = "Bottleneck severity"
)

# Apply parameter name mapping
combined_table$Parameter_Name <- param_mapping[combined_table$Parameter]

# Remove unmapped parameters and reorder
final_table <- combined_table[!is.na(combined_table$Parameter_Name), ]
final_table <- final_table[, c("Parameter_Name", "Grp2_bridge_indGrp1_5", "Grp2_step_indGrp1_5", 
                               "Grp3_bridge_indGrp1_5", "Grp3_step_indGrp1_5")]

# Rename first column
names(final_table)[1] <- "Model"

# Reorder rows to match your desired order
desired_order <- c("Group 1 Ne", "Group 2_Hotel Ne", "Group 2_WestGC Ne", 
                   "Group 3_NE Ne", "Group 3_SW Ne", "Group 4 Ne", "Group 5 Ne",
                   "Group 2 migration", "Group 3 migration", "Bottleneck severity")

final_table$Model <- factor(final_table$Model, levels = desired_order)
final_table <- final_table[order(final_table$Model), ]
final_table$Model <- as.character(final_table$Model)

# Print the table
print("Parameter Estimation Summary Table (Neural Network, 0.05 tolerance)")
print("================================================================")
print(final_table, row.names = FALSE)

# Save to CSV
write.csv(final_table, paste0("tables/tableS8_parameter_estimation_summary_table_", today, ".csv"), row.names = FALSE)

# Also create a nicely formatted version for copy-paste
cat("\n\nFormatted for copy-paste:\n")
cat("========================\n")
for(i in 1:nrow(final_table)) {
  cat(sprintf("%s\t%s\t%s\t%s\t%s\n", 
              final_table$Model[i],
              final_table$Grp2_bridge_indGrp1_5[i],
              final_table$Grp2_step_indGrp1_5[i], 
              final_table$Grp3_bridge_indGrp1_5[i],
              final_table$Grp3_step_indGrp1_5[i]))
}
