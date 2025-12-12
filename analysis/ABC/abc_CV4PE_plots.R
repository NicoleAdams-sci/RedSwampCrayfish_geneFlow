# Plot cross validation of parameter estimation
# Figure S9
# 12/16/2023

setwd("/Users/neasci/Documents/crayfish_lab/RedSwampCrayfish_geneFlow/")

today <- format(Sys.time(), format = "%Y-%m-%d")

# read in cross validation (cv4pe) files
cv4pe.files <- list.files(path="output/ABC_output", pattern = "rsc_CV4PE_new.*0\\.05\\.csv$", full.names = TRUE)

parmorder = c(8, 5, 6, 3, 4, 7, 9, 12, 11, 13);

longname = c(
  "Ne - LA", # Contemporary Ne
  "Ne - UNK",
  "Ne - Group 3_SW",
  "Ne - Group 3_NE",
  "Ne - Group 2_Hotel",
  "Ne - Group 2_WestGC",
  "Ne - Group 4",
  "Ne - Group 1",
  "Ne - Group 5",
  "Migration rate w/in native range",
  "Migration rate w/in Group 3",
  "Migration rate w/in Group 2",
  "Bottleneck severity",
  "Time of Divergence - UNK", 
  "Time of Divergence - Group 3_NE",
  "Time of Divergence - Group 3_SW",
  "Time of Divergence - Group2_Hotel",
  "Time of Divergence - Group 4",
  "Time of Divergence - Group2_WestGC",
  "Time of Divergence - Group 1",
  "Time of Bottleneck - Group 5");


# Re-format data
cv4pe.list<-list()
for (FILE in cv4pe.files) {
  cv4pe.df <- read.csv(FILE, row.names = NULL, header=TRUE)
  dfNamA <- unlist(strsplit(FILE, "[/|]+"))[c(3)]
  dfNamB <- unlist(strsplit(dfNamA, "[.|_]+"))
  dfNam <- paste(dfNamB[4], dfNamB[5], dfNamB[6], dfNamB[7], sep = "_")
  
  cv4pe.df$model <- paste(dfNamB[4], dfNamB[5], dfNamB[6], sep = "_")
  cv4pe.df$method <- dfNamB[7]
  cv4pe.df$tol <- dfNamB[9]
  
  cv4pe.list[[dfNam]] <- cv4pe.df
}

# plot 
for (i in 1:length(cv4pe.list)) {
  ind.df <- as.data.frame(cv4pe.list[i])
  true <- ind.df[,1:21] #N.LA - T.Drain1
  estim <- ind.df[,22:42] #est.N.LA - est.T.Drain1
  pltnam <- paste0("figures/figureS9", names(cv4pe.list[i]), "_CV4PE_plot_", today,".pdf")
  
  pdf(pltnam)
  par(mfrow=c(4,3), mar = c(4.1, 4.1, 2.1, 1.1));  #mar = c(2,2,2,2)
  for(x in parmorder[c(1:10)]) {
    plot(true[,x], estim[,x], pch = 19, main = longname[x], xlim = c(range(c(true[,x], na.omit(estim[,x])))), ylim = c(range(c(true[,x], na.omit(estim[,x])))), xlab = "Simulated", ylab = "Estimated")
    abline(a = 0, b = 1, lwd = 1.5, col="red", lty="dashed")
  }
  plot.new()
  #legend("center", bty = "n", legend = names(cv4pe.list[i]))
  dev.off()
}
