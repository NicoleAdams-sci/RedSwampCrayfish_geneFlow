#Neural network model selection analysis - for cluster
args <- commandArgs(TRUE)
tol <-as.numeric(args[1])

library(abc)
library(tidyverse) 

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2")

# read in simulation results tables
#res <- read.csv("rsc_resultsTest.csv") 

# remove NAs from results table
#res.nona <- res %>% na.omit()

#res.nona.100k <- read.csv("rsc_100k_subset.csv", header = TRUE)
res.nona.100k <- read.csv("rsc_100k_subset_new.csv", header = TRUE)

#Now make objects that just have 1. the parameters (parms), 2. summary statistics (ss), and 3. the model index 
parms <- res.nona.100k[,c(1:37)] # Model-totSNPs
ss <- res.nona.100k[,-c(1:37)] # S.pop-Psi.pop1.pop2
mod <- as.character(res.nona.100k$Model) 


setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC")
rscobs <- read.csv("rsc.empirical.stats.csv", header = TRUE)
rscobs <- rscobs[,-c(1)]  #rm totSNPs 

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2/rsc.modSel")

NNET <- postpr(target = rscobs, index = mod, sumstat = ss, tol = tol, method = "neuralnet", numnet = 10, sizenet = 20, MaxNWts = 50000, maxit=1000)
save(NNET, file = paste0("rsc_modsel_subset_NNET-",tol,".rda"))
