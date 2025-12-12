#Neural network model selection analysis - for cluster
#Usage: Rscript rsc_PE_subset.R Grp2_bridge_indGrp1-5 neuralnet 0.05
args <- commandArgs(TRUE)
model <- as.character(args[1])
method <- as.character(args[2])
tol <-as.numeric(args[3])

library(abc)

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2")

#rt <- read.csv("rsc_100k_subset.csv", header = TRUE)
rt <- read.csv("rsc_100k_subset_new.csv", header = TRUE)
rt <- rt[rt$Model == model,]

#Now make objects that just have 1. the parameters (parms), 2. summary statistics (ss), and 3. the model index
parms <- rt[,c(8:28)]   #N.LA to T.Grp5
ss <- rt[,-c(1:37)] # S.pop-Psi.pop1.pop2

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC")

rscobs <- read.csv("rsc.empirical.stats.csv", header = TRUE)
rscobs <- rscobs[,-c(1)]  #rm totSNPs 


setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2/rsc.modSel")

if(method == "neuralnet") {
  PE <- abc(target = rscobs, param = parms, sumstat = ss, tol = tol, method = method, MaxNWts = 50000, numnet = 10, sizenet = 20, maxit=1000, transf = "logit", logit.bounds = t(apply(parms,2,range)))
} else {
  PE <- abc(target = rscobs, param = parms, sumstat = ss, tol = tol, method = method, transf = "logit", logit.bounds = t(apply(parms,2,range)))
}


save(PE, file = paste0("rsc_PE_subset_new_", model, "_", method,"-",tol,".rda"))

