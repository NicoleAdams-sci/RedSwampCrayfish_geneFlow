#USAGE: Rscript rsc_CV4PE.R 'model' 'method' 'nreps' 'tol'
#Rscript rsc_CV4PE.R Grp2_bridge_indGrp1-5 neuralnet 10 0.05

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2")
args <- commandArgs(TRUE)
model <- as.character(args[1])
method <- as.character(args[2])
nreps <- as.numeric(args[3])
tol <- as.numeric(args[4])


fn <- paste0("rsc_CV4PE_new_",model, "_", method,"_", tol,".csv")


#Neural network model selection analysis - for cluster
library(abc)
library(tidyverse)

#rt <- read.csv("rsc_100k_subset.csv", header = TRUE)
rt <- read.csv("rsc_100k_subset_new.csv", header = TRUE)
rt <- rt %>% filter(Model == model)

#Now make objects that just have 1. the parameters (parms), 2. summary statistics (ss)
parms <- rt[,c(8:28)]   #N.LA to T.Grp5
ss <- rt[,-c(1:37)] # S.pop-Psi.pop1.pop2


#The loop for model selection...
for(x in 1:nreps) {
  #Select a random dataset simulated under "model"
  tmp.row <- sample(c(1:nrow(ss)), 1, replace = FALSE)
  tmp.obs <- ss[tmp.row,]
  tmp.true <- parms[tmp.row,]
  tmp.rt <- ss[-tmp.row,]
  tmp.parms <- parms[-tmp.row,]
  
  #Classify with ABC using "method"
  if(method == "neuralnet") {
    tmp.PE <- abc(target = tmp.obs, param = tmp.parms, sumstat = tmp.rt, tol = tol, method = method, numnet = 10, sizenet = 20, MaxNWts = 50000, transf = "logit", logit.bounds = t(apply(tmp.parms,2,range)))
  } else {
    tmp.PE <- abc(target = tmp.obs, param = tmp.parms, sumstat = tmp.rt, tol = tol, method = method, transf = "logit", logit.bounds = t(apply(tmp.parms,2,range)))
  }
  
  if(method == "neuralnet") {
    tmp.summ <- summary(tmp.PE)
    tmp.med <- tmp.summ[3,]
    tmp.low <- tmp.summ[2,]
    tmp.hi <- tmp.summ[6,]
    
    } else {
    #tmp.summ <- quantile(tmp.PE$adj.values, c(0.025, 0.5, 0.975), na.rm=T)
    tmp.summ <- apply(tmp.PE$adj.values, 2 , quantile , probs = c(0.025, 0.5, 0.975) , na.rm = TRUE )
    tmp.med <- tmp.summ[2,]
    tmp.low <- tmp.summ[1,]
    tmp.hi <- tmp.summ[3,]
    }
  
  names(tmp.med) <- paste0("est.", names(tmp.med))
  names(tmp.low) <- paste0("lowHPD.", names(tmp.low))
  names(tmp.hi) <- paste0("highHPD.", names(tmp.hi))
  
  #Summary of the result and save output
  tmp.out <- c(tmp.true, tmp.med, tmp.low, tmp.hi)
  
  if(!file.exists(fn)) {
    write.table(tmp.out, file = fn, sep = ",", quote = FALSE, row.names = FALSE)
  } else {
    write.table(tmp.out, file = fn, quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
  }
  
  #Remove temporary files
  rm(tmp.row, tmp.obs, tmp.true, tmp.rt, tmp.parms, tmp.PE, tmp.summ, tmp.med, tmp.low, tmp.hi, tmp.out)
  
}
