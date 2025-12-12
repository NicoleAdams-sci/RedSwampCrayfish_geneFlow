#USAGE: Rscript rsc_CV_subset.R 'method' 'model' 'nreps' 'tol' 'node'

args <- commandArgs(TRUE)
method <- as.character(args[1])
model <- as.character(args[2])
nreps <- as.numeric(args[3])
tol <- as.numeric(args[4])
node <- as.numeric(args[5])

#Neural network model selection analysis - for cluster
library(abc)

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2")

#rt <- read.csv("rsc_100k_subset.csv", header = TRUE)
rt <- read.csv("rsc_100k_subset_new.csv", header = TRUE)

#Now make objects that just have 1. the parameters (parms), 2. summary statistics (ss), and 3. the model index
parms <- rt[,c(1:37)] # Model-totSNPs
ss <- rt[,-c(1:37)] # S.pop-Psi.pop1.pop2
mod <- as.character(rt$Model)
rm(rt) #No need to keep this in the environment... try to save some space?

fn <- paste0("rsc_100k_CV_", method, "_", model, "-", tol, "-", node, ".csv")

ct <- 1
#The loop for model selection...
while(ct <= nreps) {
  #Select a random dataset simulated under "model"
  tmp.row <- sample(which(mod == model), 1, replace = FALSE)
  tmp.obs <- ss[tmp.row,]
  tmp.rt <- ss[-tmp.row,]
  tmp.mod <- mod[-tmp.row]
  
  #Classify with ABC using "method"
  if(method == "mnlogistic") {
    tmp.modsel <- postpr(target = tmp.obs, index = tmp.mod, sumstat = tmp.rt, tol = tol, method = "mnlogistic")
  } else if(method == "neuralnet") {
    tmp.modsel <- postpr(target = tmp.obs, index = tmp.mod, sumstat = tmp.rt, tol = tol, method = "neuralnet", numnet = 10, sizenet = 20, MaxNWts = 50000)
  }
  
  if(tmp.modsel$method != "rejection") {
    #Summary of the result and write individual lines of output to a file
    tmp.out <- c(true.model = model, tmp.modsel$pred)
    
    #Write line of output to .csv file
    if(!file.exists(fn)) {
      write(names(tmp.out), fn, ncolumns = length(tmp.out), sep = ",")
      write(tmp.out, fn, ncolumns = length(tmp.out), sep = ",", append = TRUE)
    } else {
      write(tmp.out, fn, ncolumns = length(tmp.out), sep=",", append = TRUE)
    }
    
    ct <- ct + 1
  } else {
    print("Not all models are represented in the accepted set, rejection only performed. Choosing another POD.")
  }
  
  #Remove temporary files
  rm(tmp.row, tmp.obs, tmp.rt, tmp.mod, tmp.modsel, tmp.out)
}
