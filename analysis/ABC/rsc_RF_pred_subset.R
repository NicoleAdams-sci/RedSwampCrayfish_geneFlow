#Random Forest - build the forest -- for cluster
library(abc)
library(abcrf)

setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2")
dir.create("rsc.modSel")
cores <- 6

#rt <- read.csv("rsc_100k_subset.csv", header = TRUE)
rt <- read.csv("rsc_100k_subset_new.csv", header = TRUE)
rt[,1] <- gsub(" ","-",rt[,1])

# It is better to convert Model to a factor.
# I think the function abcrf has issues with it if not.

rt$Model <- as.factor(rt$Model)

#Now make objects that just have 1. the parameters (parms), 2. summary statistics (ss), and 3. the model index
parms <- rt[,c(1:37)] # Model-totSNPs
ss <- rt[,-c(1:37)] # S.pop-Psi.pop1.pop2

mod <- rt$Model

abcrf_df <- data.frame(mod,ss)
if(is.factor(abcrf_df$mod)){
  print("Mod is a factor, proceeding...")
}

forest <- abcrf(mod~., data = abcrf_df, lda = FALSE, ntree = 1000, paral=TRUE, ncores=cores)
save(forest, file="rsc.modSel/rsc_RF_forest_100ksubset.rda")   #Added this so we still save the forest, even if prediction fails

rscobs <- read.csv("../SRC/rsc.empirical.stats.csv", header = TRUE)
rscobs <- rscobs[,-c(1)]  #rm totSNPs 

# # For some reason, predict() fails when the obs argument is a vector, so rbind the ashes object to itself.  The output is duplicated, but the analysis will run.

#setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT/rsc.modSel")
setwd("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2/rsc.modSel")

obs2 <- rbind(rscobs, rscobs)

rscpred <- predict(object = forest, obs = obs2, training = abcrf_df, ntree = 1000)
rscpred

save.image("RF_rscpred_100ksubset.Rws")
