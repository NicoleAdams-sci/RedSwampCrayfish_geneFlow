# Model ABC simulations
# 07/22/2025
# Function-ish

# caution the model and pop names were changed retroactively, so there might be bugs in this script

# Adding a command line argument to incorporate ${SLURM_ARRAY_TASK_ID} into the filename for the "node"
# Based on https://github.com/stranda/holoSimCell/blob/master/R/holoStats.R
# Then you can run this as 100 array jobs (sbatch -a 1-100 <job submissions script>.sh)
# The job submission script will then have a line that looks something like Rscript constant.R ${SLURM_ARRAY_TASK_ID} 10000 (if you were running 10000 reps per job)
args <- commandArgs(TRUE)
node <- as.numeric(args[1])
nreps <- as.numeric(args[2])
runBy <- as.character(args[3])


# print("working directory")
# getwd()
# print("args")
# print(args)
# print("event.file")
# print(event.file)


## load in libraries
library(swfscMisc, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1")
library(strataG, lib.loc="/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
#library(dplyr, lib.loc="/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1")
library(KScorrect, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1")
library(LDcorSV, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1")
library(parallel, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # for mclapply() needed in pwise.het() in helper-functions.R
library(prodlim, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # for row.match() needed to order event times
library(popgraph, lib.loc="/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # for calculate conditional genetic distance, betweenness centrality, closeness centrality
library(gstudio, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # used with popgraph
library(igraph, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # needed for betweenness in graph_thy
library(hierfstat, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # to calc summary stats
library(adegenet, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # to convert df2genind
library(tictoc, lib.loc= "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/myR/x86_64-pc-linux-gnu-library/4.1") # code timing

tic("Full job")
# Preliminary environ and file settings
options(scipen=999) #A penalty to be applied when deciding to print numeric values in fixed or exponential notation. Positive values bias towards fixed and negative towards scientific notation
num_cores <- 1
exec <- "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/fsc2709"  #fs27 on the cluster

# Specify a filename for output
fn <- paste0("paramYstats_", runBy, ".",node ,".csv")



##########
#Simulation begins from this line down

# Model options
model.set <- c("5cols_step_x", "Grp3_bridge_mROW", "Grp3_bridge_indGrp1-5", "Grp3_step_mROW", "Grp3_step_indGrp1-5", "Grp2_bridge_mROW", "Grp2_bridge_indGrp1-5", "Grp2_step_mROW", "Grp2_step_indGrp1-5", "W2E_step_x")

# While loop to iterate the simulation over a specified number of repeats
num.rep <-0

while(num.rep < nreps) {
  # When you go to run this on the cluster, it's best to point the temporary outputs to a directory in your $SCRATCH space
  simdir <- system("echo $SCRATCH", intern = TRUE)   #capture the output of the command as an R character vector.
  # You'll also need a directory for the output, and you'll want to move back and forth between these during the loop
  #outdir <- "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT" 
  outdir <- "/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2" 
  
  # Change to simdir before the simulation runs
  setwd(simdir)
  print("new working directory")
  print(getwd())
  
  ## Fixed Parameters
   label <- sample(model.set, 1, replace = FALSE)
   event.file <- paste0(label, "_sourceSink.txt")
  
  # Mutation rate from [Liu et al. 2016](https://doi.org/10.1093/molbev/msw226)
  G <- 1	#Generation time
  mu <- 3.6e-9	#Mutation rate (nuisance parameter - fix or specify a range, but don't worry about estimating) 
  
  
  #sample_size <- rep(10,9)
  #sample_size <- c(10, 10, 225, 215, 210, 210, 40, 60, 60)
  sample_size <- c(6, 10, 135, 264, 206, 189, 36, 54, 54) # matches empirical data
  sample_size[2] <- 0 # unsampled UNK pop
  popids <- c("LA","UNK", "Grp3A", "Grp3B", "Grp2A", "Grp2B", "Grp4", "Grp1", "Grp5")
  
  indexids <- rbind(popids, indexids=c(0:8))
  
  
  
  ## Create population info (demelist)
  
  N.LA = round(KScorrect::rlunif(1,100,10000))		#Ne Louisana pop  # changed from uniform dist to log uniform 9/29/23
  N.UNK = round(KScorrect::rlunif(1,100,10000))		#Ne Unsampled native pop  
  N.Grp3A = round(KScorrect::rlunif(1,10,2000))		#Ne Group 3 group 1
  N.Grp3B = round(KScorrect::rlunif(1,10,2000))		#Ne Group 3 group 2
  N.Grp2A = round(KScorrect::rlunif(1,10,2000))		#Ne Group 2 hotel group
  N.Grp2B = round(KScorrect::rlunif(1,10,2000))		#Ne Group 2 WestGC group
  N.Grp4 = round(KScorrect::rlunif(1,10,2000))		#Ne Group 4 group
  N.Grp1 = round(KScorrect::rlunif(1,10,2000))		#Ne Group 1
  N.Grp5 = round(KScorrect::rlunif(1,10,2000))		#Ne Group 5
  pop_size <- c(N.LA, N.UNK, N.Grp3A, N.Grp3B, N.Grp2A, N.Grp2B, N.Grp4, N.Grp1, N.Grp5)
  
  
  demelist <- vector("list",length(sample_size))
  for(p in 1:length(sample_size)) {
    demelist[[p]] <- fscDeme(deme.size = pop_size[p], sample.size = sample_size[p], sample.time = 0, inbreeding = 0, growth = 0)
  }
  
  #names(demelist) <- popids => rename after run sims
  
  demes <- do.call(fscSettingsDemes, c(demelist, ploidy = 2))
  
  
  ## Genetic markers
  seq_length <- 150
  nloci <- 1808
  num_markers_sim <- nloci*1.50
  
  genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = seq_length, mut.rate = mu), num.chrom = num_markers_sim)  
  
  
  
  ## Migration matrix
  
  migmat <- matrix(data = 0, nrow = 9, ncol = 9)
  
  # migration between Louisiana and Unsampled (within native range)
  m.NA <-runif(1, 0.01, 0.05) # high m rate
  # migration between EastGC_NE and EastGC_SW (within all of Group 3)
  m.FC <- KScorrect::rlunif(1, 0.0001, 0.01) # med m rate - log uniform distribution rlunif(), changed from uniform dist runif()
  # migration between Hotel group and WestGC (within Group 2)
  m.Grp2A <- KScorrect::rlunif(1, 0.0001, 0.01) # lower m rate - log uniform distribution rlunif(), changed from uniform dist runif()
  
  
  migmat[1,2] <- m.NA
  migmat[2,1] <- m.NA
  
  migmat[3,4] <- m.Grp3
  migmat[4,3] <- m.Grp3
  
  migmat[5,6] <- m.Grp2
  migmat[6,5] <- m.Grp2
  
  migration <- fscSettingsMigration(migmat)
  
  
  ## Historical events
  ## Event times (coalescence and bottlenecks)
  # always true
  T.UNK = round(runif(1,1000,5000)) #Time of coalescence to LA pop
  T.MI = 50
  
  if (grepl("indGrp1-5", label)){  # plus separate Group1 and Group 5 introduction
    T.Grp1 = round(runif(1,2, T.MI)) #Time of coalescence for Group 1
  }
  
  ### FC ###
  if (grepl("FC", label)){  
    T.Grp3A = round(runif(1,2, T.MI))		
    T.Grp2A = round(runif(1,2, T.Grp3A)) 	
    
    if (grepl("bridge", label)) {
      T.Grp4 = round(runif(1,2, T.Grp3A)) #Time of coalescence for Group 4
      if (grepl("mROW", label)){  # minus (without) separate Group1 and Group 5 introductions (aka interior)
        T.Grp1 = round(runif(1,2, T.Grp3A)) #Time of coalescence for Group 1
      }
    }
    
    if (grepl("step", label)) {
      T.Grp4 = round(runif(1,2, T.Grp2A))   #Time of coalescence for Group 4
      if (grepl("mROW", label)){  
        T.Grp1 = round(runif(1,2, T.Grp2A)) #Time of coalescence for Group 1
        T.Grp4 = round(runif(1,2, T.Grp1)) #Time of coalescence for Group 4
      }
    }
  }
  
  ### Grp2A ###
  if (grepl("Grp2A", label)){  
    T.Grp2A = round(runif(1,2, T.MI))		#Time of coalescence for Group 2 subgroup hotel
    T.Grp3A = round(runif(1,2, T.Grp2A)) 	#Time of coalescence for Group 3 subgroup 1
    
    if (grepl("bridge", label)) {
      T.Grp4 = round(runif(1,2, T.Grp2A))   #Time of coalescence for Group 4
      if (grepl("mROW", label)){  
        T.Grp1 = round(runif(1,2, T.Grp2A)) #Time of coalescence for Group 1
      }
    }
    
    if (grepl("step", label)) {
      T.Grp4 = round(runif(1,2, T.Grp3A))   #Time of coalescence for Group 4
      if (grepl("mROW", label)){  
        T.Grp1 = round(runif(1,2, T.Grp3A)) #Time of coalescence for Group 1
        T.Grp4 = round(runif(1,2, T.Grp1)) #Time of coalescence for Group 4
      }
    }
  }
  
  
  ### W2E ###
  if (grepl("W2E", label)){  
    T.Grp1 = round(runif(1,2, T.MI))   #Time of coalescence for Group 1
    T.Grp2A = round(runif(1,2, T.Grp1))		#Time of coalescence for Group 2 subgroup hotel
    T.Grp3A = round(runif(1,2, T.Grp2A)) 	#Time of coalescence for Group 3 subgroup 1 
    T.Grp4 = round(runif(1,2, T.Grp2A))   #Time of coalescence for Group 4
  }
  
  ### 5COL ###
  #if (grepl("4col", label)){  
   if (grepl("5col", label)){
    T.Grp3A = round(runif(1,2, T.MI)) 	#Time of coalescence for Group 3 subgroup 1 
    T.Grp2A = round(runif(1,2, T.MI))		#Time of coalescence for Group 2 subgroup hotel
    T.Grp1 = round(runif(1,2, T.MI))   #Time of coalescence for Group 1
    T.Grp4 = round(runif(1,2, T.MI))   #Time of coalescence for Group 4
  }
  
  # always true - constrained
  T.Grp3B = round(runif(1,2, T.Grp3A))		#Time of coalescence for Group 3 subgroup 2, constrained by T.Grp3A
  T.Grp2B = round(runif(1,2, T.Grp2A)) #Time of coalescence for Group 2 subgroup WestGC, constrained by T.Grp2A
  T.Grp5 = round(runif(1,2, T.Grp1)) #Time of coalescence for Group 5, constrained by T.Grp1
  
  # bottleneck times
  T.UNK_1 = T.UNK-1	#Time of bottleneck for UNK
  T.Grp3A_1 = T.Grp3A-1	#Time of bottleneck for Group 3 subgroup 1
  T.Grp3B_1 = T.Grp3B-1 #Time of bottleneck for Group 3 subgroup 2
  T.Grp2A_1 = T.Grp2A-1 #Time of bottleneck for Group 2
  T.Grp4_1 = T.Grp4-1 #Time of bottleneck for Group 4
  T.Grp2B_1 = T.Grp2B-1 #Time of bottleneck for Group2 WestGC
  T.Grp1_1 = T.Grp1-1 #Time of bottleneck for Group 1 into UNK
  T.Grp5_1 = T.Grp5-1 #Time of bottleneck for Group5 into UNK
  
  bott.sev <- round(KScorrect::rlunif(1,0.001,0.25), 3) # changed from unifrom dist 9/29/23
  
  
  # build dataframe for historical event info
  #eventA <- as.data.frame(read.delim(paste0(system("echo $HOME", intern = TRUE), event.file)))
  print("current directory")
  print(getwd())
  print("change directory")
  setwd(outdir)
  print("new current directory")
  print(getwd())
  
  eventA <- as.data.frame(read.delim(paste0("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SHELL/", event.file)))
  
  
  prop.migrants <- c(rep(1,8), rep(0,8))
  new.size <- c(rep(1,8), rep(bott.sev,8))
  
  event.df <- as.data.frame(cbind(eventA, prop.migrants, new.size, new.growth=rep(0,16), migr.mat=rep(0,16)))
  
  # fill in event times
  event.df$E.times <- NA
  for(x in 1:length(event.df$E.times)){
    event.df$E.times[x] <- get(event.df$time.names[x])
  }
  
  
  # loop to check each event time for instances where a population serves as both the recipient of lineages from another deme and the source of lineages moving to another deme (backward in time) and make sure event type #2 happens after event type #1
  for(t in unique(event.df$E.times)){
    tmp <- event.df[event.df$E.times == t,]
    tmp <- tmp[tmp$source != tmp$sink,]
    if(length(tmp$E.times) > 1) {
      srcsnkPOP <- tmp$source[which(tmp$source %in% tmp$sink)]
      if(length(srcsnkPOP) > 0) {
        if(length(srcsnkPOP) == 1) {
          srcline <- tmp[tmp$source == srcsnkPOP,]
          snkline <- tmp[tmp$sink == srcsnkPOP,]
          srcrow <- prodlim::row.match(srcline, event.df)
          snkrow <- prodlim::row.match(snkline, event.df)
          if(srcrow < max(snkrow)) {
            event.df[max(snkrow),] <- srcline
            event.df[srcrow,] <- snkline[nrow(snkline),]
          }
          rm(srcline,snkline,srcrow,snkrow)
        } else {
          srclinelist <- vector("list", length(srcsnkPOP))
          snklinelist <- vector("list", length(srcsnkPOP))
          srcrowlist <- vector("list", length(srcsnkPOP))
          snkrowlist <- vector("list", length(srcsnkPOP))
          issue <- vector("list", length(srcsnkPOP))
          for(z in 1:length(srcsnkPOP)) {
            srclinelist[[z]] <- tmp[tmp$source == srcsnkPOP[z],]
            snklinelist[[z]] <- tmp[tmp$sink == srcsnkPOP[z],]
            srcrowlist[[z]] <- prodlim::row.match(srclinelist[[z]], event.df)
            snkrowlist[[z]] <- prodlim::row.match(snklinelist[[z]], event.df)
            if(sum(srcrowlist[[z]] < snkrowlist[[z]]) > 0) {
              issue[[z]] <- TRUE
            } else {
              issue[[z]] <- FALSE
            }
          }
          if(sum(unlist(issue)) > 0) {
            problem = TRUE
          } else {
            problem = FALSE
          }
          
          while(problem == TRUE) {
            for(z in 1:length(srcsnkPOP)) {
              srclinelist[[z]] <- tmp[tmp$source == srcsnkPOP[z],]
              snklinelist[[z]] <- tmp[tmp$sink == srcsnkPOP[z],]
              srcrowlist[[z]] <- prodlim::row.match(srclinelist[[z]], event.df)
              snkrowlist[[z]] <- prodlim::row.match(snklinelist[[z]], event.df)
              if(sum(srcrowlist[[z]] < snkrowlist[[z]]) > 0) {
                issue[[z]] <- TRUE
                wrongrows <- c(srcrowlist[[z]],snkrowlist[[z]])
                event.df[max(wrongrows),] <- srclinelist[[z]]
                wrongrows <- wrongrows[-which(wrongrows == max(wrongrows))]
                event.df[wrongrows,] <- snklinelist[[z]]
              } else {
                issue[[z]] <- FALSE
              }
            }
            
            if(sum(unlist(issue)) > 0) {
              problem = TRUE
            } else {
              problem = FALSE
            }
          }
        }
      }
      rm(srcsnkPOP)
    }
    rm(tmp)
  }
  
  
  event.df <- event.df[, c(1, 4, 2, 3, 5, 6, 7, 8)]
  
  event.df[,2:(length(event.df))] <- sapply(event.df[,2:(length(event.df))], as.numeric)
  
  # build list for historical event info
  eventlist <- list()
  for(ev in 1:nrow(event.df)) {
    eventlist[[ev]] <- fscEvent(event.time = event.df$E.times[ev], source = event.df$source[ev], sink = event.df$sink[ev], prop.migrants = event.df$prop.migrants[ev], new.size = event.df$new.size[ev], new.growth = event.df$new.growth[ev], migr.mat = event.df$migr.mat[ev])
  }
  
  # create events obj for fsc
  events <- do.call(fscSettingsEvents, eventlist)
  
  
  # Add event times to kill ea deme (needed for fsc27)
  kill.df <- events[events$source == events$sink,] 
  kill.df$event.time <- kill.df$event.time+1
  kill.df$new.size <- 0
  events <- rbind(events, kill.df)
  
  
  # order events newest to oldest
  events <- events[order(events$event.time),]
  
  print("go back to scratch")
  setwd(simdir)
  print(getwd())
  #setwd("~/RSC/ABC")
  
  #Define functions needed below
  #Modified function from Eric Archer's strataG fscTutorial() markdown
  #Adding pair.per.loc, simulating diploid data, need to retain 2 columns for each locus
  source("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/holoSimCell_helper-functions.R") #on Cluster
  source("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/holoSimCell_additional_stats_functions_NEA.R") #on Cluster
  source("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/holoSimCell_segmentedGLM.R") #on Cluster
  
  missMat <- read.delim("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/jy22yLoA_adults_baseFilter_imiss99_dp5miss50_imiss75_hz.ab_1pRt_100kb_rmdups_orderedMissMat.txt")
  
  # Run fsc in loop to eval no of SNPs and rerun if too low
  # from John's github https://github.com/stranda/holoSimCell/blob/master/R/runFSC_step_agg3.R
  ### Run fsc in loop 
  #MAF= 3/(2*944)  #MAF - minor allele frequency filter
  #nloci = 900     # moved up in the script 6/13/23
  maxloc = 500000       #max number of marker loci to attempt in fastsimcoal
  
  
  
  varSNPs <- 0
  while(varSNPs < nloci) {
    p <- fscWrite(demes = demes, genetics = genetics, events = events, migration = migration, label = paste0(label, ".", node, ".", num.rep), use.wd = TRUE)
    p <- fscRun(p, all.sites = FALSE, inf.sites = FALSE, dna.to.snp = TRUE, quiet = FALSE, num.cores = num_cores, exec = exec)
    out <- fscReadArp(p)
    
    print(date())
    #print(paste("cleaning up fsc files:",label, "-", nreps))
    #fscCleanup(paste0(label, "-", nreps))
    print(paste("cleaning up fsc files:",label, ".", node, ".", num.rep))
    fscCleanup(paste0(label, ".", node, ".", num.rep))
    
    message(paste("Coalescent simulation with",attr(genetics,"num.chrom"), "loci resulted in", (ncol(out)-2)/2, "polymorphic sites"), appendLF=TRUE)
    
    fscout <- sampleOnePerLocus(mat = out) # 'MAF = MAF' rm 6/13/23
    
    
    varSNPs_nomask <- sum(numAlleles(df2gtypes(fscout, ploidy = 2))$num.alleles == 2)  #A new variable to see how many markers are variable in the unmasked dataset
    if(varSNPs_nomask >= nloci) {   #Only proceed with masking if the matrix is big enough
      fscextra <- fscout[,c(1:2,((nloci*2)+3):ncol(fscout))] #Take the first nloci*2+2 columns as the initial dataset of nloci genotypes + 2 info columns (individual and deme ID), the rest are extra
      fscout <- fscout[,1:((nloci*2)+2)]    #The first 1808 loci + 2 ID columns are the dataset to mask
      fscout[missMat == TRUE] <- NA
      tmp_gtype <- df2gtypes(fscout, ploidy = 2)
      SNPalleles <- numAlleles(tmp_gtype)
      varSNPs <- sum(SNPalleles == 2)   #How many markers have 2 alleles?
      if(varSNPs != nloci) {   #What to do if some loci have lost polymorphism?
        for(p in SNPalleles$locus[SNPalleles$num.alleles < 2]) {
          # print(p)
          tmp_extra <- fscextra
          tmp_extra[missMat[,grep(p, colnames(fscout))[1]] == TRUE,-c(1,2)] <- NA   
          tmp_alleles <- numAlleles(df2gtypes(tmp_extra, ploidy = 2))
          if(sum(tmp_alleles$num.alleles == 2) > 0) {
            swapin <- tmp_alleles$locus[tmp_alleles$num.alleles == 2][1]
            print(paste("Swapping in", swapin, "from fscextra for", p, "from fscout"))
            fscout[,grep(p,colnames(fscout))] <- tmp_extra[,grep(swapin, colnames(tmp_extra))]
            colnames(fscout)[grep(p, colnames(fscout))] <- colnames(tmp_extra)[grep(swapin, colnames(tmp_extra))]
            fscextra <- fscextra[,-grep(swapin, colnames(fscextra))]
            # print(paste("we now have", sum(numAlleles(df2gtypes(fscout, ploidy = 2)) == 2), "variable loci after the mask"))  #This print statement was slowing things down (b/c of the df2gtypes, only included for testing)
          } else {
            print(paste("No simulated markers from fscextra are still variable when applying the mask for", p))
          }
          rm(tmp_extra, tmp_alleles, swapin)    #remove temporary variables before moving on to the next monomorphic locus in fscout
        }
        tmp_gtype <- df2gtypes(fscout, ploidy = 2)  #Recalculate the number of variable SNPs after the mask is applied
        SNPalleles <- numAlleles(tmp_gtype)$num.alleles
        varSNPs <- sum(SNPalleles == 2)
        print(paste("Finished masking data -", varSNPs, "polymorphic markers after applying NA mask"))
      }
    } else {
      tmp_gtype <- df2gtypes(fscout, ploidy = 2)
      varSNPs <- sum(numAlleles(tmp_gtype)$num.alleles == 2)
    }
    
    message(paste("Subsampling to one SNP per locus...", varSNPs, "loci have at least one polymorphic sites"), appendLF=TRUE)
    
    if(varSNPs < nloci) {
      if(attr(genetics, "num.chrom") == maxloc) {
        #stop("Too few SNPs pass the MAF!  Redrawing parameter values for this replicate!  Setting maxSNPstried to a larger value will lead to longer simulations")
        stop("Too few SNPs pass!  Redrawing parameter values for this replicate!  Setting maxSNPstried to a larger value will lead to longer simulations")
      }
      scaleSNP <- 1.1*(attr(genetics, "num.chrom")/varSNPs)
      newSNPnum <- round(scaleSNP*attr(genetics, "num.chrom"),0)
      if(newSNPnum > maxloc) {
        newSNPnum <- maxloc
        message(paste("Coalescent simulation with", attr(genetics, "num.chrom"), "loci resulted in", varSNPs, "variable markers. Trying again with the maximum number of loci allowable -",newSNPnum, "!"), appendLF=TRUE)
      } else {
        message(paste("Coalescent simulation with", attr(genetics, "num.chrom"), "loci resulted in", varSNPs, "variable markers. Trying again with",newSNPnum, "loci!"), appendLF=TRUE)
      }
      attr(genetics, "num.chrom") <- newSNPnum
      rm(scaleSNP, newSNPnum)
    } else if(varSNPs > nloci) {
      nalleles <- numAlleles(tmp_gtype)
      varchromnames <-  nalleles$locus[nalleles$num.alleles == 2]
      keep.loci <- sample(varchromnames, nloci, replace = FALSE)
      keep.columns <- c(1,2)
      for(locus in keep.loci) {
        keep.columns <- c(keep.columns, grep(locus, colnames(fscout)))
      }
      fscout <- fscout[,keep.columns]
    }
  }
  
  # gather params data
  params <- data.frame(Model= label, date = date(), node=node, runBy=runBy,
                       GenTime=G, mu=format(mu, scientific = T), seqLength=seq_length,
                       N.LA=N.LA, N.UNK=N.UNK, N.Grp3A=N.Grp3A, N.Grp3B=N.Grp3B, N.Grp2A=N.Grp2A, N.Grp2B=N.Grp2B, N.Grp4=N.Grp4, N.Grp1=N.Grp1, N.Grp5=N.Grp5,
                       migNat=round(m.NA,4), migFC=round(m.FC, 4), migGrp2A=round(m.Grp2A,4), bottleneckSeverity=bott.sev,
                       T.UNK=T.UNK, T.Grp3A=T.Grp3A, T.Grp3B=T.Grp3B, T.Grp2A=T.Grp2A, T.Grp4=T.Grp4, T.Grp2B=T.Grp2B, T.Grp1=T.Grp1, T.Grp5=T.Grp5, 
                       T.UNK_1=T.UNK_1, T.Grp3A_1=T.Grp3A_1, T.Grp3B_1=T.Grp3B_1, T.Grp2A_1=T.Grp2A_1, T.Grp4_1=T.Grp4_1, T.Grp2B_1=T.Grp2B_1, T.Grp1_1=T.Grp1_1, T.Grp5_1=T.Grp5_1) #MAF=MAF took out 
  
  
  
  
  # S U M M A R Y  S T A T I S T I C S from simulated SNP data
  message(paste("Start summary stats"), appendLF=TRUE)
  
  # Statistics measuring levels of genetic diversity within populations:
  popids.post <- c("LA", "Grp3A", "Grp3B", "Grp2A", "Grp2B", "Grp4", "Grp1", "Grp5")
  
  # add population names
  #fscout$deme <- popids.post[as.numeric(fscout$deme)]
  for (i in 1:(length(unique(fscout$deme)))) {
    fscout$deme <- gsub(paste0("Deme.", i), popids.post[i], fscout$deme)
  }
  
  # convert fscout dataframe to other formats (gtypes, genind)
  tic("convert to gtypes")
  fsc.gtypes <- df2gtypes(fscout, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3) 
  toc() # convert to gtypes
  
  # minor allele frequency
  tic("mafreq")
  allMAF <- mafreq(fscout)
  toc() # mafreq
  
  # heterozygosity per locus
  tic("totalHe")
  totalHe <- 2*allMAF*(1-allMAF)
  toc() # "totalHe"
  
  # total number of SNPs
  SNPs <- sum(allMAF<1)
  names(SNPs) <- "tot_SNPs"
  
  
  # split fscout by deme
  split_out <- vector("list", length(unique(fscout$deme)))
  for(p in 1:length(unique(fscout$deme))) {
    split_out[[p]] <- fscout[fscout$deme == unique(fscout$deme)[p],]
    names(split_out)[p] <- unique(split_out[[p]]$deme)
  }
  
  # id minor allele for ea locus
  minor <- sapply(names(allMAF), FUN=function(x){strsplit(x, split = "[.]")[[1]][2]})
  
  # minor allele freq per locus per deme
  locMAF <- loc.mafreq(split_out, minor) 
  
  # Frequency-weighted marker values (see Schonswetter & Tribsch 2005) for each deme
  tic("DW")
  DW <- colSums(locMAF/allMAF)
  names(DW) <- paste0("DW_",names(DW))
  toc() # DW
  
  # sample size per deme
  locN <- sapply(split_out,function(o){nrow(o)})
  
  # number of SNPs per deme
  localSNP <- apply(locMAF,2,function(x){sum(x<1 & x>0)})
  names(localSNP) <- paste0("S.", names(localSNP))
  
  # allelic richness per deme
  tic("allelic richness")
  richnessSNP <- allelicRichness(fsc.gtypes, by.strata=TRUE)
  richnessSNP <- richnessSNP %>% dplyr::group_by(stratum) %>% dplyr::summarize(mean=mean(allelic.richness, na.rm=TRUE)) # added na.rm=T 7/10/23
  richnessSNP <- as.data.frame(t(richnessSNP))
  names(richnessSNP) <- paste0("rS.", richnessSNP[1,])
  richnessSNP <- richnessSNP[-1,]
  toc() # allelic richness
  
  # private alleles per deme
  tic("private alleles")
  privateSNP <- colSums(privateAlleles(fsc.gtypes))
  privateSNP <- privateSNP[match(unique(fscout$deme), names(privateSNP))]
  names(privateSNP) <- paste0("pS.", names(privateSNP))
  toc() # "private alleles"
  
  # total private alleles across all demes
  total_priv = sum(privateSNP)
  names(total_priv) <- "tot_priv"
  
  # Pairwise private SNPs and their frequencies
  tic("pw private snps")
  pops <- colnames(locMAF)
  combos <- combn(pops,2)
  pwpriv_names <- c()
  pwpriv <- c()
  for(pair in 1:length(combos[1,])) {
    pwpriv_names <- c(pwpriv_names, paste0("pwpriv_",combos[1,pair],".",combos[2,pair]), paste0("pwpriv_", combos[2,pair],".",combos[1,pair]))
    tmp_locMAF <- locMAF[,c(combos[1,pair], combos[2,pair])]
    tmp_locMAF <- tmp_locMAF[rowSums(tmp_locMAF > 0),]
    #Sites that are polymorphic in pop 1, but not 2
    pwpriv <- c(pwpriv, sum(tmp_locMAF[,1] > 0 & tmp_locMAF[,2] == 0))
    #Sites that are polymorphic in pop 2, but not 1
    pwpriv <- c(pwpriv, sum(tmp_locMAF[,2] > 0 & tmp_locMAF[,1] == 0))
    rm(tmp_locMAF)
  }
  names(pwpriv) <- pwpriv_names
  toc() # pw private snps
  
  # expected heterozygosity per deme
  tic("locHe")
  locHe <- colMeans(2*locMAF*(1-locMAF))
  names(locHe) <- paste0("he.", names(locHe)) # added 7/11/23
  toc() # locHe
  
  # variation in expected het per deme
  tic("varlocHe")
  varlocHe <- apply(2*locMAF*(1-locMAF),2,var)
  names(varlocHe) <- paste0("vHe.", names(varlocHe)) # added 7/11/23
  toc() # varlocHe
  
  # combine expected het and variation
  hedf <- data.frame(he=locHe, pop=paste0(names(locHe)),stringsAsFactors=F)
  
  # Fis (based on heterozygosity; 1 - (Ho/He))
  ho <- heterozygosity(fsc.gtypes, type="observed", by.strata=TRUE)
  he <- heterozygosity(fsc.gtypes, type="expected", by.strata=TRUE)
  hohe <- merge(ho,he)
  loc_Fis <- hohe %>% dplyr::group_by(stratum) %>% dplyr::mutate(loc_Fis=1-(obsvd.het/exptd.het))
  FisA <- loc_Fis %>% dplyr::group_by(stratum) %>% dplyr::summarize(meanFis=mean(loc_Fis, na.rm=TRUE))
  Fis <- as.vector(FisA$meanFis)
  names(Fis) <- paste0("Fis_", FisA$stratum)
  
  
  message(paste("Finished summary stats on populations"), appendLF=TRUE)
  
  
  
  message(paste("Start bt pop summary stats"), appendLF=TRUE)
  
  # Statistics measuring levels of genetic differentiation between populations:
  # pairwise expected heterozygosity
  cores=1
  tic("pwhet")
  pwhet <- pwise.het(locMAF,locN,cores)
  toc() #pwhet
  
  # pairwise Fst (based on Nei 1987)
  tic("pwFst")
  pwFst <- pairwise.neifst(gtypes2genind(fsc.gtypes))
  fst_names <- t(combn(colnames(pwFst), 2))
  pwFst.df <- data.frame(fst_names, fst=pwFst[fst_names])
  pwFst.df$label <- paste0("Fst_",pwFst.df$X1, ".", pwFst.df$X2)
  pairFst <- as.vector(pwFst.df$fst)
  names(pairFst) <- pwFst.df$label
  toc() #pwFst
  
  # format pairwise Fst for use later in ibd test
  fsts <- cbind(pwFst.df, data.frame(t(sapply(strsplit(names(pairFst),"\\."),function(nms){c(from=strsplit(nms[1],"_")[[1]][2],to=nms[2])})),stringsAsFactors=F))
  rownames(fsts) <- pwFst.df$label
  fsts <- fsts[, c("fst", "from", "to")]
  fsts <- fsts[order(fsts$to,fsts$from),]
  
  
  # Test for isolation by distance (IBD)
  rsc.latlong <- read.csv("/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/SRC/RSC_clusters_median_lat-long.csv")
  
  # geographic distance
  tic("geog dist")
  pdist=as.matrix(dist(rsc.latlong[,c("med.lat","med.long")]))
  colnames(pdist) <- rsc.latlong$cluster
  rownames(pdist) <- colnames(pdist)
  diag(pdist) <- NA
  pdist[upper.tri(pdist)] <- NA
  dsts = as.data.frame(as.table(pdist))
  dsts = dsts[complete.cases(dsts),]
  names(dsts) <- c("to","from","d")
  dsts$from <- as.character(dsts$from)
  dsts$to <- as.character(dsts$to)
  dsts <- dsts[order(dsts$to,dsts$from),]
  toc() # geog dist
  
  # merge geographic distance and fsts (Nei 1987)
  dsts <- merge(dsts,fsts)
  
  # linear regression for linearized Fst vs geographic distance
  tic("ibd")
  IBDfst <- lm((fst/(1-fst))~log(d),dsts) #changed from the previous line to implement Rousset's (1997) version
  ibdfst.slope <- c(coef(IBDfst)[2])
  ibdfst.int <- c(coef(IBDfst)[1])
  bsfst <- segmentGLM(c(dsts$d),log(c(dsts$fst+1)))
  toc() #ibd
  
  # directionality index
  tic("psi")
  psi <- psiCalc(locMAF, samplen=nsamples)
  toc() # psi
  
  
  message(paste("Finished bt pop summary stats"), appendLF=TRUE)
  
  
  # combine all stats
  stats = c(SNPs, localSNP, privateSNP, total_priv, DW, richnessSNP, locHe, varlocHe,
            ibdfst.slope=ibdfst.slope,ibdfst.int=ibdfst.int,
            pwpriv, pairFst, psi)
  
  # stats left out
  # Fis, pairnei, graphstats, ibdnei.slope=ibdnei.slope,ibdnei.int=ibdnei.int, ibdedist.slope=ibdedist.slope,ibdedist.int=ibdedist.int,bsfst.break=bsfst[1], bsfst.ll=bsfst[2], bsedist.break=bsedist[1], bsedist.ll=bsedist[2], bsnei.break=bsnei[1], bsnei.ll=bsnei[2],
  
  stats1 = matrix(data=stats, nrow = 1)
  colnames(stats1) = names(stats)
  stats2 = as.data.frame(stats1)
  
  # combine params and stats 
  params.stats <- cbind(params, stats)
  params.stats <- data.frame(lapply(params.stats, as.character), stringsAsFactors=FALSE)
  
  
  # Switch to outdir before writing the file
  setwd(outdir)
  
  
  if(!file.exists(fn)) {
    write.table(params.stats, fn, sep = ",", quote = FALSE, row.names = FALSE)
  } else {
    write.table(params.stats, fn, quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
  }
  
  # Then you can add one to num.rep after writing the row to the output table
  num.rep <- num.rep + 1
  
  # You can remove stats before the next replicate begins to avoid writing duplicate rows (e.g. if a coalescent simulation failed)
  rm(stats, params.stats, out)
  
}
toc() # full job tic
