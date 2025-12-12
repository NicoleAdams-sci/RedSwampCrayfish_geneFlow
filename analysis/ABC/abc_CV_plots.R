### Plot ABC cross validataion - Figure S7 ###
### 11/20/2023 ###

library(tidyverse)

setwd("~/Documents/crayfish_lab/RedSwampCrayfish_geneFlow/")

# read in data already reformated and combined
mn.nn <- read.csv("output/ABC_output/rsc_100k_CV_combo_new.csv")

# melt dataframes
mm <- mn.nn %>% filter(method == "mnlog")
mm2 <- mm %>% dplyr::select(-c(method, tol))

mm.m <- mm2 %>% dplyr::select(-c(node, source_file)) %>% pivot_longer(!true.model, names_to = "model", values_to = "postP")
mm.m$method <- "mnlog"

nn <- mn.nn %>% filter(method == "nnet")
nn2 <- nn %>% dplyr::select(-c(method, tol))

nn.m <- nn2 %>% dplyr::select(-c(node, source_file)) %>% pivot_longer(!true.model, names_to = "model", values_to = "postP")
nn.m$method <- "nnet"


mn.nn.m <- rbind(mm.m, nn.m)

# fix name of 5cols bc R doesn't like numbers first
mn.nn.m$true.model <- gsub("5cols_step_x", "cols5_step_x", mn.nn.m$true.model)
mn.nn.m$model <- gsub("X5cols_step_x", "cols5_step_x", mn.nn.m$model)

# clean up model names - remove "mROW"
mn.nn.m <- mn.nn.m %>%
  mutate(
    true.model = str_remove(true.model, "_mROW"),
    model = str_remove(model, "_mROW")
  )

mn.nn.m$model <- factor(mn.nn.m$model, levels = c("cols5_step_x", "Grp2_bridge", "Grp2_bridge_indGrp1_5", "Grp2_step",
                                                  "Grp2_step_indGrp1_5", "Grp3_bridge", "Grp3_bridge_indGrp1_5", "Grp3_step",
                                                  "Grp3_step_indGrp1_5", "W2E_step_x"))
mn.nn.m$true.model <- factor(mn.nn.m$true.model, levels = c("cols5_step_x", "Grp2_bridge", "Grp2_bridge_indGrp1_5", "Grp2_step",
                                                            "Grp2_step_indGrp1_5", "Grp3_bridge", "Grp3_bridge_indGrp1_5", "Grp3_step",
                                                            "Grp3_step_indGrp1_5", "W2E_step_x"))


# Add an indication of the correct model 
cv.all.p2 <-  ggplot(mn.nn.m, aes(x=model, y=postP, fill=method)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.4, end = 0.85) +
  facet_wrap(~true.model) +
  geom_rect(data= mn.nn.m %>% filter(model == as.character(true.model)), 
            aes(xmin = as.integer(true.model)-0.5, xmax = as.integer(true.model)+0.5, ymin = -0.01, ymax = 1.00),
            color = "red", fill = NA, linetype = "dashed", alpha = 0.5) +
  ylab("Posterior probability") +
  xlab("Models tested") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

today <- format(Sys.time(), format = "%Y-%m-%d")
filename <- paste0("~/Documents/crayfish_lab/RedSwampCrayfish_geneFlow/figures/figureS7_CV_Allwboxes.plots_newMods_", today, ".png")

#ggsave(cv.all.p2, file=filename, width = 11, height =8, bg = "white")

  
  
   
  