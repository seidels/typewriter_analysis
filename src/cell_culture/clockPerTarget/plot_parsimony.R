## ---------------------------
##
## Script name: calculate parsimony
##
## Purpose of script: calculate parsimony score from logged values in individual barcodes
##
## Author: Antoine Zwaans
##
## Date Created: 2024-03-08
##
## Copyright (c) Antoine Zwaans, 2024
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory for Mac 

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/parsimony")   

library(ggpubr)
library(ggplot2)

mcc_file <- "parsimony_mcc.log"
mcc <- read.table(mcc_file, header = T) 

upgma_file <- "parsimony_upgma.log"
upgma <- read.table(upgma_file, header = T) 

#37666
upgma_parsimony <- sum(upgma[,paste0(rep("treeLikelihood.",13),1:13)])

#37560
mcc_parsimony <- sum(mcc[,paste0(rep("treeLikelihood.",13),1:13)])

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/individual_trees")   

# save the trees
trees <- ape::read.nexus(file = "thinned4000000.trees")
trees <- sciphy_trees
for(i in 1:length(trees)) {
  name <- paste0("sciphy_tree.",i,".txt")
  ape::write.tree(trees[i], file=name)
}

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_OU/parsimony")   

sciphy_parsimony <-c()
for(i in 1:649) {
  file <- paste0("typewriter_13Sites_1000Cells_DataSet1_PARSIMONY.",i,".log")
  trace  <- read.table(file, header = T) 
  sciphy_parsimony <- c(sciphy_parsimony,sum(trace[,paste0(rep("treeLikelihood.",13),1:13)]))
  
}


plot(sciphy_parsimony)

parsimony_df <- data.frame(number=c(1:649,324.5,324.5),Tree = c(rep("SciPhy",649),"SciPhy MCC","UPGMA"), parsimony_score = c(sciphy_parsimony,mcc_parsimony,upgma_parsimony),size=c(rep(1,649),2,2))

plot <- ggplot(data=parsimony_df,aes(x=number,y=parsimony_score,colour=Tree,size=size)) + geom_point() + ylab("Parsimony score") + xlab("")  + theme_bw() + guides( size = "none") + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
plot
ggsave("parsimony_scores.png", plot, width = 30, height = 30, units = "cm", dpi = 300)

sciphy_trees <- ape::read.nexus(file = "~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/thinned4000000.trees")

sample_nr_tree <- as.numeric(unlist(strsplit(names(sciphy_trees),"_"))[seq(2,2*length(sciphy_trees),by=2)])

## Load log data and extract likelihood values corresponding to the matching sample nrs
log <- read.table("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells/combined.log", header = T)

step_tree <- sample_nr_tree[2] - sample_nr_tree[1]
step_log <- log$Sample[2] - log$Sample[1]

#resample the log file at the same frequency
subsampled_log <- log[seq(1,length(log$Sample),by=step_tree/step_log),]

#resample the same number
subsampled_log <- subsampled_log[1:min(length(sciphy_trees),length(subsampled_log$Sample)),]

#check that all tree sample nrs and likelihood sample nrs match
which(subsampled_log$Sample != sample_nr_tree)

#extract likelihood values
tree_likelihood <- subsampled_log$likelihood


parsimony_likelihood_df <- data.frame( parsimony_score = sciphy_parsimony,likelihood = tree_likelihood)


plot_like_pars <- ggplot(data=parsimony_likelihood_df,aes(x=parsimony_score,y=likelihood)) + geom_point() + ylab("Log-likelihood") + xlab("Parsimony score")  + theme_bw()  + sm_statCorr()
ggsave("parsimony_likelihood_corr.png", plot_like_pars, width = 30, height = 30, units = "cm", dpi = 300)
