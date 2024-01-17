## ---------------------------
##
## Script name: plot_inference_results_annotated
##
## Purpose of script: Plot estimates from SciPhy on HEK293 cell culture data, clock per target, with color annotations matching the tree.
##
## Author: Antoine Zwaans
##
## Date Created: 2023-10-13
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory where the log files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_skyline/")   

## load up the packages we will need:  

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the combined log file 

typewriter_file <- "combined_burnin10.log"
typewriter <- read.table(typewriter_file, header = T) 

# -----------------------------
#  plot insertion probabilities
# -----------------------------

#extract insertion probabilities from the full dataframe
insert_probs <- typewriter[,c(7:25)]
trinucleotides_names <- c("CAT","CCG","GCC","ATA","GAT","ACG","ACA","TCG","TAT","GCT","CTA","TGT","AGA","TAA","CAG","TAG","GAG","ACC","GCG")

#reorder by median
names(insert_probs) <- trinucleotides_names
insert_probs <- insert_probs[order(-sapply(insert_probs, median))]
insert_probs <- bind_cols(insert_probs,Prior=rdirichlet(nrow(insert_probs), rep(1.5,19))[,2])
insert_probs_medians <- sapply(insert_probs,median)
insert_probs_low <- as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.025)}))
insert_probs_up <-  as.numeric(sapply(insert_probs,function(x) {quantile(x, 0.975)}))
datafra <- data_frame(name=c(trinucleotides_names,"Prior"),median=insert_probs_medians,low=insert_probs_low,up=insert_probs_up)


datafra$name <- factor(datafra$name, levels = datafra$name)

#assign standard colors to the trinucleotides alphabetical to match the tree plot. 
colors <- c(hue_pal()(19),"#E1E1F7")
names(colors) <- c(sort(trinucleotides_names),"Prior")
p_inserts <-  ggplot(datafra) +
  geom_bar(aes(x=name,y=median,fill=name),stat="identity",colour="black") +
  xlab("Insert") + 
  theme_bw() + 
  ylab("Posterior insert probability") + 
  theme(legend.position = "none" ) + 
  scale_fill_manual(values=colors) +
  
  geom_errorbar(aes(x=name,ymin=low, ymax=up), width=.2,
                position=position_dodge(.9)) +
  
  coord_cartesian(ylim=c(0,0.175),expand = FALSE) + theme(axis.text.x = element_text(angle = 90),text = element_text(size = 22), 
                                                          panel.grid.minor = element_blank(),
                                                          panel.border = element_blank(),
                                                          panel.background = element_blank(),panel.grid.major.x = element_blank(),axis.text.x.bottom = element_text(size = 22)) 


ggsave("insert_probs.png", p_inserts, width = 40, height = 12, units = "cm", dpi = 1000)

# --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate <- typewriter[,c(33:45)]

#rename by the targetBC
names(clock_rate) <- c("ATGGTAAG","ATTTATAT",
                       "ATTTGGTT", "GCAGGGTG",
                       "GTAAAGAT", "TAGATTTT",
                       "TGCGATTT", "TGGACGAC",
                       "TGGTTTTG", "TTAGATTG",
                       "TTGAGGTG","TTTCGTGA",
                       "TTCACGTA")


#reorder by median
clock_rate <- clock_rate[order(-sapply(clock_rate, median))]
ordered_names <- names(clock_rate)

#add a prior column
clock_rate <- bind_cols(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_long <- pivot_longer(clock_rate,seq(1,ncol(clock_rate)))

#order columns 
clock_rate_long <- mutate(clock_rate_long,name = fct_relevel(name,ordered_names,"Prior"))

p_clock_pos <- ggplot(clock_rate_long,aes(x=name,value,fill=name)) +
  theme_bw() +
  geom_violin(draw_quantiles =  c(0.5)) + 
  xlab("Tape") + 
  ylab(parse(text = paste0('"Posterior editing rate "', '(~ day^-1)'))) +
  
  theme(legend.position = "none") + 
  scale_fill_manual(values=c(rep("#5CA17D",13),"#E1E1F7")) + 
  coord_cartesian(
    ylim=c(0,0.4),expand=FALSE) + theme(text = element_text(size = 22),panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90) )

ggsave("clock_rate.png", p_clock_pos, width = 30, height = 15, units = "cm", dpi = 1000)


# ---------------------------
#  plot all estimates aligned
# ---------------------------

combined <- cowplot::plot_grid(p_inserts+ theme(axis.line = element_line(colour = "black")),p_clock_pos+ theme(axis.line = element_line(colour = "black")),nrow = 1)
ggsave("combined_estimates.png", combined, width = 50, height = 15, units = "cm", dpi = 800)










