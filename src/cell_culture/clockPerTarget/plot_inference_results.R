## ---------------------------
##
## Script name: plot_inference_results
##
## Purpose of script: Generating violin plots from SciPhy estimates on HEK293 cell culture data 
##
## Author: Antoine Zwaans
##
## Date Created: 2023-07-19
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory where the log files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target/1000_cells")   

## load up the packages we will need:  

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)

#load the combined log file 

typewriter_file <- "combined.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)

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
insert_probs_long <- pivot_longer(insert_probs,seq(1,ncol(insert_probs)))
insert_probs_long$name <- factor(insert_probs_long$name, levels = names(insert_probs)) 

p_inserts <-  ggplot(insert_probs_long,aes(x=name,value,fill=name)) +
              geom_violin() + 
              xlab("Insertion probabilities") + 
              theme_bw() + 
              ylab("Value") + 
              theme(legend.position = "none" ) + 
              scale_fill_manual(values=c(rep("#4545A8",19),"#E1E1F7")) +
              ylim(c(0,0.2))

ggsave("insert_probs.png", p_inserts, width = 30, height = 15, units = "cm", dpi = 300)

# --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate <- typewriter[,c(29:41)]

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

p_clock <- ggplot(clock_rate_long,aes(x=name,value,fill=name)) +
           theme_bw() +
           geom_violin(draw_quantiles =  c(0.5)) + 
           xlab("Clock rates") + 
           ylab("Value") + 
           theme(legend.position = "none") + 
           scale_fill_manual(values=c(rep("#4545A8",13),"#E1E1F7")) + 
           ylim(c(0,0.4))

ggsave("clock_rate.png", p_clock, width = 30, height = 15, units = "cm", dpi = 300)

# --------------------
#  plot the growth rate
# -------------------- 

#extract and substract the birth and death rates 
bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1) - rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))

#order columns 
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))
p_growth <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) + 
            theme_bw() + 
            geom_violin(draw_quantiles = c(0.5)) + 
            theme(legend.position = "none" ) + 
            xlab("Growth rate") + 
            ylab("Value") + 
            ylim(c(-0.2,1.5)) +
            scale_fill_manual(values=c("#4545A8","#E1E1F7"))

ggsave("growth_rate.png", p_growth, width = 30, height = 30, units = "cm", dpi = 300)


# ---------------------------
#  plot all estimates aligned
# ---------------------------
combined <- cowplot::plot_grid(p_inserts + theme(axis.text.x = element_text(angle = 90)),p_clock + theme(axis.title.y = element_blank(),axis.text.x = element_text(angle = 90)),p_growth +theme(axis.title.y = element_blank()),nrow = 1,rel_widths = c(2,2,1))
ggsave("combined_estimates.png", combined, width = 30, height = 10, units = "cm", dpi = 300)
