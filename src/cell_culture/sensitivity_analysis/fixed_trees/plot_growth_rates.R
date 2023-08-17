## ---------------------------
##
## Script name: plot_inference_results_UPGMA_BDS
##
## Purpose of script: Generating violin plots from phylodyn estimates on fixed UPGMA/MCC tree
##
## Author: Antoine Zwaans
##
## Date Created: 2023-07-19
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

## set working directory where the log files are

setwd("~/typewriter_analysis/results/analysis_cell_culture_data/sensitivity_analysis/fixed_trees/")   

## load up the packages we will need:  

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)


# ---------------------------------------------------------
#  plot the growth rate for the run on the fixed UPGMA tree
# --------------------------------------------------------- 

typewriter_file <- "1000_UPGMA.1691754256522.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)

#extract and substract the birth and death rates 
bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1) - rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))

prior <- bind_cols(name=rep("Prior",1000),value=rlnorm(1000, meanlog = -0.6, sdlog = 1) - rlnorm(1000, meanlog = -2, sdlog = 1) )
bd_rates_long <- rbind(bd_rates_long,prior)

#order columns 
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))
p_growth <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) + 
  theme_bw() + 
  geom_violin(draw_quantiles = c(0.5)) + 
  theme(legend.position = "none" ) + 
  xlab("Growth rate on UPGMA") + 
  ylab("Value") + 
  ylim(c(0.3,0.5)) +
  scale_fill_manual(values=c("red","pink"))

ggsave("growth_rate_UPGMA.png", p_growth, width = 15, height = 30, units = "cm", dpi = 300)

# ----------------------------------------------------------------
#  plot the growth rate for the run on the fixed scaled UPGMA tree
# ---------------------------------------------------------------- 

typewriter_file <- "1000_UPGMA_medianPosteriorHeight.1691754434739.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)

#extract and substract the birth and death rates 
bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1) - rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))
#add more prior data
prior <- bind_cols(name=rep("Prior",10000),value=rlnorm(1000, meanlog = -0.6, sdlog = 1) - rlnorm(1000, meanlog = -2, sdlog = 1) )

bd_rates_long <- rbind(bd_rates_long,prior)

#order columns 
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))


p_growth_median_height <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) + 
  theme_bw() + 
  geom_violin(draw_quantiles = c(0.5)) + 
  theme(legend.position = "none" ) + 
  xlab("Growth rate on scaled UPGMA") + 
  ylab("Value") + 
  ylim(c(0.3,0.5)) +
  scale_fill_manual(values=c("red","pink"))

ggsave("growth_rate_UPGMA_median_height.png", p_growth_median_height, width = 15, height = 30, units = "cm", dpi = 300)

# -------------------------------------------------------
#  plot the growth rate for the run on the fixed MCC tree
# -------------------------------------------------------

typewriter_file <- "1000_MCC.1691759850766.log"
typewriter <- read.table(typewriter_file, header = T) %>% slice_tail(prop = 0.10)

#extract and substract the birth and death rates 
bd_rates <- typewriter[,"birthRate"] - typewriter[,"deathRate"]

#add a prior column
bd_rates <- bind_cols(Rate=bd_rates,Prior= rlnorm(length(bd_rates), meanlog = -0.6, sdlog = 1) - rlnorm(length(bd_rates), meanlog = -2, sdlog = 1) )
bd_rates_long <- pivot_longer(bd_rates,seq(1,ncol(bd_rates)))
#add more prior data
prior <- bind_cols(name=rep("Prior",1000),value=rlnorm(1000, meanlog = -0.6, sdlog = 1) - rlnorm(1000, meanlog = -2, sdlog = 1) )

bd_rates_long <- rbind(bd_rates_long,prior)

#order columns 
bd_rates_long <- mutate(bd_rates_long,name = fct_relevel(name,"Rate","Prior"))


p_growth_MCC <- ggplot(bd_rates_long, aes(x=name,value,fill=name)) + 
  theme_bw() + 
  geom_violin(draw_quantiles = c(0.5)) + 
  theme(legend.position = "none" ) + 
  xlab("Growth rate on MCC") + 
  ylab("Value") + 
  ylim(c(0.3,0.5)) +
  scale_fill_manual(values=c("#4545A8","#E1E1F7"))

ggsave("growth_rate_MCC.png", p_growth_MCC, width = 15, height = 30, units = "cm", dpi = 300)


# ---------------------------
#  plot all estimates aligned
# ---------------------------
combined <- cowplot::plot_grid(p_growth_median_height,p_growth ,p_growth_MCC,nrow = 1)
ggsave("combined_UPGMA_median.png", combined, width = 30, height = 10, units = "cm", dpi = 300)
