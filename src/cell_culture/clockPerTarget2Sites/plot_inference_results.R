## ---------------------------
##
##
## Purpose of script: Plot inference results for the analysis where only 2 sites were included
##
## Author: Sophie Seidel
##
## Date Created: 2024-03-18
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/")      # Sophie's working directory (mac)

## ---------------------------


## set working directory where the log files are

setwd("~/Projects/typewriter_analysis/results/analysis_cell_culture_data/inference_results/clock_per_target_2_sites/")

## load up the packages we will need:

library(tidyverse)
library(lubridate)
library(dirmult)
library(coda)
library(LaplacesDemon)
library(cowplot)
library(scales)

#load the combined log file

log_file_2_sites <- "analysis_2sites_DataSet1.combined_b30.log"
dat_2_sites <- read.table(log_file_2_sites, header = T)

log_file_all_sites = "../clock_per_target/1000_cells/combined.log"
dat_all_sites = read.table(log_file_all_sites, header = T)


# --------------------
#  plot the clock rates
# --------------------

#extract the clock rate from the tract
clock_rate_all_sites = dat_all_sites[, startsWith(x = colnames(dat_all_sites), prefix = "clock")]
clock_rate_2_sites <- dat_2_sites[, startsWith(x = colnames(dat_2_sites), prefix = "clock")]

#rename by the targetBC
names(clock_rate_2_sites) = names(clock_rate_all_sites) <- c("ATGGTAAG","ATTTATAT",
                       "ATTTGGTT", "GCAGGGTG",
                       "GTAAAGAT", "TAGATTTT",
                       "TGCGATTT", "TGGACGAC",
                       "TGGTTTTG", "TTAGATTG",
                       "TTGAGGTG","TTTCGTGA",
                       "TTCACGTA")


#reorder by median
ord = order(-sapply(clock_rate_all_sites, median))
clock_rate_all_sites <- clock_rate_all_sites[ord]
clock_rate_2_sites <- clock_rate_2_sites[ord]
ordered_names <- names(clock_rate_all_sites)

#add a prior column
#clock_rate <- cbind(clock_rate,Prior=rlnorm(nrow(clock_rate), meanlog = -2, sdlog = 0.5))
clock_rate_all_sites_long <- pivot_longer(clock_rate_all_sites, seq(1,ncol(clock_rate_all_sites)))
clock_rate_all_sites_long$data = "all sites"
clock_rate_2_sites_long <- pivot_longer(clock_rate_2_sites, seq(1,ncol(clock_rate_2_sites)))
clock_rate_2_sites_long$data = "2 sites"

#order columns
clock_rate_all_sites_long <- mutate(clock_rate_all_sites_long,name = fct_relevel(name,ordered_names))
clock_rate_2_sites_long <- mutate(clock_rate_2_sites_long,name = fct_relevel(name,ordered_names))

clock = rbind(clock_rate_all_sites_long, clock_rate_2_sites_long)


p_clock_pos <- ggplot(clock, aes(x=name,value,fill=name)) +
  facet_grid(data~.)+
  theme_bw() +
  geom_violin(draw_quantiles =  c(0.5)) +
  xlab("Tape") +
  ylab(parse(text = paste0('"Posterior editing rate "', '(~ day^-1)'))) +

  theme(legend.position = "none") +
  scale_fill_manual(values=c(rep("#5CA17D",13),"#E1E1F7")) +
  coord_cartesian(
    ,expand=FALSE) + theme(text = element_text(size = 22),panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                           strip.background  = element_rect(colour="black", fill="white"),
                           panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90) )

p_clock_pos


ggsave("compare_clock_2vsall_sites.png", p_clock_pos, width = 30, height = 20, units = "cm", dpi = 1000)









