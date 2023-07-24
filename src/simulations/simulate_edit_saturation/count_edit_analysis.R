setwd("~/Projects/typewriter_analysis/")

clock_rate <- 0.1 #0.1 0.25 0.5

# load libs
library(ggplot2)

# Define the file names
file_path_start <- "src/simulations/simulate_edit_saturation/typewriter_data_cut_tree_"
file_path_mid <- "_seed"
file_path_end <- paste0(".clock_", clock_rate, ".alignment.nexus")

# Function to count the number of zeros (i.e. unedited sites) in a target
count_zeros <- function(x) {
  num_zeros <- sum(strsplit(as.character(x), "")[[1]] == 0)
  return(num_zeros)
}

count_edits_from_alignment = function(file_path){
  # Read the file
  lines <- readLines(file_path)

  # Drop the first 11 lines and the last line
  lines <- lines[12:(length(lines)-1)]

  # Read the remaining lines as tsv data
  data <- read.table(text = lines, sep = " ", header = FALSE)

  # Apply the function to count non-zero digits for each line
  data$nonzero_digits <- 5 - sapply(data[, ncol(data)], count_zeros)

  return(data$nonzero_digits)
}

# count the edits for each cut time
cut_times = c(5, 10, 15, 20)

result_table <- data.frame(Time = numeric(),
                           Mean = numeric(),
                           SD = numeric(),
                           seed = numeric())

result_table_avg <- data.frame(Time = numeric(),
                               Mean = numeric(),
                               SD = numeric())

for (time in cut_times) {

  mean_edits = rep(0, 20)

  for (seed in 1:20){

    file_path <- paste0(file_path_start, time, file_path_mid, seed, file_path_end)
    edit_count <- count_edits_from_alignment(file_path)


    mean_val <- mean(edit_count)
    mean_edits[seed] = mean_val
    sd_val <- sd(edit_count)

    result_table <- rbind(result_table, data.frame(Time = time, Mean = mean_val, SD = sd_val, seed=seed))
  }

  mean_of_means = mean(mean_edits)
  sd_of_means = sd(mean_edits)
  # make sd = -1 because I do not know yet how to calculate the sd of the mean of means
  result_table_avg = rbind(result_table_avg, data.frame(Time = time, Mean = mean_of_means, SD = sd_of_means))
}

# add edit count 0 at time 0 (i.e age 25)
result_table = rbind(result_table, data.frame(Time=25, Mean=0, SD=0, seed=1:20))
#result_table[5, ] = c(25, 0, 0, 1)
result_table_avg[5, ] = c(25, 0, 0)


# Plot edit counts through time
g <- ggplot(result_table, aes(x = 25-Time, y = Mean, col=as.factor(seed))) +
  geom_line()+
  #geom_point(position = position_dodge(width = 2)) +
  #geom_errorbar(aes(ymin = Mean - SD/2, ymax = Mean + SD/2), width = 0.2,
  #              position=position_dodge(width = 2)) +
  labs(x = "Time", y = "Edit count") +
  ggtitle("Edit counts across time") +
  xlim(0,15)+
  ylim(0,4)+
  theme_minimal()+
  theme(legend.position = "none")
g

ggsave(plot = g, width = 10, height = 10, units = "cm", dpi = 330,
       filename = paste0("results/simulation/simulate_edit_saturation/edit_counts_through_time_line_clock_",
                                   clock_rate, ".jpg"))

g <- ggplot(result_table, aes(x = 25-Time, y = Mean, col=as.factor(seed))) +
  #geom_line()+
  geom_point(position = position_dodge(width = 2)) +
  geom_errorbar(aes(ymin = Mean - SD/2, ymax = Mean + SD/2), width = 0.2,
                position=position_dodge(width = 2)) +
  labs(x = "Time", y = "Edit count") +
  ggtitle("Edit counts across time") +
  xlim(0,15)+
  theme_minimal()+
  theme(legend.position = "none")
g

ggsave(plot = g, width = 10, height = 10, units = "cm", dpi = 330,
       filename = paste0("results/simulation/simulate_edit_saturation/edit_counts_through_time_clock_",
                         clock_rate, ".jpg"))



## average results over seeds
g <- ggplot(result_table_avg, aes(x = 25-Time, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SD/2, ymax = Mean + SD/2), width = 0.2,
                ) +
  labs(x = "Time", y = "Edit count") +
  ggtitle("Edit counts across time") +
  theme_minimal()

g

ggsave(plot = g, width = 10, height = 10, units = "cm", dpi = 330,
       filename = paste0("results/simulation/simulate_edit_saturation/avg_edit_count_through_time_clock_",
                                   clock_rate, ".jpg"))
