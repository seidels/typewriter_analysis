# load libs
library(ggplot2)

# Define the file names
file_path_start <- "~/frameworks/new_beasts/beast_typewriter/typewriter_data_cut_tree_"
file_path_end <- "_seed41.clock_0.1.alignment.nexus"

# Function to count the number of zeros in a number
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
                           SD = numeric())

for (time in cut_times) {
  file_path <- paste0(file_path_start, time, file_path_end)
  edit_count <- count_edits_from_alignment(file_path)
  mean_val <- mean(edit_count)
  sd_val <- sd(edit_count)
  
  result_table <- rbind(result_table, data.frame(Time = time, Mean = mean_val, SD = sd_val))
}

# add edit count 0 at time 0 (i.e age 25)
result_table[5, ] = c(25, 0, 0)

# Plot edit counts through time
g <- ggplot(result_table, aes(x = 25-Time, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  labs(x = "Time", y = "Edit count") +
  ggtitle("Edit counts across time") +
  theme_minimal()

g

ggsave(plot = g, filename = "results/simulation/simulate_edit_saturation/edit_count_through_time.jpg")
