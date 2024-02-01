library(ggplot2)

source("~/Projects/typewriter_analysis/src/cell_culture/understand_growth/functions.R")


#traj_file = "~/frameworks/beasts2.7/beast_typewriter/analyses/traj_per_state_3.traj"
traj_file = "~/frameworks/beasts2.7/beast_typewriter/analyses/traj_per_state.traj"

# traj_text = read_file(traj_file)
#
# traj_parsed = parseTrajectory(traj_text)
#
# head(traj_parsed)
# traj_parsed$N
#
# dat = data.frame(time=traj_parsed$time, N = traj_parsed$N[,1])
# library(ggplot2)
#
# p = ggplot(dat, aes(x=time, y=N))+
#   geom_point()
# #p
# #ggsave(p, filename="~/Projects/typewriter_analysis/src/cell_culture/understand_growth/plot_N_50_particles.pdf")
#
# #ggsave(p, filename="~/Projects/typewriter_analysis/src/cell_culture/understand_growth/plot_N_50_particles.jpg")
# ggsave(p, filename="~/Projects/typewriter_analysis/src/cell_culture/understand_growth/plot_N_50_particles_e7thinned.jpg")
# p
#
# head(traj_parsed
#      )

###
traj_loaded = loadTrajectories(filename = traj_file, burninFrac = 0)
states = traj_loaded$states

states = readRDS(file = "~/Projects/typewriter_analysis/src/cell_culture/understand_growth/state_files/states_bdsky.100.rds")
states = readRDS(file = "/Volumes/stadler/People/Sophie_Antoine_shared/cell_culture/trajectory_inference/given_bds_model/trajectory_data/states.100.rds")

#p = ggplot(states, aes(x=time, y=N))+
#  geom_point()
#p

#grids = gridTrajectories(trajStates = states, times = seq(0,25))
times = seq(0, 25)

grids = states %>%
  summarize(N=approx(time, N, times, method="constant", f=1, yleft=0, rule=2)$y,
            time=times)

rm(states)

N_summarised = grids %>%
  group_by(time) %>%
  summarize(
    mean_population = mean(N),
    median_population = median(N),
    lower95_hpd = HDInterval::hdi(N)[[1]],
    upper95_hpd = HDInterval::hdi(N)[[2]]
            )


p = ggplot(N_summarised, aes(x=time, y=median_population))+
  geom_line()+
  geom_line(aes(y=lower95_hpd), linetype="dashed")+
  geom_line(aes(y=upper95_hpd), linetype="dashed")+
  geom_ribbon(aes(ymin=lower95_hpd, ymax=upper95_hpd), fill="blue",
              alpha=0.1)+
  #scale_y_continuous(trans = "log10") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000),
                     labels = expression(1, 10, 10^2, 10^3, 10^4, 10^5, 10^6),
                     expand = c(0, 0))+
  geom_point(aes(x=25,y = 1200000), col="#B00020", size = 4, shape=3)+
  theme_minimal() +
  ylab(label = "Inferred number of cells")+
  xlab(label = "Time")

p
#ggsave(p, width = 10, height = 7, units = "cm",
#       filename="~/Projects/typewriter_analysis/src/cell_culture/understand_growth/plot_N_post_particles_100_thin_e7.jpg")

ggsave(p, width = 10, height = 7, units = "cm",
       filename="~/Projects/typewriter_analysis/src/cell_culture/understand_growth/plot_N_bdsky_post_particles_100_thin_e7.jpg")
