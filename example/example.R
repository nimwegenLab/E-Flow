library(dplyr)
library(ggplot2)
library(EFlow)
library(furrr)
library(stringr)
library(gridExtra)

# store processed results in the same folder as the data, 
# but with "data" substituted by "data_proc"
data2preproc <- function(.d) sub('data', 'data_proc', .d)

# create info data frame (one line per plate)
pl_index_file <- "index_plates.csv"
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#")
pl_info <- pl_index %>% select(dir, empty, condition, date)

# create the structure containing the information for processing the FCM data
f_par <- list(
  channels = c('FSC-H', 'FSC-W', 'SSC-H', 'SSC-W', 'GFP-H'),
  scattering_threshold=0.5,
  scattering_frac_cells=NA,
  file_pattern = "fcs$",
  delta_shot_noise = 12.73,
  delta_shot_noise.err = 0.63
)

# we can test whether the channels have been correctly specified.
# If a channel is not present, the following function returns an error.
check.channels(pl_info$dir, f_par)

#Process the data
plan(multiprocess, workers=5)
processed <- analyse_raw(pl_info$dir, data2preproc, .f_par=f_par, .force = FALSE) %>% 
  mutate(path = str_match(path, './data/Condition[0-9]/.+$')[,1]) %>%  #Required to correctly join pl_info
  propagate_index_info(pl_info) %>% 
  mutate(plasmid=str_match(basename(as.character(path)), '[:alnum:]+')[,1],
         condition = paste('Condition:', condition))  
plan(sequential)

#Make histograms of fitted fluo distribution
make_distribution_plots(pl_info$dir, data2preproc)

#Analyze the autofluorescence
autofluo <- processed %>% dplyr::filter(grepl('puA', plasmid)) %>% group_by(date, condition) %>% do((function(.df){
  .m <- get_autofluo_stats(.df$fl_mean_lin, .df$fl_mean_lin.err)
  .v <- get_autofluo_stats(.df$fl_var_lin, .df$fl_var_lin.err)
  return(data.frame(autofluo_mean = .m['mu'], autofluo_mean.err = .m['tau'], autofluo_var = .v['mu'], autofluo_var.err = .v['tau']))
})(.)) 
processed <- processed %>% left_join(autofluo)

#Plot autofluo estimates
p1 <- processed %>% dplyr::filter(grepl('puA', plasmid)) %>% 
  ggplot() + 
  geom_point(aes(well, fl_mean_lin, col=plasmid)) + 
  geom_errorbar(aes(x = well, ymin=fl_mean_lin-fl_mean_lin.err, ymax = fl_mean_lin+fl_mean_lin.err, col = plasmid)) + 
  geom_hline(data = autofluo, aes(yintercept = autofluo_mean)) +
  geom_rect(data = autofluo, aes(xmin = 0.5, xmax = 2.5, ymin = autofluo_mean-autofluo_mean.err, ymax = autofluo_mean+autofluo_mean.err), col=NA, alpha = 0.35) +
  theme_bw() + xlab('Well') + ylab('Mean GFP [A.U.]') + ggtitle('Mean autofluorescence') + facet_wrap(~condition)
p2 <- processed %>% dplyr::filter(grepl('puA', plasmid)) %>% 
  ggplot() + 
  geom_point(aes(well, fl_var_lin, col=plasmid)) + 
  geom_errorbar(aes(x = well, ymin=fl_var_lin-fl_var_lin.err, ymax = fl_var_lin+fl_var_lin.err, col = plasmid)) + 
  geom_hline(data = autofluo, aes(yintercept = autofluo_var)) +
  geom_rect(data = autofluo, aes(xmin = 0.5, xmax = 2.5, ymin = autofluo_var-autofluo_var.err, ymax = autofluo_var+autofluo_var.err), col=NA, alpha = 0.35) +
  theme_bw() + xlab('Well') + ylab('var GFP [A.U.]') + ggtitle('Var autofluorescence') + facet_wrap(~condition)
grid.arrange(p1, p2)

#Remove autofluo and noise
processed.noise_rm <- remove_noise(processed, f_par) 

#Plot the final results
processed.noise_rm %>% dplyr::filter(fl_mean_lin_autofluo_rm>exp(4)) %>% 
  ggplot() + 
  geom_point(aes(log(fl_mean_lin_autofluo_rm), fl_var_lin_autofluo_rm/fl_mean_lin_autofluo_rm^2, col='Only autfluo removed')) + 
  geom_point(aes(log(fl_mean_lin_autofluo_rm), fl_var_lin_noise_rm/fl_mean_lin_autofluo_rm^2, col = 'Autofluo + shot noise removed')) +
  theme_bw() + xlab("Log mean") + ylab('CV^2')
