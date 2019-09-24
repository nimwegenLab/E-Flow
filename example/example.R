library(dplyr)
library(ggplot2)
library(EFlow)

# store processed results in the same folder as the data, 
# but with "data" substituted by "data_proc"
data2preproc <- function(.d) sub('data', 'data_proc', .d)

# create info data frame (one line per plate)
pl_index_file <- "index_plates.csv"
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#")
pl_info <- pl_index %>% select(dir, empty, condition, date)

# create the structure containing the information for processing the FCM data
f_par <- list(
  channels = c('FSC-H', 'FSC-W', 'SSC-H', 'SSC-W', 'GFP-H', 'GFP-W'),
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
analyse_raw(pl_info$dir, data2preproc, .f_par=f_par, .preproc_script = 'preproc.R')
processed <- collect_raw(pl_info$dir, data2preproc, .f_par = f_par, .drop_raw=FALSE, .plot = TRUE, .use_doParallel = FALSE)
processed$stats <- propagate_index_info(.info = pl_info, .pls = processed$stats)
processed$preproc <- propagate_index_info(.info = pl_info, .pls = processed$preproc)

#Mark empty promoters
processed$stats <- processed$stats %>% mutate(is.empty = well %in% unlist(strsplit(empty, ',')))

#Analyze the autofluorescence 
autofluo <- processed$stats %>% dplyr::filter(is.empty) %>%  group_by(condition, date) %>% do((function(.df){
  .m <- get_autofluo_stats(.df$fl_mean_lin, .df$fl_mean_lin.err)
  .v <- get_autofluo_stats(.df$fl_var_lin, .df$fl_var_lin.err)
  return(data.frame(autofluo_mean = .m['mu'], autofluo_mean.err = .m['tau'], autofluo_var = .v['mu'], autofluo_var.err = .v['tau']))
})(.)) 
processed$stats <- processed$stats %>% left_join(autofluo)

#Plot of mean autofluo and its estimate
processed$stats %>% dplyr::filter(is.empty) %>% 
  ggplot(aes(well, fl_mean_lin)) + facet_wrap(condition~date, nrow=1) +
  geom_point(size=2) + geom_errorbar(aes(ymin=fl_mean_lin-fl_mean_lin.err, ymax = fl_mean_lin+fl_mean_lin.err)) + 
  geom_hline(aes(yintercept = autofluo_mean)) + geom_hline(aes(yintercept = autofluo_mean+autofluo_mean.err)) +
  geom_hline(aes(yintercept = autofluo_mean-autofluo_mean.err)) +
  xlab('Replicate') + ylab('Mean GFP [A.U.]') + ggtitle('Autofluorescence')

#Remove autofluo and noise
processed$stats <- remove_noise(processed$stats, f_par) 

#Plot CV2 vs mean
processed$stats %>% dplyr::filter(!is.empty) %>% 
  ggplot(aes(log(fl_mean_lin_autofluo_rm), fl_var_lin_autofluo_rm/fl_mean_lin_autofluo_rm^2)) + 
  geom_point() + ggtitle('CV vs mean') + xlab('Log mean expression') + ylab('CV2')
