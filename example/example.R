library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(vngFCM)

# store cache file in preproc subdir
data2preproc <- function(.d) sub('data_gwendoline', 'data_gwendoline_preproc', .d)

# create info data frame (one line per plate)
pl_index_file <- "index_plates.csv"
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#")
pl_info <- pl_index %>% select(dir, empty, cip, date)

f_par <- list(
  channels = c('FSC-H', 'FSC-W', 'SSC-H', 'SSC-W', 'GFP-H', 'GFP-W'),
  scattering_threshold=0.5,
  scattering_num_cells=NA,
  file_pattern = "*fcs$", #puA139|puA66|recA|highexp|medexp|polB|recN|ruvA|ftsK|rpsU", #Files to process
  delta_shot_noise = 14.51,
  autofluo_wells = pl_info #Tibble containing for every dir the names of the empty plasmids wells
)


# Call the function to preprocess the data. This must be called twice: the first time to generate the filtered
# data sending the jobs to the cluster and the second time to merge the results from the cluster
preproc <- preproc_facs_plates(pl_info$dir, data2preproc, .f_par=f_par, .plot=TRUE, .verbose=1, .force=FALSE, .drop_preproc=TRUE, .use_doParallel=FALSE)

# Info we want to merge
info <- dplyr::select(pl_info, dir, date, cip)
# Merge with the preproc statistics
preproc <- propagate_index_info(preproc$stats, info)
# Remove columns that we don't need anymore and transform to a tibble
preproc <-  preproc %>% select(-dir) %>% as_data_frame()
# From the file name extract the name of the promoter
preproc <- preproc %>% extract_('path', 'promoter', '[0-9]/([[:alpha:]]+)_[A-Z]')

# Do the plot var vs mean
ggplot(preproc) + geom_point(aes(fl_mean_noise_rm, fl_var_noise_rm, col=factor(promoter)))
