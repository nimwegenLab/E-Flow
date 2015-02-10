### PARAMETERS
proj_path <- "~/Documents/Biozentrum/Projects/vngWetLabR/facs/example"
pl_index_file <- "exampleIndex_plates.csv" # index file by plate
w_index_file <- "exampleIndex_wells.csv"  # index file by well
scripts_path <- "~/Documents/Biozentrum/Projects/vngWetLabR/facs/facs_toolbox.R"
data2preproc <- function(.d) sub('/data/', '/preproc/', .d) # store cache file in preproc subdir
# data2preproc <- identity # store cache files with raw data


# FACS related parameters
# can be overriden with a .r file in the data directory (see README)
f_par <- list(
  channels = c(fsc = 1, ssc = 3, # channels to analyze
                fl1 = 5), #, fl2 = 6),
  lims = c(2.5, 4.5), # limits for fsc and ssc (in log10 scale)
  min.cells = 5000,  # minimum number of cells to include in preprocesing
  facs.max = 262143, #
  file.pattern = '.fcs$')
#   od.pattern = '')


### ANALYSIS
setwd(proj_path)
source(scripts_path)

### Load data defined per plate

# create info data frame (one line per plate, with `dir` and `discard` columns)
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#") # with path, and media and/or time
pl_info <- pl_index %>%
  select(-comment) %>%
  extract_('dir', c('condition', 'date'), '.*/([^/]+)/([0-9]{8})/[^/]+$', remove=FALSE)

delete_preproc_files(pl_index$dir)#, .silent=TRUE) 

# create gates and log transform (with interactive control) 
f_utils <- set_fsc_ssc_gates(pl_index$dir[1], f_par)#, .interactive=TRUE)

# facs preprocessing
mypls <- preproc_facs_plates(pl_info$dir, data2preproc, f_par, f_utils, .min_cells=5000,
                             .plot=TRUE, .verbose=2, .force=TRUE)
mypls <- propagate_index_info(mypls, pl_info)

qplot(fsc, ssc, data=mutate(mypls$gates, fc=interaction(condition, date)), 
      col=factor(plate), alpha=I(.8), group=interaction(fc, path, well), geom='path') +
  facet_wrap(~fc)
qplot(fsc, ssc, data=mypls$gates %>% split_well_col, 
      col=interaction(condition, date), group=interaction(condition, date, plate), geom='path', facets=row~col)

mypls$preproc %>%
  group_by(condition, date, plate, well) %>%
  summarise(n_cells = n())


### Load data defined per wells

# create info data frame (one line per well, with `dir` and `well` columns)
w_index <- read.csv2(w_index_file, stringsAsFactors=FALSE, comment.char="#")
.ws <- strsplit(w_index$well, split='[,;] *')
w_info <- w_index %>%
  select(-comment) %>%
  # add date column
  extract_('dir', c('date'), '.*/([0-9]{8})/[^/]+$', remove=FALSE) %>%
  # split well column and create one line for each value
  slice(rep(1:dim(w_index)[1], times=sapply(.ws, length))) %>%
  mutate(well=unlist(.ws))

# create gates and log transform (with interactive control) 
f_utils <- set_fsc_ssc_gates(w_index$dir[1], f_par)#, .interactive=TRUE)

# facs preprocessing
mywells <- preproc_facs_plates(unique(w_index$dir), data2preproc, f_par, f_utils, .min_cells=5000,
                               .plot=TRUE, .verbose=2, .force=TRUE)
mywells <- propagate_index_info(mywells, w_info)

qplot(fsc, ssc, data=mywells$gates, col=well, alpha=I(.8), group=path, geom='path') +
  scale_colour_periodic_brewer() +
  facet_wrap(~strain)


### ANALYSIS
# mean / var plots in different envts
ggplot(mypls$stats, aes(x=gfp_mean_noise_rm, y=gfp_var_noise_rm, col=interaction(condition, date))) +
  geom_point() +
  labs(x='mean(log10(gfp)))', y='var(log10(gfp)))', col='condition')

