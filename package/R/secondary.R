#'Preprocessing of FCM files.
#'
#'\code{remove_noise} is the main function to preprocess a set of FCM
#'files contained in one single folder (usually corresponding to a single plate)
#'after they have been filtered on their scattering profile and their
#'fluorescence has been fitted with a lognormal distribution through the
#'function \code{scattering_filter}.
#'
#'
#'This function is called by \code{remove_noise} after having filtered the
#'files based on their scattering profile; usually the user shouldn't directly
#'call this function. It takes all the FCM files in a given folder whose cells
#'have already been filtered based on their scattering profile and whose results
#'have been saved in a Rdata structure with the name specified by
#'\code{.filter_preproc_namer} (this is done with the function
#'\code{scattering_filter}). Specifically, it merges the filtered files found in
#'the directory together (possibly in a parallel way) and it removes the
#'autofluorescence and the shot noise.
#'
#'@param .dir Directory to look for the files to preprocess
#'@param .use_doParallel If it is set to TRUE then the merging of the results from
#'  all the files in the folder is done parallely on the cores of the local
#'  machine. Please note that this refers only to the final step of the
#'  analysis, when the results of the filtering of the files in the folder must
#'  be merged together; the filtering step is \emph{always} peformed on the
#'  cluster, even when .use_doParallel is set to FALSE. Also note that this
#'  parallelization on the local machine requires the package \emph{doParallel}
#'  and that the user sets up the parallel environment. If the package
#'  doParallel is present and this function is called through
#'  \code{preproc_facs_plates} then the parallel environment is automatically
#'  set to half of the available cores.
#'@return A list containing \enumerate{ \item \emph{stats} the summary
#'  statistics of all the FCM analized. \item \emph{preproc} the preprocessed
#'  data, i.e., raw scattering and fluorescence values of the cells which
#'  survived the scattering filter. If \emph{.drop_raw} is FALSE, then this
#'  item is NULL.} The function also saves the final result in the file
#'  specified by \code{.cache_namer} in the folder \code{.out_dir} so that it
#'  can be restored later if required. Also the function writes the statistiscs
#'  in a CSV file.
make_distribution_plots_dir <- function(.out_dir, .pdf_dim)
{
  # Get the fcs files to be analysed
  .files <- list.files(.out_dir, '.RData', full.names = TRUE)
  
  # Load all the stats and preproc (df with scattering, fluo and posterior for the fluo) from Rdata files.
  # Store them in a single df called .total
  .preproc <- .files %>% purrr::map_dfr(function(.x){
    load(.x)
    return(.gfp_stats$preproc)
  })
  .stats <- .files %>% purrr::map_dfr(function(.x){
    load(.x)
    return(.gfp_stats$stats)
  })
  .total <- left_join(.preproc, .stats, by = c('path', 'well')) %>% mutate(base = basename(path))
  
  # Remove unecessary df
  rm(.stats, .preproc)
  
  #Create the fits lines
  .fit <- .total %>% group_by(path, well) %>% summarise(MAX = max(log(GFP.H)), MIN = min(log(GFP.H)), mu = unique(fl_mean), s = unique(fl_var), w = unique(w)) %>% 
    group_by_all %>% do(data.frame(x = seq(.$MIN, .$MAX, length.out = 100))) %>% 
    mutate(y = w*dnorm(x, mu, sqrt(s)) + (1-w)/(MAX-MIN)) %>% mutate(base = basename(path))
  
  #Make the plots in a pdf
  pdf(file = file.path(.out_dir, paste0(basename(.out_dir), '_plots.pdf')), height=.pdf_dim[2]*.pdf_dim[3], width=.pdf_dim[1]*.pdf_dim[4])
  par(mfrow=c(.pdf_dim[3],.pdf_dim[4]))
  par(mar=c(5.5, 2.5, 2, .5) + 0.1,
      mgp=c(1.5, .5, 0),
      tcl=-0.25)
  
  print(
    ggplot() + 
    geom_histogram(data = .total, aes(log(GFP.H), y=..density..), bins = 30) + 
    geom_line(data = .fit, aes(x, y), col = 'red') + 
    facet_wrap(base~well)
  )
  
  dev.off()
}

#'Filtering of FCM files based on their scattering profile.
#'
#'\code{scattering_filter} filters all the FCM files in a set of folders based
#'on their scattering profile and it fits the survived cells with a mixture of
#'lognormal and uniform distribution. It stores the final result for each FCM
#'file in the directory specified by \code{.data2preproc}, such that they can be
#'loaded for subsequent analysis (e.g. by \code{remove_noise}).
#'
#'This function loops over all the dirs provided and check for FCM files that
#'have not yet been filtered. The filtered information are stored in a path
#'given by the function .data2preproc and with a filename given by
#'.filter_preproc_namer. More specificaly, what the function does is to check
#'for every FCM files if the file storing the filtered info is present. If it is
#'not present, or if .force is TRUE, it adds the command to filter the file in
#'the file CMD in the working directory. At the end, it creates a shell job file
#'to execute an array job which is sent to the cluster. IMPORTANT: dirnames must
#'be given as absolute paths or as relative the current working directory.
#'
#'
#'@return The function doesn't return any value, but it saves the statistics of
#'  the filtered data in Rdata files in the folder \code{.data2preproc} with the
#'  name given by \code{.filter_preproc_name}.

scattering_filter <- function(.dirs, .data2preproc, .f_par, .filter_preproc_namer, .preproc_func, .force)
{
  # List of commands to send to the cluster, in SIMD mode. Each line is a different file to be processed
  .cmds <- data.frame()

  # List all the fcs files that have to be analyzed
  .files <- do.call(list.files, list(path=.dirs, pattern=.f_par$file_pattern, full.names=TRUE, recursive=TRUE))
  # If the files are named with the relative path, transform it to absolute path, required to be correctly loaded in the cluster.
  # I don't use normalizePath, because files can be symbolically linked and normalizePath will give the true path and not the linked one
  .files[which(substr(.files, 1, 1)!='/')] <- file.path(getwd(), .files[which(substr(.files, 1, 1)!='/')])
  # Check if everything went good and the paths point to existing files
  if(any(!file.exists(.files))){
    stop("Some of files to be processed do not exist. Check that the directories are given as absoulte paths or as relative to the current working directory", call.=TRUE)
  }

  .cmds <- tibble::enframe(.files, name = NULL, value='file') %>% mutate(
    out = .filter_preproc_namer(.data2preproc(file)),
    fsc1=.f_par$channels[1], fsc2=.f_par$channels[2], 
    ssc1=.f_par$channels[3], ssc2=.f_par$channels[4],
    gfp=.f_par$channels[5], 
    force=.force,
    threshold=.f_par$scattering_threshold, scattering_frac_cells=.f_par$scattering_frac_cells
  )
  
  #Create the paths if they don't exist. I use purrr because creating a dir is quick, I don't want to parallelize it
  .cmds %>% mutate(out_dir = dirname(out)) %>% group_by(out_dir) %>% summarise %>% 
    rename(path = out_dir) %>% purrr::pmap(dir.create, showWarnings=FALSE, recursive=TRUE)
  
  #Filter the fcs files
  if(is.null(.preproc_func)) .preproc_func <- preproc
  .stats <- .cmds %>% ungroup() %>% mutate(tmp = seq(n())) %>% split(.$tmp) %>% future_map_dfr(preproc, .progress = TRUE) %>% as_tibble()
  
  return(.stats)
}
