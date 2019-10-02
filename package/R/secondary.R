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
collect_dir <- function(.dir, .out_dir, .plot, .pdf_dim, .cache_namer, .f_par, .filter_preproc_namer, .drop_raw, .use_doParallel, .force)
{
  cat("\nCurrent directory:", .dir, "\n")
  
  # Check if the analysis has already been performed and saved. If so and .force==FALSE, load the saved analysis
  if (file.exists(.cache_namer(.out_dir)) && !.force) {
    load(.cache_namer(.out_dir)) # load objet .pl
    cat("Cache loaded from ", basename(.cache_namer(.out_dir)), "\n", sep='')
    return (.pl)
  }
  
  # Delete any previous saved analysis
  if(file.exists(.cache_namer(.out_dir)))
    file.remove(.cache_namer(.out_dir))
  
  # Create the directory where to store cached files
  dir.create(.out_dir, recursive=TRUE, showWarnings=FALSE)
  
  if(.plot) {
    pdf(file = file.path(.out_dir, paste0(basename(.out_dir), '_plots.pdf')), height=.pdf_dim[2]*.pdf_dim[3], width=.pdf_dim[1]*.pdf_dim[4])
    par(mfrow=c(.pdf_dim[3],.pdf_dim[4]))
    par(mar=c(5.5, 2.5, 2, .5) + 0.1,
        mgp=c(1.5, .5, 0),
        tcl=-0.25)
  }
  
  # Load the filtered data.
  # Get the fcs files to be analysed
  .files <- file.path(.out_dir, .filter_preproc_namer(list.files(.dir, pattern = .f_par$file_pattern, full.names = FALSE)))
  # Load all the Rdata files
  .preproc <- data.frame()
  .stats <- data.frame()
  
  if(.use_doParallel){
    .pp <- foreach::foreach(f=.files) %dopar% {
      load(f)
      if(.drop_raw)
        .tmp1 <- NULL
      else
        .tmp1 <- .gfp_stats$preproc
      .tmp2 <- .gfp_stats$stats
      
      ##If required, plot the histogram on the pdf file
      if(.plot){
        hist(log(.gfp_stats$preproc$GFP.H), breaks='fd', main=basename(f), prob=TRUE, xlab='Mean log fluo - A.U.', ylab='Density')
        .x <- seq(min(log(.gfp_stats$preproc$GFP.H)), max(log(.gfp_stats$preproc$GFP.H)), length.out=100)
        lines(.x, col='red', lwd=2,
              y=.gfp_stats$stats$w*dnorm(.x, .gfp_stats$stats$fl_mean, sqrt(.gfp_stats$stats$fl_var)) + 
                (1-.gfp_stats$stats$w) * dunif(.x, min(.x), max(.x)))
      }
      
      rm(.gfp_stats)
      return(list(preproc=.tmp1, stats=.tmp2))
    }
    .stats <- do.call(bind_rows, lapply(1:length(.pp), function(i) .pp[[i]]$stats))
    .preproc <- do.call(bind_rows, lapply(1:length(.pp), function(i) .pp[[i]]$preproc))
  }else{
    for(f in .files){
      load(f)
      
      #Keep preproc or not?
      if(.drop_raw) .tmp1 <- NULL
      else .tmp1 <- .gfp_stats$preproc
      
      ##Get the summary statistics
      .tmp2 <- .gfp_stats$stats
      .stats <- rbind(.stats, .tmp2)
      .preproc <- rbind(.preproc, .tmp1)
      
      ##If required, plot the histogram on the pdf file
      if(.plot){
        hist(log(.gfp_stats$preproc$GFP.H), breaks='fd', main=basename(f), prob=TRUE, xlab='Mean log fluo - A.U.', ylab='Density')
        .x <- seq(min(log(.gfp_stats$preproc$GFP.H)), max(log(.gfp_stats$preproc$GFP.H)), length.out=100)
        lines(.x, col='red', lwd=2,
              y=.gfp_stats$stats$w*dnorm(.x, .gfp_stats$stats$fl_mean, sqrt(.gfp_stats$stats$fl_var)) + 
                (1-.gfp_stats$stats$w) * dunif(.x, min(.x), max(.x)))
      }
      
      ##Remove the loaded data
      rm(.gfp_stats)
    }
  }
  if(.plot) dev.off()
  
  ### SAVE THE RESULTS ###
  # Write the final result to a directory
  .pl <- list(preproc=.preproc %>% mutate(dir=.dir), stats=.stats %>% mutate(dir=.dir), method='EFlow')
  save(.pl, file=.cache_namer(.out_dir))
  write.csv (.stats, file=file.path(.out_dir, paste0(basename(.out_dir), '_stats.csv')))
  
  return(.pl)
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

  # For every file check if it has already been filtered, i.e. if the Rdata file specified by .filter_preproc_namer is present.
  # If it is not present, add the command to filter the file to the .cmd_file file.
  for(.f in .files ) {
    # Determine how the file containing the RData of the analyzed fcs file is called
    .out_file <- .filter_preproc_namer(.data2preproc(.f))
    # If the file exist and the reanalisys is not forced, then current file doesn't need to be processed
    if(file.exists(.out_file) & !.force) next
    # Otherwhise we have to process it, so put it in the list of commands to be sent to the cluster
    .cmds <- bind_rows(.cmds,
                       data.frame(fsc1=.f_par$channels[1], fsc2=.f_par$channels[2], 
                                  ssc1=.f_par$channels[3], ssc2=.f_par$channels[4],
                                  gfp=.f_par$channels[5], 
                                  threshold=.f_par$scattering_threshold, scattering_frac_cells=.f_par$scattering_frac_cells,
                                  out=.out_file, file=.f,
                                  stringsAsFactors = FALSE))
      
      # c(.cmds, sprintf("Rscript %s --file %s --fsc1 %s --fsc2 %s --ssc1 %s --ssc2 %s --gfp %s --threshold %f --scattering_frac_cells %f --out %s",
      #                         .preproc_script, .f,
      #                         .f_par$channels[1], .f_par$channels[2], .f_par$channels[3], .f_par$channels[4], .f_par$channels[5],
      #                         .f_par$scattering_threshold, .f_par$scattering_frac_cells, .out_file))
    # If the directory for the outputfile doesn't exist, we also have to create it
    dir.create(dirname(.out_file), showWarnings = FALSE, recursive = TRUE)
  }
  .cmds <- .cmds %>% as_tibble()
  
  # If some files need to be filtered, then write the list of files to be processed on a cmd file and send the instructions to the clustes
  # (I think it is quick to store the cmd in the var cmds and then write them all at once to the command file)
  if(nrow(.cmds)){
    if(.preproc_func=='None') .preproc_func <- preproc
    .cmds %>% ungroup() %>% mutate(tmp = seq(n())) %>% split(.$tmp) %>% future_map_dfr(preproc, .progress = TRUE)
  }
  else{
    cat("All the data have already been filtered\n")
  }
}
