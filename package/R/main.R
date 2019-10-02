#'Preprocessing of FCM files.
#'
#'\code{preproc_facs_plates} is the main function to preprocess a set of FCM
#'files contained in several folders.
#'
#'This functions is the \emph{main function} that should be called to process
#'FCM data. It takes all the FCM files in a set of folders and applies the
#'filtering in the scattering space to remove non viable cells; it fits a
#'mixture of a lognormal and a uniform distribution to extract the mean and
#'variance from the fluorescence of the remaining population; it removes the
#'autofluorescence and the shot noise. The scattering filtering and the
#'fluorescence fit are performed parallel as an array job. Moreover if the
#'package \emph{doParallel} is present, the code is further parallelized
#'exploiting half of the cores available on the local machine, as further
#'explained in the documentation of \code{preproc_facs_plate}.
#'
#'@param .dirs List of directories containing the files to be preprocessed.
#'@param .data2preproc Function that transforms the directory name of the FCM
#'  files to the directory name where to store preprocessed data.
#'@param .cache_namer Function that takes the directory name where to store the
#'  preprocessed files and returns the name and path of a Rdata file where to
#'  store the preprocessed data of all the files in the directory. By default,
#'  if \emph{path} is the directory of the preprocessed data, the file is stored
#'  in \emph{path/path_preproc_luca.Rdata}.
#'@param .filter_preproc_namer Function that takes the directory name where to
#'  store the preprocessed files and returns the name and path of a Rdata file
#'  where to store the preprocessed data of the single files. By default, if
#'  \emph{path} is the directory of the preprocessed data and the file is called
#'  abc.fcs, then the preprocessed data for that file is stored in
#'  \emph{path/abc.Rdata}.
#'@param .f_par List containing the parameters used to preprocess the data. It
#'  must contain the elements \enumerate{ \item \emph{channels}: the names of
#'  the channels containing the signals to be analysed, in the order: forward
#'  scattering first signal, forward scattering second signal, side scattering
#'  first signal, side scattering second signal, fluorescence. \item
#'  \emph{scattering_threshold}: a cell is considered viable if its probability
#'  to come from the main scattering distribution is bigger or equal to this
#'  threshold. \item \emph{scattering_frac_cells}: retain only the fraction
#'  \emph{scattering_frac_cells} of cells with the highest poterior distribution; this is
#'  a way to specify a specific number of cells to retain from the filtering. If
#'  \emph{scattering_threshold} is specified, \emph{scattering_frac_cells} will
#'  be ignored. \item \emph{file_pattern}: the pattern fo the FCM files to be
#'  analysed. Useful if the directories contain files that must not be analysed.
#'  \item \emph{delta_shot_noise}: the amplitude of the shot noise. 
#'  \item \emph{delta_shot_noise_err}: the error associated to the shot noise. It can be zero\item
#'  \emph{autofluo}: Tibble containing for every dir information about the autofluorescence to be subtracted. 
#'  These are, in order: mean, error of the mean, standard deviation, error of the standard deviation. 
#'  These can be computed with the function get_autofluo_stats()}
#'@param .plot Should some plots summarizing the filtering be produced and
#'  stored to a pdf? Default to TRUE
#'@param .force_noise_rm Forces the removal of the autofluorescence and of the shot noise
#'@param .pdf_dim Dimensions of the pdf file where to print the plots: c(width,
#'  height, rows, columns)
#'@param .verbose Should some information about the several stages the function
#'  goes through be displayed? Default to FALSE
#'@param .drop_raw Should the raw scattering and fluorescence values of the
#'  cells which survived the scattering filter (i.e. viable cells) be dropped
#'  and only the summary statistics be returned? Dropping them can save space
#'  and make the algorithm faster. Default to TRUE.
#'@param .preproc_script Name of the script that perform the filtering of the
#'  scattering signals and the lognormal fit of the fluorescence. Deafult to
#'  \emph{preproc.R}.
#'@param .jobs_cmd_name Name of the file where to store the list of commands to
#'  send to the cluster. Default to \emph{filtering_cmd}
#'@param .use_doParallel If TRUE take advantage of the package doParallel to
#'  spped up the code by exploiting half of the cores of the local machine. The
#'  filtering of the scattering profile is in any case sent to the cluster.
#'@return A list containing \enumerate{ \item \emph{stats} the summary
#'  statistics of all the FCM analized. \item \emph{raw} the raw
#'  data, i.e., raw scattering and fluorescence values of the cells which
#'  survived the scattering filter. If \emph{.drop_raw} is FALSE, then this
#'  item is NULL.}
#'@export
analyse_raw <- function(.dirs, .data2preproc,
                        .filter_preproc_namer=(function(.d) sub("fcs$", 'RData', .d)),
                        .f_par,
                        .force = FALSE,
                        .preproc_func=NULL) {
  
  ## The first thing to do is to filter all the datasets based on their scattering profile.
  ## This step requires some time, so I parallelize it through an array job and I store the results in Rdata files named acording to .filter_preproc_name
  ##in the preproc directory specified by .data2preproc
  cat("Filtering and log normal fitting the fcs files\n")
  stats <- scattering_filter(.dirs, .data2preproc, .f_par, .filter_preproc_namer, .preproc_func, .force)

  return(stats)
}


#'@export
make_distribution_plots <- function(.dirs, .data2preproc, .pdf_dim=c(2, 2.4, 6, 8)){
  .data2preproc(.dirs) %>% purrr::map(make_distribution_plots_dir, .pdf_dim=.pdf_dim)
}

#'@export
remove_noise <-  function(.pls, .f_par){
  #Add shot noise info
  .pls <- .pls %>% mutate(delta_shot_noise = .f_par$delta_shot_noise, delta_shot_noise.err = .f_par$delta_shot_noise.err)
  
  .pls <- .pls %>% 
    mutate(fl_mean_lin_autofluo_rm = fl_mean_lin-autofluo_mean, 
           fl_mean_lin_autofluo_rm.err = sqrt(fl_mean_lin.err^2 + autofluo_mean.err^2),
           fl_var_lin_autofluo_rm = fl_var_lin-autofluo_var,
           fl_var_lin_autofluo_rm.err = sqrt(fl_var_lin.err^2+autofluo_var.err^2),
           fl_var_lin_noise_rm = fl_var_lin_autofluo_rm-delta_shot_noise^2*(fl_mean_lin_autofluo_rm),
           fl_var_lin_noise_rm.err = sqrt(fl_var_lin_autofluo_rm.err^2 +
                                            (2*delta_shot_noise*(fl_mean_lin_autofluo_rm)*delta_shot_noise.err)^2 + 
                                            (delta_shot_noise^2*fl_mean_lin_autofluo_rm.err)^2))
  return(.pls)
}
