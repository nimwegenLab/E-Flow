#'Preprocessing of FCM files.
#'
#'\code{analyse_raw} reads all FCM files in a group of directories, filters them based on the scattering signals, and 
#'returns the mean and variance of the populations of cells.
#'
#'This functions is the \emph{main function} that should be called to process
#'FCM data. It takes all the FCM files in a set of folders and applies the
#'filtering in the scattering space to remove non viable cells. Then it fits a
#'mixture of a lognormal and a uniform distribution to extract the mean and
#'variance from the fluorescence of the remaining population. Notice that to speed up the 
#'processing, the function \emph{plan} from the package \emph{future} can be called
#'just before this function.
#'
#'@param .dirs Character vector of directories containing the files to be processed.
#'@param .data2preproc Function that transforms the directory names of the FCM
#'  files to the directory names where to store the processed data. For example, if the files are in a directory
#'  called \emph{data/condition0/}, you want to store the processed data in \emph{data_processed/condition0/}. In 
#'  order to do that, you can set \code{data2preproc <- function(.d) sub('data', 'data_proc', .d)}
#'@param .filter_preproc_namer Function that takes the directory name where to
#'  store the processed files and returns the name and path of a Rdata file
#'  where to store the processed data of the single FCM files. By default, if
#'  \emph{path} is the directory of the processed data and the file is called
#'  abc.fcs, then the processed data for that file is stored in
#'  \emph{path/abc.Rdata}.
#'@param .f_par List containing the parameters used to process the FCM files. It
#'  must contain the elements \enumerate{ \item \emph{channels}: the names of
#'  the channels containing the signals to be analysed, in the order: forward
#'  scattering first signal, forward scattering second signal, side scattering
#'  first signal, side scattering second signal, fluorescence. Usually the first 
#'  signal is the height and the second one is the width. \item
#'  \emph{scattering_threshold}: a cell is considered viable if its probability
#'  to come from the main scattering distribution is bigger or equal to this
#'  threshold. \item \emph{scattering_frac_cells}: retain only the fraction
#'  \emph{scattering_frac_cells} of cells with the highest poterior distribution; this is
#'  a way to specify a specific number of cells to retain from the filtering. If
#'  \emph{scattering_threshold} is specified, \emph{scattering_frac_cells} will
#'  be ignored. \item \emph{file_pattern}: the pattern of the FCM files to be
#'  analysed. Useful if the directories contain files that must not be analysed.
#'  \item \emph{delta_shot_noise}: the amplitude of the shot noise. 
#'  \item \emph{delta_shot_noise_err}: the error associated to the shot noise. It can be zero.}
#'@param .force Forces the removal of the autofluorescence and of the shot noise, even if the FCM have already been 
#'processed.
#'@param .preproc_func Name of the function that performs the filtering of the
#'  scattering signals and the lognormal fit of the fluorescence. If NULL, the 
#'  default function is used.
#'  
#'@return A dataframe containing the summary
#'  statistics of all the FCM analized \enumerate{\item \emph{n}: the number of cells after the scattering filtering.
#'  \item \emph{fl_mean, fl_var, fl_mean.err, fl_var.err, fl_mean_lin, fl_mean_lin.err, fl_var_lin, fl_var_lin.err}: 
#'  estimated means and variances of the populations of cells in log and real space, with their error.
#'  \item \emph{scattering_threshold, scattering_frac_cells}: parameters used to perform the scattering filtering. They are
#'  the ones passed to the function through the \emph{.f_par} structure}.
#'@export
analyse_raw <- function(.dirs, .data2preproc,
                        .filter_preproc_namer=(function(.d) sub("fcs$", 'RData', .d)),
                        .f_par,
                        .force = FALSE,
                        .preproc_func=NULL) {
  
  ## The first thing to do is to filter all the datasets based on their scattering profile.
  ## This step requires some time, so I parallelize it through an array job and I store the results in Rdata files named according to .filter_preproc_name
  ##in the preproc directory specified by .data2preproc
  cat("Filtering and log normal fitting the fcs files\n")
  stats <- scattering_filter(.dirs, .data2preproc, .f_par, .filter_preproc_namer, .preproc_func, .force)

  return(stats)
}

#'Plot the fitted fluorescence distributions
#'
#'\code{make_distribution_plots} goes through all the directories and creates a pdf file containing the fluorescence 
#'distribution of the cells survived to the filtering based on the scattering profile. With each distribution it is also shown 
#'the fitted mixture of the gaussian and the uniform, from which the mean and variance are extracted. This can be useful to visually
#'inspect the goodness of the datasets.
#'
#'@param .dirs Character vector of directories where to look for the processed files created by the function 
#'\code{analyse_raw}
#'@param .data2preproc Function that transforms the directory name of the FCM
#'  files to the directory name where to store preprocessed data.
#'@param .pdf_dim Width and height of the pdf to be generated. In inches.
#'@export
make_distribution_plots <- function(.dirs, .data2preproc, .pdf_dim=c(16, 14)){
  .data2preproc(.dirs) %>% purrr::map(make_distribution_plots_dir, .pdf_dim=.pdf_dim)
}


#'Remove the autofluorescence and the shot noise
#'
#'\code{remove_noise} is responsible to remove the autofluorescence and the shot noise from the processed data returned by 
#'\code{analyse_raw}. 
#'
#'@param .pls The results from \code{analyse_raw} with added information about the autofluorescence. The autofluorescence
#'information must be stored in the columns: 
#'\enumerate{\item \emph{autofluo_mean, autofluo_mean.err, autofluo_var, autofluo_var.err}: the mean and variance of the 
#'autofluorescence and their errors}
#'See the example file to see how to compute this information
#'
#' @return The same dataframe as the input \emph{.pls}, but with the extra columns
#' \enumerate{\item \emph{fl_mean_lin_autofluo_rm, fl_var_lin_autofluo_rm}: The mean fluorescence of the population after 
#' autofluorescence removal. \item \emph{fl_var_lin_noise_rm}: The autofluorescence of the population after removal of 
#' autofluorescence and shot noise.}
#' Notice that the mean is not affected by the shot noise removal, therefore we only return the mean corrected by the 
#' autofluorescence.
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
