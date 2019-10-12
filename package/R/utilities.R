#'Check if the channels are correctly selected
#'
#'The function \code{check_channels} checks that the channel numbers or names
#'specified in \code{f_par} matches the desired channels in the FCM files.
#'
#'In order to call \code{preproc_facs_plates} the user must supply a structure
#'\code{f_par} which contains the names or the numbers of the channels that must
#'be analyzed. This function is useful to check whether the specified channels
#'match the ones in the FCM files.
#'
#'@export
check.channels <- function(.dirs, .f_par){
  for (.dir in .dirs) {
    .fs <- flowCore::read.flowSet(path=.dir, full.names=TRUE, pattern=.f_par$file_pattern, phenoData=list(Filename="$FIL", Time="$BTIM"))
    .fs <- .fs[[1]]@exprs[, f_par$channels]
    cat("\n", .dir, "Columns:", colnames(.fs), sep=" - ")
  }
}


#' Mapping of mean and standard deviation between log and linear space.
#'
#' This function maps the mean and standard deviation of a lognormal distribution between log and linear space.
#'
#' @param mean Mean of the distribution either in linear or log space
#' @param std Standard deviation of the distribution either in linear or log space
#' @param to.linear Should the mean and standard deviation provided be transformed from linear to log space or viceversa?
#'
#' @return A list containing the transformed mean and standard deviation.
#' @export
log_stats <- function(mean, std, to.linear=TRUE){
  if(to.linear){
    mu <- exp(mean+std**2/2)
    sigma <- sqrt((exp(std**2)-1)*mu**2)
  }
  else{
    sigma <- sqrt( log((std/mean)**2+1) )
    mu <- log(mean)-sigma**2/2
  }
  return(c(mu, sigma))
}


#' Creation of the bash file containing the instruction for the array job to
#' send to che cluster.
#'
#' This function creates the bash file \emph{job.sh} containing the instructions
#' to be sent to the cluster in order to filter the data based on their
#' scattering profile. The bash file contains the number of datasets to be
#' processed, that is why we cannot have a universal job.sh to be called for any
#' analysis, but we have to create it every time. The log of the jobs are stored
#' in the folder 'qlog'.
#'
#' @param .num_jobs Number of jobs to be sent to the cluster. Every job is a fcs
#'   file to e processed.
#' @param .cmd_file Name of the file containing the commands to filter the data.
#'   It must correspond to the number of lines in \emph{.num_jobs}.
#' @return The function writes the file \emph{job.sh} and it doesn't return any
#'   value..
create_job_script <- function(.num_jobs, .cmd_file){
  #Create directory where to store jobs results
  suppressWarnings(dir.create(file.path(getwd(), 'qlog')))

  #Create the text of the job
  .text <- sprintf(
"#!/bin/bash
#SBATCH --job-name=FCM
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --mem=1G
#SBATCH --output=%s
#SBATCH --error=%s

## here is the important flag for array jobs
#SBATCH --array=1-%i

# load your required modules
#######################
module purge
module load R
SEEDFILE=%s
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

", file.path(getwd(), 'qlog', 'output_%A_%a.out'), file.path(getwd(), 'qlog', 'error_%A_%a.out'), .num_jobs, .cmd_file)
  write.table(.text, quote=FALSE, row.names = FALSE, col.names = FALSE, file="./job.sh", append = FALSE)
}


##Clean any loaded dynamic library
.onUnload <- function (libpath) {
  library.dynam.unload("EFlow", libpath)
}


#' Add info to the preprocessed statistics.
#'
#' Takes information sotred in a tibble \emph{.info} and adds them to the tibble
#' of preprocessed statistics. The incorporation is done by matching the
#' folders.
#'
#' @param .pls a tibble containing the preprocessed statistics.
#' @param .info a tibble containing the column \emph{dir} and the information to
#'   be added. The matching between the information in \emph{.info} and the
#'   datasets in \emph{.pls} is done according to the directory names.
#' @return A tibble containing the original statistics and the added
#'   information.
#' @export
propagate_index_info <- function(.pls, .info) {
  left_join(.pls %>% mutate(dir = dirname(path)), .info, by='dir') %>% select(-dir)
}




#' Merging different measures supposing that they come form one single gaussian plus some variations 
#'
#' This function merges the values of different replicates under the following generative model:
#' It exists a true mean M, nonetheless different replicates have some variation tau around this mean, so that the measured
#' value on the dataset d is Md=M+tau. This value is measured experimentally with some error bar D, so the final measured value is
#' M.meas = M + tau + D.
#' This function estimates M with the error sigma of the estimate, and the width tau of the variations among datasets.
#'
#' @param measures The experimental values
#' @param errors The error bars on each measures
#'
#' @return A vector containing the estimated mean mu, the standard error of the mean estimate sigma and the width of the variations tau
#' @export
get_autofluo_stats <- function(measures, errors){
  tau0 <- as.numeric(diff(quantile(measures, c(.25, .75))))
  opt <- maxLik(tau_LL, start = c(tau=tau0), measures=measures, errors=errors, method = 'BFGS')
  
  if(!opt$code %in% c(0, 1, 2, 8)) stop("The likelihood optimization has not converged")
  
  tau <- abs(opt$estimate['tau']) #The likelihood is invariant for tau -> -tau
        
  alpha = 1. / (errors^2 + tau^2)
  m0 = sum(alpha)
  m1 = sum(alpha * measures)
  mu = m1 / m0
  new_sig = sqrt(1 / m0)
  return (c(tau, mu=mu, sigma=new_sig, return_code = opt$code))
}

#This is the likelihood that the previous function maximises
tau_LL <- function(tau, measures, errors){
  alpha = 1. / (errors^2 + tau^2)
  mSqrt = sum(sqrt(alpha))
  m0 = sum(alpha)
  m1 = sum(alpha * measures)
  m2 = sum(alpha * measures^2)
  mu = m1 / m0

  return(-(0.5 * (m0 * mu * mu - 2 * mu * m1 + m2) - 0.5*sum(log(alpha))))
}
