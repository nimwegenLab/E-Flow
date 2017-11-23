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
#'@inheritParams preproc_facs_plates
#'@export
check.channels <- function(.dirs, .f_par){
  for (.dir in .dirs) {
    .fs <- flowCore::read.flowSet(path=.dir, full.names=TRUE, pattern=.f_par$file_pattern, phenoData=list(Filename="$FIL", Time="$BTIM"))
    .fs <- .fs[, .f_par$channels]
    cat("\nColumns selected to analyze\n[order should be: fsc, ssc, fl1]:", colnames(.fs), sep="\n")
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
  .text <- sprintf("#!/bin/bash
#$ -N facs_preproc
##$ -pe smp 1
##$ -l membycore=2G
##$ -l runtime=04:00:00
#$ -o %s
#$ -e %s
#$ -cwd

## here is the important flag for array jobs
## 50 is the number of lines in your cmds file. To get this number just 'wc -l commands.cmd'
#$ -t 1-%i
## you can limit the number of simultaneous running jobs
##$ -tc 15

# load your required modules
#######################
module load R
SEEDFILE=%s
SEED=$(sed -n ${SGE_TASK_ID}p $SEEDFILE)
eval $SEED

", file.path(getwd(), 'qlog', '$JOB_NAME.o$JOB_ID'), file.path(getwd(), 'qlog', '$JOB_NAME.e$JOB_ID'), .num_jobs, .cmd_file)
  write.table(.text, quote=FALSE, row.names = FALSE, col.names = FALSE, file="./job.sh", append = FALSE)
}


##Clean any loaded dynamic library (e.g. EM_Luca)
.onUnload <- function (libpath) {
  library.dynam.unload("vngFCM", libpath)
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
  # uniformize directory names
  .pls <- mutate(.pls, tmp=basename(dirname(as.character(path))))
  .dir_lv <- .pls$tmp %>% c(basename(.info$dir)) %>% unique
  .info <- mutate(.info, dir=basename(as.character(dir)))
  .pls <- left_join(.info, .pls, by=c("dir"="tmp"))

  return(.pls)
}
