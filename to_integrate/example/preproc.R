##Libraries and utility functions
library(flowCore)
library(vngFCM)

log.stats <- function(mean, std, to.real=TRUE){
  if(to.real){
    mu <- exp(mean+std**2/2)
    sigma <- sqrt((exp(std**2)-1)*mu**2)
  }
  else{
    sigma <- sqrt( log((std/mean)**2+1) )
    mu <- log(mean)-sigma**2/2
  }
  return(c(mu, sigma))
}

##Parse command line arguments
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="fsc file to be filtered", metavar="character"),
  make_option(c("-f1", "--fsc1"), type="character", default=NULL,
              help="name of the channel containing the first forward scattering signal", metavar="character"),
  make_option(c("-f2", "--fsc2"), type="character", default=NULL,
              help="name of the channel containing the second forward scattering signal", metavar="character"),
  make_option(c("-s1", "--ssc1"), type="character", default=NULL,
              help="name of the channel containing the first side scattering signal", metavar="character"),
  make_option(c("-s2", "--ssc2"), type="character", default=NULL,
              help="name of the channel containing the second side scattering signal", metavar="character"),
  make_option(c("-g", "--gfp"), type="character", default=NULL,
              help="name of the channel containing the fluorescence signal", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=NULL,
              help="threshold for the posterior probability", metavar="character"),
  make_option(c("-n", "--scattering_num_cells"), type="numeric", default=NULL,
              help="number of cells to retain after the filtering. It is ignored if -t is specified", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file where to store the filtered result", metavar="character")
);

#since -t or -tp can be NA, the following commnd can give a warning, but I'm not concerned with that
opt_parser = suppressWarnings(OptionParser(option_list=option_list))
opt = suppressWarnings(parse_args(opt_parser))
.file <- opt$file
.fsc1 <- opt$fsc1
.fsc2 <- opt$fsc2
.ssc1 <- opt$ssc1
.ssc2 <- opt$ssc2
.gfph <- opt$gfp
.out_file <- opt$out
.scattering_threshold <- opt$threshold
.scattering_num_cells <- opt$scattering_num_cells


.ff <- read.FCS(.file)
X = exprs(.ff)
## Remove negative data
positive <- which(apply(X, 1, function(.x) !any(.x<=0)))
X <- X[positive, ]

## Transform to log
fsc.h <- log(X[, .fsc1])
fsc.w <- log(X[, .fsc2])
ssc.h <- log(X[, .ssc1])
ssc.w <- log(X[, .ssc2])
gfp.h <- log(X[, .gfph])

#Filter the scattering
data.scat <- cbind(fsc.h, fsc.w, ssc.h, ssc.w)
EM.data.scat <- EM_Luca(
  data.scat,
  apply(data.scat, 2, mean),
  c(0.9, 0.1),
  list(var(data.scat)),
  0, FALSE)

if(is.na(.scattering_threshold)){
  data.use <- order(EM.data.scat$post[,1], decreasing = TRUE)[1:.scattering_num_cells]
}else{
  data.use <- which(EM.data.scat$post[,1]>.scattering_threshold)
}

data <- cbind(gfp.h[data.use])
EM.data.gfp <- EM_Luca(
  data,
  apply(data, 2, mean),
  c(0.9, 0.1),
  list(var(data)),
  0,
  FALSE)

## For every cell, keep track if it is discarded or not
.survived <- rep(FALSE, nrow(data.scat))
.survived[data.use] <- TRUE

## Get info about flow frame
.path <- unlist(keyword(.ff,"FILENAME"), use.names=FALSE)
.well <- unlist(keyword(.ff,"TUBE NAME"), use.names=FALSE)
## Get number of cells after the filtering of the scattering
.n <- length(data.use)
## Get mean and var from fluo lognormal fit
.fl_mean <- EM.data.gfp$mean[1,1]
.fl_var <- EM.data.gfp$covariances[[1]][1,1]
## Transform to real
.tmp <- log.stats(.fl_mean, sqrt(.fl_var))
.fl_mean_lin <- .tmp[1]
.fl_var_lin <- .tmp[2]**2

## Collect everything in a dataframe
.gfp_stats <- list(stats=data.frame(path=.path, well=.well, n=.n, n_before_filter=nrow(data.scat),
                                    fl_mean=.fl_mean, fl_var=.fl_var,
                                    fl_mean_lin=.fl_mean_lin, fl_var_lin=.fl_var_lin),
                   preproc=data.frame(path=.path, well=.well, post_fluo=EM.data.scat$post[,1], survived=.survived,
                                      fsc.h, fsc.w, ssc.h, ssc.w, gfp.h)
)

save(.gfp_stats, file=.out_file)

