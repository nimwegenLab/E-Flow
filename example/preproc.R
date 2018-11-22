##Libraries and utility functions
#library(dplyr)
#library(tidyr)
library(flowCore)
library(vngFCM)
#dyn.load("./EM.so", sep='')

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
              help="fsc file to be filtered"),
  make_option(c("-f1", "--fsc1"), type="character", default=NULL, 
              help="name of the channel containing the first forward scattering signal"),
  make_option(c("-f2", "--fsc2"), type="character", default=NULL, 
              help="name of the channel containing the second forward scattering signal"),
  make_option(c("-s1", "--ssc1"), type="character", default=NULL, 
              help="name of the channel containing the first side scattering signal"),
  make_option(c("-s2", "--ssc2"), type="character", default=NULL, 
              help="name of the channel containing the second side scattering signal"),
  make_option(c("-g", "--gfp"), type="character", default=NULL, 
              help="name of the channel containing the fluorescence signal"),
  make_option(c("-t", "--threshold"), type="double", default=NULL, 
              help="threshold for the posterior probability", metavar="character"),
  make_option(c("-tn", "--scattering_frac_cells"), type="double", default=NULL, 
              help="fraction of cells to keep. It is ignored if -t is specified"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file where to store the filtered result")
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
.scattering_frac_cells <- opt$scattering_frac_cells

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

if(is.na(.scattering_threshold) | .scattering_threshold=='NA'){
  .em.thres <-  quantile(EM.data.scat$post[,1], 1-.scattering_frac_cells)
}else{
  .em.thres <- .scattering_threshold
}

data.use <- which(EM.data.scat$post[,1]>.em.thres)

data <- cbind(gfp.h[data.use])
EM.data.gfp <- EM_Luca(
  data,
  apply(data, 2, mean),
  c(0.9, 0.1),
  list(var(data)),
  0,
  FALSE)


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
.fl_var_lin <- .tmp[2]**2

## Get the error on the estimation of the mean using the hessian matrix of the likelihood
d2ll <- function(x_i, mu, s, rho){
delta <- diff(range(x_i))
rho_U <- 1-rho
g_i <- dnorm(x_i, mean = mu, sd = sqrt(s))
P_i <-  rho*g_i+rho_U/delta

dg_mu <- (x_i-mu)/s * g_i
d2g_mu <- 1/s * ( (x_i-mu)*dg_mu-g_i )

dg_s <- g_i/(2*s) * ( (x_i-mu)^2/s -1 )
d2g_s <- dg_s/(2*s) * ( (x_i-mu)^2/s -1 ) - g_i/(2*s^2) * ( (x_i-mu)^2/s -1 ) - g_i/(2*s) * ( (x_i-mu)^2/s^2 )

hess_mu <- rho/P_i*(d2g_mu-rho/P_i*dg_mu^2)
hess_s <- rho/P_i*(d2g_s-rho/P_i*dg_s^2)

return(data.frame(hess_mu=sum(hess_mu), hess_s=sum(hess_s)))
}


.hessian <- d2ll(data[,1], .fl_mean, .fl_var, EM.data.gfp$weight[1])
.fl_mean.err <- sqrt(-1/.hessian$hess_mu)
.fl_var.err <- sqrt(-1/.hessian$hess_s)

## Get the estimates in real space
.fl_mean_lin <- exp(.fl_mean + .fl_var/2. + .fl_mean.err^2/2. + .fl_var.err^2/8)
.fl_mean_lin.err <- sqrt( (exp(.fl_mean.err^2+.fl_var.err^2/4.)-1)*.fl_mean_lin^2 )

.fl_var_lin.err <- sqrt( (2*.fl_var_lin*.fl_mean.err)^2 + .fl_var.err^2*(.fl_var_lin + exp(2*(.fl_mean+.fl_var))) )



## Collect everything in a dataframe
.gfp_stats <- list(stats=data.frame(path=.path, well=.well, n=.n, w=EM.data.gfp$weight[1],
                                    fl_mean=.fl_mean, fl_mean.err=.fl_mean.err, 
                                    fl_var=.fl_var, fl_var.err=.fl_var.err, 
                                    fl_mean_lin=.fl_mean_lin, fl_mean_lin.err=.fl_mean_lin.err, 
                                    fl_var_lin=.fl_var_lin, fl_var_lin.err = .fl_var_lin.err),
                   preproc=data.frame(path=.path, well=.well, data.frame(X[data.use,]), post_fluo=EM.data.gfp$post[,1])
)

save(.gfp_stats, file=.out_file) 

