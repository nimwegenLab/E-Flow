#'Filtering of the scattering profile
#'
#'\code{preproc} Is the functiont hat actually performs the fitting of the
#'Gaussian mixture to the scattering profile and extracts the mean and variance
#'estimates of the fluorescence of the population.
#'
#'@param .params Dataframe containing 
#'\enumerate{
#'\item file: the fsc file to analyze
#'\item fsc1, fsc2, ssc1, ssc2: the names of the scattering channels
#'\item gfp: the name of the channel with the GFP intensities
#'\item out: the name of the output file where to store the statistics
#'\item scattering_threshold, scattering_frac_cells: the threshold on the posterior or the fraction 
#'of cells to retain. 
#'}
preproc_func <- function(.params){
  .file <- .params$file
  .fsc1 <- .params$fsc1
  .fsc2 <- .params$fsc2
  .ssc1 <- .params$ssc1
  .ssc2 <- .params$ssc2
  .gfph <- .params$gfp
  .out_file <- .params$out
  .force <- .params$force
  .scattering_threshold <- .params$threshold
  .scattering_frac_cells <- .params$scattering_frac_cells
  
  #If the file already exists and I don't want to force the filtering, just load it and return it
  if(file.exists(.out_file) & !.force){
    print(paste("Loading from cache", .out_file))
    load(.out_file)
    return(.gfp_stats$stats)
  }
  
  .ff <- flowCore::read.FCS(.file)
  X = flowCore::exprs(.ff)
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
  EM.data.scat <- EM_mixture(
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
  EM.data.gfp <- EM_mixture(
    data,
    apply(data, 2, mean),
    c(0.9, 0.1),
    list(var(data)),
    0,
    FALSE)
  
  
  ## Get info about flow frame
  .path <- unlist(flowCore::keyword(.ff,"FILENAME"), use.names=FALSE)
  .well <- unlist(flowCore::keyword(.ff,"TUBE NAME"), use.names=FALSE)
  ## Get number of cells after the filtering of the scattering
  .n <- length(data.use)
  ## Get mean and var from fluo lognormal fit
  .fl_mean <- EM.data.gfp$mean[1,1] 
  .fl_var <- EM.data.gfp$covariances[[1]][1,1]
  ## Transform to real
  .tmp <- log_stats(.fl_mean, sqrt(.fl_var))
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
  .gfp_stats <- list(
    stats=data.frame(path=.path, well=.well, n=.n, w=EM.data.gfp$weight[1],
                     fl_mean=.fl_mean, fl_mean.err=.fl_mean.err, 
                     fl_var=.fl_var, fl_var.err=.fl_var.err, 
                     fl_mean_lin=.fl_mean_lin, fl_mean_lin.err=.fl_mean_lin.err, 
                     fl_var_lin=.fl_var_lin, fl_var_lin.err = .fl_var_lin.err, 
                     scattering_threshold = .scattering_threshold, scattering_frac_cells=.scattering_frac_cells,
                     stringsAsFactors = FALSE) %>% as_tibble(),
    preproc=data.frame(path=.path, well=.well, data.frame(X[data.use,]), post_fluo=EM.data.gfp$post[,1], stringsAsFactors = FALSE) %>% as_tibble()
  )
  
  save(.gfp_stats, file=.out_file) 
  
  return(.gfp_stats$stats)
}
