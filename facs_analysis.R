
mean_log10_to_ln <- function(.m)
  return(.m / log10(exp(1)))
var_log10_to_ln <- function(.v)
  return(.v / log10(exp(1))^2)
mean_log10_to_n_ln <- function(.m, .gfp_per_FCMunit=2.884)
  # can be written .m / log10(e) + log(.gfp_per_FCMunit)
  # or log( .gfp_per_FCMunit * exp(.m / log10(e)))
  return(.m / log10(exp(1)) + log(.gfp_per_FCMunit) )

noise_vs_n <- function(.sigma2, .beta, .n_bg) {
  # minimal noise relationship as defined in Wolf, et al. 2014 (eqn 66, http://dx.doi.org/10.1101/007237).
  function(.n) ifelse(.n > 2*.n_bg, .sigma2 * (1 - .n_bg/.n)^2 + .beta/.n * (1 - .n_bg/.n), NA)
}

noise_vs_n_ln <- function(.sigma2, .beta, .n_bg) {
  # minimal noise relationship derived from Wolf, et al. 2014 (eqn 66, http://dx.doi.org/10.1101/007237).
  function(.n_ln) (noise_vs_n(.sigma2, .beta, .n_bg))(exp(.n_ln))
}

d_sum_norm <- function(.x, .means, .sds) {
  # d_sum_norm returns the density of the sum of normal RV with parameters (.means, .sds)
  if (length(.means) != length(.sds)) 
    stop('.means and .sds must have the same length.')
  lapply(seq_along(.means), function(.i) dnorm(.x, .means[.i], .sds[.i])) %>%
    do.call(rbind, .) %>%
    apply(., 2, mean)
}


min_noise_manual_fit <- function(.mean_ln, .var_ln, .m_bg=NA, .bsize=NA, .extr_noise=NA,
                                 .gfp_per_FCMunit=2.884) {
  # min_noise_manual_fit helps the user choosing parameter values to fit manually
  # the minimal noise relationship as defined in Wolf, et al. 2014 (eqn 66, http://dx.doi.org/10.1101/007237).
  # mn_fit <- min_noise_manual_fit(mean_log10_to_ln(pr_stats$gfp_mean_noise_rm),
  #                                var_log10_to_ln(pr_stats$gfp_var_noise_rm))
  stopifnot(length(.mean_ln) == length(.var_ln))
  cat("\nHint: make sure that your data are transformed with ln and not log10...\n")
  
  if (is.na(.m_bg)) {
    cat("\nPlease enter your estimation of background fluo mean (in linear scale; e.g. from empty plasmid measurements):")
    .m_bg <- scan(n=1, what=numeric())
  }
  .n_bg <- .gfp_per_FCMunit * .m_bg
  .n_ln <- log( .gfp_per_FCMunit * exp(.mean_ln)) # - .n_bg 
  plp <- qplot(.n_ln, .var_ln, size=I(.5)) + expand_limits(y=0)
  print(plp)  
  
  # if burst size is given, adjust extrinsic noise
  if (is.na(.extr_noise) && !is.na(.bsize)) {
    value_set <- FALSE
    while (!value_set) {
      cat("\nPlease enter your estimation of extrinsic noise (i.e. var at high mean levels):")
      cat("\nHint: you can enter several values to see them plotted all at once...")
      .extr_noise <- scan(what=numeric())
      
      .xs <- seq(6.5, 13, length.out=500)
      myfns <- data.frame(extr_noise=.extr_noise) %>%
        group_by(extr_noise) %>%
        do((function(.df) data.frame(n_ln=.xs, var_ln=(noise_vs_n_ln(.df$extr_noise, .bsize, .n_bg))(.xs)))(.))
      
      pl <- plp +
        geom_line(aes(x=n_ln, y=var_ln, col=factor(extr_noise)), data=myfns, size=1) +
        labs(col='extrinsic\nnoise')
      print(pl)  
      
      cat("\nIs this plot correct (y or n)?\n")
      try(y <- scan(n=1, what=character()))
      if (length(y) == 1) {
        if (y == "y") {
          if (length(.extr_noise) > 1) {
            cat("\nPlease choose a single value of extrinsic noise before validating.")
          } else { 
            value_set <- TRUE
          }
        }
      }
    }
  }
  
  # if no params are given, adjust extrinsic noise first, then burst size
  if (is.na(.extr_noise) && is.na(.bsize)) {
    value_set <- FALSE
    while (!value_set) {
      cat("\nPlease enter your estimation of extrinsic noise (i.e. var at high mean levels):")
      .extr_noise <- scan(n=1, what=numeric())
      
      pl <- plp +
        geom_hline(yintercept=.extr_noise, col='red')
      print(pl)
      
      cat("\nIs this plot correct (y or n)?\n")
      try(y <- scan(n=1, what=character()))
      if (length(y) == 1) {
        if (y == "y") value_set <- TRUE
      }
    }
  }
  
  if (is.na(.bsize)) {
    value_set <- FALSE
    while (!value_set) {
      cat("\nPlease enter your estimation of transcription burst size:")
      cat("\n(yes, you must guess! 400 might be a good start :)")
      cat("\nHint: you can enter several values to see them plotted all at once...")
      .bsize <- scan(what=numeric())
      
      #       min_noise_b <- function(.b) noise_vs_n_ln(.extr_noise, .b, .n_bg)
      #       pl <- qplot(.n_ln, .var_ln) + expand_limits(y=0)
      #       for (.b in .bsize) 
      #         pl <- pl + stat_function(aes(col=factor(.b)), fun=min_noise_b(.b))
      #       print(pl + labs(col="burst size"))
      .xs <- seq(6.5, 13, length.out=500)
      myfns <- data.frame(bsize=.bsize) %>%
        group_by(bsize) %>%
        do((function(.df) data.frame(n_ln=.xs, var_ln=(noise_vs_n_ln(.extr_noise, .df$bsize, .n_bg))(.xs)))(.))
      
      pl <- plp +
        geom_line(aes(x=n_ln, y=var_ln, col=factor(bsize)), data=myfns, size=1) +
        labs(col='burst size')
      print(pl)
      
      cat("\nIs this plot correct (y or n)?\n")
      try(y <- scan(n=1, what=character()))
      if (length(y) == 1) {
        if (y == "y") {
          if (length(.bsize) > 1) {
            cat("\nPlease choose a single value of transcription burst size before validating.")
          } else { 
            value_set <- TRUE
          }
        }
      }
    }
  }
  
  .excess_noise <- .var_ln - (noise_vs_n_ln(.extr_noise, .bsize, .n_bg))(.n_ln)
  return(list(n_bg=.n_bg, bsize=.bsize, extr_noise=.extr_noise, n_ln=.n_ln, excess_noise=.excess_noise))
}

### PARAMETERS INFERENCE FOR ERIK'S MODEL ####

epsilon_opt <- function(.vp1, .vp2) {
  # eq. 38
  stopifnot(length(.vp1)==length(.vp2)) # "different number of measurements in the two replicates"
  .n <- length(.vp1)
  .epsilon <- sum((.vp1-.vp2)^2) / (2*.n)
  return(.epsilon)
}

sig2_coef <- function(.np, .nbg) {
  # eq. 4
  .out <- (1-(.nbg/.np))^2
  return(.out)
}
beta_coef <- function(.np, .nbg) {
  # eq. 5
  .out <- (1-(.nbg/.np)) / .np
  return(.out)
}

w_prom <- function(.ap,.bp,.epsilon,.gamma,.lambda){
# eq. 17
  .out <- 1 / (.epsilon^2 + (.bp*.gamma)^2 + (.ap*.lambda)^2)
  return(.out)
}

avg_prod <- function(.a, .b, .c){
# average product (eq. 21)
  stopifnot(length(.a)==length(.b) && length(.a)==length(.c))
  .out <- sum(.a * .b * .c) / length(.a)
  return(.out)
}

sig2_opt <- function(.w,.a,.b,.v){
# eq. 19
  .out <- ( avg_prod(.w,.b,.v)*avg_prod(.w,.a,.b) - avg_prod(.w,.b,.b)*avg_prod(.w,.a,.v) ) /
    ( avg_prod(.w,.a,.b)^2 - avg_prod(.w,.a,.a)*avg_prod(.w,.b,.b) )
  return(.out)
}
beta_opt <- function(.w,.a,.b,.v){
  # eq. 20
  .out <- ( avg_prod(.w,.a,.v)*avg_prod(.w,.a,.b) - avg_prod(.w,.a,.a)*avg_prod(.w,.b,.v) ) /
    ( avg_prod(.w,.a,.b)^2 - avg_prod(.w,.a,.a)*avg_prod(.w,.b,.b) )
  return(.out)
}

logL_factory <- function(.np,.vp,.nbg,.epsilon){
  .ap <- sig2_coef(.np, .nbg)
  .bp <- beta_coef(.np, .nbg)
  #   return(
  function(.lambda, .gamma, .print=FALSE, .plot=FALSE){
    stopifnot(length(.lambda)==1 && length(.gamma)==1)
    # L is defined only for positive lambda and gamma
    if (.lambda<0 || .gamma<0) {
      return(list(L=NA, sig2=NA, beta=NA))
    }
    # eq. 17
    .wp <- w_prom(.ap,.bp,.epsilon,.gamma,.lambda)
    .sig2_opt <- sig2_opt(.wp,.ap,.bp,.vp)
    .beta_opt <- beta_opt(.wp,.ap,.bp,.vp)
    if (.print) {
      print(sprintf('%f %f : %f %f', .lambda, .gamma, .sig2_opt, .beta_opt))
    }
    # plot
    if (.plot) {
      .pl <- ggplot(data=data.frame(np=.np, vp=.vp, ve=.ap*.sig2_opt+.bp*.beta_opt)) +
        geom_point(aes(x=np, y=vp, col="data")) +
        geom_point(aes(x=np, y=ve, col="estimation")) +
        scale_x_continuous(trans="log10", breaks=10^(1:7), labels=sprintf('log(%d)', 10^(1:7))) +
        expand_limits(y=0) +
        guides(col='none') + labs(title=sprintf('sig2=%.4f beta=%.1f', .sig2_opt, .beta_opt))
      print(.pl)
    }
    # eq. 16 (do not use eq. 18 to estimate lambda and gamma)
    .Ls <- -sum( (.vp-.ap*.sig2_opt-.bp*.beta_opt)^2 * .wp + log(1/.wp) )
    return(list(L=.Ls, sig2=.sig2_opt, beta=.beta_opt))
  }#)
}

logL_factory_optimWrapper <- function(.fun) {
  return( function(.pars) {
    .lambda <- .pars["lambda"]
    .gamma <- .pars["gamma"]
    .fun(.lambda, .gamma)$L
  } )
}

# lambda_gamma_optimum <- optim(logL_factory)
# test <- logL_factory(a,b,c,d)
# optim(c(lambda=0,gamma=0),test)

variance_per_optimum <- function(.wp,.ap,.bp){
  .denom <- length(.wp) * (avg_prod(.wp,.ap,.ap)*avg_prod(.wp,.bp,.bp) - avg_prod(.wp,.ap,.bp)^2)
  # eq. 26
  .var_sig2 <- avg_prod(.wp,.bp,.bp)/.denom
  # eq. 27
  .var_beta <- avg_prod(.wp,.ap,.ap)/.denom
  # eq. 28
  .covar <- avg_prod(.wp,.ap,.bp)/.denom
  return(list(sig2=.var_sig2, beta=.var_beta, covar=.covar))
}

variance_promoter_opt <- function(.vp,.ap,.bp,.sig2_opt,.beta_opt){
  # eq. 32
  .out <- .vp - .ap*.sig2_opt - .bp*.beta_opt
  return(.out)
}

sig2_promoter <- function(.epsilon,.ap,.bp,.wp){
  # eq. 33
  .out <- .epsilon^2 + 
           (.ap^2*avg_prod(.wp,.bp,.bp) + .bp^2*avg_prod(.wp,.ap,.ap) - 2*.ap*.bp*avg_prod(.wp,.ap,.bp)) /
              (length(.ap) * (avg_prod(.wp,.ap,.ap)*avg_prod(.wp,.bp,.bp) - avg_prod(.wp,.ap,.bp)^2))
  return(.out)
}

prom_contrib_opt <- function(.vp_opt,.ap,.bp,.sig2p,.lambda,.gamma){
  .factor <- .vp_opt / (.sig2p + .ap^2*.lambda^2 + .bp^2*.gamma^2)
  # eq. 30
  .eta_opt <- .factor * .ap*.lambda^2
  # eq. 31
  .delta_opt <- .factor * .bp*.gamma^2
  return(list(eta=.eta_opt, delta=.delta_opt))
}

prom_contrib_var <- function(.ap,.bp,.lambda,.gamma,.sig2p){
  .denom <- .sig2p + .ap^2*.lambda^2 + .bp^2*.gamma^2
  # eq. 34
  .etap_var= .lambda^2 * (.sig2p+.bp^2*.gamma^2) / .denom
  # eq. 35
  .deltap_var= .gamma^2 * (.sig2p+.ap^2*.lambda^2) / .denom
  # eq. 36
  .covar= - .ap*.bp*.lambda^2*.gamma^2 / .denom
  return(list(eta=.etap_var,delta=.deltap_var,covar=.covar))
}

beta_promoter <- function(.beta,.deltap){
# eq. 3
  return(.beta + .deltap)
}

sig2_promoter <- function(.var,.etap){
# eq. 2
  return(.var + .etap)
}

variance_promoter <- function(.np,.sig2p,.betap,.nbg){
  #eq 1
  .out <- .sig2p * (1-(.nbg/.np))^2 + .betap*(1-(.nbg/.np))/.np
  return(.out)
}

# data simulation for inference validation
simul_inference_data <- function(.n_p, # promoters' mean expression (in nb of GFP)
                                 .sig2=.02, .beta=500, .n_bg=600,
                                 .lambda=.01, .gamma=100,
                                 .epsilon=0.001 ) {
# simulate variance from means, according to Erik's model. Example:
# p_sim <- simul_inference_data(rlnorm(500, 7, 1.5), .lambda=.01)
# qplot(n, v_exp, data=p_sim, log="x") + expand_limits(y=0)
# llf <- logL_factory(p_sim$n, p_sim$v_exp, p_sim$n_bg, 0.001)
# lg_opt <- optim(c(lambda=.01,gamma=100), logL_factory_optimWrapper(llf),
#                 control=list(fnscale=-1)) # fnscale=-1 to find maximum
# 
# X <- 10^seq(-4, -1, length=110)
# Y <- 10^seq(1.8, 2.2, length=110)
# Z <- outer(X, Y, (function(.l, .g) llf(.l, .g)$L) %>% Vectorize())
# contour(log10(X), log10(Y), Z)
# points(log10(lg_opt$par["lambda"]),log10(lg_opt$par["gamma"]), col="red")
  
  .n <- length(.n_p)
  .a <- sig2_coef(.n_p, rep(.n_bg, .n))
  .b <- beta_coef(.n_p, rep(.n_bg, .n))
  .out <- data.frame(id=1:.n,
                     n=.n_p,
                     eta=rnorm(.n, 0, .lambda),
                     delta=rnorm(.n, 0, .gamma))
  
  .out <- mutate(.out, v = .a*.sig2 + .b*.beta + .a*eta + .b*delta,
                 v_exp = v + rnorm(.n, 0, .epsilon), n_bg=.n_bg) %>%
    filter(n > 2*.n_bg)
}
