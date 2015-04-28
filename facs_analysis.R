
mean_log10_to_ln <- function(.m)
  return(.m / log10(exp(1)))
var_log10_to_ln <- function(.v)
  return(.v / log10(exp(1))^2)
mean_log10_to_n_ln <- function(.m, .gfp_per_FCMunit=2.884)
  # can be written .m / log10(e) + log(.gfp_per_FCMunit)
  # or log( .gfp_per_FCMunit * exp(.m / log10(e)))
  return(.m / log10(exp(1)) + log(.gfp_per_FCMunit) )

min_noise_vs_n <- function(.sigma2_ab, .beta, .n_bg) {
  # minimal noise relationship as defined in Wolf, et al. 2014 (eqn 66, http://dx.doi.org/10.1101/007237).
  function(.n) ifelse(.n > 2*.n_bg, .sigma2_ab * (1 - .n_bg/.n)^2 + .beta/.n * (1 - .n_bg/.n), NA)
}

min_noise_vs_n_ln <- function(.sigma2_ab, .beta, .n_bg) {
  # minimal noise relationship derived from Wolf, et al. 2014 (eqn 66, http://dx.doi.org/10.1101/007237).
  function(.n_ln) (min_noise_vs_n(.sigma2_ab, .beta, .n_bg))(exp(.n_ln))
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
        do((function(.df) data.frame(n_ln=.xs, var_ln=(min_noise_vs_n_ln(.df$extr_noise, .bsize, .n_bg))(.xs)))(.))
      
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
      
      #       min_noise_b <- function(.b) min_noise_vs_n_ln(.extr_noise, .b, .n_bg)
      #       pl <- qplot(.n_ln, .var_ln) + expand_limits(y=0)
      #       for (.b in .bsize) 
      #         pl <- pl + stat_function(aes(col=factor(.b)), fun=min_noise_b(.b))
      #       print(pl + labs(col="burst size"))
      .xs <- seq(6.5, 13, length.out=500)
      myfns <- data.frame(bsize=.bsize) %>%
        group_by(bsize) %>%
        do((function(.df) data.frame(n_ln=.xs, var_ln=(min_noise_vs_n_ln(.extr_noise, .df$bsize, .n_bg))(.xs)))(.))
      
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
  
  .excess_noise <- .var_ln - (min_noise_vs_n_ln(.extr_noise, .bsize, .n_bg))(.n_ln)
  return(list(n_bg=.n_bg, bsize=.bsize, extr_noise=.extr_noise, n_ln=.n_ln, excess_noise=.excess_noise))
}
