# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# biocLite(c("flowCore", "flowViz"))
# install.packages(c('Hmisc', 'tools', 'tidyr', 'dplyr', 'ggplot2'))

suppressPackageStartupMessages( {
  library(flowCore)
  library(flowViz)
  # library(geneplotter)
  # library(MASS)
  
  library(Hmisc)
  library(tools)
  library(tidyr)
  library(dplyr)
  
  library(ggplot2)
  library(RColorBrewer)
})

# ggplot2 config
theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")
scale_colour_periodic_brewer <- 
  function(...) scale_colour_manual(..., values = rep(c(brewer.pal(4, 'Set1'), 'gray42'), 100))


### FUNCTIONS

split_well_col <- function(.df, .col='well', .override=FALSE) {
# split_well_col splits the column named after .col to `row` and `col`
  if ('row' %in% names(.df) || 'col' %in% names(.df)) {
    if (.override) { warning("Column(s) named 'row' and/or 'col' have been overriden by split_well_col().")
    } else {stop("split_well_col(): column(s) named 'row' and/or 'col' already exist and cannot be overriden.")}
  }
  .df <- separate_(.df, .col, c('row', 'col'), sep=1, remove=FALSE)
  .df <- mutate_(.df, .dots=list(row=~factor(row), col=~factor(as.numeric(col))))
  return(.df)
}

read_od_file <- function(.filename, .format='wide') {
# read_od_file reads raw data files produced by xxx spectrophotometer.
# data is output to wide format by default, but is reshaped if .format == 'long' 
# do no use file header since the last column does not have a name. rename to 1:12 instead.
  .od <- read.table(.filename, header=FALSE, row.names=1, fill=TRUE)
  .od <- .od[3:10,1:12]
  names(.od) <- 1:12
  
  if (.format == 'long') {
    .od$row <- row.names(.od)
    .od <- gather_(.od, 'col', 'od', as.character(1:12))
    .od <- mutate_(.od, .dots=list(well=~paste(row, col, sep='')))
  }
  return(.od)
}

get_flowset_timestamp <- function(.fs) {
  .t_locale <- Sys.getlocale(category = "LC_TIME")
  .en_locale <- Sys.setlocale("LC_TIME", "en_US.UTF-8")
  .timestamp <- strptime(paste(keyword(.fs,"$DATE"), keyword(.fs,"$BTIM")),
                         '%d-%b-%Y %T')
  Sys.setlocale("LC_TIME", .t_locale)
  return(.timestamp)
}

facs_contour_with_polyg <- function(.ff, .chs, .lims, .polyg=NULL, .highlight_no_polyg=TRUE, .fill=topo.colors(4, alpha=0.1)) {
# draw a contour plot given 2 channels .ch of a flowframe, within .lims for both axis
# if polygon is a 2 column matrix a coordinates, the corrseponding gate is overlaid in red
  contour(.ff, y=.chs, xlim=.lims, ylim=.lims, nlevels=8, fill=.fill, col='gray70', lwd=0.05)
  polygon(.polyg, lwd=.2, border='red', col=rgb(1, 1, 1, .5))
  if (is.null(.polyg) && .highlight_no_polyg)  # draw a red background if no polygon provided
    rect(.lims[1]-10, .lims[1]-10, .lims[2]+10, .lims[2]+10, col=rgb(1, 0, 0, .1))
}

density_with_normal <- function(.xs, .normal, .name, .sub, .lab, .lims=c(0.0,6.0)) {
# plot the density of a vector .xs, and overlays normal densities as definied by .normal
# (a list of (mu, sigma, color) lists)
  .dens <- density(.xs)
  plot(.dens, ty="n", lwd=1, cex.main=.9, xlim=.lims, xlab=.lab, ylab="kernel density", main=.name, sub=.sub)
  polygon(.dens, col="gray60", border=NA)
  curve(dnorm(x, mean(.xs), sd(.xs)), col="gray30", lwd=.3, add=TRUE)
  for (.n in .normal) curve(dnorm(x, .n$mu, .n$sigma), col=.n$col, lwd=.3, add=TRUE)
}

set_fsc_ssc_gates <- function(.dir, .f_par, .pattern='A1', .interactive=FALSE) { 
# set_fsc_ssc_gates reads in a single flowFrame in order to create the gates and define log transform.
# gates can be controlled visually with .interactive=TRUE
  .fs <- read.flowSet(path=.dir, full.names=TRUE, pattern=.pattern, phenoData = list(Filename = "$FIL"))
  if (.interactive) cat("All colnames:", colnames(.fs), sep="\n")
  .fs <- .fs[, .f_par$channels]
  
  # print out the selected columns to be sure
  if (.interactive) {
    cat("\nColumns selected to analyze\n[order should be: fsc, ssc, fl1, fl2]:", colnames(.fs), sep="\n")
    cat("\nNum cells to gate:", .f_par$min.cells, "\nLimits of FSC-SSC density plots:", .f_par$lims, "\n")
    cat("File pattern: ", .f_par$file.pattern, "\n")
    cat("\nIs this all correct (y or n)?\n")
    y <- scan(n=1, what=character())
    if(y != "y") { stop("\nincorrect parameter values or \n you didn't type 'y'\n") }
  }
    
  # DEFINE log TRANSFORMATION
  logT <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
  ListlogT <- transformList(from=colnames(.fs), tfun=logT)
  
  # DEFINE debris FILTERS and GATES
  .n_channels <- length(colnames(.fs))
  .m <- matrix( rep(c(1.01, .f_par$facs.max), .n_channels), ncol=.n_channels)
  colnames(.m) <- c(colnames(.fs))
  not.debris.gate <- rectangleGate(filterId="NonDebris", .m)
  ### scale.factor is the factor of stdev; n is number of events used to calculate the bivariate normal
  n2f.gate <- norm2Filter(x=c(colnames(.fs)[1], colnames(.fs)[2]),
                           method="covMcd", scale.factor=1, n=20000, filterId="Norm2Filter")
  
  if (.interactive) {
    .fst <- ListlogT %on% .fs[[1]]
    facs_contour_with_polyg(.fst, 1:2, .f_par$lims, .highlight_no_polyg=FALSE)
    cat("\nIs this plot okay (y or n)?\n")
    y <- scan(n=1, what=character())
    if(y != "y") { stop("\nincorrect parameter values or \n you didn't type 'y'\n") }
    dev.off()
  }
  
  return(list(ListlogT=ListlogT, not.debris.gate=not.debris.gate, n2f.gate=n2f.gate))
}

preproc_facs_plate <- function(.dir, .out_dir, .f_par, .f_utils, .plot=TRUE, 
                               .verbose=0,              # console output (0: silent, 1: minimalist, 2: detailled)
                               .write_format=NULL,      # formats in which to write subseted data to .out_dir (vector of strings in 'fcs', 'tab' or both)
                               .min_cells=5000,         # fraction of data to remove beofre calculating trimmed stats
                               .pdf_dim=c(2, 2, 6, 8), 	# width, height, rows, columns
                               .cache_namer = (function(.d) file.path(.d, paste(basename(.d), '_preproc.Rdata', sep=''))),
                               .force=FALSE) {          # flag to force analysing raw fcs files even if a cache file exists
# preproc_facs_plate loads .fcs files in .dir as a flowSet, apply FSC/SSC gating to select .min_cells
# subsequently this subset is filtered on channel fl1 to remove background (using a normal + uniform mixture model)
# optionnaly OD files can be loaded in order to append od values to the statistics
  if (.verbose) cat("\nCurrent directory:", .dir, "\n")
  
  # if a cache file exists, load it and skip the analysis
  if (file.exists(.cache_namer(.out_dir)) && !.force) {
    load(.cache_namer(.out_dir)) # load objet .pl
    if (.verbose) cat("Cache loaded from ", basename(.cache_namer(.out_dir)), "\n", sep='')
    return(.pl)
  }
  
  # if a local FACS params file exists, load it
  .fpar_file <- list.files(.dir, pattern=paste(basename(.dir), '_fpar.[rR]', sep=''), full.names=TRUE)
  if (length(.fpar_file) > 0) {
    source(.fpar_file[1], local=TRUE)
    .f_par <- f_par
    .f_utils <- set_fsc_ssc_gates(.dir, .f_par)
    if (.verbose) cat("Local FACS parameters loaded...\n")
  }
  
  # Prepare output (directory, variables, plot)
  dir.create(.out_dir, recursive=TRUE, showWarnings=FALSE)  # here we will store the gated facs data and cache files
  .gates <- data.frame()
  .preproc <- data.frame()
  .stats <- data.frame()
  if(.plot) { # open a new pdf for plots (in raw data dir)
    pdf(file = file.path(.dir, paste(basename(.dir), 'pdf', sep='.')), height=.pdf_dim[2]*.pdf_dim[3], width=.pdf_dim[1]*.pdf_dim[4])
    par(mfrow=c(.pdf_dim[3],.pdf_dim[4]))
    par(mar=c(3.5, 2.5, 2, .5) + 0.1,
        mgp=c(1.5, .5, 0),
        tcl=-0.25)
  }
  
  # Load OD file if a pattern is provided
  .f_od = NA
  try({.od_file <- list.files(.dir, pattern=.f_par$od.pattern)
  .f_od <- read_od_file(file.path(.dir, .od_file), .format='long')
  if (.verbose>1) cat("OD file read\n\n")}, silent=TRUE)
  #   if (exists('od.pattern', where=.f_par)) {
  #     .od_file <- list.files(.dir, pattern=.f_par$od.pattern)
  #     .f_od <- read_od_file(file.path(.dir, .od_file), .format='long')
  #   } else {
  #     .f_od = NA
  #   }
  
  .fs <- read.flowSet(path=.dir, full.names=TRUE, pattern=.f_par$file.pattern, phenoData=list(Filename="$FIL", Time="$BTIM"))
  .fs <- .fs[, .f_par$channels]
  if (.verbose>1) cat("number of files read: ", length(.fs), "\n\n")  
  
  # Do all the primary gating on the entire flowSet
  fs.cells  <-  Subset(.fs, .f_utils$not.debris.gate)
  if (.verbose>1) cat("Debris gated\n\n")
  fs.cellsT <- .f_utils$ListlogT %on% fs.cells
  if (.verbose>1) cat("Data Log transformed\n\n")
  fs.normT <- Subset(fs.cellsT, .f_utils$n2f.gate)
  if (.verbose>1) cat("Bv Normal filter applied\n\n")
  
  # Do the secondary gating on each flowFrame (in chronological order)
  if (.verbose>1) cat("Subsetting and filtering files:\n")
  .ch <- c('fsc' = which(names(.f_par$channels) == 'fsc'),
           'ssc' = which(names(.f_par$channels) == 'ssc'),
           'gfp' = which(names(.f_par$channels) == 'fl1'))
  
  for(.file in order(get_flowset_timestamp(.fs)) ) { # process by increasing acquisition time
    .ff <- fs.normT[[.file]]
    if (.verbose>1) cat(.file, ":", identifier(.ff), "\t", sep="")
    
    .od <- try(.f_od$od[.f_od$well==keyword(.ff, "TUBE NAME")], silent=TRUE) # fetch OD value if od has been loaded
    if (class(.od) == "try-error") .od <- NA # set OD to NA otherwise
    plot.name <- sprintf("%d: %s", .file, keyword(.ff, "ORIGINALGUID"))
    
    .spl <- preproc_facs_sample_secondary(.ff, .out_dir, .ch, .f_par, .od, .min_cells, 
                                          .write_format, .plot, plot.name)
    .gates <- rbind(.gates, .spl$gate)
    .preproc <- rbind(.preproc, .spl$preproc)
    .stats <- rbind(.stats, as.data.frame(.spl$stats))
  }
  if (.verbose>1) cat("\n\n")
  if (.plot) dev.off()
  # output stats
  write.csv (.stats, file=file.path(.out_dir, "stats.csv"))
  
  .pl <- list(gates=.gates, preproc=.preproc, stats=.stats)
  save('.pl', file=.cache_namer(.out_dir))
  return(.pl)
}

preproc_facs_sample_secondary <- function(.ff, .out_dir, .ch, .f_par, .od=NA, .min_cells=5000, 
                                .write_format=NULL, .plot=TRUE, plot.name ) {
# preproc_facs_sample_secondary performs gating and filtering operation that cannot be handled in flowSets
# it is primarily used by preproc_facs_plate and is not intended to be called directly.
  .g <- gate_fsc_ssc(.ff, .ch['fsc'], .ch['ssc'], .min_cells, .write_format)
  
  # writes gated values to file(s) using format(s) given in .write (values in NULL, "tab", "fcs") 
  if ('tab' %in% .write_format)
    write.FCS.tab(.g$ff, file.path(.out_dir, paste(file_path_sans_ext(identifier(.g$ff)), 'sub.txt', sep=".")))
  if ('fcs' %in% .write_format)
    write.FCS(.g$ff, file.path(.out_dir, paste(file_path_sans_ext(identifier(.g$ff)), 'sub.fcs', sep=".")))
  
  # Filter noisy values using a mixture model (gaussian + uniform)
  ff.g <- data.frame(exprs(.g$ff))
  names(ff.g) <- names(.ch)
  gfp_normal <- EM_normal_unif_mixture(ff.g$gfp, .p=0.99, .min=1, .max=log10(.f_par$facs.max)) #, .plot=TRUE)
  
  # prepare preprocessed data for output
  .path <- unlist(keyword(.g$ff,"FILENAME"), use.names=FALSE)
  .well <- unlist(keyword(.g$ff,"TUBE NAME"), use.names=FALSE)
  .ff_preproc <- data.frame(path=.path, well=unlist(keyword(.g$ff,"TUBE NAME"), use.names=FALSE),
                            ff.g, weight=gfp_normal$weights)
  
  # Compute stats summary
  .ff_stats <- list(path=.path, well=.well, n=dim(.g$ff)[1], od=.od)
  .t <- keyword(.g$ff, c("$BTIM","$ETIM"))
  names(.t) <- c('start', 'stop')
  .ff_stats <- c(.ff_stats, .t, stats_df_summary(data.frame(ff.g), .ch['gfp']))
  # stats <- c(stats, stats_var_summary(ff.g[, gfp_ch], gfp_normal$weights, trim, dimnames(ff.g)[[2]][gfp_ch]))
  .ff_stats <- c(.ff_stats, stats_var_summary(ff.g$gfp, gfp_normal$weights, .xname='gfp'))
  
  if (.plot) {
    facs_contour_with_polyg(.ff, .ch[c('fsc', 'ssc')], .f_par$lims, .polyg=.g$polygon)
    
    if (missing(plot.name)) plot.name <- keyword(.g$ff, "ORIGINALGUID")
    plot.label <- sprintf("m=%.2f  sd=%.2f  n=%d", gfp_normal$theta$mu1,  gfp_normal$theta$sigma1, dim(.g$ff)[1])     
    density_with_normal(exprs(.g$ff[, .ch['gfp']]), list(list(mu=gfp_normal$theta$mu1, sigma=gfp_normal$theta$sigma1, col='red')),
                             plot.name, plot.label, .lab="log10 GFP (AU)")
  }
  .gate <- .g$polygon
  try(names(.gate) <- c('fsc', 'ssc'), silent=TRUE)
  if (!is.null(.gate)) .gate <- data.frame(path=.path, well=.well, .gate)
  
  return(list(gate=.gate, preproc=.ff_preproc, stats=.ff_stats))
}

find_densest_area <- function(.x, .y, .prop, n.bins=50) {
# given a set of points (.x, .y), find_densest_area returns the coordinates of the polygon that 
# includes .prop of all points in the densest region of the (.x, .y) space:
# 1. compute the 2D density of points using the KernSmooth library on a n.bins x n.bins lattice
# 2. compute the exptal cumulative of the density over the grid
# 3. find the density value within which .prop of all points are included (here called .level)
  .dens <- KernSmooth::bkde2D(cbind(.x, .y), bandwidth=c(MASS::bandwidth.nrd(.x)/2, MASS::bandwidth.nrd(.y)/2), # need to reduce bandwidth to mimick kde2d
                              gridsize = c(n.bins, n.bins), range.x=list(range(.x), range(.y)))
  .dx <- diff(.dens$x1)[1]
  .dy <- diff(.dens$x2)[1]
  .da <- .dx * .dy
  
  .tot.dens <- sum(.dens$fhat) * .da
  .ord <- order(.dens$fhat, decreasing=TRUE)
  .ecdf <- cumsum(.dens$fhat[.ord])*.da / .tot.dens # sanity check: max must ≈ 1
  .level.idx <- max( which(.ecdf < .prop))
  .level <- .dens$fhat[.ord][.level.idx]
  .contours <- contourLines(.dens$x1, .dens$x2, .dens$fhat, levels=.level)
  .polyg <- matrix(c(.contours[[1]]$x, .contours[[1]]$y), ncol=2)
  return(.polyg)
  
  # plot(.ecdf)
  
  # .keep <- .ord[ which(.dens$fhat[.ord]>.level) ]    
  # .dens2 <- .dens
  # .dens2$fhat[.keep] <- max(.dens2$fhat)
  # image(.dens2$x1, .dens2$x2, .dens2$fhat)
  # lines(.polyg[[1]]$x, .polyg[[1]]$y)
  # points(.polyg[[1]]$x, .polyg[[1]]$y)
  
  # # find the convex hull of the polygon
  # .keep.x <- .keep %% n.bins
  # .keep.y <- floor(.keep/n.bins) + 1
  # .polyg <- chull(cells.dens$x1[.keep.x], cells.dens$x2[.keep.y])
  # points(cells.dens$x1[.keep.x][.polyg], cells.dens$x2[.keep.y][.polyg])
}

gate_fsc_ssc <- function(.ff, .fsc_ch, .ssc_ch, .min_cells=5000, .write_format='tab') {
# Apply a gate selecting .min_cells in the densest area of the (.fsc, .ssc) space.  
  .n_cells <- dim(.ff)[1]
  
  # Gate population of cells with similar fsc and ssc values (called .ff.S)
  if(.n_cells > .min_cells) {
    .polyg <- find_densest_area(exprs(.ff[, .fsc_ch]), exprs(.ff[, .ssc_ch]), .prop=.min_cells / .n_cells)
    colnames(.polyg) <- c(colnames(.ff)[.fsc_ch], colnames(.ff)[.ssc_ch]) # required for gating
    poly.gate <- polygonGate(filterId="similar.cells", .gate=.polyg)
    .ff_s <- Subset(.ff, poly.gate)
    colnames(.polyg) <- c('fsc', 'ssc') # convenient for export
  }  else { 
    .ff_s <- .ff
    .polyg <- NULL # matrix(c(contours[[1]]$x, contours[[1]]$y), ncol=2)
  }  
  return(list(polygon=.polyg, ff=.ff_s))
}

write.FCS.tab <- function (.x, .filename, .channel, .digits=NULL) {
# write.FCS.tab writes the raw values of the flowframe .x to a text file,
# possibly subsetting channels and specifying the number of digits used for the output
  if (!is(.x, "flowFrame")) 
    stop("Argument '.x' has to be a 'flowFrame'.")
  if (missing(.channel))
    .channel <- 1:dim(.x)[2] # if no channel specified, output all channels
  .m <- exprs(.x[, .channel])
  write.table(format(.m, digits=.digits), file=.filename, quote=FALSE, row.names=FALSE)
  # write.table is more flexible and not slower than write.matrix. format controls the number of digits
}

EM_normal_unif_mixture <- function(.xs, .p=0.5, .m, .sd, .min, .max, .plot=FALSE, .plot_title='iter #%d', .thresh=1e-3) {
# EM_normal_unif_mixture fits a mixture model with normal(.m, .sd) + uniform(.min, .max)
# to a univariate sample .xs using iterative EM.
# adapted from example code in http://commons.wikimedia.org/wiki/File:Em_old_faithful.gif
  
  E.step <- function(.xs, .theta)
    t(apply(cbind(
      .theta$tau[1] * dnorm(.xs, mean=.theta$mu1, sd=.theta$sigma1),
      .theta$tau[2] * dunif(.xs, min=.min, max=.max)
    ), 1, function(.x) .x/sum(.x)))
  
  M.step <- function(.xs, .T) { list(
    tau= apply(.T, 2, mean),
    mu1= weighted.mean(.xs, .T[, 1]),
    sigma1= sqrt(Hmisc::wtd.var(.xs, .T[, 1]))) }
  
  theta.dif <- function(.theta, .theta_old) {
    .dif <- abs(.theta$tau[1]-.theta_old$tau[1]) / (.theta$tau[1]+.theta_old$tau[1])
    .dif <- .dif + abs(.theta$mu1-.theta_old$mu1) / (.theta$mu1+.theta_old$mu1)
    .dif <- .dif + abs(.theta$sigma1-.theta_old$sigma1) / (.theta$sigma1+.theta_old$sigma1)
    .dif <- 2 * .dif 
    return(.dif)
  }
  
  plot.em <- function(.xs, .theta, .xpts, .main){
    hist(.xs, breaks=60, freq=FALSE, xlim=range(.xpts), main=.main)
    lines(.xpts, .theta$tau[1] * dnorm(.xpts, .theta$mu1, .theta$sigma1), col='red')
    lines(.xpts, .theta$tau[2] * dunif(.xpts, min=.min, max=.max), col='blue')
  }
  
  # if init params are missing, use statistics of the empirical distribution
  if (missing(.m)) .m <- mean(.xs) # min(.xs) + diff(range(.xs)) / 2
  if (missing(.sd)) .sd <- sd(.xs) # diff(range(.xs)) / 4
  if (missing(.min)) .min <- min(.xs)
  if (missing(.max)) .max <- max(.xs)
  .xpts <- seq(from=.min, to=.max, length.out=1000)
  
  theta <- list(
    tau=c(.p, 1-.p),
    mu1=.m,
    sigma1=.sd
  )
  
  # .t0 <- Sys.time()
  .dif <- 1
  .iter <- 0
  while (.dif > .thresh) {
    .iter <- .iter + 1
    if (.plot) plot.em(.xs, theta, .xpts,
                       .main=sprintf(.plot_title, .iter))
    
    T <- E.step(.xs, theta)
    theta_old <- theta
    theta <- M.step(.xs, T)
    .dif <- theta.dif(theta, theta_old)
  }
  # .dt <- Sys.time() - .t0
  
  return(list(theta=theta, weights=T[, 1]))
}

stats_df_summary <- function(.df, .lin_cols) {
  .stats <- summarise_each(.df, funs(mean, var))
  
  .s_lin <- .df %>% select(.lin_cols) %>% 
    mutate_each_(funs(lin=10^.), names(.lin_cols)) %>%
#     (function(.x) {names(.x) <- paste(names(.x), 'lin', sep='_'); .x}) %>%
    summarise_each_(funs(lin_mean=mean, lin_var=var), names(.lin_cols))
  if (length(.lin_cols) == 1) names(.s_lin) <- paste(names(.lin_cols), names(.s_lin), sep="_")
  .stats <- c(.stats, .s_lin)

  return(unlist(.stats))
}

stats_var_summary <- function(.xs, .w, .xname, .trim) {
  .df <- data.frame(.xs)
  # Compute weighted stats
  .stats <- summarise_each_(.df, funs(mean_noise_rm=wtd.mean(., .w),
                                      var_noise_rm=wtd.var(., .w),
                                      lin_mean_noise_rm=wtd.mean(10^., .w), 
                                      lin_var_noise_rm=wtd.var(10^., .w)), '1')
  
  # Not used: compute stats from trimmed data (Olin's legacy)
  if (!missing(.trim)) {
    .qt <- quantile(.xs, c(.trim, 1-.trim))
    .stats <- c(.stats, .df %>%
                  filter(.xs>.qt[1] & .xs<.qt[2]) %>% # xxx nse
                  #       filter(.dots=~.xs>.qt[1] & .xs<.qt[2]) %>% 
                  summarise_each_(funs(mean_trim=mean, var_trim=var), '1'))
  }
  
  # Rename stats variables
  if (!missing(.xname)) names(.stats) <- paste(.xname, names(.stats), sep="_")  
  return(unlist(.stats))
}

propagate_index_info <- function(.pls, .info) {
# propagate_index_info joins each dataframes included in the .pls list (output
# of preproc_facs_plate) with .info, a dataframe describing experimental
# conditions, etc (must have a column `dir` and either `well` or `discard`; 
# beware the fact that all other columns are appended to facs data, so better
# keep only relevant ones).
# in case where .info has a `well` column (1 row per well expected),
# inner_join() guaranties that wells that are not found in .info are discarded
# from all .pls dataframes.
# in case where .info has a `discard` column (1 row per plate expected), it
# is parsed (using , and ; to split) and corresponding wells are filtered out
# from all .pls dataframes.

  if (all(c('well', 'discard') %in% names(.info)))
    stop("propagate_index_info: .info data frame should have column `discard` or column `well`, not both.")
  
  # CASE: DISCARD WELLS
  if ('discard' %in% names(.info)) {
    # join dataframes with all wells
    .pls <- lapply(.pls, function(.df) {
      .df %>% 
        mutate(dir=dirname(as.character(path))) %>%
        inner_join(.info, by="dir")
    })
    
    # discard wells specified in .info
    .ws <- strsplit(.info$discard, split='[,;] *')
    .pl_discard <- .info %>%
      select(dir, discard) %>%
      # split well column and create one line for each value
      slice(rep(1:dim(.info)[1], times=sapply(.ws, length))) %>%
      mutate(discard=unlist(.ws))
    
    .pls <- lapply(.pls, function(.df)
      .df %>% 
        mutate(dir=dirname(as.character(path))) %>%
        filter(!interaction(dir, well) %in% interaction(.pl_discard$dir, .pl_discard$discard)) )
    return(.pls)
  }
  
  # CASE: SELECT WELLS
  if ('well' %in% names(.info))
    # join dataframes with well specified in .info only
    lapply(.pls, function(.df) {
      .df %>% 
        mutate(dir=dirname(as.character(path))) %>%
        mutate(well=as.character(well)) %>%
        inner_join(.info, by=c("dir", "well")) %>%
        mutate(well=factor(well))
    })
}

