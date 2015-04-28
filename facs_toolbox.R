# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# biocLite(c("flowCore", "flowViz"))
# install.packages(c('Hmisc', 'tidyr', 'dplyr', 'ggplot2')) # 'tools'

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

# flowViz config
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

# ggplot2 config
theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")
# to revert to the default ggplot2 discrete colour scale, use: + ggplot2::scale_colour_discrete()
scale_colour_periodic_brewer <-
  function(...) scale_colour_manual(..., values = rep(c(brewer.pal(4, 'Set1'), 'gray42'), 100))
scale_shape_periodic <- 
  function(...) scale_shape_manual(..., values = rep(15:18, 5))


