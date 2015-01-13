f_par <- list(
  channels = c(fsc=2, ssc=5, # channels to analyze
               fl1=8), #, fl2=6),
  lims = c(2.5, 4.5), # limits for fsc and ssc (in log10 scale)
  min.cells = 5000,  # minimum number of cells to include in preprocesing
  facs.max = 262143, #
  file.pattern = '.fcs$')
