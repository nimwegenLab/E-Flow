library(dplyr)
#To read into the file downloaded directly from Uri Alon's webpage that contains all the information about the plate. 
#The file is located in the facs_toolbox and the .join variable refers to to which dataframe of the final list you
#want to add the information about the genes. 
read_combine_gene_name_file <- function(file ="~/projects/facs_toolbox/alon_all_strains_forweb.csv", .join){
  alon_pwg <- read.csv(alonpwg_path) %>%
    select (Plate_Number, Well, Gene_Name) %>%
    rename(plate_upper=Plate_Number, well=Well, gene=Gene_Name)   
  
  plate <- data.frame(plate = tolower(alon_pwg$plate_upper))
  
  alon_pwg <- cbind(plate, alon_pwg)
  alon_pwg$plate_upper <- NULL
  left_join(.join, alon_pwg, by = c("plate", "well"))
}
