This package has been designed to preprocess FCM data by removing non viable cells, subtracting autofluorescence and reduce the effect of electronic noise coming from the measurement machine.

#Installing the package
As for any R package, just download the source folder and install it using the command
```{r, eval=FALSE}
install.packages("/path/to/vngFCM", repos=NULL, type="source")
```
Please notice that in order to compile the C++ sources GSL is required. In the cluster it must be loaded with the command
```
module load GSL
```

The package depends from 

* dplyr
* tidyr
* flowCore
* Rcpp

Optionally it takes advantage of *doParallel*,  if present, to parallelize the merging of all the results of the analysis in a single final data frame.

##Examples
In the folder example you can find a working source and data to show you how to use the package.

#How does the preprocessing work
##Remove the debris
The first step is to eliminate outbris from viable cells. To do this we suppose that viable cells must share similar geometrical properties and hence their scattering profile should be similar. We thus identify debris as cells which appear as outliers in the scattering profile of the population. 

Every signal measured by the FACS machine is just an electrical signal that creates an analog impulse in some detectors. This impulse is then represented by a number by considering the value of the highest point, the area under the curve of the impulse or the width of the impulse. Of these, only two are independent and since a priori for the scattering we don't have any reason to privilege one or the other, we consider both of them. The choice of what pair of signals is the indipendent one is arbitrary and we usually chose the height and the width of the electric impulse. 

In the end the scattering profile is measured by four signals, two coming from the forward scattering (e.g. the height and the width) and the corresponding two from the side scattering. To identify outliers we perform a 4D mixture of a gaussian (representing the main cloud of scattering points from viable cells) and a uniform distribution (representing the outliers). The mixture gives us the posterior probability for each cell to come from the gaussian or from the uniform and we discard all the cells whose posterior is lower than a cretain threshold that has to be set by the user.

##Extract mean and variance
To extract the mean and variance fluorescence from the cell population, we decided to focus on the height of the fluorescence impulse. The distribution of fluorescence is fit in the logarithmic space (to reduce the effect of outliers) with a mixture of a uniform and a normal distribution. Once we have the mean and variance from the normal fit, we can compute the mean and variance in the linear space by using the lognormal transformation formula.

##Remove autofluorescence and noise from the device
The means and variances we have are corrupted by the autofluorescence of the cells and by noise introduced by the electronics of the measuring apparatus. We remove the autofluorescence by subtracting the mean and variance of empty plasmids; if in the experiment there are more empty plasmids, we consider their average mean and variance.

The electronic noise is modeled as shot noise
$$I_m = I_t + \sqrt{I_t}\epsilon(0;\delta^2)$$
where $I_m$ is the measured signal, $I_t$ is the true signal not contaminated by the shot noise,and $epsilon$ is a random gaussian variable of mean 0 and variance $\delta^2$. 
The amplitude $\delta^2$ is machine dependent and must be set manually by the user. We have estimated it using artificial beads and we found that for our machine it is $\delta=14.51$.


#Setting up the environment
##Setting the location where to store the final preprocessed files
Before starting the preprocessing of the files, we need to setup the environment.
First of all we need to create a function that takes the name of the directory containing the FCM files and outputs the name of the directory where to store the preprocessed files. We can choose to store the preprocessed files in the same directory of the FCM files; in this case the function would be
```{r, eval=FALSE}
data2preproc <- function(.d) .d
```
Or we can decide to store the preprocessed files in a different directory; for example, if the directory containing the FCM is called *./data_facs/20160824/Cip0_3/*, we might want to store the final files a directory whose name is the same, except that *data_facs* is substituted by *data_facs_preprocessed*. In this case the function would be
```{r, eval=FALSE}
data2preproc <- function(.d) sub('data_facs', 'data_facs_preproc', .d) 
```

##Setting the parameters of the preprocessing algorithm
The preprocessing algorithm depends on some parameters that must be set manually. For this we need to create a list of parameters of 
the form
```{r, eval=FALSE}
f_par <- list(
  channels = c('FSC-H', 'FSC-W', 'SSC-H', 'SSC-W', 'GFP-H'),
  scattering_threshold=0.5,
  scattering_num_cells=NA,
  file_pattern = "puA139|puA66|recA|highexp|medexp|polB|recN|ruvA|ftsK|rpsU",
  delta_shot_noise = 14.51,
  autofluo_wells = pl_info
)
```
The meaning of the fields are

* *channels*: The name, or number, of the channels containing the relevant signals. They must be in the order:
  * The two signals from the forward scattering (typically the hight and the width).
  * The two signals from the side scattering.
  * The signal of the fluorescence (typically the height of the GFP channel).
* *scattering_threshold*: the minimum value of the posterior distribution of the scattering filter to consider a point to be a viable cell and not a debris. Typically is *0.5*.
* *scattering_num_cells*: number of cells to retain after the filtering. Only the \emph{scattering_num_cells} with the highest posterior probability are kept. If scattering_threshold is also specified, scattering_num_cells is ignored.
* *file_pattern*: the pattern of the FCM files to preprocess. Set '*' to analyse all the files.
* *delta_shot_noise*: the amplitude of the shot noise. This is machine dependent.
* *autofluo_wells*: a dataframe containing for every row the name of a directory of FCM files and the names of the wells with empty plasmids, separated by a comma.

To check if the names of the channels are correct, you can use the function $check.channels*.


##Using an external file to hold information about the directories
We suggest to use an external file to hold information about directories. For example you can have a CSV file that for every directory holds the names of the wells with empty plasmids and some comments which are ignored by the preprocessing algorithm, but can be useful to provide a description of the directory. For example we can have a file that looks like

```
#Example of a CSV file containing info about the directories
#Always remember to end the file with an empty line!
dir;empty;cip;date
./data_gwendoline/20160122_Cip_induction/20160122_Cip0;F4,F5,F6,E4,E5,E6;0;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_125;F4,F5,F6,E4,E5,E6;0.125;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_25;F4,F5,F6,E4,E5,E6;0.25;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_5;F4,F5,F6,E4,E5,E6;0.5;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_75;F4,F5,F6,E4,E5,E6;0.75;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip1;F4,F5,F6,E4,E5,E6;1;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_25;F4,F5,F6,E4,E5,E6;1.25;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_5;F4,F5,F6,E4,E5,E6;1.5;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_75;F4,F5,F6,E4,E5,E6;1.75;20160122
./data_gwendoline/20160122_Cip_induction/20160122_Cip2;F4,F5,F6,E4,E5,E6;2;20160122
```

We can then load the information using 
```{r, eval=FALSE}
pl_index_file <- "index_plates.csv"
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#")
pl_info <- pl_index %>% select(dir, empty, cip, date)
```

Please notice that the structure of the file is arbitrary and it can be adapted to the user needs. The preprocessing algorithm doesn't care how you obtain the list of directories to preprocess and the names of the wells with the empty plasmids! But be aware that if you want to merge this information with the preprocessed statistics, you need the field \empth{dir}.

##Preprocessing script
The filtering on the scattering is the most intensive part and it requires quite long time to be executed. To increase the speed the process has been parallelized using an array job. To do this the package needs to send the preprocessing script to each node in a cluster and thus this script cannot be inside the package. The external file must be placed somewhere in the local file system and its location must be provided later on when we will call the main preprocessing function.

#Starting the preprocessing
To start the preprocessing it is enough to call the function 
```{r, eval=FALSE}
preproc_facs_plates(.dirs, .data2preproc,
                    .cache_namer=(function(.d) file.path(.d, paste0(basename(.d), '_preproc_luca.Rdata'))),
                    .filter_preproc_namer=(function(.d) sub("fcs$", 'RData', .d)),
                    .f_par,
                    .plot=TRUE,
                    .pdf_dim=c(2, 2.4, 6, 8),
                    .verbose=FALSE,
                    .drop_preproc=TRUE,       
                    .force=FALSE,   
                    .preproc_script='preproc.R',
                    .use_doParallel=TRUE,
                    .jobs_cmd_name='filtering_cmd')
```
You are encouraged to look at the function documentation to understand all the parameters it depends on. The most important ones are 

* *.dirs*: a vector containing all the directories with the FCS files. using our CSV file, this is just *pl_nfo$dir*.
* *data2preproc* a function that takes the name of a directory containing the FCM files and returns the name of the directory where to store the preprocessed files.
* *.cache_namer*: a function that takes the name of a directory and returns the name of a file where to save a Rdata object containing the results of the preprocessing of the FCM files in that directory.
* *.f_par*: the list containing the parameters of the algorithm.
* *.drop_preproc*: should the algorithm stores the fluorescence and scattering values of all the cells which survived the scattering filtering? If TRUE this information is discarded, with an increase in speed and a reduction in disk space.
* *.force*: if some files have already been preprocessed, should we force the algorithm to preprocess them again?
* *.preproc_script*: location of the script containing the commands to be sent to the cluster. This is provided aong with the package.

What the function does is to loop over all the directories and to check for the FCM files that must be filtered according to the scattering profile. If a file has already been filtered, the result is written in a Rdata file in the directory specified by *.data2preproc* and it is not filtered again (unless *.force* is set). If some file need to be filtered, a job is sent to the cluster and the execution stops. When the cluster has finished the job the user must run again the function *preproc_facs_plates* to proceed with the autofluorescence and shot noise removal. Once the preprocesing is finish the results are saved in a Rdata file whose name is specified by *.cache_namer*. Specifically, the result is a list containing the field *stats* with the summary statistics of every FCM file, and *preproc*  with the raw scattering and fluorescence values of the cells which survived the scattering filter. If *.drop_preproc* is FALSE, then this field is NULL.

#After the preprocessing
Once you have called the preprocessing function with e.g.
```{r, eval=FALSE}
preproc <- preproc_facs_plates(pl_info$dir, data2preproc, .f_par=f_par, .plot=FALSE, .verbose=1, .force=FALSE, .drop_preproc=TRUE, .use_doParallel=FALSE)
```
you have a list \emph{preproc} which contains the preprocessed statistics and the raw data with the outliers removed.

##Add additional information
Now you may want to add additional information to the statistics, for example the date and the cip concentration. If the information applies to whole directories, the easiest way to do this is to use the function \emph{propagate_index_info} which takes information contained in a tibble and merge it to the statistics by matching the column \emph{dir}. In our example we have dates and cip concentrations stored in the csv file and loaded in the tibble \emph{pl_info} and their values are the same inside a specific director; we merge them using
```{r, eval=FALSE}
# Info we want to merge
info <- dplyr::select(pl_info, dir, date, cip)
# Merge with the preproc statistics
preproc <- propagate_index_info(preproc$stats, info) 
# Remove columns that we don't need anymore and transform to a tibble
preproc <-  preproc %>% select(-dir) %>% as_data_frame()
```

Information which doesn't apply globally to whole directories must be inserted by hand. For example if we want to add the promoter name we need to manually extract it from the directory name
```{r, eval=FALSE}
preproc <- preproc %>% extract_('path', 'promoter', '[0-9]/([[:alpha:]]+)_[A-Z]')
```

##Postprocessing
Now \emph{preproc} stores a tibble containing the preprocessed statistics for each .fcs file and additional information like date and cip concentration and the promoter name. We are then ready to start the analysis of our data; as an example, you can plot the variance of each promoter as a function of its mean expression
```{r, eval=FALSE}
ggplot(preproc) + geom_point(aes(fl_mean_noise_rm, fl_var_noise_rm, col=factor(promoter)))
```
