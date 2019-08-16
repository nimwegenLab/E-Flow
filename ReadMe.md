# Fluocytometer analysis of bacterial cells

## Purpose of the package
This package has been designed to preprocess FCM data by removing non viable cells, subtracting autofluorescence and reduce the effect of electronic noise coming from the measurement machine.

## Installing the package
As for any R package, just download the source folder and install it using the command
```
install.packages("/path/to/vngFCM", repos=NULL, type="source")
```
Please notice that in order to compile the C++ sources GSL is required. 

The package depends from 

* dplyr
* tidyr
* flowCore
* Rcpp

Optionally it takes advantage of *doParallel*,  if present, to parallelize the merging of all the results of the analysis in a single final data frame.

## Examples
In the folder example you can find a working source and data to show you how to use the package.

## Brief theoretical introduction
### Remove the debris (optional)
Although in the companion paper we have shown that for bacterial cells is difficult to separate good and bad data points, the algorithm offers the possibility to filter the particles based on their scattering profile. The idea is that viable cells must share similar geometrical properties and hence their scattering profile should be similar. We thus identify debris as cells which appear as outliers in the scattering profile of the population. 

Every signal measured by the FACS machine is just an electrical signal that creates an analog impulse in some detectors. This impulse is then represented by a number by considering the value of the highest point, the area under the curve of the impulse or the width of the impulse. Of these, only the width and the height are independent and since a priori for the scattering we don't have any reason to privilege one or the other, we consider both of them.

In the end, the scattering profile is measured by four signals, two coming from the forward scattering (height and width) and the corresponding two from the side scattering. To identify outliers, we perform a 4D mixture of a Gaussian (representing the main cloud of scattering points from viable cells) and a uniform distribution (representing the outliers). The mixture gives us the posterior probability for each cell to come from the gaussian or from the uniform and we discard all the cells whose posterior is lower than a threshold that has to be set by the user (e.g. 0.5).

### Extract mean and variance
To extract the mean and variance of the fluorescence from the cell population, we showed in the companion paper that the height is the most informative component. The distribution of fluorescence is fit in the logarithmic space with a mixture of a uniform and a normal distribution. Once we have the mean and variance from the normal fit, we can compute the mean and variance in the linear space.

### Remove autofluorescence and noise from the device
The means and variances are corrupted by the autofluorescence of the cells and by noise introduced by the electronics of the measuring apparatus. We remove the autofluorescence by subtracting the mean and variance of empty plasmids; if in the experiment there are more empty plasmids, we merge their means and variances as explained in the companion paper.

The electronic noise is modeled as shot noise with an offset $O$
$$I_m = I_t + B_T \sqrt{I_T+B_T}\epsilon(0;\delta^2)+O$$
where $I_m$ is the measured signal, $I_t$ is the true signal not contaminated by the shot noise,and $epsilon$ is a random gaussian variable of mean 0 and variance $\delta^2$. $B_T$ is the true autofluorescence signal, which the algorithm infers looking at the fluorescence of non-expressing cells.
The amplitude $\delta^2$ and the offset $O$ are machine dependent and must be set manually by the user. We have estimated them using artificial beads and we found that for our machine $\delta=12.7\pm 0.6$ and $O=97/pm 29$.

Notice that the shift $O$ affects also the observed autofluorescence. In the companion paper we show that for this reason $O$ simplifies from the equations once we provide the observed values of non-expressing cells.


## Setting up the environment
### Setting the location where to store the final preprocessed files
Before starting the preprocessing of the files, we need to setup the environment.
First of all, we need to create a function that takes as argument the name of the directory containing the FCM files and outputs the name of the directory where to store the preprocessed files. We can choose to store the preprocessed files in the same directory of the FCM files; in this case the function would be
```
data2preproc <- function(.d) .d
```
Or we can decide to store the preprocessed files in a different directory; for example, if the directory containing the FCM is called *./data_facs/20160824/Cip0_3/*, we might want to store the final files a directory whose name is the same, except that *data_facs* is substituted by *data_facs_preprocessed*. In this case the function would be
```
data2preproc <- function(.d) sub('data_facs', 'data_facs_preproc', .d) 
```

### Setting the parameters of the preprocessing algorithm
The preprocessing algorithm depends on some parameters that must be set manually. For this we need to create a list of parameters of 
the form
```
f_par <- list(
  channels = c('FSC-H', 'FSC-W', 'SSC-H', 'SSC-W', 'GFP-H', 'GFP-W'),
  scattering_threshold=0.5,
  scattering_frac_cells=NA,
  file_pattern = "fcs$", 
  delta_shot_noise = 12.73,
  delta_shot_noise.err = 0.63
)
```
The meaning of the fields are

* *channels*: The name, or number, of the channels containing the relevant signals. They must be in the order:
  * The two signals from the forward scattering (in our example the hight and the width).
  * The two signals from the side scattering.
  * The signal for the fluorescence (here, the height of the GFP channel).
* *scattering_threshold*: the minimum value of the posterior distribution of the scattering filter to consider a point to be a viable cell and not a debris. Typically is *0.5*. Set it to $0$ to not perform any filtering).
* *scattering_frac_cells*: fraction of cells to retain after the filtering. Only the fraction \emph{scattering_frac_cells} with the highest posterior probability is kept. If scattering_threshold is also specified, scattering_frac_cells is ignored.
* *file_pattern*: the pattern of the FCM files to preprocess. In our example we analyse all files with extension *fcs*. Set '*' to analyse all the files.
* *delta_shot_noise*: the amplitude of the shot noise. This is machine dependent and can be determined e.g. using calibration beads.
* *delta_shot_noise.err*: the error on the amplitude. Set it to $0$ if it is not known. 

To check if the names of the channels are correct and corresponds to measurements in the fcs files, you can use the function *check.channels*.


### Using an external file to hold information about the directories
We suggest to use an external file to store information about directories. For example you can have a CSV file with a list of directories containing the files to be analysed and for every directory stores optional comments which are ignored by the preprocessing algorithm, but can be useful to provide a description of the directory. For example we can have a file that looks like

```
#Lines starting with # are ignored by the script
#Example of a CSV file containing info about the directories
#Always remember to end the file with an empty line!
dir;antibiotic;date;author
./data_gwendoline/20160122_Cip_induction/20160122_Cip0;author1
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_125;author1
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_25;author2
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_5;author1
./data_gwendoline/20160122_Cip_induction/20160122_Cip0_75;author3
./data_gwendoline/20160122_Cip_induction/20160122_Cip1;author2
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_25;author1
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_5;author2
./data_gwendoline/20160122_Cip_induction/20160122_Cip1_75;author1
./data_gwendoline/20160122_Cip_induction/20160122_Cip2;author3
```

We can then load the information using 
```
pl_index_file <- "index_plates.csv"
pl_index <- read.csv2(pl_index_file, stringsAsFactors=FALSE, comment.char="#")
pl_info <- pl_index %>% select(dir, antibiotic, date, author)
```

Please notice that the structure of the file is arbitrary and it can be adapted to the user needs. But be aware that if at the end you want to merge this information with the final output of the algorithm, you need the field \empth{dir}.

## Preprocessing step
The filtering on the scattering is the most intensive part and it requires quite long time to be executed. To increase the speed, the process has been parallelized using an array job. To do this the package needs to send the preprocessing script to each node in a cluster and thus this script cannot be inside the package. The external file must be placed somewhere in the local file system and its location must be provided later on when we will call the main preprocessing function.

### Starting the preprocessing
To start the preprocessing it is enough to call the function 
```
analyse_raw(.dirs, 
            .data2preproc, 
            .filter_preproc_namer = (function(.d) sub("fcs$", "RData", .d)), 
            .f_par, 
            .force = FALSE, 
            .preproc_script = "preproc.R", 
            .jobs_cmd_name = "filtering_cmd") 
```
You are encouraged to look at the function documentation to understand all the parameters it depends on. The most important ones are 

* *.dirs*: a vector containing all the directories with the FCS files. using our CSV file, this is just *pl_nfo$dir*.
* *data2preproc* a function that takes the name of a directory containing the FCM files and returns the name of the directory where to store the preprocessed files.
* *.f_par*: the list containing the parameters of the algorithm.
* *.force*: if some files have already been preprocessed, should we force the algorithm to preprocess them again?
* *.preproc_script*: location of the script containing the commands to be sent to the cluster. This is provided aong with the package.

What the function does is to loop over all the directories and to check for the FCM files that must be filtered according to the scattering profile. For each file, the result is written in a Rdata file in the directory specified by *.data2preproc*. To speedup the process, the scattering commands for each file are collected in an array job that is sent to the cluster. In case the machine doesn't have access to a cluster, the user can manually execute the commands in the file generated by the function. The name of the file is given by the parameter *.filter_preproc_namer*.

For our example we would write
```
analyse_raw(pl_info$dir, data2preproc, .f_par=f_par)
```

If the function sends the commands to a cluster, please wait for the job to be completed before proceeding.
When the preprocessing is terminated, the results of all the single files can be pooled together using the function *collect_raw*.
```
preproc <- collect_raw(pl_info$dir, data2preproc, .f_par = f_par)
```
Preproc is a list containing different fields

* *stats*: contains the summary statistics of the populations of cells.
* *preproc*: contains for each directory the scattering and fluorescence information of the single cells. It also contains a field *post* that gives the posterior probability for a given cell to come from the main cloud of the scattering profile.

Notice that to reduce the size, the user can choose to discard the *preproc* field by passing *.drop_preproc=TRUE* when calling the function.

### Add additional information
Now you may want to add additional information to the statistics, for example the date, author and antibiotic used. If the information applies to whole directories, the easiest way to do this is to use the function \emph{propagate_index_info} which merges the information with the output of the algorithm by matching the column \emph{dir}. In our example we have dates, antbiotic and author stored in the tibble \emph{pl_info} and their values are the same inside a specific director; we merge them using
```
# Info we want to merge
info <- dplyr::select(pl_info, dir, date, cip)
# Merge with the preproc statistics
preproc <- propagate_index_info(preproc$stats, info) 
# Remove columns that we don't need anymore and transform to a tibble
preproc <-  preproc %>% select(-dir) %>% as_data_frame()
```

Information which doesn't apply globally to whole directories must be inserted by hand. For example, if we want to add the promoter name and it is contained in the name of the fcs file, we can manually extract it from the fileanme using regular expressions
```
preproc <- preproc %>% extract_('path', 'promoter', '[0-9]/([[:alpha:]]+)_[A-Z]')
```

## Autofluorescence and noise correction
### Autofluorescence removal
Now we can start analysing the autofluorescence by looking at the fluorescence of non-expressing cells. In our example, these cells are the ones with the promoter *puA139* and we have two different dates. The algorithm allows to merge the fluorescence from different wells and dates as explained in the companion paper. The function to merge them is *get_autofluo_stats* and takes as arguments a vector of means and a vector of variances to merge.
```
autofluo <- preproc$stats %>% filter(promoter=='puA139') %>% do((function(.df){
  .m <- get_autofluo_stats(.df$fl_mean_lin, .df$fl_mean_lin.err)
  .v <- get_autofluo_stats(.df$fl_var_lin, .df$fl_var_lin.err)
  return(data.frame(autofluo_mean = .m['mu'], autofluo_mean.err = .m['tau'], autofluo_var = .v['mu'], autofluo_var.err = .v['tau']))
})(.)) 
```
Notice that if we only have one measurement of autofluorescence, we don't need to run this function.

Once we have the mean, variance and their corresponding errors for the autofluorescence, we need to add them to the output of the algorithm
```
preproc$stats <- preproc.all$stats %>% mutate(autofluo_mean = autofluo$autofluo_mean, 
                                              autofluo_mean.err = autofluo$autofluo_mean.err,
                                              autofluo_var = autofluo$autofluo_var, 
                                              autofluo_var.err = autofluo$autofluo_var.err)
```

### Shot noise removal
At this point we have all we need to remove the autofluorescence and the shot noise. Please verify that in the output of the algorithm, in our example *preproc*, there are the following fields

* *fl_mean_lin, fl_mean_lin.err*: the mean fluorescence and its error in linear space (i.e. not in log space).
* *fl_var_lin, fl_var_lin.err*: the variance of the fluorescence and its error in linear space (i.e. not in log space).
* *autofluo_mean, autofluo_mean.err*: the mean and its error for the autofluorescence in the linear space.
* *autofluo_var, autofluo_var.err*: the variance and its error of the autofluorescence in the linear space.

If some of these parameters are missing, some error eccor in the previous steps.

To remove the autofluorescence and the shot noise you have to call the function *remove_noise*
```
preproc$stats <- remove_noise(preproc$stats, f_par) 
```
The result is a tibble with the following fields

* *fl_mean*, *fl_mean.err*, *fl_var*, *fl_var.err*: mean, variance and their standard errors in log space. Without any autofluo or noise removal.
* *fl_mean_lin*, *fl_mean_lin.err*, *fl_var_lin*, *fl_var_lin.err*: same as above, but in linear space. Without any autofluo or noise removal.
* *fl_mean_lin_autofluo_rm*, *fl_mean_lin_autofluo_rm.err*, *fl_var_lin_autofluo_rm*, *fl_var_lin_autofluo_rm.err*: statistics in linear space with autofluo removed. The shot noise is still present. Notice that the noise removal doesn't affect the mean.
* *fl_var_lin_noise_rm*, *fl_var_lin_noise_rm.err*: variance of the fluorescence with autofluo and shot noise removed. Since the shot noise doesn't affect the mean, the autofluo removal is the only correction to the mean expression.

For later reference, the function returns also the value of the $\delta$ used to remove the shot noise.
