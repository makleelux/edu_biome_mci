# R Manuscript

This repo serves as data documentation for the submitted manuscript Education as risk factor of mild cognitive impairment: the link to the gut microbiome.

To set up your session, run 20230908_packages.R which will download and install packages required for analyses that are not yet installed on your machine.

The notebooks folder contains the rendered html and the Rmd file it is dervied from. It is supposed to guide you through the creation of outputs presented in the submitted manuscript.

The workflow was as follows:
- the analysis pipeline is defined in the _targets.R file
- the markdown file loads required packages and initializes the targets pipeline, using functions defined in R/20230908_func.R
- results of targets defined in the pipeline will be stored in a dedicated folder
- said results are consecutively loaded into the markdown for neat plots and tables

