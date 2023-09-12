# R Manuscript

This repo serves as data documentation for the submitted manuscript Education as risk factor of mild cognitive impairment: the link to the gut microbiome.

To set up your session, run packages.R which will download and install packages required for analyses but not yet installed on your machine.

The notebooks folder contains the markdown guiding you through the creation of outputs contained in the submitted manuscript. 
The workflow was as follows:
- the analysis pipeline is defined in the _targets.R file
- the markdown file loads required packages and initializes the pipeline, therby using functions defined in R/20230909_func.R
- completion of targets defined in the pipeline will store analysis outputs
- said outputs are consecutively loaded into the markdown for neat plots and tables

