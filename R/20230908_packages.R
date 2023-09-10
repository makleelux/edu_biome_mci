# packages
# adapted from: https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# note that some of the bioconductor packages will cause issues with 
# R-4.3.1-arm64.pkg for Apple silicon (M1/M2) Macs
# You may have to install R-4.3.1-x86_64.pkg

# check and install packages available via CRAN/Github
list.of.packages <- c('devtools',
                      'targets', 
                      'tarchetypes', 
                      'readxl', 
                      'tidyverse', 
                      'Hmisc', 
                      'reshape2', 
                      'vegan', 
                      'CMAverse',
                      'LDM',
                      'cowplot',
                      'ggpubr', 
                      'kableExtra',
                      'ggbeeswarm',
                      'gridExtra')

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BiocParallel")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if('devtools' %in% new.packages) install.packages("devtools")
if('CMAverse' %in% new.packages) devtools::install_github("BS1125/CMAverse")
new.packages = setdiff(new.packages, c('devtools', 'CMAverse'))
if(length(new.packages)) utils::install.packages(new.packages)

# check and install packages available via Bioconductor
list.of.bioc.packages <- c('phyloseq',
                           'DESeq2',
                           'densvis',
                           'scater',
                           'ANCOMBC')

new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages) 

