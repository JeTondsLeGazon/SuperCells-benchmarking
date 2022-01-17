# README

This repo contains the necessary tools for a differential expression analysis and benchmarking
for the package ['SuperCell'](https://github.com/GfellerLab/SuperCell) in R.

The datasets used were found in a [paper](https://www.nature.com/articles/s41467-021-25960-2) analysing the false discoveries of single-cell differential expression.


## Folder structure

.  
|----- data  
|----- configs  
|----- figs  
|----- tests  
|----- src  

## Dependencies

... in progress ...


## Data

All data are available [here](https://doi.org/10.5281/zenodo.5048449). Run:
```
bash fetchData.sh
```
in the project folder in order to fetch the data.

Alternatively, you can download the files 'bulk_rnaseq.tar.gz' and 'sc_rnaseq.tar.gz', unzip them and place them in the data folder.


## Run files

The pipeline is runnable via the run files: 'runProcessing.R', 'runDE.R' and 'runAnalysis.R', respectively.

Any file can be run with:
```
Rscript runFile.R configs/myConfigFile 
```
where both 'runFile.R' and 'configs/myConfigFile' should be replaced with the corresponding files. The configuration file contains all the meta-data that should be used for a pipeline.

Alternatively, one can source the file in Rstudio directly but should replace under HEADER the line:
```
args <- commandArgs(trailingOnly = TRUE)
```
by:
```
args <- "myConfigFilePath"
```

**Note**: Each run files includes manually some library paths, which may be needed in some cases. Please comment or replace under HEADER the line:
```
.libPaths("C:/Users/miche/OneDrive/Documents/R/win-library/4.1")
```
with what is needed accordingly.
