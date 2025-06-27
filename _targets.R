library(targets)
library(Seurat)
library(methods)
library(dplyr)
source("R/sc_workflow.R")

# StandaSeurat# Standard start of a targets files

# Load packages needed for the workflow and set options for them.
tar_option_set(
  packages = c("tidyverse", "rmarkdown", "knitr", "kableExtra", "tinytex", 
               "here", "fs", "janitor", "readxl", "readr", "stringr", "lubridate", "ggplot2",
               "ggpubr", "ggrepel", "ggthemes", "ggsci", "ggforce", "cowplot", "patchwork",
               "scales", "gridExtra", "grid", "ggridges", "gghighlight", "ggbeeswarm",
               "ggrepel", "ggtext", "gganimate", "magick", "ggdark",
               "ggpointdensity", "ggalluvial", "ggforce"))

# Load data
dir <- "~/tank/RProjects/Kashkar"

samples <- paste0(dir,list("/Casp8_CS", 
                           "/Casp8_WT_REP1",
                           "/Casp8_WT_REP2", 
                           "/Casp8_CS_REP2", 
                           "/Casp8_WT_REP3", 
                           "/Casp8_WT_REP4", 
                           "/Casp8_CS_REP3", 
                           "/Casp8_KO_REP1", 
                           "/Casp8_KO_REP2")) %>% as.list()

print(samples)

sample_names <- list("CSR1", "WTR1", "WTR2", "CSR2", "WTR3", "WTR4", "CSR3", "KOR1", "KOR2")

# Create an instance of Seurat_Strategy
seurat_strategy <- new("ConcreteSeuratWorkflow")

target_list <- vector("list", length(sample_names))

resolutions <- c(.75, .75,.75,.5,.5, .5, .5, .5, .5)
for(sn in 1:length(sample_names)){ 
  # Generate the targets using the run_analysis method
  target_list[[sn]] <- run_analysis(seurat_strategy, 
                                    directory=samples[[sn]], 
                                    resolution = resolutions[sn], 
                                    dims = 1:15,
                                    id=sample_names[[sn]]
  )
}

target_list
