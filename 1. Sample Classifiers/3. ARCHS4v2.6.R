rm(list = ls(all.names = TRUE)) # clears global environ

library(rhdf5)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(GEOquery)

#################################################

setwd("./Input")

ARCHS_file_location = 'human_gene_v2.6.h5'

destination_file = ARCHS_file_location # Originally accessed on 8/4/2025

# Retrieve information from compressed data
Samples = as.data.frame(h5read(destination_file, '/meta/samples'))

h5closeAll()

Samples <- Samples[, c('sample', 'series_id', 'characteristics_ch1', 'source_name_ch1', 'title')]

# Unite methods metadata together for BERT prediction
Samples <- Samples %>% unite('Combined_Methods', characteristics_ch1:title, sep = ' ', remove = FALSE)

setwd("./Output")

save(Samples, file = 'ARCHS4v2.6_Samples.RData')