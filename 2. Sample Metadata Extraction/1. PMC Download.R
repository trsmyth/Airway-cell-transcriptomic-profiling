rm(list = ls(all.names = TRUE)) # clears global environ

library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)

`%ni%` <- negate(`%in%`)

setwd("./Input")

Sample_Classification_Results <- read.csv('Sample_Classification_Results.csv')

# PMC open access file list containing: File, Article Citation, Accession ID, Last Updated, PMID, License
PMC_dataset <- read.csv('oa_file_list.csv')

setwd("./Output")

# Get unique PMIDs
PMID <- unique(Sample_Classification_Results$PMID)
PMID <- PMID[which(PMID != ' ')]
PMID <- unique(unlist(strsplit(PMID, ' ')))

# Subset PMC_dataset to match unique PMIDs
Paper_database_locations <- PMC_dataset[which(PMC_dataset$PMID %in% PMID), ]

# Get PMIDs not in the oa database
PMID_of_papers_not_in_database <- PMID[PMID %ni% Paper_database_locations$PMID]

#################

url <- Paper_database_locations$File
destination <- Paper_database_locations$PMID

Paper_database_locations <- Paper_database_locations[, c('PMID', 'Accession.ID')] %>% setNames(c('PMID', 'PMCID'))

# Download the paper associated with each OA PMID
for(i in seq_along(url)){

  download.file(paste0('https://ftp.ncbi.nlm.nih.gov/pub/pmc/', url[i]),
                paste0('./Output/PMC Files/', destination[i]),
                mode = "wb")

}

#################

library(RSelenium)
library(rvest)
library(xml2)

# Set up firefox as driver
rD <- rsDriver(browser = "firefox",
               chromever = NULL,
               phantomver = NULL)

# Set up client
remDr <- rD$client

Master_dir <- './Output/Papers'
setwd(Master_dir)

# For papers which are not OA and cannot be directly downloaded, automatically navigate to the
# PMC page for that PMID and manually download the manuscript and supplemental data if available
PMCID <- lapply(PMID_of_papers_not_in_database, function(single_paper){
  
  # Print which paper is being located
  print(paste0('Locating paper ',
               match(single_paper, PMID_of_papers_not_in_database), 
               ' of ', 
               length(PMID_of_papers_not_in_database)))
  
  # Navigate to the pmc page for the current PMID
  remDr$navigate(paste0("https://pubmed.ncbi.nlm.nih.gov/", single_paper))

  # Manually add the PMCID
  PMCID <- c((readline(prompt = "Enter PMCID: ")))
  
  # If a PMCID is found, create a folder to store the downloaded paper
  if(PMCID != ''){

    print(PMCID)
    dir.create(file.path(paste0(Master_dir, '/', PMCID)))
    readline(prompt = "Press [enter] to continue")

  }

  # Else make a folder named after the PMID
  else{

    dir.create(file.path(paste0(Master_dir, '/', single_paper)))
    readline(prompt = "Press [enter] to continue")

  }

  PMCID
  
})

# Create a df storing PMCID data from manual downloads
PMCID <- data.frame('PMID' = PMID_of_papers_not_in_database, 
                    'PMCID' = do.call(rbind, PMCID))

# Add all data together and save
Paper_IDs <- bind_rows(lapply(Paper_database_locations, as.character), lapply(PMCID, as.character)) %>% data.frame()
rownames(Paper_IDs) <- seq(1, nrow(Paper_IDs), by = 1)

save(Paper_IDs, file = 'Paper_IDs.RData')
