rm(list = ls(all.names = TRUE)) # clears global environ

library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)

setwd("./Output")

######################################################

load('Media_Checkpoint.RData')

# Isolate cell culture and tissue metadata
Isolated_Methods_Data <- Isolated_Methods_Data[which(Isolated_Methods_Data$RNAseq_target == 'Cell_culture' |
                                                       Isolated_Methods_Data$RNAseq_target == 'Tissue'), ]

# Remove 'No details provided' examples
Isolated_Methods_Data <- Isolated_Methods_Data[which(Isolated_Methods_Data$Specific_media != 'No details provided' &
                                                       Isolated_Methods_Data$Culture_Type != 'No details provided'), ]

# Isolate tissue samples and unique ALI cultures by [series + cell type]
Tissue_samples <- Isolated_Methods_Data[which(Isolated_Methods_Data$Culture_Type == 'Tissue'), ]
Assign_media <- Isolated_Methods_Data[which(Isolated_Methods_Data$Culture_Type == 'ALI_Culture'), ]
Unique_Assign_media <- Assign_media[which(!duplicated(paste(Assign_media$series_id, Assign_media$Cell_Type))), ]

# Save the raw media details for manual assignment
# This is done based on a hierarchy for listed media based on known culturing procedures
# I.e. listing DMEM and UNC media consistently means a sample was collected or processed in DMEM
# before being grown/differentiated in UNC media
write.csv(Unique_Assign_media, file = 'Media_Assignment.csv')

# Load manually assigned media data
Unique_Assign_media <- read.csv('Assigned_Media.csv')
Unique_Assign_media <- Unique_Assign_media %>% select(-c(X))
Unique_Assign_media <- Unique_Assign_media[, c('series_id', 'Cell_Type', 'Media')]

Assigned_Methods <- merge(Isolated_Methods_Data, Unique_Assign_media, by = c('series_id', 'Cell_Type'))

####################################################

# Series removed due to unclear methods and high prevalence of treated/infected samples following filtering
Assigned_Methods <- Assigned_Methods[which(Assigned_Methods$series_id != 'GSE189613'), ]

# ALI time course data. Removed due to clear difference in sample population between samples (doi:10.1242/dev.177428 Supplemental Figure S3)
# Also, this is the only series which utilizes BEGM_DMEM for HNEC culturing (~1/2 of samples)
Assigned_Methods <- Assigned_Methods[which(Assigned_Methods$series_id != 'GSE121600'), ]

# Expanded tracheal aspirates from preterm and full term infants. Children and adults are known to have 
# different cell behaviors and this is the only series using UNC_ALI media for HNEC cultures
Assigned_Methods <- Assigned_Methods[which(Assigned_Methods$series_id != 'GSE164358'), ]

####################################################

# For each row
for(i in 1:nrow(Assigned_Methods)){
  
  # If a media is selected
  if(Assigned_Methods[i, 'Media'] == ''){
    
    # Isolate the current sample metadata
    GEO_id <- Assigned_Methods[i, 'geo_accession']
    Sample_title <- Control_GEO_Data[which(Control_GEO_Data$geo_accession == GEO_id), 'title'][1]
    
    # Change 'SC' to pneumacult and 'UNC' to unc_ali
    Assigned_Methods[i, ] <- Assigned_Methods[i, ] %>% mutate(Media = case_when(str_detect(Sample_title, 'SC') ~ 'pneumacult',
                                                                                str_detect(Sample_title, 'UNC') ~ 'unc_ali'))
    
  }
  
}

# Bind data together, remove duplicates, and count results before saving
GEO_Data <- bind_rows(Assigned_Methods, Tissue_samples)
GEO_Data[is.na(GEO_Data)] <- ''
GEO_Data_by_Series <- GEO_Data[!duplicated(paste(GEO_Data$series_id, GEO_Data$Cell_Type)), ]

Initial_Count <- GEO_Data %>% group_by(Cell_Type, Media) %>% summarize(n = n(), series_n = n_distinct(series_id))

save(GEO_Data, GEO_Data_by_Series, Control_GEO_Data, file = 'Media_Metadata.RData')
