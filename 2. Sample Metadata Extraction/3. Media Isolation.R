rm(list = ls(all.names = TRUE)) # clears global environ

options(java.parameters = "-Xmx8000m")

library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(pdftools)
library(tabulapdf)

setwd("./Functions")

source('Extract_methods_section.R')
source('Isolate_text_around_keyword.R')
source('Isolate_Single_Series_Methods.R')

`%ni%` <- negate(`%in%`)

setwd("./Output")

###################################################################################################

load('PDF_Checkpoint.RData')

# Define filtering terms
Cell_culture_terms <- list('\\bmedia\\b', '\\bmedium\\b', 'culture', 'culturing', 'propagate', 'grown',
                           'confluency', 'seeded', 'plated',  'expanded', 'monolayer', '\\bali\\b', 'air-liquid', 'air liquid')

ALI_terms <- list('air-liquid', 'air liquid', '\\bali\\b', 'semipermeable', 'transwell', 'snapwell', 'apical chamber', 'basal chamber', 
                  'basolateral chamber', 'two-chamber', 'two chamber', 'bio-one', 'bio one')

Submerged_terms <- list('submerged', 'monolayer', 'dish', '\\bplastic\\b', 'phenoplate')

RNA_terms <- list('rna seq', 'rna-seq', 'rnaseq', 'trizol', 'lysis buffer', 'lysing buffer', 'lysed', 'omic', 'sequencing', 'library generation', '\\brna\\b')

Single_cell_rnaseq <- list('single cell rna-seq', 'single cell rna seq', 'single cell rnaseq', 
                           'single-cell rna-seq', 'single-cell rna seq', 'single-cell rnaseq', 
                           'scrna-seq', 'scrna seq', 'sc rnaseq', 'sc-rnaseq', 'scrnaseq')

Tissue_terms <- list('brush', 'scrape', 'scraping', 'bronchoscopy', 'curettage', 'surgery', 'surgical', 
                     'polypectomy', 'lobectomy', 'biopsy', 'biopsies', 'turbinate', 'turbinoplasty', 'whole lung', 
                     'snap-frozen', 'snap frozen', 'frozen tissue', 'polyp', 'swab', 'autopsy', 'in vivo', 'in-vivo')

Cell_terms <- list('hbec', 'hbtec', 'hnec', 'saec', '\\bali\\b', 'air-liquid', 'air liquid', 'culture', 'grown in', 
                   'semipermeable', 'transwell', 'differentiated', 'seeded', 'matrigel', 'stimulation', 'stimulated', 'confluen',
                   'exposed', 'exposure', 'treated', 'transfect', 'transfected', 'in vitro', 'in-vitro')

Transwell_terms <- list('semipermeable', 'transwell', 'snapwell', 'two-chamber', 'two chamber', 'bio-one', 'bio one', 'thincert', 'millicell')

Growth_surface <- list(pet = list('transwell-clear', 'transwell clear', '\\bpet\\b', 'polyester', 'polyethylene', 
                                  '\\b3801\\b', '\\b3450\\b', '\\b3460\\b', '\\b3470\\b', 'thincert', 'pe membrane'), 
                       
                       polystyrene = list('polystyrene'), 
                       
                       polycarbonate = list('polycarbonate', '\\bpc\\b', 'pihp01250', 
                                            '\\b3401\\b', '\\b3412\\b', '\\b3413\\b', 
                                            '\\b7910\\b', '\\b3407\\b', '\\b3801\\b'), 
                       
                       ptfe = list('ptfe', 'polytetrafluoroethylene', 'picm01250'), 
                       
                       bio_one = list('bio one', 'bio-one'), 
                       millicell = list('millicell', 'millipore'),
                       transwell = list('corning', 'costar', 'falcon', 'transwell'))

Media_terms <- list(ck_dci = list('ck\\+dci', 'ck[[:punct:]]\\+dci', 'ck_dci'), 
                    dccm = list('dccm'), 
                    dmem = list('dmem', 'dulbecco'), 
                    f12 = list('f12', 'f-12'), 
                    rpmi = list('rpmi'),
                    f_media = list('f-media', 'f-medium', 'f_media', 'f_medium'),
                    lhc_8 = list('lhc-8', 'lhc_8'),
                    mucilair = list('mucilair'), 
                    gray = list('gray'),
                    b_ali = list('b-ali', 'lonza ali', 'b_ali'),
                    lm_0050 = list('lm-0050', 'hbtec ali differentiation medium'),
                    bronchialife = list('bronchialife', 'bleam'),
                    keratinocyte = list('keratinocyte'),
                    pneumacult = list('pnuemacult', 'pneumacult'), 
                    aec_medium = list('airway epithelial cell basal', 'airway epithelial cell growth', 
                                      'airway epithelial culture media', 'aec growth medium', 'aec medium', 
                                      'aec_medium', 'basal epithelial growth medium'),
                    begm = list('bronchial epithelial growth', 'bronchial epithelium growth', 'bulletkit', 'singlequot', 'begm'),
                    bebm = list('bronchial epithelial basal', 'bronchial epithelium basal', 'bebm'),
                    sagm = list('small airway epithelial cell growth', 'small airway growth basal', 'sagm'),
                    unc_ali = list('lhc:basal/dmem', 'lhc basal medium and dmem', 'lhc basal medium/dmem', 
                                   'lhc basal/dmem', 'unc ali', 'unc_ali'), 
                    
                    Confirm_ALI_Media = list('ali media', 'ali medium', 'air liquid interface media', 'air liquid interface medium',
                                             'air liquid interface (ali) media', 'air liquid interface (ali) medium',
                                             'air-liquid interface media', 'air-liquid interface medium',
                                             'air-liquid interface (ali) media', 'air-liquid interface (ali) medium'))

# Create a progress bar
pb = txtProgressBar(min = 0, max = length(unique(United_Metadata[, 'geo_accession'])), initial = 0, style = 3) 

# Isolate metadata information regarding sample type and culture conditions
Isolated_Methods_Data <- lapply(unique(United_Metadata[, 'geo_accession']), function(Single_sample){
  
  Single_series <- United_Metadata[which(United_Metadata$geo_accession == Single_sample), ] # Isolate single sample
  Single_method <- tolower(unique(unlist(Single_series$Method))) # Isolate methods details of sample
  
  ###################
  
  # Replace predetermined key words with 'Cell_culture' or 'Tissue' to ID method describes cell culturing or tissue biopsies
  RNAseq_method <- Single_method
  
  RNAseq_target <- gsub(pattern = paste(unlist(Cell_terms), collapse = '|'),
                        replacement = 'Cell_culture',
                        RNAseq_method,
                        perl = TRUE)
  
  RNAseq_target <- gsub(pattern = paste(unlist(Tissue_terms), collapse = '|'),
                        replacement = 'Tissue',
                        RNAseq_target,
                        perl = TRUE)
  
  if(str_detect(RNAseq_target, 'Cell_culture')){RNAseq_target <- 'Cell_culture'}
  else{if(str_detect(RNAseq_target, 'Tissue')){RNAseq_target <- 'Tissue'}
    else{
      
      RNAseq_method <- list()
      
      for(i in 1:length(unique(Single_series[, 'PMCID']))){
        
        if(Single_series[i, 'PMCID'] != 'No PMCID'){
          
          RNAseq_method[[i]] <- Isolate_Single_Series_Methods(Single_method =  Single_series[i, 'PMCID'],
                                                              Search_terms = unlist(RNA_terms),
                                                              Paper_list = PDFs,
                                                              GEO_Methods = FALSE, 
                                                              Secondary_search_terms = NULL)
        }
        
        else{RNAseq_method[[i]] <- ''}
        
      }
      
      # Turn RNAseq_method list to a str and replace cell culture and tissue terms
      RNAseq_method <- paste(unlist(RNAseq_method), collapse = ', ')
      
      RNAseq_target <- gsub(pattern = paste(unlist(Cell_terms), collapse = '|'),
                            replacement = 'Cell_culture',
                            RNAseq_method,
                            perl = TRUE)
      
      RNAseq_target <- gsub(pattern = paste(unlist(Tissue_terms), collapse = '|'),
                            replacement = 'Tissue',
                            RNAseq_target,
                            perl = TRUE)
      
      if(str_detect(RNAseq_target, 'Cell_culture')){RNAseq_target <- 'Cell_culture'}
      else{if(str_detect(RNAseq_target, 'Tissue')){RNAseq_target <- 'Tissue'}
        else{RNAseq_target <- 'Manually Assign'}
      }
    }
    }
  
  ###################
  
  # Check for single cell RNAseq details
  RNAseq_method <- Single_method
  
  scRNAseq_check <- gsub(pattern = paste(unlist(Single_cell_rnaseq), collapse = '|'),
                         replacement = 'single cell rnaseq',
                         RNAseq_method,
                         perl = TRUE)
  
  if(str_detect(scRNAseq_check, 'single cell rnaseq')){scRNAseq <- 'scRNA-seq detected'}
  else{
    
    scRNAseq_method <- list()
      
    for(i in 1:length(unique(Single_series[, 'PMCID']))){
      
      if(Single_series[i, 'PMCID'] != 'No PMCID'){
        
        scRNAseq_method[[i]] <- Isolate_Single_Series_Methods(Single_method =  Single_series[i, 'PMCID'],
                                                              Search_terms = unlist(RNA_terms),
                                                              Paper_list = PDFs,
                                                              GEO_Methods = FALSE, 
                                                              Secondary_search_terms = NULL)
      }
      
      else{scRNAseq_method[[i]] <- ''}
      
    }
    
    scRNAseq_method <- paste(unlist(scRNAseq_method), collapse = ', ')
    
    scRNAseq_method <- gsub(pattern = paste(unlist(Single_cell_rnaseq), collapse = '|'),
                            replacement = 'single cell rnaseq',
                            scRNAseq_method,
                            perl = TRUE)
    
    if(str_detect(scRNAseq_method, 'single cell rnaseq')){scRNAseq <- 'scRNA-seq detected'}
    else{scRNAseq <- ''}
    }
  
  ###################
  
  # Search methods information for text around search terms
  # See Isolate_Single_Series_Methods -> Isolate_text_around_keyword
  Media <- Isolate_Single_Series_Methods(Single_method = Single_method,
                                         Search_terms = unlist(Media_terms),
                                         Paper_list = NULL,
                                         GEO_Methods = TRUE, # Isolate_text_around_keyword function used
                                         Secondary_search_terms = NULL)
  
  Single_media <- list()
  
  for(i in 1:length(unique(Single_series[, 'PMCID']))){
    
    if(Single_series[i, 'PMCID'] != 'No PMCID'){
      
      # For each section of the extracted paper, search each section for text around search terms
      Single_media[[i]] <- Isolate_Single_Series_Methods(Single_method = Single_series[i, 'PMCID'],
                                                         Search_terms = unlist(Cell_culture_terms),
                                                         Paper_list = PDFs,
                                                         GEO_Methods = FALSE, 
                                                         Secondary_search_terms = NULL, 
                                                         Check_subsection_body = TRUE)
    }
    
    else{Single_media[[i]] <- ''}
    
  }
  
  Media <- paste(Media, unlist(Single_media), collapse = ', ')
  
  Specific_media <- lapply(names(Media_terms), function(Single_media_type){
    
    # Convert media names shared by same media type to set designation for that type
    gsub(pattern = paste(unlist(Media_terms[[Single_media_type]]), collapse = '|'),
         replacement = Single_media_type,
         Media,
         perl = TRUE)
    
  })
  
  Specific_media <- names(Media_terms)[which(str_detect(Specific_media, names(Media_terms)))]
  
  ###################
  
  Full_search_terms <- list(unlist(ALI_terms), unlist(Submerged_terms))
  
  # Search for ALI or submerged culture details
  ALI <- Isolate_Single_Series_Methods(Single_method = Single_method,
                                       Search_terms = unlist(Full_search_terms),
                                       Paper_list = PDFs,
                                       GEO_Methods = TRUE, 
                                       Secondary_search_terms = NULL)
  
  ALI <- paste(ALI, collapse = ', ')
  
  # Check for ALI or submerged culture details
  ALI_check <- gsub(pattern = paste(unlist(ALI_terms), collapse = '|'),
                    replacement = 'ALI_info_found',
                    ALI,
                    perl = TRUE)
  
  ALI_check <- gsub(pattern = paste(unlist(Submerged_terms), collapse = '|'),
                    replacement = 'Submerged_culture_found',
                    ALI_check,
                    perl = TRUE)
  
  if(str_detect(ALI_check, 'ALI_info_found')){Culture_Type <- 'ALI_Culture'}
  else{if(str_detect(ALI_check, 'Submerged_culture_found')){Culture_Type <- 'Submerged_culture'}
    else{
      
      ALI <- list()
      
      for(i in 1:length(unique(Single_series[, 'PMCID']))){
        
        if(Single_series[i, 'PMCID'] != 'No PMCID'){
          
          # For each section of the extracted paper, search for ALI terms
          ALI[[i]] <- Isolate_Single_Series_Methods(Single_method = Single_series[i, 'PMCID'],
                                                    Search_terms = unlist(Full_search_terms),
                                                    Paper_list = PDFs,
                                                    GEO_Methods = FALSE, 
                                                    Secondary_search_terms = NULL, 
                                                    Check_subsection_body = TRUE)
        }
        
        else{ALI[[i]] <- ''}
        
      }
      
      ALI <- paste(unlist(ALI), collapse = ', ')
      
      ALI_check <- gsub(pattern = paste(unlist(ALI_terms), collapse = '|'),
                        replacement = 'ALI_info_found',
                        ALI,
                        perl = TRUE)
      
      ALI_check <- gsub(pattern = paste(unlist(Submerged_terms), collapse = '|'),
                        replacement = 'Submerged_culture_found',
                        ALI_check,
                        perl = TRUE)
      
      if(str_detect(ALI_check, 'ALI_info_found')){Culture_Type <- 'ALI_Culture'}
      else{if(str_detect(ALI_check, 'Submerged_culture_found')){Culture_Type <- 'Submerged_culture'}
        else{Culture_Type <- 'Manually Assign'}
      }
    }
    }
  
  ###################
  
  # Search for transwell surface details
  Transwell <- Isolate_Single_Series_Methods(Single_method = Single_method,
                                             Search_terms = unlist(Transwell_terms),
                                             Paper_list = PDFs,
                                             GEO_Methods = TRUE, 
                                             Secondary_search_terms = NULL)
  
  Transwell <- paste(ALI, collapse = ', ')
  
  Transwell_data <- list()
  
  for(i in 1:length(unique(Single_series[, 'PMCID']))){
    
    if(Single_series[i, 'PMCID'] != 'No PMCID'){
      
      # For each section of the extracted paper, search for transwell surface details
      Transwell_data[[i]] <- Isolate_Single_Series_Methods(Single_method = Single_series[i, 'PMCID'],
                                                           Search_terms = unlist(Transwell_terms),
                                                           Paper_list = PDFs,
                                                           GEO_Methods = FALSE, 
                                                           Secondary_search_terms = NULL, 
                                                           Check_subsection_body = TRUE)
    }
    
    else{Transwell_data[[i]] <- ''}
    
  }
  
  Transwell <- paste(unlist(Transwell_data), Transwell, collapse = ', ')
  
  # Check for growth surface information
  Growth_surface_details <- lapply(names(Growth_surface), function(Single_surface_type){
    
    gsub(pattern = paste(unlist(Growth_surface[[Single_surface_type]]), collapse = '|'),
         replacement = Single_surface_type,
         Transwell,
         perl = TRUE)
    
  })
  
  Growth_surface_details <- names(Growth_surface)[which(str_detect(Growth_surface_details, names(Growth_surface)))]
  
  ###################
  
  setTxtProgressBar(pb, 
                    match(Single_series[1, 'geo_accession'], 
                          unique(United_Metadata[, 'geo_accession'])))
  
  # Combine all extracted and found metadata into an output df
  Output <- data.frame(Single_series[1, c('geo_accession', 'series_id', 'PMID', 'PMCID', 'Cell_Type')], 
                       RNAseq_target = RNAseq_target, 
                       scRNAseq = scRNAseq,
                       Specific_media = paste(Specific_media, collapse = ', '),
                       Culture_Type = Culture_Type, 
                       Growth_surface = paste(Growth_surface_details, collapse = ', '),
                       GEO_Method = Single_method,
                       Media = Media, 
                       ALI = ALI,
                       Growth_surface_info = Transwell)
})

# Bind output metadata
Isolated_Methods_Data <- do.call(bind_rows, Isolated_Methods_Data)

# Remove duplicate (series + cell type) metadata
Unique_Methods_Data <- Isolated_Methods_Data[which(!duplicated(paste(Isolated_Methods_Data$series_id, Isolated_Methods_Data$Cell_Type))), ]

write.csv(Unique_Methods_Data, file = 'Unique_Methods_Data.csv')

save(Isolated_Methods_Data, Unique_Methods_Data, Control_GEO_Data, file = 'Methods_Data_Checkpoint.RData')

####################

# Series which required manual assignment of any method details were manually assigned
Unique_Methods_Data <- read.csv('Assigned_Unique_Methods_Data.csv')

Unique_Methods_Data <- Unique_Methods_Data[, c("series_id", "PMID", "PMCID", "Cell_Type",
                                               "RNAseq_target", "Specific_media", "Culture_Type", "Growth_surface" )]

Isolated_Methods_Data <- Isolated_Methods_Data[, c("geo_accession", "series_id", 'Cell_Type')]

Isolated_Methods_Data <- merge(Unique_Methods_Data, Isolated_Methods_Data, by = c("series_id", 'Cell_Type'), relationship = 'many-to-many')

Isolated_Methods_Data <- Isolated_Methods_Data[, c("geo_accession", "series_id", "PMID", "PMCID", 
                                                   "Cell_Type", "RNAseq_target", "Specific_media", 
                                                   "Culture_Type", "Growth_surface")]

# This particular experimental series has both tissue and cultured cells derived from those tissues
Isolated_Methods_Data[which(Isolated_Methods_Data$series_id == 'GSE172232' & str_detect(tolower(Isolated_Methods_Data$title), 'naive')), 
                      c("RNAseq_target", "Specific_media", "Culture_Type", "Growth_surface")] <- 'Tissue'

save(Control_GEO_Data, Isolated_Methods_Data, file = 'Media_Checkpoint.RData')
