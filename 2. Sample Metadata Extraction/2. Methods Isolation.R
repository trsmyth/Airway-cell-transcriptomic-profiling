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

Control_GEO_Data <- read.csv("./Input/Sample_Classification_Results.csv")
load("./Output/Paper_IDs.RData")

###################################################################################################

# Split and duplicate rows which have multiple PMIDs
Control_GEO_Data <- Control_GEO_Data %>% separate_rows(PMID, sep = ' ')
Control_GEO_Data <- Control_GEO_Data %>% select(-c(X))

# Remove series with single samples
Control_GEO_Data <- Control_GEO_Data %>% group_by(series_id, Cell_Type) %>% filter(n() > 1)

# Isolate unsorted samples and export one example per series id to assign tissue type
Unsorted_metadata <- Control_GEO_Data[which(Control_GEO_Data$Cell_Type == 'Unsorted Pulmonary Epithelium'), ]
Unique_unsorted_metadata <- Unsorted_metadata[!duplicated(paste(Unsorted_metadata$series_id, Unsorted_metadata$source_name_ch1)), 
                                              c('series_id', 'source_name_ch1', 'Cell_Type')]

write.csv(Unique_unsorted_metadata, file = 'Unsorted_metadata.csv')

# Remove unsorted samples
Control_GEO_Data <- Control_GEO_Data[which(Control_GEO_Data$Cell_Type != 'Unsorted Pulmonary Epithelium'), ]

# Import manually sorted samples and reassign
Reassigned_metadata <- read.csv('Reassigned_metadata.csv')
Reassigned_metadata <- Reassigned_metadata %>% select(-c(X))
Reassigned_metadata[is.na(Reassigned_metadata)] <- ''

# Add sorted labels back to dataset
Unsorted_metadata <- Unsorted_metadata[, colnames(Unsorted_metadata) %ni% 'Cell_Type']
Unsorted_metadata <- merge(Unsorted_metadata, Reassigned_metadata, by = c('series_id', 'source_name_ch1'))
Unsorted_metadata <- Unsorted_metadata[, colnames(Control_GEO_Data)]

# Add samples back to total dataset
Control_GEO_Data <- bind_rows(Control_GEO_Data, Unsorted_metadata)

# Set blank PMIDs to 'No PMID'
Control_GEO_Data[which(Control_GEO_Data$PMID == ''), 'PMID'] <- 'No PMID'

# Add a no PMCID row to the data frame containing PMIDs and corresponding PMCIDs
Paper_IDs <- add_row(Paper_IDs, PMID = 'No PMID', PMCID = 'No PMCID')

# Remove basal and alveolar cells
Control_GEO_Data <- Control_GEO_Data[which(Control_GEO_Data$Cell_Type == 'HBEC' | 
                                             Control_GEO_Data$Cell_Type == 'HNEC' | 
                                             Control_GEO_Data$Cell_Type == 'SAEC'), ]

# Merge PMCID data to metadata data frame
Control_GEO_Data <- merge(Paper_IDs, Control_GEO_Data, by = 'PMID')

###################################################################################################

# Create a list of file names in directory containing each paper
All_folders <- list.files(path = './Output/Papers')

###########

# Create function to uniformly process any search terms
# By default, add a start and end key to each term and its plural form
# If Add_start and/or Add_end are set to TRUE, also create terms with just the start or end tag
process_terms <- function(single_search_term, Add_start = FALSE, Add_end = FALSE){
  
  start_key <- '\r\n\\s*\\d*[.]*\\s*' # \r\n is newline, \\s* is any number of spaces including none, \\d* is any number of digits including none, [.]* is a period or none
  end_key <- '\\s*\\d*[.]*\\s*\r\n'
  
  # Make the term plural so both can be checked
  plural_term <- paste0(single_search_term, 's')
  
  # For the singular and plural version of the search term, process according to optional parameters above
  modified_terms <- lapply(c(single_search_term, plural_term), function(one_term){
    
    term_both <- paste0(start_key, one_term, end_key) # Add start and end key to term
    modified_terms <- list(term_both) # Create storage list for terms and potential start/end key tags
    
    if(Add_start == TRUE){
      
      term_start <- paste0(start_key, one_term)
      modified_terms <- append(modified_terms, term_start)
      
    }
    
    if(Add_end == TRUE){
      
      term_end <- paste0(one_term, end_key)
      modified_terms <- append(modified_terms, term_end)
      
    }
    
    modified_terms
    
  })
  
  modified_terms <- unlist(modified_terms)
  
}

Start_terms_preprocess <- c('material', 'method', 
                            'materials and method', 'materials & method', 
                            'procedure', 'experimental approach', 
                            'experimental section', 'experimental procedure',
                            'experimental model and subject detail', 
                            'experimental model and subject availability')

Start_terms_initial <- process_terms(Start_terms_preprocess)
Start_terms_secondary <- process_terms(Start_terms_preprocess, Add_start = TRUE, Add_end = TRUE)

End_terms_preprocess <- c('result', 'discussion', 'funding', 'acknowledgment', 'acknowledgement', 'reference')
End_terms <- process_terms(End_terms_preprocess)

Extraction_terms <- c('\r\nmethods:', '\r\nresults:', '\r\ndiscussion:')
Modified_terms <- c('star\\+methods')
Removal_terms <- c('star methods', 'key resources table', 'lead contact and materials availability', 'experimental model and subject details', 'figure')

dash_patterns <- c("\\\u002d", # Hyphen-Minus
                   "\\\u2010", # Hyphen
                   "\\\u2011", # Non-Breaking Hyphen 
                   "\\\u2012", # En Dash 
                   "\\\u2013", # Em Dash 
                   "\\\u2014", # Minus Sign 
                   "\\\u2212") # Figure Dash 

###########

# For each folder
PDFs <- lapply(All_folders, function(single_folder){
  
  # print(single_folder)
  
  # Define the current paper directory
  path <- paste0('./Output/Papers/', single_folder)
  
  # List files within the current directory that are pdfs
  files <- list.files(path = path, pattern = "pdf$")
  
  # For each pdf within the current paper directory
  PDFs <- lapply(files, function(single_file){
    
    # Define the pdf file information
    Current_file <- paste0(path, '/', single_file)
    
    #############################################################
    
    # Extract the text from the pdf. This function orients data within columns to maintain sentence structure.
    PDFs <- extract_text(Current_file)
    PDFs <- stri_trans_tolower(PDFs)
    
    Encoding <- Encoding(PDFs)
    PDFs <- iconv(PDFs, "", "ASCII", sub = " ")
    PDFs <- gsub("[^\x20-\x7E]", " ", PDFs)
    PDFs <- gsub(paste(dash_patterns, collapse = '|'), " ", PDFs)
    
    # If a list of terms to be removed from the pdf is present, loop through it to remove the terms
    if(is.null(Extraction_terms) == FALSE){for(i in 1:length(Extraction_terms)){PDFs <- gsub(pattern = Extraction_terms[[i]], '', PDFs)}}
    
    # If a list of terms to be modified from the pdf is present, loop through it to modify the terms
    if(is.null(Modified_terms) == FALSE){for(i in 1:length(Modified_terms)){
      
      Current_term <- Modified_terms[[i]]
      Replacement_term <- gsub('\\s+', ' ', gsub('[[:punct:]]', ' ', Current_term))
      
      PDFs <- gsub(pattern = Current_term,
                   replacement = Replacement_term,
                   PDFs,
                   perl = TRUE)}}
    
    PDFs <- gsub(pattern = paste(Removal_terms, collapse = '|'),
                 replacement = '',
                 PDFs)
    
    #############################################################
    
    # Isolate the pdf table of contents
    toc <- pdf_toc(Current_file)

    # If a TOC is present
    if(length(toc) > 0){
      
      # Isolate elements of TOC (Introduction, Results, etc)
      Individual_sections <- unlist(toc)
      
      # Isolate the toc objects which relate to methods sections
      Methods <- Individual_sections[which(str_detect(tolower(Individual_sections), 'materials|methods|procedure|experimental'))]
      
      # If at least one methods section title found somewhere in the TOC
      if(length(Methods) > 0){
        
        # Determine the number of times the word 'children' appear in the name of the methods section location
        # Find the minimum, as that will correspond to the beginning of the overall methods section
        Methods_length <- min(unlist(lapply(names(Methods), function(x){str_count(x, 'children')})))
        
        # Duplicate the TOC
        Section <- toc
        
        # Isolate the methods subsection of the TOC and its corresponding children by isolating
        # the child
        for(i in 1:Methods_length){
          
          if(i != Methods_length){Section <- Section$children[[1]]}
          else{Section <- Section$children}
          
        }
        
        
        Section
        
        Individual_section_names <- list()
        
        for(i in 1:length(Section)){
          
          Individual_section_names[[i]] <- Section[[i]]['title']
          
        }
        
        Individual_section_names <- unlist(Individual_section_names)
        Individual_section_names <- stri_trans_tolower(Individual_section_names)
        Individual_section_names <- gsub("[^\x20-\x7E]", " ", Individual_section_names)
        
        Individual_sections <- as.list(Individual_section_names)
        names(Individual_sections) <- Individual_section_names
        
        # Isolate the Methods/Procedure section
        # A LOT of pre-processing has to be done here based on the MANY ways methods sections are set up across journals
        Methods_section_index <- which(str_detect(names(Individual_sections), 'materials|methods|procedure|experimental')) # Isolate the methods section
        Methods_section_titles <- stri_trans_tolower(unlist(Section[Methods_section_index])) # Set methods subsections to lower case and unlist
        
        Methods_section_titles <- gsub(paste(dash_patterns, collapse = '|'), " ", Methods_section_titles)
        Methods_section_titles <- gsub('^[^[:alpha:]]*', '', Methods_section_titles) # Remove everything up to first letter
        Methods_section_titles <- gsub("[^\x20-\x7E]", " ", Methods_section_titles) # Remove non-ascii characters
        Methods_section_titles <- gsub("\\.*", "", Methods_section_titles) # Remove periods from titles
        Methods_section_titles <- gsub('(?<= [A-Z]) (?=[a-z])', '. ', Methods_section_titles, perl = TRUE) # Add . between capital letter and lowercase for species names, i.e. M. tuberculosis
        Methods_section_titles <- str_trim(Methods_section_titles) # Trim white space from front/back of title
        
        # Papers which have a 'star★methods' section cause massive problems. In particular, a methods TOC sometimes appears within the text which breaks the start/end search.
        # Remove the below section titles to skip the in text TOC and remove annoying sections such as the key resources table.
        # Also, remove figures. Some figures which are placed in text within a methods section can cause issues.
        Methods_section_titles <- Methods_section_titles[which(!str_detect(Methods_section_titles, paste(Removal_terms, collapse = '|')))]
        
        Search_TOC <- ifelse(length(Methods_section_titles) > 0, 
                             TRUE,
                             FALSE)
      }
      
      else{Search_TOC <- FALSE}
      
    }
    
    else{Search_TOC <- FALSE}
    
    # If a TOC is present
    if((length(toc) > 0 & Search_TOC == TRUE) == TRUE){
      
      # Remove excess white space from PDF
      squished_pdf <- PDFs %>% str_squish()
      
      ##############
      
      # If the first method section title is the methods section header, remove it.
      # Start/end locations will be selected with a minimum character count between start and stop which will break if the section title is included
      if(str_detect(Methods_section_titles[[1]], 'materials|methods|procedure|experimental')){
        Methods_section_titles <- Methods_section_titles[-1]
      }
      
      names(Methods_section_titles) <- Methods_section_titles # Reset section titles to match search terms
      
      #############
      
      # Create empty list to store extracted strings
      String <- list()
      
      # If a methods section is found and there are methods subsection titles
      if(length(Methods_section_titles) > 0){
        
        # For each method subsection title
        for(i in 1:length(Methods_section_titles)){
          
          ################################
          
          Split_title <- unlist(str_split(Methods_section_titles[[i]], "\\s+"))
          
          Start_index <- str_locate_all(squished_pdf, fixed(paste(Split_title, collapse = " ")))
          
          for(j in 1:length(Split_title) - 1){
            
            if(length(as.data.frame(Start_index)[, 2]) == 0){
              
              Start_index <- ifelse(j == 1,
                                    str_locate_all(squished_pdf, fixed(paste(Split_title[-j], collapse = " "))),
                                    str_locate_all(squished_pdf, fixed(paste(Split_title[-1:-j], collapse = " "))))
              
            }
            
            else{break}
            
          }
          
          Start_index <- as.data.frame(Start_index)[, 2] # Turn 2D-int into a df and isolate the second column (end of start term)
          
          ################################
          
          # If the current subsection is not the final subsection
          if(i != length(Methods_section_titles)){
            
            Split_title <- unlist(str_split(Methods_section_titles[[i + 1]], "\\s+"))
            
            End_index <- str_locate_all(squished_pdf, fixed(paste(Split_title, collapse = " ")))
            
            for(k in 1:length(Split_title) - 1){
              
              if(length(as.data.frame(End_index)[, 2]) == 0){
                
                End_index <- ifelse(k == 1,
                                    str_locate_all(squished_pdf, fixed(paste(Split_title[-k], collapse = " "))),
                                    str_locate_all(squished_pdf, fixed(paste(Split_title[-1:-k], collapse = " "))))
                
              }
              
              else{break}
              
            }
          }
          
          # If the current subsection is the final subsection
          else{
            
            # If the methods section is not the final section listed in the TOC:
            # Set the end index as the index of the first element of the next section
            # Else set the end index to the end of the document
            End_index <- ifelse(Methods_section_index < length(Individual_sections),
                                str_locate_all(squished_pdf, fixed(names(Individual_sections[Methods_section_index + 1])[1])),
                                str_locate_all(squished_pdf, "$"))
            
            # If no end index is found or the end index is before the start index, set the end index to the
            # data availability, acknowledgements, author contributions, or references
            if(any(length(as.data.frame(End_index)[, 2]) == 0 | (max(unlist(End_index)) < max(Start_index)))){
              
              End_index <- str_locate_all(squished_pdf, 'data availability|acknowledgment|acknowledgement|author contributions|references|bibliography|work cited|works cited')
              
              # If no end index is found or the end index is before the start index, set end index to the end of the file
              # This occurs when the methods section is at the end of the document or if this is a supplemental file
              if(any(length(as.data.frame(End_index)[, 2]) == 0 | (max(unlist(End_index)) < max(Start_index)))){
                
                End_index <- str_locate_all(squished_pdf, "$")
                
              }
            }
          }
          
          End_index <- as.data.frame(End_index)[, 1] # Turn 2D-int into a df and isolate the first column (start of end term)
          
          # Create an expanded grid with each start with each end index
          Index_matrix <- expand.grid(unlist(Start_index), unlist(End_index)) %>% setNames(c('Start_index', 'End_index'))
          
          # Isolate indexes which are at least 50 characters between start and stop
          # Then find the index with the minimum distance between start and stop
          Index_matrix <- Index_matrix[which(Index_matrix[, 'End_index'] - Index_matrix[, 'Start_index'] > 50), ]
          Differences <- Index_matrix[which.min(Index_matrix[, 'End_index'] - Index_matrix[, 'Start_index']), ]
          
          
          # Isolate text between the beginning and end index
          Subsection_text <- paste('[Title start]', Methods_section_titles[[i]], '[Title end]',
                                   '[Section start]', substr(squished_pdf, Differences['Start_index'], Differences['End_index']), '[Section end]')
          
          String[[i]] <- Subsection_text %>% str_squish()
          
          # If no text data is found, print the current paper folder and section title to manually check 
          if(String[[i]] == paste('[start]', Methods_section_titles[[i]], 'NA [end]')){
            
            print(paste(single_folder, Methods_section_titles[[i]]))
            
          }
        }
        
        String
        
      }
      
      # If there are not methods subsections, isolate the entire methods section
      else{
        
        # Set the start index to the start of the section
        Start_index <- str_locate_all(squished_pdf, fixed(names(Individual_sections[Methods_section_index])[1]))
        Start_index <- as.data.frame(Start_index)[, 1]
        
        # Set the end index to the start of the next section if it is acknowledgements, contributions, or references
        End_index <- str_locate_all(squished_pdf, 'acknowledgment|acknowledgement|author contributions|references|bibliography|work cited|works cited')
        End_index <- as.data.frame(End_index)[, 1]
        
        # Create an expanded grid with each start with each end index
        Index_matrix <- expand.grid(unlist(Start_index), unlist(End_index)) %>% setNames(c('Start_index', 'End_index'))
        
        # Isolate indexes which are at least 50 characters between start and stop
        # Then find the index with the minimum distance between start and stop
        Index_matrix <- Index_matrix[which(Index_matrix[, 'End_index'] - Index_matrix[, 'Start_index'] > 50), ]
        Differences <- Index_matrix[which.min(Index_matrix[, 'End_index'] - Index_matrix[, 'Start_index']), ]
        
        # Isolate the methods section
        String[[i]] <- substr(squished_pdf,  Differences['Start_index'],  Differences['End_index'])
        
      }
    }
    
    else{
      
      # Extract the text from the pdf. This function orients data within columns to maintain sentence structure.
      PDFs <- extract_text(Current_file)
      PDFs <- stri_trans_tolower(PDFs)
      
      Encoding <- Encoding(PDFs)
      PDFs <- iconv(PDFs, "", "ASCII", sub = " ")
      PDFs <- gsub("[^\x20-\x7E]", " ", PDFs)
      
      if(!str_detect(PDFs, 'peer review file|peer review document|reviewer comment|reviewer #')){
      
        String <- Extract_methods_section(path = path,
                                          single_file = single_file,
                                          Start_terms = Start_terms_initial,
                                          End_terms = End_terms,
                                          Extraction_list = Extraction_terms)
        
        if(is.null(String)){
          
          String <- Extract_methods_section(path = path,
                                            single_file = single_file,
                                            Start_terms = Start_terms_secondary,
                                            End_terms = End_terms,
                                            Extraction_list = Extraction_terms)
          
        }
        
        String <- String %>% str_squish()
      }
      
      else{
        print(paste0(single_folder, ' ', single_file, ' appears to be a peer review document'))
        String <- ''
        }
      
    }
  })
  
  # Remove elements without a length
  PDFs <- PDFs[lapply(PDFs, length) > 0]
  
  # Collapse manuscript and any supplement data together if both present
  PDFs <- as.list(unlist(PDFs))
  
  PDFs
  
})

names(PDFs) <- All_folders

###################################################################################################

# Isolate text from papers with PMCIDs
PMCID <- data.frame('PMCID' = names(PDFs[which(str_detect(names(PDFs), 'PMC'))]))

rownames(PMCID) <- seq(1, nrow(PMCID), by = 1)

# Combine text data with control sample metadata
PMCID <- merge(PMCID, Control_GEO_Data, by = 'PMCID')

# Rearrange metadata
PMCID <- PMCID[, c('geo_accession', 'series_id', 'PMID', 'PMCID', 'contact_name', 'contact_email', 'Cell_Type',
                   'title', 'experimental_series_title', 'summary',  'overall_design', 'description', 'source_name_ch1', 
                   'characteristics_ch1', 'treatment_protocol_ch1', 'growth_protocol_ch1', 'extract_protocol_ch1')]

################

# Isolate text from papers without PMCIDs
PMID <- data.frame('PMID' = names(PDFs[which(!str_detect(names(PDFs), 'PMC'))]))

rownames(PMID) <- seq(1, nrow(PMID), by = 1)

# Combine text data with control sample metadata
PMID <- merge(PMID, Control_GEO_Data, by = 'PMID')

# Rearrange metadata
PMID <- PMID[, colnames(PMCID)]

################

# Isolate metadata from samples without PMIDs. Will attempt to find culturing data from GEO entries
No_PMID <- Control_GEO_Data[which(Control_GEO_Data$PMID == 'No PMID'), colnames(PMCID)]

# Add all metadata together
Methods_data <- bind_rows(PMCID, PMID, No_PMID) 

# Add period at end of each cell in indicated columns
Methods_data[, 9:ncol(Methods_data)] <- lapply(Methods_data[, 9:ncol(Methods_data)], 
                                               function(x) paste(x, '.', sep = '')) 

Methods_data <- Methods_data[, c("geo_accession", "series_id", "PMID", "PMCID", "contact_name", 
                                 "contact_email", "Cell_Type", "experimental_series_title", 
                                 "summary", "overall_design", "title", "description", "source_name_ch1", 
                                 "characteristics_ch1", "treatment_protocol_ch1", "growth_protocol_ch1", 
                                 "extract_protocol_ch1")]

# Unite all text metadata together into a single column for text detection
United_Metadata <- Methods_data %>% unite('Method', title:extract_protocol_ch1, remove = FALSE)
United_Metadata <- United_Metadata[, c("geo_accession", "series_id", "PMID", "PMCID", "contact_name", "contact_email", "Cell_Type", "Method")]
United_Metadata$Method <- United_Metadata$Method %>% str_squish()
United_Metadata$Method <- gsub(pattern = '_', replacement = ' ', United_Metadata$Method)

save(Control_GEO_Data, United_Metadata, PDFs, file = 'PDF_Checkpoint.RData')
