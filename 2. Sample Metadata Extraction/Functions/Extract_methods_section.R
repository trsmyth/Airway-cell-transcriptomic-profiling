options(java.parameters = "-Xmx8000m")

library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(tabulapdf)

Extract_methods_section <- function(path, 
                                    single_file, 
                                    Start_terms, 
                                    End_terms, 
                                    Extraction_list = NULL, 
                                    Modify_list = NULL){
  
  # Define the pdf file information
  Current_file <- paste0(path, '/', single_file)
  
  # Extract the text from the pdf. This function orients data within columns to maintain sentence structure.
  PDFs <- extract_text(Current_file)
  PDFs <- stri_trans_tolower(PDFs)
  File_encoding <- Encoding(PDFs)
  
  # If file encoding is not unknown or UTF-8, convert text to UTF-8
  if(File_encoding != 'unknown'){PDFs <- iconv(PDFs, Encoding(PDFs), "UTF-8", sub = " ")}
  
  # If a list of terms to be removed from the pdf is present, loop through it to remove the terms
  if(is.null(Extraction_list) == FALSE){for(i in 1:length(Extraction_list)){PDFs <- gsub(pattern = Extraction_list[[i]], '', PDFs)}}
  if(is.null(Modify_list) == FALSE){for(i in 1:length(Modify_list)){PDFs <- gsub(pattern = Modify_list[[i]], '', PDFs)}}
  
  # Find the bounding indices for any occurrence of the following patterns
  # Each are some form of the materials/procedures/approach which may appear in a text
  Start_index <- str_locate_all(pattern = Start_terms, PDFs)
  
  # Retain results with at least 1 in the first dimension (not empty integer)
  Start_index <- Start_index[sapply(Start_index, function(x) dim(x)[1]) > 0]
  
  # If at least one starting index is found
  if(length(Start_index) > 0){
    
    # Isolate the lowest index (first occurrence) of first index (start of word)
    Start_index <- as.numeric(lapply(Start_index, `[[`, 1))
    
    # Check for an abstract section. Abstracts commonly have methods/results/discussion section
    # which can be detected causing incomplete method extraction
    Abstract_index <- str_locate_all(pattern = c('\r\n\\s*\\d*[.]*\\s*abstract\\s*\\d*[.]*\\s*\r\n'), PDFs)
    
    # Isolate indices with at least one dimension (not empty)
    Abstract_index <- Abstract_index[sapply(Abstract_index, function(x) dim(x)[1]) > 0]
    
    # Isolate the first index of non-empty indices
    Abstract_index <- as.numeric(lapply(Abstract_index, `[[`, 1))
    
    # If at least one abstract index found, set the start index to the first index at least 1000 
    # characters from the abstract start. If there is not an abstract index, start at the minimum start index.
    Start_index <- ifelse(length(Abstract_index) > 0,
                          Start_index[which.min(sapply(Start_index, function(i) i[(i - min(Abstract_index)) > 1000]))],
                          min(Start_index))
    
    # Find the bounding indices for any occurrence of the following patterns
    # Each are some form of the results, discussion, or references which may appear in a text after the methods
    End_index <- str_locate_all(pattern = End_terms, PDFs)
    
    # Retain results with at least 1 in the first dimension (not empty integer)
    End_index <- End_index[sapply(End_index, function(x) dim(x)[1]) > 0]
    End_index <- unlist(End_index) # Unlist results so each are in a single vector for testing
    
    # If at least one end index is found
    if(length(End_index) > 0){
      
      # Find the index of the next closest end word that is at least 100 characters away. 100 characters chosen in case of
      # listed sections. I.e. if a table of contents or such lists methods, results, discussions, references, ignore those
      # end terms and find the index of the results section to end the methods instead of the TOC results.
      End_index <- End_index[which.min(sapply(End_index, function(i) i[(i - as.numeric(Start_index)) > 100]))]
      
      # If end index is greater than 0, the end filtering term is after the method section and the string can be isolated
      if(length(End_index) > 0){
        
        # Isolate text between the beginning and end index
        String <- substr(PDFs, Start_index, End_index)
        
      }
      
      # If the end index is not after the method section, check if the current document is a supplemental methods file
      else{
        
        # If it is and is not a peer review report, isolate to the final character. If not, do not extract from the file
        if(!str_detect(PDFs, paste('peer', 'reviewer', 'comment', sep = '|')) == TRUE){
          
          # Find index of the last word
          End_index <- str_locate(PDFs, "$")
          
          # Isolate text between the beginning and end index
          String <- substr(PDFs, Start_index, End_index[2])
          
        }
      }
    }
    
    else{
      
      # Find index of the last word
      End_index <- str_locate(PDFs, "$")
      
      # Isolate text between the beginning and end index
      String <- substr(PDFs, Start_index, End_index[2])

    } 
  }
}