Isolate_Single_Series_Methods <- function(Single_method,
                                          Search_terms, 
                                          Paper_list = NULL, 
                                          GEO_Methods = FALSE, 
                                          Secondary_search_terms = NULL, 
                                          Check_subsection_body = FALSE){
  
  if(GEO_Methods == TRUE){
    
    Sentence_of_interest <- Isolate_text_around_keyword(paste(Search_terms, collapse = '|'), Single_method)
    Section_Data <- Sentence_of_interest
    
  }
  
  else{
  
    # If there is a paper ID at the current index
    if(Single_method != 'No PMCID'){
      
      # Isolate the list of extracted methods
      Extracted_methods <- Paper_list[[Single_method]]
  
      # Create list to store results
      Section_Data <- list()
      
      #################
      
      # If the length of the extracted methods is greater than one, search each
      # extracted section for the indicated key words and surrounding text
      if(length(Extracted_methods) > 1){
      
        # For each extracted method, whether it is a whole methods section or method subsection
        for(i in 1:length(Extracted_methods)){
          
          # Isolate the current section or subsection
          Current_Section <- tolower(Extracted_methods[i])
          
          # If the current section is NA, change it to a blank
          Current_Section[is.na(Current_Section)] <- ''
          
          # Check if the current section/subsection has a title marker
          # If yes, extract the section title and search for key terms
          if(str_detect(Current_Section, 'title end')){
          
            # Isolate the section title
            Current_Section_Title <- gsub(pattern = '(?<=\\[title end\\]).*', replacement = '\\1', Current_Section, perl = TRUE)
            Current_Section_Title <- gsub(pattern = '\\[title end\\]*', replacement = '', Current_Section_Title)
            Current_Section_Title <- gsub(pattern = '\\[title start\\]', replacement = '', Current_Section_Title)
            Current_Section_Title <- str_trim(Current_Section_Title)
          
            # If the current section contains any of the search terms, isolate the sentence with the term
            if(str_detect(replace_na(Current_Section_Title, ''), paste(Search_terms, collapse = '|'))){
              
              Sentence_of_interest <- Isolate_text_around_keyword(paste(Search_terms, collapse = '|'), Current_Section)
              Section_Data[[i]] <- Sentence_of_interest
              
            }
            
            else{
              if(Check_subsection_body == TRUE){
              Sentence_of_interest <- Isolate_text_around_keyword(paste(Search_terms, collapse = '|'), Current_Section)
              Section_Data[[i]] <- Sentence_of_interest
              }
            }
          }
          
          # If not, search for the keyword and surrounding sentence
          else{
            
            Sentence_of_interest <- Isolate_text_around_keyword(paste(Search_terms, collapse = '|'), Current_Section)
            Section_Data[[i]] <- Sentence_of_interest
            
          }
        }
      }
      
      # If not, search for the keyword and surrounding sentence
      else{
        
        Sentence_of_interest <- Isolate_text_around_keyword(paste(Search_terms, collapse = '|'), tolower(Extracted_methods))
        Section_Data <- Sentence_of_interest
        
      }
      
      Section_Data
      
    }
  }
  
  Section_Data <- as.list(unique(unlist(Section_Data)))
  
  if(all(!is.null(Secondary_search_terms) & length(Section_Data) > 0)){
    
    Section_Data <- Section_Data[which(lapply(Section_Data, function(x) str_detect(string = x, pattern = paste(Secondary_search_terms, collapse = '|'))) == TRUE)]
    
  }
  
  return(Section_Data)
  
}