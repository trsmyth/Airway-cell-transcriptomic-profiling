library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(stringi)

`%ni%` <- negate(`%in%`)

# Isolate text around the chosen keywords
Isolate_text_around_keyword <- function(Single_keyword, Str_to_search){
  
  # To account for situations where non-ASCII characters are present, use stri_trans_tolower()
  Str_to_search <- stri_trans_tolower(Str_to_search)
  
  Str_to_search <- gsub(pattern = '- ', replacement = '', Str_to_search) # Remove hyphenations at end of lines
  Str_to_search <- gsub(pattern = ' inc\\.', replacement = ' inc', Str_to_search) # Replace . in incorporated name
  Str_to_search <- gsub(pattern = ' co\\.', replacement = ' co', Str_to_search) # Replace . in company name
  Str_to_search <- gsub(pattern = ' ltd\\.', replacement = ' ltd', Str_to_search) # Replace . in ltd name
  Str_to_search <- gsub(pattern = '(\\([^)]*)\\.([^)]*\\))', replacement = ' ', Str_to_search, perl = TRUE) # Remove period within parentheses
  Str_to_search <- gsub(pattern = '(?<=\\D)tm ', replacement = ' ', Str_to_search, perl = TRUE) # Remove trademark symbol
  Str_to_search <- gsub(pattern = '(?<=\\d) h', replacement = 'h', Str_to_search, perl = TRUE) # Replace space between hour and number in time
  Str_to_search <- gsub(pattern = '(?<=\\d)h.', replacement = 'h', Str_to_search, perl = TRUE) # Replace . in time
  Str_to_search <- gsub(pattern = "(?<=\\d)\\,(?=\\d)", replacement = "", Str_to_search, perl = TRUE) # Replace comma in number. Is an issue with some product numbers.
  Str_to_search <- gsub(pattern = "(?<=\\d)\\.(?=\\d)", replacement = ",", Str_to_search, perl = TRUE) # Replace decimal place with comma
  Str_to_search <- gsub(pattern = "\\s{2,}", replacement = " ", Str_to_search, perl = TRUE) # Replace multiple spaces with single space
  
  # Create a data frame with the start and end index location for the current key term
  Character_index <- data.frame(str_locate_all(pattern = Single_keyword, Str_to_search))

  # Find all periods in string to determine sentence bounding
  Sentence_index <- data.frame(str_locate_all(pattern = '\\.', Str_to_search))

  # Create a list to store results
  Text <- list()

  # If at least one row is present (one example of key term found)
  if(nrow(Character_index) > 0){

    # For each row corresponding to a start and stop index
    for(i in 1:nrow(Character_index)){

      start <- Sentence_index$start[which(Sentence_index$start < Character_index[i, 'start'])]
      start <- start[length(start)]

      if(length(start) == 0){start = 0}

      end <- Sentence_index$start[which(Sentence_index$start > Character_index[i, 'end'])]
      end <- end[1]

      if(i == 1){

        # Isolate the up and downstream of the key term
        String <- substr(Str_to_search, start + 1, end)
        String <- sub(pattern = '^[^[:alnum:]]+', replacement = '', String) # Remove leading characters up to first letter or number
        Text[[i]] <- String

      }

      else{

        if(Character_index[i, 'start'] > Character_index[i-1, 'end']){

          String <- substr(Str_to_search, start + 1, end) # Isolate the up and downstream of the key term
          String <- sub(pattern = '^[^[:alnum:]]+', replacement = '', String) # Remove leading characters up to first letter or number

        }
      }

      # If current string is not already in the list, add it to the list
      if(String %ni% Text){

        Text[[i]] <- String

      }
    }
  }

  Text

}