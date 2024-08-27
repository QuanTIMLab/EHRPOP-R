#' @import dplyr
#' @import tidyr
#' @import networkD3
#' @import stringr
#' @import jsonlite

library(dplyr)
library(tidyr)
library(networkD3)
library(stringr)
library(jsonlite)



# Read the JSON file
json_file_path <- system.file("extdata", "all_codes.json", package = "EHRPOP")
data <- fromJSON(json_file_path)


treatment_data <- data$Treatment

SnakeyDiagram <- function(df) {

    # Function to map code to category
    map_code <- function(code,treatment_data) {
        # All surgery codes (CCAM, ICD10, ATC)
        surgeryCodes <- c(
          treatment_data$Surgery$CCAM,
          treatment_data$Surgery$ICD10,
          treatment_data$Surgery$ATC
        )

        # All Chemotherapy codes (ICD10, CCAM, ATC)
        chemoCodes <- c(
          treatment_data$CT$CCAM,
          treatment_data$CT$ICD10,
          treatment_data$CT$ATC
        )

        # Codes Targeted Therapy (TT)
        TT_codes <- c(
          treatment_data$TT$CCAM,
          treatment_data$TT$ICD10,
          treatment_data$TT$ATC
        )

        # Codes Radio Therapy (RT)
        RtCodes <- c(
          treatment_data$RT$CCAM,
          treatment_data$RT$ICD10,
          treatment_data$RT$ATC
        )

        # ET Codes
        ET_codes <- c(
          treatment_data$ET$CCAM,
          treatment_data$ET$ICD10,
          treatment_data$ET$ATC
        )

        # Other codes
        Other_codes <- c(
          treatment_data$Other$CCAM,
          treatment_data$Other$ICD10,
          treatment_data$Other$ATC
        )




      if (code %in% RtCodes) {
        return('RT')
      } else if (code %in% TT_codes) {
        return('TT')
      } else if (code %in% chemoCodes) {
        return('CT')
      } else if (code %in% surgeryCodes) {
        return('Surgery')
      } else if (code %in% ET_codes) {
        return('ET')
      } else if (code %in% Other_codes) {
        return('OTHER')
      } else {
        return(NA)
      }
    }
    
    
    # Create CODE_ACT column
df <- df %>%
  mutate(CODE_ACT = case_when(
    CODE_ATC != "" ~ CODE_ATC,
    CODE_CCAM != "" ~ CODE_CCAM,
    CODE_ICD10 != "" ~ CODE_ICD10,
    TRUE ~ NA_character_
  ))

    # Drop the original columns CODE_ICD10, CODE_CCAM, and CODE_ATC
df <- df %>%
  select(-CODE_ICD10, -CODE_CCAM, -CODE_ATC,-CODE_BIO,-CODE_LPP,-DATE_SORTIE)

# Drop rows where CODE_ACT is NA
df <- df %>%
  drop_na(CODE_ACT)

    
    # Apply the map_code function to CODE_ACT column
df <- df %>%
  mutate(CODE_ACT = sapply(CODE_ACT, map_code , treatment_data = treatment_data))


# Convert DATE to numeric and replace NA with a large number (999999)
df <- df %>%
  mutate(DATE = as.numeric(as.character(DATE))) %>%
  mutate(DATE = ifelse(is.na(DATE), 999999, DATE))

# Sort the DataFrame by ID_PATIENT and DATE
df <- df %>%
  arrange(ID_PATIENT, DATE)
    
    
    
df <- df %>%
  drop_na(CODE_ACT)
    
        
# Define a function to remove consecutive duplicates
remove_consecutive_duplicates <- function(treatments) {
  # Initialize a vector to store the result
  result <- c()
  
  # Iterate through the treatments and remove consecutive duplicates
  for (i in seq_along(treatments)) {
    if (i == 1 || treatments[i] != treatments[i - 1]) {
      result <- c(result, treatments[i])
    }
  }
  
  # Collapse the result into a single string with "->" separator
  return(str_c(result, collapse = "->"))
}


result <- df %>%
  group_by(ID_PATIENT) %>%
  summarize(Traitements = remove_consecutive_duplicates(CODE_ACT)) %>%
  ungroup()
    
    

label_treatments <- function(treatment_string) {
  treatments <- unlist(strsplit(treatment_string, "->"))
  labeled_treatments <- paste0(treatments, seq_along(treatments))
  return(labeled_treatments)
}
    
    

# Calculate treatment counts and filter sequences
treatment_counts <- table(result$Traitements)
frequent_treatments <- names(treatment_counts[treatment_counts >= 10])
filtered_df <- result %>% filter(Traitements %in% frequent_treatments)

    

# Create a sequence counter
sequence_counter <- list()
for (treatments in filtered_df$Traitements) {
  sequence_list <- label_treatments(treatments)
  for (i in seq_along(sequence_list)[-length(sequence_list)]) {
    pair <- c(sequence_list[i], sequence_list[i + 1])
    key <- paste(pair, collapse = "->")
    if (key %in% names(sequence_counter)) {
      sequence_counter[[key]] <- sequence_counter[[key]] + 1
    } else {
      sequence_counter[[key]] <- 1
    }
  }
}

# Create data for Sankey diagram
pairs <- do.call(rbind, strsplit(names(sequence_counter), "->"))
df_links <- data.frame(
  source = pairs[,1],
  target = pairs[,2],
  value = unlist(sequence_counter)
)

    # Create unique nodes
all_nodes <- unique(c(df_links$source, df_links$target))
node_indices <- setNames(seq_along(all_nodes) - 1, all_nodes)

# Map source and target to their indices
df_links$source <- node_indices[df_links$source]
df_links$target <- node_indices[df_links$target]

# Create Sankey diagram
sankey <- sankeyNetwork(
  Links = df_links,
  Nodes = data.frame(name = all_nodes),
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "Ties",
  fontSize = 12,
  nodeWidth = 30
)
    
    
sankey

    

}
