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

# Define the path to the package's JSON file
pkg_json_file_path <- system.file("extdata", "all_codes.json", package = "EHRPOP")
treatment_data <- fromJSON(pkg_json_file_path)$Treatment


#' @importFrom jsonlite fromJSON
#' @importFrom dplyr mutate case_when select arrange group_by summarize ungroup
#' @importFrom tidyr drop_na
#' @importFrom stringr str_c
#' @importFrom networkD3 sankeyNetwork
#' @importFrom magrittr %>%
SnakeyDiagram <- function(df) {

    treatment_data <- fromJSON(pkg_json_file_path)$Treatment

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

# Function to add a code to the Surgery section of the JSON data
addCodeSurgery <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(pkg_json_file_path, simplifyVector = FALSE)

  print(data)
  print(treatment_data)

  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$Surgery)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$Surgery[[code_type]])) {
    data$Treatment$Surgery[[code_type]] <- c(data$Treatment$Surgery[[code_type]], code)
    message(sprintf("Code '%s' added to Surgery '%s'.", code, code_type))
    
    # Save the updated data back to the file
    write_json(data, path = pkg_json_file_path, pretty = TRUE, auto_unbox = TRUE)
    assign("treatment_data", data$Treatment, envir = .GlobalEnv)
  } else {
    message(sprintf("Code '%s' already exists in Surgery '%s'.", code, code_type))
  }
}


deleteCodeSurgery <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(pkg_json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the Surgery section
  for (code_type in names(data$Treatment$Surgery)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$Surgery[[code_type]]) {
      # Remove the code from the list
      data$Treatment$Surgery[[code_type]] <- setdiff(data$Treatment$Surgery[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from Surgery '%s'.", code, code_type))
      
      assign("treatment_data", data$Treatment, envir = .GlobalEnv)
      # Save the updated data back to the file
      write_json(data, path = pkg_json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in Surgery section.", code))
  }
}


# Function to add a code to the RT section of the JSON data
addCodeRT <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$RT)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$RT[[code_type]])) {
    data$Treatment$RT[[code_type]] <- c(data$Treatment$RT[[code_type]], code)
    message(sprintf("Code '%s' added to RT '%s'.", code, code_type))
          assign("treatment_data", data$Treatment, envir = .GlobalEnv)

    # Save the updated data back to the file
    write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    message(sprintf("Code '%s' already exists in RT '%s'.", code, code_type))
  }
}


deleteCodeRT <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the RT section
  for (code_type in names(data$Treatment$RT)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$RT[[code_type]]) {
      # Remove the code from the list
      data$Treatment$RT[[code_type]] <- setdiff(data$Treatment$RT[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from RT '%s'.", code, code_type))
            assign("treatment_data", data$Treatment, envir = .GlobalEnv)

      # Save the updated data back to the file
      write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in RT section.", code))
  }
}


# Function to add a code to the CT section of the JSON data
addCodeCT <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$CT)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$CT[[code_type]])) {
    data$Treatment$CT[[code_type]] <- c(data$Treatment$CT[[code_type]], code)
    message(sprintf("Code '%s' added to CT '%s'.", code, code_type))
                assign("treatment_data", data$Treatment, envir = .GlobalEnv)

    # Save the updated data back to the file
    write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    message(sprintf("Code '%s' already exists in CT '%s'.", code, code_type))
  }
}


deleteCodeCT <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the CT section
  for (code_type in names(data$Treatment$CT)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$CT[[code_type]]) {
      # Remove the code from the list
      data$Treatment$CT[[code_type]] <- setdiff(data$Treatment$CT[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from CT '%s'.", code, code_type))
                  assign("treatment_data", data$Treatment, envir = .GlobalEnv)

      # Save the updated data back to the file
      write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in CT section.", code))
  }
}



# Function to add a code to the ET section of the JSON data
addCodeET <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$ET)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$ET[[code_type]])) {
    data$Treatment$ET[[code_type]] <- c(data$Treatment$ET[[code_type]], code)
    message(sprintf("Code '%s' added to ET '%s'.", code, code_type))
                      assign("treatment_data", data$Treatment, envir = .GlobalEnv)

    # Save the updated data back to the file
    write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    message(sprintf("Code '%s' already exists in ET '%s'.", code, code_type))
  }
}


deleteCodeET <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the ET section
  for (code_type in names(data$Treatment$ET)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$ET[[code_type]]) {
      # Remove the code from the list
      data$Treatment$ET[[code_type]] <- setdiff(data$Treatment$ET[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from ET '%s'.", code, code_type))
      
      assign("treatment_data", data$Treatment, envir = .GlobalEnv)

      # Save the updated data back to the file
      write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in ET section.", code))
  }
}


# Function to add a code to the TT section of the JSON data
addCodeTT <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$TT)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$TT[[code_type]])) {
    data$Treatment$TT[[code_type]] <- c(data$Treatment$TT[[code_type]], code)
    message(sprintf("Code '%s' added to TT '%s'.", code, code_type))
          assign("treatment_data", data$Treatment, envir = .GlobalEnv)

    # Save the updated data back to the file
    write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    message(sprintf("Code '%s' already exists in TT '%s'.", code, code_type))
  }
}


deleteCodeTT <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the TT section
  for (code_type in names(data$Treatment$TT)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$TT[[code_type]]) {
      # Remove the code from the list
      data$Treatment$TT[[code_type]] <- setdiff(data$Treatment$TT[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from TT '%s'.", code, code_type))
            assign("treatment_data", data$Treatment, envir = .GlobalEnv)

      # Save the updated data back to the file
      write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in TT section.", code))
  }
}


# Function to add a code to the Other section of the JSON data
addCodeOther <- function(code, code_type) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  # Check if the code_type is valid
  if (!code_type %in% names(data$Treatment$Other)) {
    message(sprintf("Error: Invalid code type '%s'", code_type))
    return(NULL)
  }
  
  # Add the code if it does not already exist
  if (!(code %in% data$Treatment$Other[[code_type]])) {
    data$Treatment$Other[[code_type]] <- c(data$Treatment$Other[[code_type]], code)
    message(sprintf("Code '%s' added to Other '%s'.", code, code_type))
                assign("treatment_data", data$Treatment, envir = .GlobalEnv)

    # Save the updated data back to the file
    write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    message(sprintf("Code '%s' already exists in Other '%s'.", code, code_type))
  }
}


deleteCodeOther <- function(code) {
  # Read the JSON data from the file
  data <- fromJSON(json_file_path, simplifyVector = FALSE)
  
  found <- FALSE
  
  # Iterate over each code type in the Other section
  for (code_type in names(data$Treatment$Other)) {
    # Check if the code is in the list for the current code type
    if (code %in% data$Treatment$Other[[code_type]]) {
      # Remove the code from the list
      data$Treatment$Other[[code_type]] <- setdiff(data$Treatment$Other[[code_type]], code)
      found <- TRUE
      message(sprintf("Code '%s' removed from Other '%s'.", code, code_type))
                  assign("treatment_data", data$Treatment, envir = .GlobalEnv)

      # Save the updated data back to the file
      write_json(data, path = json_file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  }
  
  if (!found) {
    message(sprintf("Code '%s' not found in Other section.", code))
  }
}


yesOrNo <- function(df, CCAM_codes, ATC_codes, ICD_Codes, columnName,
                    indexDate = "DATE", indexCodeCCAM = "CODE_CCAM", indexCodeATC = "CODE_ATC",
                    indexCodeICD = "CODE_ICD10", indexID = "ID_PATIENT") {
  
  # Create a logical vector for matches within the specified timeframe and code lists
  df <- df %>%
    mutate(MATCH = ((df[[indexCodeCCAM]] %in% CCAM_codes |
                     df[[indexCodeATC]] %in% ATC_codes |
                     df[[indexCodeICD]] %in% ICD_Codes) 
                    ))
  
  # Group by patient ID and check if any matches exist for each patient
  result <- df %>%
    group_by_at(indexID) %>%
    summarize(!!columnName := any(MATCH)) %>%
    ungroup()
  
  return(result)
}


isTreatedByIt <- function(df, CCAM_codes, ATC_codes, ICD_Codes, columnName,
                          indexCodeCCAM = "CODE_CCAM", indexCodeATC = "CODE_ATC",
                          indexCodeICD = "CODE_ICD10", indexID = "ID_PATIENT") {
  # Create a logical vector for matches within the specified code lists
  df <- df %>%
    mutate(MATCH = (df[[indexCodeCCAM]] %in% CCAM_codes |
                    df[[indexCodeATC]] %in% ATC_codes |
                    df[[indexCodeICD]] %in% ICD_Codes))
  
  # Group by patient ID and check if any matches exist for each patient
  result <- df %>%
    group_by_at(indexID) %>%
    summarize(!!columnName := any(MATCH)) %>%
    ungroup()
  
  return(result)
}


isTreatedByItWithQte <- function(df, CCAM_codes, ATC_codes, ICD_Codes, columnName) {
  # Helper function to check if code starts with any prefix in codes_list
  starts_with_any <- function(code, codes_list) {
    code <- as.character(code)
    any(sapply(codes_list, function(prefix) startsWith(code, prefix)))
  }
  
  # Fill NA values
  df <- df %>%
    mutate(
      DATE = ifelse(is.na(DATE), DATE_ENTREE, DATE),
      QUANTITE = ifelse(is.na(QUANTITE), 1, QUANTITE),
      SESSION = 1,
      CODE_CCAM = as.character(CODE_CCAM),
      CODE_ATC = as.character(CODE_ATC),
      CODE_ICD10 = as.character(CODE_ICD10)
    )
  
  # Check for relevant codes
  df <- df %>%
    mutate(
      Is_Relevant = mapply(starts_with_any, CODE_CCAM, MoreArgs = list(codes_list = CCAM_codes)) |
                    mapply(starts_with_any, CODE_ATC, MoreArgs = list(codes_list = ATC_codes)) |
                    mapply(starts_with_any, CODE_ICD10, MoreArgs = list(codes_list = ICD_Codes))
    )
  
  # Summarize the results
  result <- df %>%
    filter(Is_Relevant) %>%
    group_by(ID_PATIENT) %>%
    summarise(!!columnName := sum(SESSION, na.rm = TRUE)) %>%
    ungroup()
  
  # Include patients with no relevant treatments
  all_patients <- df %>%
    select(ID_PATIENT) %>%
    distinct()
  
  result <- all_patients %>%
    left_join(result, by = "ID_PATIENT") %>%
    mutate(!!columnName := ifelse(is.na(!!sym(columnName)), 0, !!sym(columnName)))
  
  return(result)
}


# Define tableValues function in R
tableValues <- function(df, listColumns, dateStart = NULL, dateEnd = NULL) {
  if (length(listColumns) == 0) {
    stop("listColumns cannot be empty.")
  }

  missing_columns <- setdiff(listColumns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(paste("The following columns are missing in the DataFrame:", paste(missing_columns, collapse = ", ")))
  }

  if (!is.null(dateStart)) {
    base_date <- as.Date("1960-01-01")
    dateStart <- as.integer(as.Date(dateStart) - base_date)
  }

  if (!is.null(dateEnd)) {
    base_date <- as.Date("1960-01-01")
    dateEnd <- as.integer(as.Date(dateEnd) - base_date)
  }

  results <- NULL

  for (coln in listColumns) {
    unique_values <- unique(na.omit(df[[coln]]))
    count <- 0
    for (value in unique_values) {
      if (startsWith(coln, "CODE_")) {
        r <- isTreatedByItWithQte(df[df[[coln]] == value, ], list(value), list(value), list(value), paste0(coln, "_", value))
      } else {
        r <- isTreatedByItWithQte(df[df[[coln]] == value, ], list(), list(), list(), paste0(coln, "_", value))
      }
      
      if (count == 0) {
        results <- r
      } else {
        results <- merge(results, r, by = "ID_PATIENT", all = TRUE)
      }
      count <- count + 1
    }
  }

  return(results)
}

# Define the function to combine codes from specified columns
combine_codes <- function(row, listColumns) {
  codes <- ""
  for (col in listColumns) {
    if (!is.na(row[[col]])) {
      codes <- row[[col]]
    }
  }
  return(codes)
}

# Define the tableSequances function
tableSequances <- function(df, listColumns) {
  if (length(listColumns) == 0) {
    stop("listColumns cannot be empty.")
  }

  missing_columns <- setdiff(listColumns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(paste("The following columns are missing in the DataFrame:", paste(missing_columns, collapse = ", ")))
  }

  # Apply the function to create CODE_ACTS
  df <- df %>%
    rowwise() %>%
    mutate(CODE_ACTS = combine_codes(cur_data(), listColumns)) %>%
    ungroup()

  df <- df %>%
    filter(CODE_ACTS != "") %>%
    select(ID_PATIENT, DATE, CODE_ACTS)

  df <- df %>%
    arrange(ID_PATIENT, DATE)

  grouped <- df %>%
    group_by(ID_PATIENT) %>%
    summarise(DATES = list(DATE), CODE_ACTS = list(CODE_ACTS)) %>%
    ungroup()

  # Sort acts by date within each group
  grouped <- grouped %>%
    rowwise() %>%
    mutate(ACTES = list(CODE_ACTS[order(DATES)])) %>%
    ungroup()

  # Create the final DataFrame with ID_PATIENT, DATES, and ACTES columns
  final_df <- grouped %>%
    select(ID_PATIENT, DATES, ACTES)

  return(final_df)
}

# Define the tableSequancesTwo function
tableSequancesTwo <- function(df, listColumns) {
  if (length(listColumns) == 0) {
    stop("listColumns cannot be empty.")
  }

  missing_columns <- setdiff(listColumns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(paste("The following columns are missing in the DataFrame:", paste(missing_columns, collapse = ", ")))
  }

  # Apply the function to create CODE_ACTS
  df <- df %>%
    rowwise() %>%
    mutate(CODE_ACTS = combine_codes(cur_data(), listColumns)) %>%
    ungroup()

  df <- df %>%
    filter(CODE_ACTS != "") %>%
    select(ID_PATIENT, DATE, CODE_ACTS)

  # Sort the DataFrame by ID_PATIENT and DATE
  df <- df %>%
    arrange(ID_PATIENT, DATE)

  # Function to create the sequence for each patient
  create_sequence <- function(group) {
    sequence <- c()
    prev_date <- NA
    for (i in 1:nrow(group)) {
      current_date <- group$DATE[i]
      current_code <- group$CODE_ACTS[i]
      if (is.na(prev_date)) {
        sequence <- c(sequence, paste0("[0,", current_code, "]"))
      } else {
        interval <- current_date - prev_date
        sequence <- c(sequence, paste0("[", interval, ",", current_code, "]"))
      }
      prev_date <- current_date
    }
    return(sequence)
  }

  # Group by ID_PATIENT and apply the function
  result <- df %>%
    group_by(ID_PATIENT) %>%
    summarise(Sequence = list(create_sequence(cur_data()))) %>%
    ungroup()

  return(result)
}

library(jsonlite)

readJSON <- function(input_filepath, json_path) {
  # Utility function to load JSON data from a given file path
  load_json_data <- function(filepath) {
    data <- fromJSON(filepath)
    return(data)
  }
  
  # Load the input JSON data
  input_data <- load_json_data(input_filepath)
  
  # Check if the JSON structure starts with "Treatment"
  if (!"Treatment" %in% names(input_data)) {
    stop("The JSON structure must start with 'Treatment'.")
  }
  
  # Load the existing data from the JSON file where the output will be stored
  if (file.exists(json_path)) {
    data <- load_json_data(json_path)
  } else {
    data <- list()
  }
  
  # Ensure the "Treatment" key exists in the existing data
  if (!"Treatment" %in% names(data)) {
    data$Treatment <- list()
  }
  
  # Concatenate the input data with the existing output data
  for (treatment_type in names(input_data$Treatment)) {
    if (!treatment_type %in% names(data$Treatment)) {
      data$Treatment[[treatment_type]] <- input_data$Treatment[[treatment_type]]
    } else {
      for (category in names(input_data$Treatment[[treatment_type]])) {
        if (!category %in% names(data$Treatment[[treatment_type]])) {
          data$Treatment[[treatment_type]][[category]] <- input_data$Treatment[[treatment_type]][[category]]
        } else {
          # Combine the existing and new codes, removing duplicates
          combined_codes <- unique(c(
            data$Treatment[[treatment_type]][[category]],
            input_data$Treatment[[treatment_type]][[category]]
          ))
          data$Treatment[[treatment_type]][[category]] <- combined_codes
        }
      }
    }
  }

  assign("treatment_data", data$Treatment, envir = .GlobalEnv)
  # Save the updated data back to the output JSON file
  write_json(data, json_path, pretty = TRUE)
}


library(dplyr)

neoadjuvantOrAdjuvantOrBoth <- function(df, CCAM_codes, ATC_codes, ICD_Codes, columnName, daysBefore, daysAfter,
                                        indexDate = "DATE", indexCodeCCAM = "CODE_CCAM", 
                                        indexCodeATC = "CODE_ATC", indexCodeICD = "CODE_ICD10",
                                        indexID = "ID_PATIENT", BC_index_surgery = "BC_index_surgery") {
    
    
    # Calculate DATE_DIFF as the difference between DATE and BC_index_surgery
    df$DATE_DIFF <- df_merged$DATE - df$BC_index_surgery
    df <- df[df$DATE_DIFF >= -daysBefore & df$DATE_DIFF <= daysAfter, ]
    
    df$Is_Relevant <- df[[indexCodeCCAM]] %in% CCAM_codes | 
                  df[[indexCodeATC]] %in% ATC_codes | 
                  df[[indexCodeICD]] %in% ICD_Codes
    
    
    # Initialize an empty list to store the patient classifications
    patient_classification <- list()

    # Get unique patient IDs
    unique_patients <- unique(df[[indexID]])

    # Loop through each patient ID
    for (patient_id in unique_patients) {
      # Filter the patient's relevant treatments
      patient_df <- df[df[[indexID]] == patient_id & df$Is_Relevant, ]

      # Skip patients with no relevant treatments
      if (nrow(patient_df) == 0) {
        patient_classification[[as.character(patient_id)]] <- 'False'
        next
      }

      # Count positive and negative DATE_DIFF values
      positive_count <- sum(patient_df[['DATE_DIFF']] > 0)
      negative_count <- sum(patient_df[['DATE_DIFF']] < 0)

      # Classify based on the counts
      if (positive_count > 0 & negative_count > 0) {
        patient_classification[[as.character(patient_id)]] <- 'Both'
      } else if (positive_count > 0) {
        patient_classification[[as.character(patient_id)]] <- 'Adjuvant'
      } else if (negative_count > 0) {
        patient_classification[[as.character(patient_id)]] <- 'Neoadjuvant'
      } else {
        patient_classification[[as.character(patient_id)]] <- 'False'
      }
    }

    # Convert the list to a data frame
    result <- setNames(
      data.frame(
        ID_PATIENT = names(patient_classification),
        unlist(patient_classification),
        stringsAsFactors = FALSE
      ),
      c("ID_PATIENT", columnName)
    )



    
  return(result)
}

chemotherapyIntervals <- function(df, CCAM_codes, ATC_codes, ICD_Codes, columnName, daysBefore, daysAfter,
                                        indexDate = "DATE", indexCodeCCAM = "CODE_CCAM", 
                                        indexCodeATC = "CODE_ATC", indexCodeICD = "CODE_ICD10",
                                        indexID = "ID_PATIENT", BC_index_surgery = "BC_index_surgery") {
    
    
    # Calculate DATE_DIFF as the difference between DATE and BC_index_surgery
    df$DATE_DIFF <- df_merged$DATE - df$BC_index_surgery
    df <- df[df$DATE_DIFF >= -daysBefore & df$DATE_DIFF <= daysAfter, ]
    
    df$Is_Relevant <- df[[indexCodeCCAM]] %in% CCAM_codes | 
                  df[[indexCodeATC]] %in% ATC_codes | 
                  df[[indexCodeICD]] %in% ICD_Codes
    
    
    # Initialize an empty list to store the patient classifications
    patient_classification <- list()

    # Get unique patient IDs
    unique_patients <- unique(df[[indexID]])
    
    # Initialize an empty list to store the patient classifications
    patient_classification <- list()

    # Get unique patient IDs
    unique_patients <- unique(df[[indexID]])

    # Loop through each patient ID
    for (patient_id in unique_patients) {
      # Filter the patient's relevant treatments
      patient_df <- df[df[[indexID]] == patient_id & df$Is_Relevant, ]

      # Skip patients with no relevant treatments
      if (nrow(patient_df) == 0) {
        patient_classification[[as.character(patient_id)]] <- 'False'
        next
      }

      Text <- ""
      temp <- 0
      c <- 0

      # Loop through each row in the filtered data
      for (i in 1:nrow(patient_df)) {
        data <- patient_df[i, ]

        if (c == 0) {
          Text <- "0"
          temp <- data[[indexDate]]
        } else {
          Text <- paste0(Text, " -> ", as.character(data[[indexDate]] - temp))
          temp <- data[[indexDate]]
        }

        c <- c + 1
      }

      patient_classification[[as.character(patient_id)]] <- Text
    }

    # Convert the list to a data frame with a dynamic column name
    result <- setNames(
      data.frame(
        ID_PATIENT = names(patient_classification),
        unlist(patient_classification),
        stringsAsFactors = FALSE
      ),
      c("ID_PATIENT", columnName)
    )

    # return the result
    return(result)

    
    }


classify_regimen_chemo <- function(text) {
  # Load necessary libraries
  library(stringr)
  
  # Never done a chemotherapy
  if (text == "False") {
    return("False")
  }
  
  # Extract numeric values from the text
  numbers <- as.numeric(str_extract_all(text, "\\b\\d+\\.?\\d*\\b")[[1]])
  numbers <- numbers[numbers != 0]
  
  # They have done only one chemotherapy session
  if (length(numbers) == 0) {
    return("ONE TREATMENT")
  }
  
  # Normalize specific values to standard intervals
  numbers[numbers %in% c(20, 22)] <- 21
  numbers[numbers %in% c(6, 8)] <- 7
  numbers[numbers %in% c(13, 15)] <- 14
  
  # Define helper functions to classify the regimens
  is_paclitaxel <- function(seq) {
    all(seq == 7)
  }
  
  is_anthracyclines_docetaxel <- function(seq) {
    len <- length(seq)
    len >= 4 && all(seq[1:3] == 14) && all(seq[4:len] == 21)
  }
  
  is_anthracyclines_paclitaxel <- function(seq) {
    len <- length(seq)
    len >= 4 && all(seq[1:3] == 14) && all(seq[4:len] == 7)
  }
  
  is_anthracyclines_paclitaxel2 <- function(seq) {
    transition_index <- which(seq == 7)[1]
    if (is.na(transition_index)) {
      return(FALSE)  # No 7-day interval found
    }
    before_transition <- all(seq[1:(transition_index - 1)] == 21)
    after_transition <- all(seq[transition_index:length(seq)] == 7)
    before_transition && after_transition
  }
  
  is_anthracyclines <- function(seq) {
    all(seq == 21)  # Unknown after March 2012
  }
  
  # Apply classification rules
  if (is_paclitaxel(numbers)) return('Paclitaxel')
  if (is_anthracyclines_docetaxel(numbers)) return('Anthracyclines/docetaxel')
  if (is_anthracyclines_paclitaxel(numbers)) return('Anthracyclines/paclitaxel')
  if (is_anthracyclines_paclitaxel2(numbers)) return('Anthracyclines/paclitaxel')
  if (is_anthracyclines(numbers)) return('Anthracyclines')
  
  return('Other')
}


isDementia <- function(df, columnName) {
    # Define the dementia criteria
    determine_dementia <- function(row) {
        if (row['IS_Dementia']) {
            return(TRUE)
        } else if (is.na(row['Qte_Dementia']) || row['Qte_Dementia'] <= 2) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
    
    # Apply the treatments and conditions
    DementiaATC <- isTreatedByItWithQte(df, CCAM_codes = NULL, ATC_codes = data$Disease$Dementia$ATC,ICD_Codes=NULL,columnName="Qte Dementia")
    DementiaICD10 <- isTreatedByIt(df, CCAM_codes = NULL, ICD_Codes = data$Disease$Dementia$ICD10,ATC_codes=NULL,columnName="IS_Dementia")
    
    # Merge results
    Dementia <- merge(DementiaATC, DementiaICD10, by = 'ID_PATIENT')
    
    # Apply the dementia determination function
    Dementia[[columnName]] <- apply(Dementia, 1, determine_dementia)
    
    # Drop unnecessary columns
    Dementia <- Dementia[, !names(Dementia) %in% c('Qte_Dementia', 'IS_Dementia')]
    
    return(Dementia)
}


isCOPD <- function(df, columnName) {
    
COPD <- isTreatedByIt(df, CCAM_codes = NULL, ICD_Codes = data$Disease$Dementia$ICD10,ATC_codes=data$Disease$Dementia$ATC,columnName=columnName)

return (COPD)

}

isHypertension <- function(df, columnName) {
    # Identifies patients with hypertension based on antihypertensive drug dispensation criteria.
    # Criteria:
    # - Antihypertensive drugs dispensed at least 3 times during the previous 12 months.

    # Determine hypertension based on the quantity of hypertension treatments
    determine_hypertension <- function(row) {
        if (is.na(row$Qte_Hypertension)) {
            return(FALSE)
        } else if (row$Qte_Hypertension > 2) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    # Assuming isTreatedByItWithQte is defined somewhere else
    Hypertension <- isTreatedByItWithQte(df, CCAM_codes = NULL, ATC_codes = data$Disease$Hypertension, ICD_Codes = NULL, columnName = 'Qte Hypertension')

    # Apply the determine_hypertension function to each row
    Hypertension[[columnName]] <- apply(Hypertension, 1, determine_hypertension)

    # Drop the Qte Hypertension column
    Hypertension <- Hypertension[, !names(Hypertension) %in% c('Qte_Hypertension')]

    return(Hypertension)
}

isDiabetes <- function(df, columnName) {
    
Diabetes <- isTreatedByIt(df, CCAM_codes = NULL, ICD_Codes = data$Disease$Diabetes$ICD10,ATC_codes=data$Disease$Diabetes$ATC,columnName=columnName)

return (Diabetes)

}

isCerebrovascular <- function(df, columnName) {
    
Cerebrovascular <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Cerebrovascular Disease`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Cerebrovascular)

}

isHeart_failure <- function(df, columnName) {
    
Heart_failure <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Heart Failure`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Heart_failure)

}


isMyocardial_infarction <- function(df, columnName) {
    
Myocardial_infarction <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Myocardial_infarction`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Myocardial_infarction)

}

isChronic_ischaemic <- function(df, columnName) {
    # Function to identify patients with chronic ischemic heart disease
  
    # Define ICD-10 codes to exclude
    without_chronic <- c('I21', 'I24')
    
    # Define ICD-10 codes for chronic ischemic heart disease
    chronic_ischaemic_codes <- data$Disease$`Chronic Ischaemic Heart Disease`$ICD10
    
    # Check for chronic ischemic heart disease in the data
    chronic_ischaemic <- df %>%
        dplyr::filter(CODE_ICD10 %in% chronic_ischaemic_codes) %>%
        dplyr::group_by(ID_PATIENT) %>%
        dplyr::summarize(Chronic_Ischaemic = TRUE)
    
    # Check for codes that exclude chronic ischemic heart disease
    without_chronic_df <- df %>%
        dplyr::filter(CODE_ICD10 %in% without_chronic) %>%
        dplyr::group_by(ID_PATIENT) %>%
        dplyr::summarize(Without_chronic = TRUE)
    
    # Merge both datasets
    merged_df <- dplyr::left_join(chronic_ischaemic, without_chronic_df, by = "ID_PATIENT")
    
    # Determine chronic ischemic heart disease, excluding specific codes
    merged_df[[columnName]] <- ifelse(is.na(merged_df$Without_chronic) & !is.na(merged_df$Chronic_Ischaemic), TRUE, FALSE)
    
    # Drop temporary columns
    merged_df <- merged_df %>%
        dplyr::select(ID_PATIENT, !!columnName)
    
    return(merged_df)
}


isStroke <- function(df, columnName) {
    
    
Acute_stroke <- isTreatedByIt(df, CCAM_codes = NULL, 
                                     ICD_Codes = data$Disease$`Acute Stroke`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Acute_stroke)

}

isRenal_disease <- function(df, columnName) {
    
Renal_disease <- isTreatedByIt(df, CCAM_codes = data$Disease$`Renal Disease`$CCAM, 
                                 ICD_Codes = data$Disease$`Renal Disease`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Renal_disease)

}




isLiver_and_Pancreas <- function(df, columnName) {
    
Liver_and_Pancreas <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Liver and Pancreas`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Liver_and_Pancreas)

}


isUndernutrition <- function(df, columnName) {
    
Undernutrition <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Undernutrition`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Undernutrition)

}


isParkinson <- function(df, columnName) {
    
Parkinson <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Parkinson`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Parkinson)

}

isEpilepsy <- function(df, columnName) {
    
Epilepsy <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Epilepsy`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Epilepsy)

}
isPsychiatric_Disease <- function(df, columnName) {
    
    
Psychiatric_Disease_ICD10 <- c(
  data$Disease$`Psychiatric Disease`$`Schizophrenia and Delusional Diseases`,
  data$Disease$`Psychiatric Disease`$`Depression and Mood Diseases`,
  data$Disease$`Psychiatric Disease`$`Mental Deficiency`,
  data$Disease$`Psychiatric Disease`$`Substance Abuse Disorders`,
  data$Disease$`Psychiatric Disease`$`Disorders of Psychological Development`
)
    
Psychiatric_Disease <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = Psychiatric_Disease_ICD10,ATC_codes=NULL,columnName=columnName)

return (Psychiatric_Disease)

}

isPeripheral_vascular <- function(df, columnName) {
    
    
    
Peripheral_vascular <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Peripheral Vascular Disease`$ICD10,ATC_codes=NULL,columnName=columnName)

return (Peripheral_vascular)

}

isDyslipidemia <- function(df, columnName) {
  #' Identifies patients with dyslipidemia based on specific treatment criteria and absence of associated pathologies.
  #'
  #' Criteria:
  #' - Delivery on at least 3 occasions in the year:
  #'   - Statins: C10AA, C10BA, C10BX
  #'   - Fibrates: C10AB
  #'   - Other lipid-lowering agents: C10AC, C10AD, C10AX
  #'
  #' Without associated pathology:
  #' - No code for coronary heart disease (I20, I21, I22, I23, I24, I25)
  #' - No code for stroke (I60, I61, I62, I63, I64)
  #' - No code for heart failure (I50, I11.0, I13.0, I13.2, I13.9, K76.1, J81)
  #' - No code for atherosclerosis of arteries of extremities (I70.2)
  #' - No code for chronic endstage kidney failure (N184, N185)
  #' - No code for diabetes mellitus ("E10", "E11", "E12", "E13", "E14", "A10A", "A10B") or complications
  #'
  #' @param df The input data frame containing patient treatment records.
  #' @param columnName The name for the output column in the resulting data frame, indicating whether each patient has dyslipidemia.
  #'
  #' @return A data frame with patient IDs and a boolean indicator of whether they have dyslipidemia.
  
  determine_dyslipidemia <- function(row) {
    # Handle NA values by treating them as FALSE in logical operations
    chronic_ischaemic <- ifelse(is.na(row["Chronic_ischaemic"]), FALSE, row["Chronic_ischaemic"])
    acute_stroke <- ifelse(is.na(row["Acute_stroke"]), FALSE, row["Acute_stroke"])
    heart_failure <- ifelse(is.na(row["Heart_failure"]), FALSE, row["Heart_failure"])
    aaoaoe <- ifelse(is.na(row["AOAOE"]), FALSE, row["AOAOE"])
    ekf <- ifelse(is.na(row["EKF"]), FALSE, row["EKF"])

    if (chronic_ischaemic | acute_stroke | heart_failure | aaoaoe | ekf) {
      return(FALSE)
    } else if (is.na(row["Qte_Dyslipidemia"]) | row["Qte_Dyslipidemia"] <= 2) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  # Define ATC codes for Dyslipidemia
  Dyslipidemia_ATC <- c(data$Disease$Dyslipidemia$Statins, data$Disease$Dyslipidemia$Fibrates, data$Disease$Dyslipidemia$`Other Lipid Lowering Agents`)

  # Identify patients treated for Dyslipidemia and count treatments
  Dyslipidemia <- isTreatedByItWithQte(df, character(0), Dyslipidemia_ATC, character(0), 'Qte_Dyslipidemia')

  # Identify patients with associated pathologies
  atherosclerosis_of_arteries_of_extremities <- isTreatedByIt(df, character(0), c('I702'), character(0), 'AOAOE')
  endstage_kidney_failure <- isTreatedByIt(df, character(0), c('N184', 'N185'), character(0), 'EKF')
  Chronic_ischaemic <- isChronic_ischaemic(df, 'Chronic_ischaemic')
  Acute_stroke <- isStroke(df, 'Acute_stroke')
  Heart_failure <- isHeart_failure(df, 'Heart_failure')

  # Merge all pathology data into the Dyslipidemia data frame
  Dyslipidemia <- merge(Dyslipidemia, Chronic_ischaemic[, c('ID_PATIENT', 'Chronic_ischaemic')], by = 'ID_PATIENT', all.x = TRUE)
  Dyslipidemia <- merge(Dyslipidemia, Acute_stroke[, c('ID_PATIENT', 'Acute_stroke')], by = 'ID_PATIENT', all.x = TRUE)
  Dyslipidemia <- merge(Dyslipidemia, Heart_failure[, c('ID_PATIENT', 'Heart_failure')], by = 'ID_PATIENT', all.x = TRUE)
  Dyslipidemia <- merge(Dyslipidemia, atherosclerosis_of_arteries_of_extremities[, c('ID_PATIENT', 'AOAOE')], by = 'ID_PATIENT', all.x = TRUE)
  Dyslipidemia <- merge(Dyslipidemia, endstage_kidney_failure[, c('ID_PATIENT', 'EKF')], by = 'ID_PATIENT', all.x = TRUE)

  # Apply the classification
  Dyslipidemia[[columnName]] <- apply(Dyslipidemia, 1, determine_dyslipidemia)

  return(Dyslipidemia[, c('ID_PATIENT', columnName)])
}

isTobacco <- function(df, columnName) {
    
    
    
Tobacco <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Tobacco`$ICD10,ATC_codes=data$Disease$`Tobacco`$ATC,columnName=columnName)

return (Tobacco)

}

isAlcohol <- function(df, columnName) {
   
Alcohol <- isTreatedByIt(df, CCAM_codes = NULL, 
                                 ICD_Codes = data$Disease$`Alcohol`$ICD10,ATC_codes=data$Disease$`Alcohol`$ATC,columnName=columnName)

return (Alcohol)

}