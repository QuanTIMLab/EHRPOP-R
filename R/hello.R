#' @import R6
#' @import dplyr
#' @import tidyr
#' @import networkD3
#' @import stringr
#' @import jsonlite


library(R6)
library(dplyr)
library(tidyr)
library(networkD3)
library(stringr)
library(jsonlite)

TreatmentProcessor <- R6Class("TreatmentProcessor",
  
  public = list(
    
    json_path = NULL,
    index_date = "DATE",
    index_code_ccam = "CODE_CCAM",
    index_code_atc = "CODE_ATC",
    index_code_icd = "CODE_ICD10",
    index_id = "ID_PATIENT",
    index_quantity = "QUANTITE",
    data = NULL,
    
    initialize = function(json_path, 
                          index_date = "DATE", 
                          index_code_ccam = "CODE_CCAM",
                          index_code_atc = "CODE_ATC",
                          index_code_icd = "CODE_ICD10",
                          index_id = "ID_PATIENT",
                          index_quantity = "QUANTITE") {
      self$json_path <- json_path
      self$index_date <- index_date
      self$index_code_ccam <- index_code_ccam
      self$index_code_atc <- index_code_atc
      self$index_code_icd <- index_code_icd
      self$index_id <- index_id
      self$index_quantity <- index_quantity
      
      self$read_traitements_data()  
    },
    
    read_traitements_data = function() {
      self$data <- fromJSON(self$json_path)  
    },
    
    print_data = function() {
      print(self$data)
    },
      
    find_code = function(code) {
        
      if ("Treatment" %in% names(self$data) && is.list(self$data$Treatment)) {

        for (treatment_type in names(self$data$Treatment)) {
          treatment_data <- self$data$Treatment[[treatment_type]]

          if (is.list(treatment_data)) {

            for (code_list in treatment_data) {
              if (is.vector(code_list) && code %in% code_list) {
                return(treatment_type)
              }
            }
          }
        }
      }

      return(NA)
    },
    
    sankey_diagram = function(df) {


        map_code <- function(code) {

            if ("Treatment" %in% names(self$data) && is.list(self$data$Treatment)) {

            for (treatment_type in names(self$data$Treatment)) {
              treatment_data <- self$data$Treatment[[treatment_type]]

              if (is.list(treatment_data)) {

                for (code_list in treatment_data) {

                  if (is.vector(code_list) && code %in% code_list) {
                    return(treatment_type)
                  }
                }
              }
            }
          }

          return(NA)
        }



        df <- df %>%
          mutate(CODE_ACT = case_when(
            !!sym(self$index_code_atc) != "" ~ !!sym(self$index_code_atc),
            !!sym(self$index_code_ccam) != "" ~ !!sym(self$index_code_ccam),
            !!sym(self$index_code_icd) != "" ~ !!sym(self$index_code_icd),
            TRUE ~ NA_character_
          ))

        df <- df %>%
          select(-!!sym(self$index_code_icd), -!!sym(self$index_code_ccam), -!!sym(self$index_code_atc))

        df <- df %>%
          drop_na(CODE_ACT)

        df <- df %>%
          mutate(CODE_ACT = sapply(CODE_ACT, map_code))

        df <- df %>%
          mutate(!!sym(self$index_date) := as.numeric(as.character(!!sym(self$index_date)))) %>%
          mutate(!!sym(self$index_date) := ifelse(is.na(!!sym(self$index_date)), 999999, !!sym(self$index_date)))


        df <- df %>%
          arrange(!!sym(self$index_id), !!sym(self$index_date))



        df <- df %>%
          drop_na(CODE_ACT)


        remove_consecutive_duplicates <- function(treatments) {
          result <- c()

          for (i in seq_along(treatments)) {
            if (i == 1 || treatments[i] != treatments[i - 1]) {
              result <- c(result, treatments[i])
            }
          }

          return(str_c(result, collapse = "->"))
        }
        

        result <- df %>%
          group_by(!!sym(self$index_id)) %>%
          summarize(Traitements = remove_consecutive_duplicates(CODE_ACT)) %>%
          ungroup()


        label_treatments <- function(treatment_string) {
          treatments <- unlist(strsplit(treatment_string, "->"))
          labeled_treatments <- paste0(treatments, seq_along(treatments))
          return(labeled_treatments)
        }


        treatment_counts <- table(result$Traitements)
        frequent_treatments <- names(treatment_counts[treatment_counts >= 10])
        filtered_df <- result %>% filter(Traitements %in% frequent_treatments)
        
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
        
        pairs <- do.call(rbind, strsplit(names(sequence_counter), "->"))
        df_links <- data.frame(
          source = pairs[,1],
          target = pairs[,2],
          value = unlist(sequence_counter)
        )

        all_nodes <- unique(c(df_links$source, df_links$target))
        node_indices <- setNames(seq_along(all_nodes) - 1, all_nodes)

        df_links$source <- node_indices[df_links$source]
        df_links$target <- node_indices[df_links$target]

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


    },
      
      
    yes_or_no = function(df, ccam_codes, atc_codes, icd_codes, column_name) {

      df <- df %>%
        mutate(MATCH = ((!!sym(self$index_code_ccam) %in% ccam_codes |
                         !!sym(self$index_code_atc) %in% atc_codes |
                         !!sym(self$index_code_icd) %in% icd_codes))) 


        result <- df %>%
        group_by_at(self$index_id) %>%
        summarize(!!column_name := any(MATCH)) %>% 
        ungroup()

      return(result)
    },
      
    is_treated_by_it <- function(df, ccam_codes, atc_codes, icd_codes, column_name) {
  
      # Handle cases where the code lists are empty
      ccam_match <- if (length(ccam_codes) > 0 && ccam_codes != "") {
        df[[self$index_code_ccam]] %in% ccam_codes
      } else {
        FALSE
      }

      atc_match <- if (length(atc_codes) > 0 && atc_codes != "") {
        df[[self$index_code_atc]] %in% atc_codes
      } else {
        FALSE
      }

      icd_match <- if (length(icd_codes) > 0 && icd_codes != "") {
        df[[self$index_code_icd]] %in% icd_codes
      } else {
        FALSE
      }

      # Combine matches
      df <- df %>%
        mutate(MATCH = ccam_match | atc_match | icd_match)

      result <- df %>%
        group_by_at(self$index_id) %>%
        summarize(!!column_name := any(MATCH)) %>%
        ungroup()

      return(result)
    },
      
    is_treated_by_it_percentage = function(df, ccam_codes, atc_codes, icd_codes) {
  
      df <- df %>%
        mutate(MATCH = (!!sym(self$index_code_ccam) %in% ccam_codes |
                        !!sym(self$index_code_atc) %in% atc_codes |
                        !!sym(self$index_code_icd) %in% icd_codes))

      treated_df <- df %>%
        group_by_at(self$index_id) %>%
        summarize(Treated = any(MATCH)) %>%
        ungroup()

      treated_patients <- treated_df %>%
        filter(Treated == TRUE) %>%
        nrow()

      total_patients <- treated_df %>%
        nrow()

      treated_percentage <- (treated_patients / total_patients) * 100

      return(paste0(treated_patients, " (", round(treated_percentage, 2), "%)"))
    },


    is_treated_by_it_with_qte = function(df, ccam_codes, atc_codes, icd_codes, column_name, quantite="QUANTITE") {

      starts_with_any <- function(code, codes_list) {
        code <- as.character(code)
        any(sapply(codes_list, function(prefix) startsWith(code, prefix)))
      }

      df <- df %>%
        mutate(
          !!sym(quantite) := ifelse(is.na(!!sym(quantite)), 1, !!sym(quantite)),
          SESSION = 1,
          !!sym(self$index_code_icd) := as.character(!!sym(self$index_code_ccam)),
          !!sym(self$index_code_icd) := as.character(!!sym(self$index_code_atc)),
          !!sym(self$index_code_icd) := as.character(!!sym(self$index_code_icd))
        )

      df <- df %>%
        mutate(
          Is_Relevant = case_when(
            starts_with_any(!!sym(self$index_code_ccam), ccam_codes) ~ TRUE,
            starts_with_any(!!sym(self$index_code_atc), atc_codes) ~ TRUE,
            starts_with_any(!!sym(self$index_code_icd), icd_codes) ~ TRUE,
            TRUE ~ FALSE
          )
        )

      result <- df %>%
        filter(Is_Relevant) %>%
        group_by(!!sym(self$index_id)) %>%
        summarise(!!sym(column_name) := sum(SESSION, na.rm = TRUE)) %>%
        ungroup()

      all_patients <- df %>%
        select(!!sym(self$index_id)) %>%
        distinct()

      result <- all_patients %>%
        left_join(result, by = self$index_id) %>%
        mutate(!!sym(column_name) := ifelse(is.na(!!sym(column_name)), 0, !!sym(column_name)))

      return(result)
    },
    traitement_characterization = function(df) {
  
      if (!exists("data") || !"Treatment" %in% names(self$data) || !is.list(self$data$Treatment)) {
        stop("The 'Treatment' section does not exist or is not a list.")
      }

      results <- list()

      for (treatment_type in names(self$data$Treatment)) {
        treatment_data <- self$data$Treatment[[treatment_type]]

        CODES_CCAM <- treatment_data[["CCAM"]]
        CODES_ATC <- treatment_data[["ATC"]]
        CODES_ICD10 <- treatment_data[["ICD10"]]

        main_result <- self$is_treated_by_it_percentage(df, CODES_CCAM, CODES_ATC, CODES_ICD10)

        results[[length(results) + 1]] <- data.frame(
          "Treatment_Type" = treatment_type,
          "Subcategory" = NA,
          "Percentage_Treated" = main_result
        )

        if (is.list(treatment_data)) {
          for (category in names(treatment_data)) {
            codes <- treatment_data[[category]]

            if (is.list(codes)) {
              for (sub_category in names(codes)) {
                sub_codes <- codes[[sub_category]]

                sub_CCM_codes <- sub_codes[["CCAM"]]
                sub_ATC_codes <- sub_codes[["ATC"]]
                sub_ICD_codes <- sub_codes[["ICD10"]]

                sub_result <- self$is_treated_by_it_percentage(df, sub_CCM_codes, sub_ATC_codes, sub_ICD_codes)

                results[[length(results) + 1]] <- data.frame(
                  "Treatment_Type" = treatment_type,
                  "Subcategory" = sub_category,
                  "Percentage_Treated" = sub_result
                )
              }
            }
          }
        }
      }

      df_results <- do.call(rbind, results)
      return(df_results)
    },
                   
    generate_patient_treatment_summary = function(df, list_columns, date_start = NULL, date_end = NULL) {
      if (length(list_columns) == 0) {
        stop("list_columns cannot be empty.")
      }

      missing_columns <- setdiff(list_columns, colnames(df))
      if (length(missing_columns) > 0) {
        stop(paste("The following columns are missing in the DataFrame:", paste(missing_columns, collapse = ", ")))
      }

      if (!is.null(date_start)) {
        base_date <- as.Date("1960-01-01")
        dateStart <- as.integer(as.Date(date_start) - base_date)
      }

      if (!is.null(date_end)) {
        base_date <- as.Date("1960-01-01")
        dateEnd <- as.integer(as.Date(date_end) - base_date)
      }

      results <- NULL

      for (coln in list_columns) {
        unique_values <- unique(na.omit(df[[coln]]))
        count <- 0
        for (value in unique_values) {

            r <- self$is_treated_by_it_with_qte(df[df[[coln]] == value, ], list(value), list(value), list(value), paste0( " ", value))

          if (count == 0) {
            results <- r
          } else {
            results <- merge(results, r, by = self$index_id, all = TRUE)
          }
          count <- count + 1
        }
      }

      return(results)
    },
                   

    # Define the tableSequances function
    generate_act_sequences = function(df, list_columns) {
      if (length(list_columns) == 0) {
        stop("list_columns cannot be empty.")
      }

      missing_columns <- setdiff(list_columns, colnames(df))
      if (length(missing_columns) > 0) {
        stop(paste("The following columns are missing in the DataFrame:", paste(missing_columns, collapse = ", ")))
      }

        combine_codes <- function(row, list_columns) {
          codes <- ""
          for (col in list_columns) {
            if (!is.na(row[[col]])) {
              codes <- row[[col]]
            }
          }
          return(codes)
        }

      df <- df %>%
        rowwise() %>%
        mutate(CODE_ACTS = combine_codes(cur_data(), list_columns)) %>%
        ungroup()

      df <- df %>%
        filter(CODE_ACTS != "") %>%
        select(!!sym(self$index_id), !!sym(self$index_date), CODE_ACTS)

      df <- df %>%
        arrange(!!sym(self$index_id), !!sym(self$index_date))

      grouped <- df %>%
        group_by(!!sym(self$index_id)) %>%
        summarise(DATES = list(!!sym(self$index_date)), CODE_ACTS = list(CODE_ACTS)) %>%
        ungroup()

      grouped <- grouped %>%
        rowwise() %>%
        mutate(ACTES = list(CODE_ACTS[order(DATES)])) %>%
        ungroup()

      final_df <- grouped %>%
        select(!!sym(self$index_id), DATES, ACTES)

      return(final_df)
    },
                   
    neoadjuvant_or_adjuvant_or_both = function(df, ccam_codes, atc_codes, 
                                               icd_codes, column_name, days_before, days_after,
                                               bc_index_surgery="BC_index_surgery") {

        df$DATE_DIFF <- df[[self$index_date]] - df[[bc_index_surgery]]
        df <- df[df$DATE_DIFF >= -days_before & df$DATE_DIFF <= days_after, ]

        df$Is_Relevant <- df[[self$index_code_ccam]] %in% ccam_codes | 
                          df[[self$index_code_atc]] %in% atc_codes   | 
                          df[[self$index_code_icd]] %in% icd_codes

        
        patient_classification <- list()

        unique_patients <- unique(df[[self$index_id]])

        for (patient_id in unique_patients) {
            if (is.na(patient_id)) {
                next
            }

            patient_df <- df[df[[self$index_id]] == patient_id & df$Is_Relevant, ]

            if (nrow(patient_df) == 0) {
                patient_classification[[as.character(patient_id)]] <- 'False'
                next
            }

            positive_count <- sum(patient_df[['DATE_DIFF']] > 0, na.rm = TRUE)
            negative_count <- sum(patient_df[['DATE_DIFF']] < 0, na.rm = TRUE)

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

        result <- setNames(
            data.frame(
                ID_PATIENT = names(patient_classification),
                Classification = unlist(patient_classification),
                stringsAsFactors = FALSE
            ),
            c(self$index_id, column_name)
        )


        return(result)
    },

    chemotherapy_intervals = function(df, ccam_codes, atc_codes, icd_codes,
                                       column_name, days_before, days_after, bc_index_surgery="BC_index_surgery") {


        df$DATE_DIFF <- df[[self$index_date]] - df[[bc_index_surgery]]
        df <- df[df$DATE_DIFF >= -days_before & df$DATE_DIFF <= days_after ,]

        df$Is_Relevant <- df[[self$index_code_ccam]] %in% ccam_codes | 
                          df[[self$index_code_atc]] %in% atc_codes   |
                          df[[self$index_code_icd]] %in% icd_codes
        
        patient_classification <- list()

        unique_patients <- unique(df[[self$index_id]])

        for (patient_id in unique_patients) {

          patient_df <- df[df[[self$index_id]] == patient_id & df$Is_Relevant, ]

          if (nrow(patient_df) == 0) {
            patient_classification[[as.character(patient_id)]] <- 'False'
            next
          }

          Text <- ""
          temp <- 0
          c <- 0


        for (i in 1:nrow(patient_df)) {
            data <- patient_df[i, ]

            if (c == 0) {
              Text <- "0"
              temp <- data[[self$index_date]]
            } else {
              Text <- paste0(Text, " -> ", as.character(data[[self$index_date]] - temp))
              temp <- data[[self$index_date]]
            }

            c <- c + 1
          }
        

          patient_classification[[as.character(patient_id)]] <- Text
        }


        valid_patients <- !is.na(names(patient_classification))

        result <- setNames(
          data.frame(
            ID_PATIENT = names(patient_classification)[valid_patients],
            unlist(patient_classification)[valid_patients],
            stringsAsFactors = FALSE
          ),
          c(self$index_id, column_name)
        )

        return(result)


        },
                   
    classify_regimen_chemo = function(text) {

      if (text == "False") {
        return(text)
      }

      numbers <- as.numeric(str_extract_all(text, "\\b\\d+\\.?\\d*\\b")[[1]])
      numbers <- numbers[numbers != 0]

      if (length(numbers) == 0) {
        return("ONE TREATMENT")
      }

      numbers[numbers %in% c(20, 22)] <- 21
      numbers[numbers %in% c(6, 8)] <- 7
      numbers[numbers %in% c(13, 15)] <- 14

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
          return(FALSE) 
        }
        before_transition <- all(seq[1:(transition_index - 1)] == 21)
        after_transition <- all(seq[transition_index:length(seq)] == 7)
        before_transition && after_transition
      }

      is_anthracyclines <- function(seq) {
        all(seq == 21) 
      }


      if (is_paclitaxel(numbers)) return('Paclitaxel')
      if (is_anthracyclines_docetaxel(numbers)) return('Anthracyclines/docetaxel')
      if (is_anthracyclines_paclitaxel(numbers)) return('Anthracyclines/paclitaxel')
      if (is_anthracyclines_paclitaxel2(numbers)) return('Anthracyclines/paclitaxel')
      if (is_anthracyclines(numbers)) return('Anthracyclines')

      return('Other')
    }











    
  )
)

