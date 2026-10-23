#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R

library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(TFBSTools)

getJASPAR2020Core <- function() {
    opts <- list()
    opts[["species"]] <- 9606
    opts[["all_versions"]] <- FALSE
    opts[["collection"]] <- "CORE"
    jaspar_motifs <- getMatrixSet(JASPAR2020, opts)
    names = sapply(jaspar_motifs, function(m) m@name )
    family = sapply(jaspar_motifs, function(m) m@tags$family)
    df = DataFrame(name = as.character(names), family = as.character(family))
    table(df$family)
                    
  # Create a new dataframe to store the expanded entries
  new_rows <- data.frame(name = character(), family = character(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(df)) {
    name_val <- df$name[i]
    family_val <- df$family[i]
    
    # Check if the family contains 'c(' indicating multiple families
    if (grepl("^c\\(", family_val)) {
      # Extract individual family names
      families <- eval(parse(text = family_val))
      
      # Add a row for each family
      for (fam in families) {
        new_rows <- rbind(new_rows, data.frame(name = name_val, family = fam, stringsAsFactors = FALSE))
      }
    } else {
      # Just add the single family row
      new_rows <- rbind(new_rows, data.frame(name = name_val, family = family_val, stringsAsFactors = FALSE))
    }
  }
  
  # Fix common typos
  new_rows$family <- gsub("Kruppel", "Krüppel", new_rows$family)
  new_rows$family <- gsub("Three-zinc finger Krüppel-related factors", "Three-zinc finger Krüppel-related factors", new_rows$family)
  new_rows$family <- gsub("Interferon-regulatory factors   ", "Interferon-regulatory factors", new_rows$family)

  new_rows$name <- gsub("::(.*)$", "..\\1", new_rows$name)
  new_rows$name <- gsub("::(.*)$", "..\\1", new_rows$name)
  new_rows$name <- gsub("\\((.*)\\)", ".\\1", new_rows$name)
  new_rows$name <- gsub("-(.*)$", ".\\1", new_rows$name)
  # Trim whitespace
  new_rows$family <- trimws(new_rows$family)
  
  return(new_rows)
}

extractTFperConditionDF <- function(se_obj, assay_name, threshold, tfdatabase = "JASPAR2020") {
   # Extract the assay matrix
  if (!assay_name %in% SummarizedExperiment::assayNames(se_obj)) {
    stop(paste("Assay", assay_name, "not found in the SummarizedExperiment object."))
  }
  assay_data <- SummarizedExperiment::assay(se_obj, assay_name)
  
  # Ensure rownames are present (TF names)
  if (is.null(rownames(assay_data))) {
    stop("Row names (transcription factors) are missing in the SummarizedExperiment object.")
  }
  
  # Long format
  df <- as.data.frame(assay_data) %>%
    tibble::rownames_to_column(var = "transcription_factor") %>%  # Preserve TF names
    pivot_longer(cols = -transcription_factor, names_to = "Condition", values_to = assay_name)  # Melt
  
  # Filter by threshold
  df <- df %>% filter(.data[[assay_name]] >= threshold)
  df <- df %>% mutate(name_tf = gsub("_.*", "", transcription_factor))

  if (tfdatabase != "JASPAR2020") {
       print("TF annotation not supported, returning df without family name")
   } else {
      annotation_df <- getJASPAR2020Core()
   }
    df <- df %>% left_join(annotation_df, by = join_by(name_tf == name), relationship = "many-to-many")
  return(df)
}


                       


