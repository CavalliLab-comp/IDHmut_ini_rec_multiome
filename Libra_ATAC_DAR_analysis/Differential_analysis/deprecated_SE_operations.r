#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

library(ArchR)
library(dplyr)
CODE_PATH <- '/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas'

# These functions were used in the methods with paired differential analysis
# with ArchR (that eventually did not work)

source(file.path(CODE_PATH, "Miscellaneous", "SE_operations.r"))

#extract list of SE from output FUllDIff
extractSElistfromFullDiff <- function(output,
                                      comparison, #name of the diff (ex : inivsrec1)
                                      se_type){ #Genes, MotifsUp, MotifsDown
    list_out <- lapply(names(output), function(cond) output[[cond]][[comparison]][[se_type]])
    names(list_out) <- names(output)
    return(list_out)
}


# give list of rownames passing a threshold of number or
# proportion passing the significance threshold
extractRowsNSig <- function(se_list,
                            se_type = "Motifs", 
                            threshold_FDR = 0.05, # converted to mlog10Padj for "Motifs"
                            threshold_Log2FC = 1,
                            frac_significant_min = NULL,
                            n_significant_min = 1) {

  if (!is.null(frac_significant_min)){
    n_significant_min <- ceiling(frac_significant_min*length(se_list))
  }
  
  list_df <- lapply(se_list, function(se) {getDFfromSESample(se)})

  if (se_type == "Genes") {
    list_filtered_names <- lapply(list_df, function(df) {
      filtered_df <- df[(df[["FDR"]] < threshold_FDR) & (!is.na(df[["FDR"]])) & (abs(df[["Log2FC"]]) > threshold_Log2FC), , drop = FALSE]
      rownames(filtered_df)
    })
  } else if (se_type == "Motifs") {
    list_filtered_names <- lapply(list_df, function(df) {
      threshold_mlogpval <- -log10(threshold_FDR)
      filtered_df <- df[(df[["mlog10Padj"]] > threshold_mlogpval) & (!is.na(df[["mlog10Padj"]])), , drop = FALSE]
      rownames(filtered_df)
    })
  }
  
  condition <- table(unlist(list_filtered_names)) >= n_significant_min
  return(names(condition[condition]))
}

# give list of rownames that are in topN
# for each column
extractRowstopN <- function(se_list,
                               se_type = "Motifs",
                               threshold_FDR = 0.05, # converted to mlog10Padj for "Motifs"
                               n_top = 10
                              ) {
  # Extract values into dataframes and get significant features
  list_df <- lapply(se_list, function(se) {getDFfromSESample(se)})

  if (se_type == "Genes") {
    list_top_names <- lapply(list_df, function(df) {
      filtered_df <- df[(df[["FDR"]] < threshold_FDR) & (!is.na(df[["FDR"]])), , drop = FALSE]
      filtered_df <- filtered_df %>% arrange(FDR, desc(Log2FC))
      nrows <- min(n_top, length(rownames(filtered_df)))
      rownames(filtered_df[seq_len(nrows),])
    })
  } else if (se_type == "Motifs") {
    list_top_names <- lapply(list_df, function(df) {
      threshold_mlogpval <- -log10(threshold_FDR)
      filtered_df <- df[(df[["mlog10Padj"]] > threshold_mlogpval) & (!is.na(df[["mlog10Padj"]])), , drop = FALSE]
      filtered_df <- filtered_df %>% arrange(desc(mlog10Padj))
      
      nrows <- min(n_top, length(rownames(filtered_df)))
      if (nrows > 0){
        rownames(filtered_df[seq_len(nrows), ])
      }
    })
  }
  # Find genes that pass thresholds in at least n_significant_min datasets
  condition <- unique(unlist(list_top_names))
  return(condition)
}


extractMatrixNSig <- function(se_list,
                              se_type = "Motifs",
                              value_in_matrix = NULL, 
                              threshold_FDR = 0.05, # converted to mlog10Padj for "Motifs"
                              threshold_Log2FC = 1,
                              frac_significant_min = NULL,
                              n_significant_min = 1) {
  names <- extractRowsNSig(se_list, se_type = se_type, threshold_FDR = threshold_FDR, 
                              threshold_Log2FC = threshold_Log2FC, frac_significant_min = frac_significant_min, n_significant_min = n_significant_min)
    
  matrix <- extractMatrixSE(se_list, names, value_in_matrix)
  return(matrix)
}

extractMatrixtopN <- function(se_list,
                              se_type = "Motifs",
                              value_in_matrix = NULL, # takes Log2FC by default if se_type = "Genes" and mlog10Padj for "Motifs"
                              threshold_FDR = 0.05,# converted to mlog10Padj for "Motifs"
                              n_top = 10) {

  names <- extractRowstopN(se_list,
                               se_type = se_type,
                               threshold_FDR = threshold_FDR,
                               n_top = n_top
                              )
  matrix <- extractMatrixSE(se_list, names, value_in_matrix)
  return(matrix)
}