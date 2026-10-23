#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

library(ArchR)
library(SummarizedExperiment)

# These functions provide tools for combining, extracting, converting, 
# and summarizing data from SummarizedExperiment objects.


#' Combine SummarizedExperiment List into Single Object
#'
#' Combines a list of SummarizedExperiment objects with common row structure
#' into a single SE object by column binding.
#'
#' @param se_list Named list of SummarizedExperiment objects with identical row dimensions
#'
#' @return Single SummarizedExperiment object with combined columns
#'
#' @examples
#' combined_se <- combineSElist(list(sample1 = se1, sample2 = se2))
#'
#' @export
combineSElist <- function(se_list) {
  se_list_out <- lapply(names(se_list), function(name) {
    colnames(se_list[[name]]) <- name
    se_list[[name]]
  })
  
  output <- do.call(cbind, se_list_out)
  if ("name" %in% colnames(rowData(output))) {
    rownames(output) <- rowData(output)$name
  }
  output
}

#' Extract Matrix from SummarizedExperiment List
#'
#' Combines SE list and extracts specified assay values for given features.
#' Handles special case of signed -log10 FDR calculation.
#'
#' @param se_list Named list of SummarizedExperiment objects
#' @param list_rownames Character vector of row names to extract
#' @param value_in_matrix Character. Assay name to extract, or "signed.neglog10FDR" for special calculation
#'
#' @return Matrix with extracted values for specified features
#'
#' @examples
#' mat <- extractMatrixSE(se_list, c("gene1", "gene2"), "Log2FC")
#' mat <- extractMatrixSE(se_list, features, "signed.neglog10FDR")
#'
#' @export
extractMatrixSE <- function(se_list, list_rownames, value_in_matrix) {
  combined_se <- combineSElist(se_list)

  # extract matrix value for selected features
  if (value_in_matrix == "signed.neglog10FDR") {
    combined_df <- -log10(as.matrix(assays(combined_se[list_rownames, ])[["FDR"]])) *
      sign(as.matrix(assays(combined_se[list_rownames, ])[["Log2FC"]]))
  } else {
    combined_df <- data.frame(assays(combined_se[list_rownames, ])[[value_in_matrix]])
  }
  # Combine the data frames into a single matrix
  return(as.matrix(combined_df))
}


#' Convert SummarizedExperiment to Data Frame
#'
#' Extracts all assays from a single SummarizedExperiment sample and combines
#' with rowData into a comprehensive data frame.
#'
#' @param se SummarizedExperiment object
#'
#' @return Data frame with rowData and all assay columns combined
#'
#' @examples
#' df <- getDFfromSESample(se_sample)
#'
#' @export
getDFfromSESample <- function(se) {
  df_out <- do.call(cbind, assays(se))
  colnames(df_out) <- names(assays(se))
  
  df_out <- cbind(rowData(se), df_out)
  rownames(df_out) <- rownames(se)

  if ("name" %in% colnames(df_out)) {
      rownames(df_out) <- df_out$name
  }
  as.data.frame(df_out)
}


#' Calculate Group Means Matrix from RangedSummarizedExperiment
#'
#' Computes group-wise means for each feature across cells for a given grouping
#' with optional logTP10K normalization. 
#'
#' @param rse RangedSummarizedExperiment object
#' @param assay_name Character. Name of assay to use (default: "PeakMatrix")
#' @param group_by Character. Column name in colData for grouping (default: "Sample")
#' @param peak_rownames Logical. Use genomic coordinates as row names (default: TRUE)
#' @param normalize Character. Normalization method: "none" or "logTP10K" (default: "none")
#'
#' @return Matrix with group means, rows as features, columns as groups
#'
#' @examples
#' means_mat <- getGroupedMeansMatFromRSE(rse, group_by = "CellType")
#' means_mat <- getGroupedMeansMatFromRSE(rse, normalize = "logTP10K")
#'
#' @export
getGroupedMeansMatFromRSE <- function(
  rse,
  assay_name = "PeakMatrix",
  group_by = "Sample",
  peak_rownames = TRUE,
  normalize = c("none", "logTP10K")
) {
  normalize <- match.arg(normalize)

  # Check inputs
  if (!assay_name %in% assayNames(rse)) {
    stop(paste("Assay", assay_name, "not found in the RangedSummarizedExperiment"))
  }
  if (!group_by %in% colnames(colData(rse))) {
    stop(paste("Column", group_by, "not found in colData"))
  }

  mat <- assay(rse, assay_name)

  # Normalization
  if (normalize == "logTP10K") {
    if (!"nFrags" %in% colnames(colData(rse))) {
      stop("Column 'nFrags' not found in colData, required for logTP10K normalization.")
    }
    nFrags <- colData(rse)$nFrags
    if (any(nFrags <= 0)) {
      stop("All values in 'nFrags' must be positive for logTP10K normalization.")
    }
    scaling_factors <- nFrags / 10000

    if (inherits(mat, "dgCMatrix")) {
      # Sparse matrix: use Matrix::Diagonal for efficient column scaling
      mat <- mat %*% Matrix::Diagonal(x = 1 / scaling_factors)
    } else {
      # Dense matrix: safe transpose scaling
      mat <- t(t(mat) / scaling_factors)
    }
    
    mat <- log1p(mat)
  }

  grouping <- colData(rse)[[group_by]]
  groups <- unique(grouping)

  group_means <- lapply(groups, function(g) {
    cols <- which(grouping == g)
    if (inherits(mat, "dgCMatrix")) {
      Matrix::rowMeans(mat[, cols, drop = FALSE])
    } else {
      rowMeans(mat[, cols, drop = FALSE])
    }
  })

  result <- do.call(cbind, group_means)
  colnames(result) <- groups

  if (peak_rownames) {
    peaks <- rowRanges(rse)
    peak_labels <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
    rownames(result) <- peak_labels
  } else if (!is.null(rownames(rse))) {
    rownames(result) <- rownames(rse)
  }

  return(result)
}