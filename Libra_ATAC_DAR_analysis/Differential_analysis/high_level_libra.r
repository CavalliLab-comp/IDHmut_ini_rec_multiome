#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

# This script contains functions for the Libra diff analysis pipeline
# Author : Thomas Berthomier



library(parallelly)
library(ArchR)
library(ggplot2)
library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)

CODE_PATH <- '/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas'
source(file.path(CODE_PATH, "Differential_analysis", "custom_libra_functions.r"))


#' Convert Gene Names to GRanges Object
#'
#' Parses genomic coordinate strings in the format "chr1:1000-2000" and converts
#' them to a GRanges object.
#'
#' @param gene_names Character vector of genomic coordinates in format "chr1:1000-2000"
#'
#' @return GRanges object with parsed genomic coordinates
#'
#' @details The function expects gene names in the format "chrN:start-end" where
#' N can be a chromosome number (1-22), sex chromosome (X, Y), or mitochondrial (M).
#' If parsing fails, the function will throw an error with a descriptive message.
#'
#' @examples
#' gene_coords <- c("chr1:1000-2000", "chr2:5000-6000", "chrX:10000-11000")
#' gr <- gr_from_gene(gene_coords)
gr_from_gene <- function(gene_names) {
    # Try to parse expected format (chr1:1000-2000)
    parsed <- tryCatch({
      strcapture("chr([0-9XYM]+):([0-9]+)-([0-9]+)",
                 gene_names,
                 proto = list(chr = character(), start = integer(), end = integer()))
    }, error = function(e) {
      stop("Gene names not in expected 'chr1:1000-2000' format")
    })
    
    # Check if parsing was successful
    if(all(is.na(parsed$chr))) {
      stop("Failed to parse genomic coordinates from gene names")
    }
    
    GenomicRanges::GRanges(seqnames = paste0("chr", parsed$chr),
            ranges = IRanges::IRanges(start = parsed$start, end = parsed$end))
}

#' Perform Pseudobulk Differential Accessibility Analysis using Modified Libra
#'
#' Conducts pseudobulk differential accessibility analysis on an ArchR project using
#' a modified version of the Libra framework. Designed specifically for scATAC-seq
#' data analysis with support for peak matrices and gene score matrices.
#'
#' @param proj_archr ArchR project object containing single-cell ATAC-seq data
#' @param diff_annotation String specifying the column name in cellColData for the
#'   condition/treatment variable to test
#' @param group_forward Tested group level for categorical comparisons (required if
#'   diff_annotation is not numeric)
#' @param group_background Control group level for categorical comparisons (required if
#'   diff_annotation is not numeric)
#' @param cell_type String specifying column name for cell type annotations.
#'   Default: "cell.type"
#' @param replicate_col String specifying column name for biological replicates. It is 
#'   the parameter used to pair samples together
#'   Default: "Patient.number"
#' @param control_col Optional string specifying column name for additional control
#'   variables. Default: NULL
#' @param sample_col String specifying column name for sample identifiers.
#'   It is used in the case where two samples have the same 'replicate_col' label but the 
#'   same 'diff_annotation' label (It is a by product of modifying the Libra package instead of 
#'   rewriting from the beggining...)
#'   Default: "Sample.number"
#' @param paired_by_replicate Logical indicating whether to perform paired analysis
#'   by replicate. Default: FALSE
#' @param useMatrix String specifying which matrix to use ("PeakMatrix" for scATAC-seq
#'   peaks or "GeneScoreMatrix" for gene activity scores). Default: "PeakMatrix"
#' @param de_method String specifying differential expression method. Default: "DESeq2"
#' @param de_type String specifying test type for differential expression. Default: "Wald"
#' @param minCell Integer specifying minimum number of cells required per group.
#'   Default: 150
#' @param impute_matrix Logical indicating whether to impute the matrix using ArchR's
#'   imputation weights with MAGIC. Default: FALSE (never actually tried it)
#' @param n_threads Integer specifying number of threads for parallel processing.
#'   Default: half of available cores
#'
#' @return List containing differential analysis results (data frame with statistics)
#'   and the pseudobulk matrices used for the analysis (output from custom_run_de)
#'
#' @details This function performs several key steps:
#' \itemize{
#'   \item Subsets the ArchR project to cells with valid annotations
#'   \item Extracts the specified matrix (peaks or gene scores) from the ArchR project
#'   \item Formats peak names as genomic coordinates (chr:start-end)
#'   \item Optionally applies matrix imputation using ArchR's MAGIC imputation weights
#'   \item Runs pseudobulk differential analysis using a modified Libra framework
#' }
#'
#' The function uses a modified version of Libra's run_de function (custom_run_de)
#' that returns both statistical results and the pseudobulk matrices used in analysis.
#'
#' For categorical comparisons, both group_forward and group_background must be specified.
#' For continuous variables, the diff_annotation column should be numeric.
#'
#' @examples
#' \dontrun{
#' # Differential accessibility between treatment groups
#' da_results <- getLibraPseudoBulkAnalysis(
#'   proj_archr = atac_proj,
#'   diff_annotation = "Treatment",
#'   group_forward = "treated",
#'   group_background = "untreated",
#'   cell_type = "cell_type",
#'   useMatrix = "PeakMatrix",
#'   minCell = 100
#' )
#'
#' }
#'
#' @export
getLibraPseudoBulkAnalysis <- function(proj_archr, diff_annotation, cell_type,
                   group_forward = NULL,
                   group_background = NULL,
                   replicate_col = "Patient.number",
                   control_col = NULL,
                   sample_col = "Sample.number",
                   paired_by_replicate = FALSE,
                   useMatrix = "PeakMatrix",
                   de_method = "DESeq2",
                   de_type = "Wald",
                   minCell = 150,
                   impute_matrix = FALSE,
                   n_threads = ceiling(parallelly::availableCores() / 2)
                   ) {

  if (!is.numeric(proj_archr@cellColData[,diff_annotation])) {
    if (is.null(group_forward) || is.null(group_background)) {
      stop("Control group or condition group is not stated")
    }
    idx_condition <- BiocGenerics::which(proj_archr@cellColData[,diff_annotation] %in% c(group_forward, group_background))
    cellsSample <- proj_archr$cellNames[idx_condition]
    proj_subset <- proj_archr[cellsSample, ]
  } else {
    idx_condition <- BiocGenerics::which(!is.na(proj_archr@cellColData[,diff_annotation]))
    cellsSample <- proj_archr$cellNames[idx_condition]
    proj_subset <- proj_archr[cellsSample, ]
  }
  if (!is.null(control_col)){
    meta <- proj_subset@cellColData[, c(replicate_col, cell_type, diff_annotation, control_col, sample_col)]
    meta[, control_col] <- as.vector(meta[, control_col])
  } else {
    meta <- proj_subset@cellColData[, c(replicate_col, cell_type, diff_annotation, sample_col)]
  }
  meta[, replicate_col] <- as.vector(meta[, replicate_col]) 
  meta[, sample_col] <- as.vector(meta[, sample_col]) 
  meta[, cell_type] <- as.vector(meta[, cell_type])
  if (!is.numeric(proj_archr@cellColData[,diff_annotation])) {
    meta[, diff_annotation] <- relevel(factor(meta[, diff_annotation]), ref = group_forward)
  } else {
    meta[, diff_annotation] <- as.vector(meta[, diff_annotation])
  }

  if (useMatrix == "PeakMatrix") {
    type_analysis <- 'scATAC'
  } else if (useMatrix == "GeneScoreMatrix") {
    type_analysis <- 'scRNA' 
  } else {
    message(paste("Matrix type", useMatrix, "unknown. Defaulting type of analysis to scRNA.", sep = " "))
    type_analysis <- 'scRNA'
  }
  
  peakmatrix <- getMatrixFromProject(
    ArchRProj = proj_subset,
    useMatrix = useMatrix,
    binarize = FALSE,
    verbose = TRUE
  )

  if (useMatrix == "PeakMatrix") {
    peaks <- rowRanges(peakmatrix)
    peak_labels <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
    rownames(peakmatrix) <- peak_labels
  } else if (useMatrix == "GeneScoreMatrix") {
    rownames(peakmatrix) <- rowData(peakmatrix)$name
  }

  if (impute_matrix) {
    iW <- getImputeWeights(proj_archr)
    # Impute
    assay(peakmatrix) <- imputeMatrix(assay(peakmatrix), iW) # TO CHECK !!!!!!!!!
  }

  DE <- custom_run_de(assay(peakmatrix),
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type,
    label_col = diff_annotation,
    control_col = control_col,
    sample_col = sample_col,
    paired_by_replicate = paired_by_replicate,
    min_cells = minCell,
    de_family = "pseudobulk",
    de_method = de_method,
    de_type = de_type,
    input_type = type_analysis,
    n_threads = n_threads
  )
  return(DE)
}

#' Convert Libra Output to ArchR-Compatible SummarizedExperiment
#'
#' Converts differential accessibility results from Libra format to a SummarizedExperiment
#' object that mimics ArchR's differential analysis output format. This enables
#' downstream integration with ArchR's analysis pipeline, including motif enrichment analysis.
#'
#' @param df Data frame containing Libra differential analysis results with required
#'   columns: cell_type, gene, avg_logFC, p_val, p_val_adj
#' @param peak_gr GRanges object containing peak coordinates for overlap matching.
#'   Typically obtained from getPeakSet() in ArchR
#'
#' @return SummarizedExperiment object formatted to match ArchR's differential
#'   analysis output, containing:
#' \itemize{
#'   \item assays: Log2FC, FDR, and Pval matrices with cell types as columns
#'   \item rowData: Genomic coordinates (seqnames, start, end) matching ArchR format
#'   \item colData: Cell type information
#'   \item metadata: Analysis parameters including useMatrix="PeakMatrix" for
#'     compatibility with ArchR's PWM scanning functions
#' }
#'
#' @details The function performs several key operations to ensure ArchR compatibility:
#' \itemize{
#'   \item Parses genomic coordinates from the gene column using gr_from_gene()
#'   \item Maps results to provided peak coordinates via findOverlaps() for quality
#'     control and validation, ensuring differential peaks correspond to actual
#'     called peaks in the dataset
#'   \item Reshapes data from long to wide format with cell types as columns
#'   \item Creates separate assay matrices for fold changes, FDR, and p-values
#'   \item Sets appropriate metadata for ArchR pipeline integration
#'   \item Removes ranges that don't overlap with the provided peak set to maintain
#'     data integrity and prevent downstream analysis errors. It should not have any
#'     effect but it is just to be sure.
#' }
#'
#' This SummarizedExperiment can be directly used with ArchR functions like:
#' \itemize{
#'   \item addMotifAnnotations() for motif scanning
#'   \item peakAnnoEnrichment() for peak annotation enrichment
#'   \item plotBrowserTrack() for genomic visualization
#' }
#'
#' Ranges that don't overlap with peak_gr are removed with a warning message
#' indicating how many ranges were filtered out.
#'
#' @examples
#' \dontrun{
#' # Convert Libra results for ArchR pipeline
#' peak_set <- getPeakSet(archr_proj)
#' libra_se <- convertLibraOutToSE(
#'   df = libra_results,
#'   peak_gr = peak_set
#' )
#' }
#'
#' @export
convertLibraOutToSE <- function(df, peak_gr) {
  # Check for required columns first
  required_cols <- c("cell_type", "gene", "avg_logFC", "p_val", "p_val_adj")
  stopifnot(all(required_cols %in% colnames(df)))
  
  # Extract cell type names before transformation
  col_names <- unique(df$cell_type)
  
  # Add index if not already present
  df <- df %>% mutate(idx = row_number())
  
  # Create assay matrices: wide format with cell_type as columns
  df_wide <- df %>%
    select(cell_type, gene, Log2FC = avg_logFC, FDR = p_val_adj, Pval = p_val) %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = c(Log2FC, FDR, Pval))
  
  # Reconstruct GRanges from gene column
  
  gene_gr <- gr_from_gene(df_wide$gene)
  
  # Find overlaps with peak_gr
  overlaps <- findOverlaps(gene_gr, peak_gr, select = "first")
  
  # Keep only valid matches
  valid_idx <- which(!is.na(overlaps))
  if (length(valid_idx) < length(overlaps)) {
    warning(sprintf("%d ranges did not overlap and were removed.", 
                  length(overlaps) - length(valid_idx)))
  }
  
  # Filter df_wide and gene_gr accordingly
  df_wide <- df_wide[valid_idx, ]
  row_gr <- peak_gr[overlaps[valid_idx]]
  
  # Consistent row names across all objects
  row_names <- as.character(seq_along(row_gr))
  names(row_gr) <- row_names
  
  # Create assay matrices
  make_matrix <- function(prefix) {
    mat <- df_wide %>%
      select(starts_with(prefix)) %>%
      as.matrix()
    colnames(mat) <- sub(paste0("^", prefix, "_"), "", colnames(mat)) # Strip prefix consistently
    rownames(mat) <- row_names  # Use consistent row names
    mat
  }
  
  assays_list <- list(
    Log2FC = make_matrix("Log2FC"),
    FDR    = make_matrix("FDR"),
    Pval   = make_matrix("Pval")
  )
  
  # Construct rowData with consistent row names
  row_df <- as.data.frame(row_gr)[, c("seqnames", "start", "end")]
  rownames(row_df) <- row_names
  
  # Final check: all parts must have the same number of rows
  stopifnot(
    length(row_gr) == nrow(df_wide),
    nrow(assays_list$Log2FC) == nrow(df_wide),
    nrow(row_df) == nrow(df_wide),
    all(dim(assays_list$Log2FC)[2] == sapply(assays_list, ncol)) # Check column counts too
  )
  
  # Build SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays_list,
    rowData = row_df,
    colData = S4Vectors::DataFrame(row.names = col_names)
  )
  
  # Set metadata for ArchR
  metadata(se)$Params$useMatrix <- "PeakMatrix"
  
  return(se)
}

