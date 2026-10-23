#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R

library(GenomicRanges)
library(readr)
library(rtracklayer)

# This is a collection of functions to extract formatted
# objects using files. 

# I agree that this first one is lazy
#' Extract and format cell names for ArchR compatibility
#'
#' Reads cell filtering file and formats cell names using ArchR naming convention
#' (sample_name#cell_barcode).
#'
#' @param name_sample Character string of sample name for prefixing
#' @param cells_path Path to tab-delimited cell filtering file (no header)
#' @return Character vector of ArchR-formatted cell names
extractCells <- function(name_sample, cells_path){ 
    cells <- read.table(cells_path,sep="\t",header=F)
    cells_formatted <- sapply(cells, function(x) paste(name_sample, x, sep = "#"))
    return(cells_formatted)
}

#' Extract and format cell annotations for ArchR compatibility
#'
#' Reads cell annotation file and formats with ArchR naming convention, adding
#' malignancy classification based on provided cell types.
#'
#' @param name_sample Character string of sample name for prefixing cell names
#' @param anno_path Path to space-delimited cell annotation file (with header)
#' @param malignant_anno_vec Character vector of cell types considered malignant
#' @return Data frame with ArchR-formatted cell names and malignancy annotations
extractCellAnnotation <- function(name_sample, anno_path, malignant_anno_vec) {
    cell_annotation <- read.table(anno_path, sep=" ", header=T)
    cell_annotation$cell.malignant <- cell_annotation$cell.type %in% malignant_anno_vec
    cell_annotation$cell.name <- paste(name_sample, cell_annotation$cell.name, sep = "#")
    return(cell_annotation)
}                             
                              
#' Extract and process Locus-specific Methylation Changes (LMC) data
#'
#' Reads LMC CSV files, converts to GRanges objects, optionally performs genome
#' coordinate conversion (hg19 to hg38), and simplifies annotations.
#'
#' @param vec_LMC Character vector of LMC file names to process
#' @param LMC_path Path to directory containing LMC files and chain file
#' @param convert Logical, whether to convert coordinates from hg19 to hg38 (default: TRUE)
#' @return Named list of GRanges objects with processed LMC data
extractLMC <- function(vec_LMC, LMC_path, convert = TRUE) {
    LMC_path_list = lapply(vec_LMC, function(x) file.path(LMC_path, x)) 
    names(LMC_path_list) = gsub("(hyper|hypo)_LMC(\\d+)_annotated.csv", "LMC\\2_\\1", vec_LMC)

    df_list <- Map(function(x) readr::read_csv(x, col_select = c(1,2,3, 4, 26), show_col_types = FALSE), LMC_path_list)
    GR_list <- Map(function(x) GRanges(x$Chromosome, IRanges(x$Start, x$End, names = x$...1), UniversalAnnotation = x$UniversalAnnotation), df_list)
    
    if (convert) {
        chain_file <- file.path(LMC_path, "hg19ToHg38.over.chain")
        chain <- import.chain(chain_file)

        GR_list <- Map(function(x) unlist(liftOver(x, chain)), GR_list)
    }

    GR_list_final <- Map(function(x) {
      elementMetadata(x)[,"Annotation_simplified"] <- gsub("[0-9]*", "", gsub("^[0-9]*_", "", elementMetadata(x)[,"UniversalAnnotation"]))
      x
    }, GR_list)
    
    return(GR_list_final)
}
