#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

# This script contains functions to build an archr project from cellranger output + sample metadata tsv
# some higher level functions are definitly overkill from when i was in my
# "overfunctionalisation" phase

library(ArchR)

CODE_PATH <- "/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas"
source(file.path(CODE_PATH, "Miscellaneous", "extractData.r"))
source(file.path(CODE_PATH, "Differential_analysis", "deprecated_SE_operations.r"))



# can't use merge directly as it breaks the archrproject object
# this function takes in a project and a metadata df to return a project with added Metadata at the sample level
addSampleMetadataDf <- function(project, metadata_df, merge_proj = "Sample", merge_meta = "Sample") {
    rownames(metadata_df) <- metadata_df[[merge_meta]]
    
    # Ensure metadata is aligned with project cell data
    if (merge_proj == "cellNames"){
        coldata_to_use <- metadata_df[as.character(as.vector(project$cellNames)), , drop = FALSE]
    } else {
        coldata_to_use <- metadata_df[as.character(as.vector(project@cellColData[[merge_proj]])), , drop = FALSE]
    }
    
    # Add metadata to project
    metadata_cols <- setdiff(colnames(metadata_df), merge_meta)
    lapply(metadata_cols, function(col) {
        project <- addCellColData(
            ArchRProj = project,
            data = coldata_to_use[, col],
            cells = project$cellNames,
            name = col, 
            force = TRUE
        )
    })
    
    return(project)
}

# sample level : 
# add cell annotation to a project sample
addCellAnnotation <- function(project, name_sample, anno_path, malignant_anno_vec){
    cell_annotation <- extractCellMetadata(name_sample, anno_path, malignant_anno_vec)
    project <- addSampleMetadataDf(project, cell_annotation, merge_proj = "cellNames", merge_meta = "cell.name")
    return(project)
}

#add sample annotation and cell annotation to a project sample
addAnnotationSample <- function(project, Sample_metadata_df, name_sample, cell_anno_path, malignant_anno_vec){
    project <- addSampleMetadataDf(project, Sample_metadata_df)
    project <- addCellAnnotation(project, name_sample, cell_anno_path, malignant_anno_vec)
    return(project)
}

#large project level 
#filter a multisample project for valid cells and annotate it 
filterAndAnnotateLargeProject  <- function(fused_project, path_pseudoraw, Sample_metadata_df, malignant_anno_vec, outputDirectory_large, logDirectory_large) {
    sample_names <- getSampleNames(ArchRProj = fused_project)
    #get paths
    cells_path_vec <- file.path(path_pseudoraw, sample_names, "cells_pass_filters.txt.gz")
    anno_path_vec <- file.path(path_pseudoraw, sample_names, "cell_annotation.txt")

    print("Loading necessary files for filtering and annotation...")
    cells_project <- unlist(mapply(extractCells, sample_names, cells_path_vec)) #gives a vector with valid cells for all samples
    names(cells_project) <- NULL
    
    anno_project <- do.call(rbind, 
                            mapply(extractCellAnnotation,
                                          sample_names, 
                                          anno_path_vec, 
                                          replicate(length(sample_names), malignant_anno_vec, simplify = FALSE),
                                          SIMPLIFY = FALSE)
                           )

    #subset large 
    print("Subsetting the project...")
    fused_project_filtered <- subsetProject_byCellVec(fused_project, "fused_filtered", cells_project, outputDirectory_large, logDirectory_large)

    #annotate large
    print("Adding sample annotations...")
    fused_project_filtered <- addSampleMetadataDf(fused_project_filtered, Sample_metadata_df)
    fused_project_annotated <- addSampleMetadataDf(fused_project_filtered, anno_project, merge_proj = "cellNames", merge_meta = "cell.name")
    return(fused_project_annotated)
}
