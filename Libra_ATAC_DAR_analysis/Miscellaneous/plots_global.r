#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R
library(ggplot2)
library(ArchR)
# Two functions to plot all cells on a UMAP + Harmony


#' Update ArchR Project for Harmony Analysis
#'
#' Adds necessary matrices and embeddings for Harmony batch correction visualization.
#' Sequentially adds IterativeLSI, UMAP, Harmony correction, and Harmony UMAP if missing.
#'
#' @param project ArchR project object
#' @param sample Character. Sample/batch column name in cellColData (default: "Sample")
#' @param force Logical. Force recreation of existing matrices (default: FALSE)
#'
#' @return Updated ArchR project with all matrices needed for plotHarmony
#'
#' @export
updateProjectforplotHarmony <- function(project, sample = "Sample", force = FALSE) {
    available_matrices <- c(names(project@reducedDims), names(project@embeddings))
    print(available_matrices)
    harmony_name <- paste("Harmony", sample, sep = "_")
    if (!("IterativeLSI" %in% available_matrices) | force) {
        project <- addIterativeLSI(project, force = force)
    }
    if (!("UMAP" %in% available_matrices) | force) {
        project <- addUMAP(ArchRProj = project, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30,
            minDist = 0.5, metric = "cosine", force = force)
    }
    if (!(harmony_name %in% available_matrices) | force) {
        project <- addHarmony(ArchRProj = project, reducedDims = "IterativeLSI", name = harmony_name,
            groupBy = sample, force = force)
    }
    if (!(paste("UMAP", harmony_name, sep = "_") %in% available_matrices) | force) {
        project <- addUMAP(ArchRProj = project, reducedDims = harmony_name, name = paste("UMAP", harmony_name, sep = "_"),
                           nNeighbors = 30, minDist = 0.5, metric = "cosine", force = force)
    }
    return(project)
}


#' Plot Harmony Batch Correction Results
#'
#' Creates four-panel comparison of standard vs Harmony-corrected UMAP embeddings,
#' showing both sample/batch effects and biological annotations.
#'
#' @param project ArchR project object with necessary matrices (run updateProjectforplotHarmony first)
#' @param annotation Character. Biological annotation column name in cellColData
#' @param sample Character. Sample/batch column name in cellColData (default: "Sample.number")
#'
#' @return List of four ggplot objects: standard UMAP (sample, annotation), 
#'   Harmony UMAP (sample, annotation)
#'
#' @export
plotHarmony <- function(project, annotation, sample = "Sample.number") {
    
    subset_idx <- BiocGenerics::which(!is.na(project@cellColData[, annotation]))
    temp <- subsetCells(project, cell = project$cellNames[subset_idx])
    UMAP_df <- temp@embeddings[["UMAP"]][["df"]][temp$cellNames,]
    temp@embeddings[["UMAP"]][["df"]] <- UMAP_df[!grepl("^NA", rownames(UMAP_df)), ]

    HARMONY_df <- temp@embeddings[[paste("UMAP_Harmony", sample, sep = "_")]][["df"]][temp$cellNames,]
    temp@embeddings[[paste("UMAP_Harmony", sample, sep = "_")]][["df"]] <- HARMONY_df[!grepl("^NA", rownames(HARMONY_df)), ]
    plotEmbedding(ArchRProj = temp, colorBy = "cellColData", name = annotation, embedding = "UMAP") + theme(plot.title = element_text(size = 10), plot.margin = margin(0, 0, 0, 0, "pt"))
    
    p1 <- plotEmbedding(ArchRProj = temp, colorBy = "cellColData", name = sample, embedding = "UMAP") + theme(plot.title = element_text(size = 10), plot.margin = margin(0, 0, 0, 0, "pt"))
    p2 <- plotEmbedding(ArchRProj = temp, colorBy = "cellColData", name = annotation, embedding = "UMAP") + theme(plot.title = element_text(size = 10), plot.margin = margin(0, 0, 0, 0, "pt"))
    p3 <- plotEmbedding(ArchRProj = temp, colorBy = "cellColData", name = sample, embedding = paste("UMAP_Harmony", sample, sep = "_")) + theme(plot.title = element_text(size = 10), plot.margin = margin(0, 0, 0, 0, "pt"))
    p4 <- plotEmbedding(ArchRProj = temp, colorBy = "cellColData", name = annotation, embedding = paste("UMAP_Harmony", sample, sep = "_")) + theme(plot.title = element_text(size = 10), plot.margin = margin(0, 0, 0, 0, "pt"))
    plot_list <- list(p1, p2, p3, p4)
    return(plot_list)
}

