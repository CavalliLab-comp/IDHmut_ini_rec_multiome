#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

# This script contains former methods used to perform differential/marker analysis using the
# ArchR pipeline. They are deprecated for reasons discussed in my report.
# Furthermore, they are VERY badly written as I did not fully understand ArchR's code at the
# time. 
# In contrary to ArchR doc, you *can* subset an ArchR project "by hand" (not using the ArchR 
# function creating a new ArchR project and copying the files on the hard drive) if you are only
# interested in the end results of your analysis (ex diff analysis).
# in this code, each time that i subset a project, I copy everything...
# My second error is that i used a very strict file architecture without thinking about the 
# reusability of the code... 
# Also, i did not care about the peak calling at the time and only used default params

library(ArchR)
library(dplyr)
CODE_PATH <- "/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas"
source(file.path(CODE_PATH, "Miscellaneous", "importantbutnowheretoputthem.r"))

# able to do differential analysis
singleDifferentialAnalysisATAC <- function(project, groupBy, Group1, Group2, threshold_FDR = 0.1, threshold_LogFC = 0.5, normBy =  NULL,  motif_database = "JASPAR2020") {
    print(paste("Computing", Group1, "VS", Group2, sep = " "))
    markerPeaks <- getMarkerFeatures(
        ArchRProj = project, useMatrix = "PeakMatrix", groupBy = groupBy,
        testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
        normBy = normBy,
        useGroups = Group1, bgdGroups = Group2, binarize = FALSE
    )

    markerGenes <- getMarkerFeatures(
        ArchRProj = project, useMatrix = "GeneScoreMatrix", groupBy = groupBy,
        testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
        normBy = normBy,
        useGroups = Group1, bgdGroups = Group2
    )

    threshold_up = paste('FDR <=', threshold_FDR, '& Log2FC >=', threshold_LogFC, sep = " ")

    motifsUp <- peakAnnoEnrichment(
        seMarker = markerPeaks,
        ArchRProj = project,
        peakAnnotation = "Motif",
        cutOff = threshold_up
    )

    threshold_down = paste('FDR <= ', threshold_FDR, ' & Log2FC <= -', threshold_LogFC, sep = "")

    motifsDown <- peakAnnoEnrichment(
        seMarker = markerPeaks,
        ArchRProj = project,
        peakAnnotation = "Motif",
        cutOff = threshold_down
    )

    return(list(Regions = markerPeaks, MotifsUp = motifsUp, MotifsDown = motifsDown, Genes = markerGenes))
}

# return a list of GRanges with all pairwise differential analysis
pairwiseDifferentialAnalysis <- function(project, groupBy, group_list = NULL,  minCells = 200, 
    minReplicates = 3, useLabels = FALSE, force_groups = TRUE,
    MACS2_PATH = findMacs2(), threshold_FDR = 0.05, threshold_LogFC = 1, normBy = NULL, motif_database = "JASPAR2020",
    library_genome = "BSgenome.Hsapiens.UCSC.hg38") {
    if (is.null(group_list)){
        group_list <- as.vector(unique(project@cellColData[, groupBy]))
    }

    if (!("GeneScoreMatrix" %in% getAvailableMatrices(project))){
        project <- addGeneScoreMatrix(project)
    }

    #group Coverages
    suppressMessages(library(library_genome, character.only = TRUE))

    project <- addGroupCoverages(
        ArchRProj = project, groupBy = groupBy, minCells = minCells,
        minReplicates = minReplicates, useLabels = useLabels, force = force_groups, maxReplicates = 10, maxCells = 3000
    )

    # peak calling
    project <- addReproduciblePeakSet(
        ArchRProj = project, groupBy = groupBy, pathToMacs2 = MACS2_PATH, reproducibility = "2", force = TRUE
    )
    project <- addPeakMatrix(project)

    project <- addMotifAnnotations(ArchRProj = project, motifSet = motif_database, annoName = "Motif", force = TRUE)

    pairs <- split(combn(group_list, 2), col(combn(group_list, 2)))
    names(pairs) <- sapply(pairs, function(x) paste(x, collapse = "vs"))

    project_name <- basename(getOutputDirectory())

    list_diffs <- list()

    for (pair in pairs) {
        diff_result <- tryCatch({
            singleDifferentialAnalysisATAC(project, groupBy, pair[1], pair[2], threshold_FDR, threshold_LogFC, normBy, motif_database)
        }, error = function(e) {
            warning(paste("Error in differential analysis for", pair[1], "vs", pair[2],"in", project_name, ":", e$message))
            warning(sys.calls())
            return(NULL)
        })

        if (!is.null(diff_result)) {
            list_diffs[[paste(pair, collapse = "vs")]] <- diff_result
        }

        # Cleanup = better memory management than lapply
        rm(diff_result)
        gc()
    }
    return(list_diffs)
}


fullDifferentialAnalysis <- function(project,
                                     annotation_differential,
                                     annotation_control = NULL,
                                     method = "subprojects",
                                     groups_in_diff = NULL, # should give a vector stating the labels I want to use : ex "astro_new" "oligo_new" and "undifferentiated" to have all paiwise combinations
                                     groups_in_control = NULL,
                                     minCells = 200,
                                     minReplicates = 3,
                                     useLabels = FALSE, # to use if I have multiple sample per condition
                                     force_subprojects = FALSE,
                                     force_group_coverage = TRUE,
                                     threshold_FDR = 0.05,
                                     threshold_LogFC = 1,
                                     normBy = NULL,
                                     motif_database = "JASPAR2020",
                                     library_genome = "BSgenome.Hsapiens.UCSC.hg38",
                                     MACS2_PATH = findMacs2(),
                                     datetime = strftime(Sys.time(), format = "%Y/%m/%d_%H:%M")) {
    message(rep("#", 100))
    message(paste("Starting fullDifferentialAnalysis of", annotation_differential, "vs", as.character(annotation_control), "using", method, sep = " "))
    message(rep("#", 100))

    if (!is.null(annotation_control)) {
        if (method == "subprojects") {
            # create subprojects
            subprojects <- subsetProject_usingAnnotation(project, annotation_control, anno_labels = groups_in_control, force = force_subprojects)
            annotation_differential_2 <- annotation_differential # I need to create new names so that i can both use a modified annotation in the analysis and have a nice title at the end
        } else if (method == "combined_annotation") { # should not use !!!!!!
            if (!is.null(groups_in_control)) {
                project_1 <- subsetProject_usingAnnotation(project, annotation_control, anno_labels = list(groups_in_control), force = force_subprojects)
            } else {
                project_1 <- list(project)
            }
            # create combined annotation
            subprojects <- list(combined_annotation = addCombinedAnnotation(project_1[[1]], annotation_differential, annotation_control))
            annotation_differential_2 <- paste(annotation_differential, "by", annotation_control, sep = "_")
        }
        annotation_control <- paste("by", annotation_control, sep = "_")
    } else {
        # format it so that it can be applied to list operations
        subprojects <- list(simple_analysis = project)
        method <- "no control"
        annotation_differential_2 <- annotation_differential
        annotation_control <- ""
    }
    # apply single analysis on them
    # this loop allow better control of the ram = delete subproj as we process and gc the intermediary results
    enrich_mat_list <- list()
    for (name in names(subprojects)) {
        subproj <- subprojects[[name]]

        message(rep("#", 100))
        message(paste("Starting Pairwise Differential analysis for"), name, sep = " ")
        message(rep("#", 100))

        if (as.character(minCells) == "Auto") {
            subproj_df <- as.data.frame(subproj@cellColData) %>%
                group_by(.data[[annotation_differential]]) %>%
                summarize(n = n())

            if (all(subproj_df$n > 200 * minReplicates)) {
                effective_minCells <- min(floor(subproj_df$n / minReplicates))
            } else {
                effective_minCells <- 200
            }
        } else {
            effective_minCells <- minCells
        }

        # Process one subproject at a time
        result <- pairwiseDifferentialAnalysis(subproj, annotation_differential,
            group_list = groups_in_diff, minCells = effective_minCells,
            minReplicates = minReplicates, useLabels = useLabels, force_groups = force_group_coverage,
            MACS2_PATH = MACS2_PATH, threshold_FDR = threshold_FDR, threshold_LogFC = threshold_LogFC, normBy = normBy, 
            motif_database = motif_database,
            library_genome = library_genome
        )

        # Store result
        enrich_mat_list[[name]] <- result
        # Remove the processed subproject from memory
        subprojects[[name]] <- NULL
        rm(subproj, result)
        gc() # Trigger garbage collection to free memory
    }
    return(enrich_mat_list)
}

# not used
# count the number of cells for 1 or 2 combined given annotations
countCellperGroup <- function(project, annotation1, annotation2) {
    project <- addCombinedAnnotation(project, annotation1, annotation2)
    return(
        table(
            project@cellColData[[paste(annotation1, "by", annotation2, sep = "_")]]
        )
    )
}


# marker analysis for a given annotation and a given project
singleMarkerAnalysis <- function(
    project, groupBy, minCells = 200, minReplicates = 3, force_groups = TRUE,
    threshold = "FDR <= 0.05 & Log2FC >= 1", motif_database = "JASPAR2020",
    library_genome = "BSgenome.Hsapiens.UCSC.hg38", MACS2_PATH = findMacs2()) {
    # creating pseudoreplicates
    suppressMessages(library(library_genome, character.only = TRUE))

    project <- addGroupCoverages(
        ArchRProj = project, groupBy = groupBy, minCells = minCells,
        minReplicates = minReplicates, force = force_groups, maxReplicates = 10, maxCells = 3000
    )

    # peak calling
    project <- addReproduciblePeakSet(
        ArchRProj = project, groupBy = groupBy, pathToMacs2 = MACS2_PATH, reproducibility = "max(ceiling((n+1)/3),2)", force = TRUE
    )
    project <- addPeakMatrix(project)

    # get marker peaks
    markerPeaks <- getMarkerFeatures(
        ArchRProj = project, useMatrix = "PeakMatrix", groupBy = groupBy,
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    #gene score markers 
    if (!("GeneScoreMatrix" %in% getAvailableMatrices(project))){
        project <- addGeneScoreMatrix(project)
    }

    markerGenes <- getMarkerFeatures(
        ArchRProj = project, useMatrix = "GeneScoreMatrix", groupBy = groupBy,
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    # add motif
    project <- addMotifAnnotations(ArchRProj = project, motifSet = motif_database, annoName = "Motif", force = TRUE)

    saveArchRProject(
        ArchRProj = project,
        outputDirectory = getOutputDirectory(project),
        overwrite = TRUE,
        load = FALSE
    )
    tryCatch(
        {
            enrichMotifs <- peakAnnoEnrichment(
                seMarker = markerPeaks, ArchRProj = project, peakAnnotation = "Motif",
                cutOff = threshold
            )
        },
        error = function(e) {
            enrichMotifs <- e$message
        }
    )

    result <- list(Regions = markerPeaks, Motifs = enrichMotifs, Genes = markerGenes)

    return(result)
}

# run a full marker analysis on regions and motifs for a given annotation while taking into account another annotation;
# possible to specify groups in annotation control to be included in the analysis
# in this case : two methods : subproject creation or use of a combined annotation
fullMarkerAnalysis <- function(project,
                               annotation_differential,
                               annotation_control = NULL,
                               method = "subprojects", # also available : "combined_annotation"
                               groups_in_diff = NULL, # only available if annotation_control = NULL
                               groups_in_control = NULL,
                               minCells = 200,
                               minReplicates = 3,
                               force_subprojects = FALSE,
                               force_group_coverage = TRUE,
                               threshold = "FDR <= 0.05 & Log2FC >= 1",
                               max_motifs = 10,
                               motif_database = "JASPAR2020",
                               library_genome = "BSgenome.Hsapiens.UCSC.hg38",
                               MACS2_PATH = findMacs2(),
                               datetime = strftime(Sys.time(), format = "%Y/%m/%d_%H:%M")) {
    if (!is.null(annotation_control)) {
        if (method == "subprojects") {
            # create subprojects
            subprojects <- subsetProject_usingAnnotation(project, annotation_control, anno_labels = groups_in_control, force = force_subprojects)
            annotation_differential_2 <- annotation_differential # I need to create new names so that i can both use a modified annotation in the analysis and have a nice title at the end
        } else if (method == "combined_annotation") {
            if (!is.null(groups_in_control)) {
                project_1 <- subsetProject_usingAnnotation(project, annotation_control, anno_labels = list(groups_in_control), force = force_subprojects)
            } else {
                project_1 <- list(project)
            }
            # create combined annotation
            subprojects <- list(combined_annotation = addCombinedAnnotation(project_1[[1]], annotation_differential, annotation_control))
            annotation_differential_2 <- paste(annotation_differential, "by", annotation_control, sep = "_")
        }
        annotation_control <- paste("by", annotation_control, sep = "_")
    } else {
        # format it so that it can be applied to list operations
        subprojects <- list(simple_analysis = project)
        method <- "no control"
        annotation_differential_2 <- annotation_differential
        annotation_control <- ""
    }
    # apply single analysis on them
    enrich_mat_list <- lapply(subprojects, function(subproj) singleMarkerAnalysis(subproj, annotation_differential_2, minCells = minCells, minReplicates = minReplicates, force_groups = force_group_coverage, threshold = threshold, motif_database = motif_database, library_genome = library_genome, MACS2_PATH = MACS2_PATH))
    names(enrich_mat_list) <- names(subprojects)

    # plot
    plots_regions <- lapply(seq_along(enrich_mat_list), function(i) {
        tryCatch(
            {
                heatmap <- plotMarkerHeatmap(seMarker = enrich_mat_list[[i]][["Regions"]], cutOff = threshold, transpose = TRUE, binaryClusterRows = TRUE, clusterCols = TRUE, plotLog2FC = (ncol(enrich_mat_list[[i]][["Regions"]]) <= 2))
                draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot", column_title = paste("Enriched regions in", names(enrich_mat_list)[[i]]))
            },
            error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Error for", names(enrich_mat_list)[[i]], ":", e$message, sep = " "), col = "red", cex = 1.5)
            }
        )
    })

    plots_genes <- lapply(seq_along(enrich_mat_list), function(i) {
        tryCatch(
            {
                heatmap <- plotMarkerHeatmap(seMarker = enrich_mat_list[[i]][["Genes"]], cutOff = threshold, transpose = TRUE, binaryClusterRows = TRUE, clusterCols = TRUE, plotLog2FC = (ncol(enrich_mat_list[[i]][["Genes"]]) <= 2))
                draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot", column_title = paste("Enriched regions in", names(enrich_mat_list)[[i]]))
            },
            error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Error for", names(enrich_mat_list)[[i]], ":", e$message, sep = " "), col = "red", cex = 1.5)
            }
        )
    })

    plots_motifs <- lapply(seq_along(enrich_mat_list), function(i) {
        tryCatch(
            {
                heatmap <- plotEnrichHeatmap(enrich_mat_list[[i]][["Motifs"]], n = max_motifs, transpose = TRUE, returnMatrix = FALSE)
                draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot", column_title = paste("Enriched motifs in", names(enrich_mat_list)[[i]]))
            },
            error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Error for", names(enrich_mat_list)[[i]], ":", e$message, sep = " "), col = "red", cex = 1.5)
            }
        )
    })
    plotPDF(plotList = plots_regions, ArchRProj = project, name = paste(datetime, "Enriched_regions_per", annotation_differential, annotation_control, "using", method, sep = "_"))
    plotPDF(plotList = plots_genes, ArchRProj = project, name = paste(datetime, "Enriched_genes_accessibility_per", annotation_differential, annotation_control, "using", method, sep = "_"))
    plotPDF(plotList = plots_motifs, ArchRProj = project, name = paste(datetime, "Enriched_motifs_per", annotation_differential, annotation_control, "using", method, sep = "_"))
    return(enrich_mat_list)
}
