
library(ArchR)
library(dplyr)

# These scripts are necessary to run some outdated analysis...
CODE_PATH <- "/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas"
source(file.path(CODE_PATH, "Miscellaneous", "importantbutnowheretoputthem.r"))

#specific function to filter by a list of cell names, badly written but yolo
subsetProject_byCellVec <- function(project, name_sample, cells_formatted, outputDirectory_individual, logDirectory_individual){
    idxSample <- BiocGenerics::which(project$cellNames %in% cells_formatted)
                             
    lapply(c(outputDirectory_individual, logDirectory_individual), check_a_path)
    outputDirectory_subproject <- file.path(outputDirectory_individual, name_sample)
    logDirectory_subproject <- file.path(logDirectory_individual, name_sample)
    lapply(c(outputDirectory_subproject), check_a_path)
                       
    projSubset <- subsetArchRProject(
      ArchRProj = project,
      cells = project$cellNames[idxSample],
      outputDirectory = outputDirectory_subproject,
      dropCells = TRUE,
      logFile = logDirectory_subproject,
      force = TRUE
    )
    saveArchRProject(
          ArchRProj = projSubset,
          outputDirectory = getOutputDirectory(projSubset),
          overwrite = TRUE,
          load = FALSE)
    return(projSubset)
}


# from a project, a given annotation in cellColData and a label (or a vector of label), will return a project for the cells with one of the specified labels stored in path_subproject_annotation/label
subsetProject_byAnno_labelvector <- function(project, annotation, label, path_subproject_annotation, force = FALSE) {
    if (!annotation %in% colnames(project@cellColData)) {
        stop("Annotation column not found in project@cellColData")
    }
    if ("UMAP" %in% names(project@embeddings)) {
    idxSample <- BiocGenerics::which(Reduce("&", list(
        unlist(project@cellColData[annotation]) %in% label,
        unlist(project$cellNames) %in% rownames(project@embeddings$UMAP$params$projectID) # I needed to add that because of a subscript error, some cells were listed in cellColData while being absent in UMAP rows...
    )))
    } else {
        idxSample <- BiocGenerics::which(project@cellColData[annotation] %in% label)
    }
    if (any(idxSample > length(project$cellNames))) {
        stop("Index out of bounds in project$cellNames")
    }

    path_sub_anno_label = file.path(path_subproject_annotation, paste(label, collapse = "_"))
    path_sub_anno_label_proj = file.path(path_sub_anno_label, "project")
    path_sub_anno_label_logs = file.path(path_sub_anno_label, "logs")
    exist = lapply(c(path_sub_anno_label,path_sub_anno_label_proj), check_a_path) #check a path returns true if the file exists or return false and creates it

    if ((!(exist[[2]]) || force)) {
        print(paste("Subsetting",paste(label, collapse = "_"), "as subproject do not exist or force = TRUE",sep = " "))
        projSubset <- subsetArchRProject(
          ArchRProj = project,
          cells = project$cellNames[idxSample],
          outputDirectory = path_sub_anno_label_proj,
          dropCells = TRUE,
          logFile = path_sub_anno_label_logs,
          force = TRUE
        )
        saveArchRProject(
          ArchRProj = projSubset,
          outputDirectory = getOutputDirectory(projSubset),
          overwrite = TRUE,
          load = FALSE)
    } else {
        projSubset <- loadArchRProject(path = path_sub_anno_label_proj)
    }
    return(projSubset)
}  

# from a project and a given annotation, returns all possible subprojects, possible to specify a list of combination of wanted labels
# the default comportement is to o subset for each single label and not overwrite existing subprojects 
subsetProject_usingAnnotation <- function(project, annotation, anno_labels = NULL, force = FALSE) {
    path_out_large = getOutputDirectory(project)

    path_sub = file.path(dirname(dirname(path_out_large)), paste(basename(path_out_large), "subprojects", sep = "_"))
    path_subproject_annotation = file.path(path_sub, annotation)
    wedonotneedtheoutput = lapply(c(path_sub,path_subproject_annotation), check_a_path)

    if (is.null(anno_labels)) { 
        anno_labels <- unique(project@cellColData[, annotation])
    }

    subprojects<- lapply(anno_labels, function(x) subsetProject_byAnno_labelvector(project, annotation, x, path_subproject_annotation, force))
    names(subprojects) = anno_labels
    return(subprojects)
}