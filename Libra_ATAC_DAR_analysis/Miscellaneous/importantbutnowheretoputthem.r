#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(SingleCellExperiment)

# Two little functions that i use a lot 

#check if a dir exists; if not creates it
check_a_path <- function(single_path){
    if (! file.exists(single_path)) {
        dir.create(single_path)
        return(FALSE)
    } else {
       return(TRUE)
    }
}

# create a new column in cellColData combining two given column
# naming it {A}_by_{B}
addCombinedAnnotation <- function(project, annotation1, annotation2){
    combined_anno = paste(annotation1, "by", annotation2, sep = "_")
    project@cellColData[[combined_anno]] <- mapply(paste, 
                                                     as.vector(project@cellColData[[annotation1]]), 
                                                     as.vector(project@cellColData[[annotation2]]), 
                                                     sep = "_")
    unname(project@cellColData[[combined_anno]])
    return(project)
}
