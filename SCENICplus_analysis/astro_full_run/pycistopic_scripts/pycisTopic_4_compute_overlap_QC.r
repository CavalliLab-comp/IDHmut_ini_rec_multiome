library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)

cancer = "astro"

base.dir = paste0("/mnt/beegfs/home/alyne")
setwd(base.dir)

#read in cell annotations csv 
#which has cells passing Seurat QC plus their annotation

cell.annot = fread(paste0("stagein/scenicplus_",cancer,"/cell_annotations.tsv.gz"), data.table=F, header=T)

#add column for cell name
cell.annot$cell.name = sapply(cell.annot$V1, function(x) {
  paste(unlist(strsplit(x,"-"))[1:2],collapse="-")
})

#read in list of cells which passed ATAC QC
cells.ATAC = fread(paste0("stageout/pycisTopic_output_",cancer,"/cells_passing_qc_ATAC.csv"), data.table=F, 
                   header=T, na.strings="NaN")

# #check, do all cells passing ATAC QC have annotations?
# # NO
# sapply(seq(ncol(cells.ATAC)), function(x) {
#   all(cells.ATAC[,x] %in% cell.annot$cell.name[which(cell.annot$sample==colnames(cells.ATAC)[x])])
# })

#make cell.annot with only cells passing both
cell.annot.list = lapply(seq(ncol(cells.ATAC)), function(x) {
  cell.annot[which(cell.annot$sample==colnames(cells.ATAC)[x]),]
})
cell.annot.filt.list = lapply(seq(ncol(cells.ATAC)), function(x) {
  cell.annot.rel = cell.annot[which(cell.annot$sample==colnames(cells.ATAC)[x]),]
  cell.annot.rel[which(cell.annot.rel$cell.name %in% cells.ATAC[,x]),]
})

#make table as input to ggplot of cell types
cell.types.all = unique(cell.annot$cell_type)
dfs.Seurat = do.call("rbind",lapply(seq(ncol(cells.ATAC)), function(x) {
  cell.type.counts = sapply(cell.types.all, function(y) {
    sum(cell.annot.list[[x]]$cell_type==y)
  })
  data.frame(sample = colnames(cells.ATAC)[x], cell.type = cell.types.all, count = cell.type.counts,
             filter = "after_Seurat")
}))
dfs.pycisTopic = do.call("rbind",lapply(seq(ncol(cells.ATAC)), function(x) {
  cell.type.counts = sapply(cell.types.all, function(y) {
    sum(cell.annot.filt.list[[x]]$cell_type==y)
  })
  data.frame(sample = colnames(cells.ATAC)[x], cell.type = cell.types.all, count = cell.type.counts,
             filter = "after_pycisTopic")
}))
dfs = rbind(dfs.Seurat, dfs.pycisTopic)
dfs$filter = factor(dfs$filter, levels = c("after_Seurat","after_pycisTopic"), ordered=T)

#make barplot showing cell types which pass both QC for each sample

ggplot(dfs) +
  geom_bar(aes(x = filter, y = count, fill = cell.type),
           position = "stack", stat = "identity") +
  facet_grid(~ sample, switch = "x")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("stageout/pycisTopic_output_",cancer,"/figures_scenicplus_",cancer,"/barplot_ct_filtering_",cancer,".pdf"), width = 20)






