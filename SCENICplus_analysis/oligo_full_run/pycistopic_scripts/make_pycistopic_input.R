rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
script.dir = as.character(args[3])
malig.types = as.character(args[4])
cancer = as.character(args[5])

cancer = "oligo"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
samp.file = paste0(data.dir,"/sample_list.txt")
malig.types = "score"


library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
library(patchwork)
library(viridis)


setwd(data.dir)

#get list of all samples
samps = read.table(samp.file, header = F, stringsAsFactors = F, colClasses = c("character"))[,1]

oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_all_",cancer,"_doub_rem.rds"))

md.col.list = which(colnames(oligo.ss@meta.data)==paste0("ICA_",malig.types,"_states"))
cell.states = oligo.ss$cell.type.sr.cl.imm.postmerge
inds.malig = which(oligo.ss$iCNV.malig.ident=="malignant")
cell.states[inds.malig] = oligo.ss@meta.data[inds.malig,md.col.list]

cell.samples = sapply(names(oligo.ss$orig.ident), function(x) {unlist(strsplit(x,"-"))[3]})
cell.names = names(oligo.ss$orig.ident)

#make dataframe with relevant meta data
df.all = data.frame(cell_type = cell.states, sample = cell.samples)
rownames(df.all) = cell.names

#filter to take only cell types for which we have a decent amount of cells
df.all = df.all[which(df.all$cell_type %in% c("Macrophage","GABAergic","Glutamatergic","Oligo","Astro",
                                              "OPC-like","astro-like","intermediate","neuro-like")),]

fwrite(df.all, file = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer,"/cell_annotations.tsv.gz"), 
       sep = "\t", col.names = T, row.names = T)












# pal = viridis(n=50)
# 
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# # setwd("/Users/alyne/scenicplus_test_data")
# data.dir = "/data/users/alyne/gliomult/svn/analyse/results"
# cancer = "astro"
# setwd(data.dir)
# samp.file = "/data/users/alyne/scrnaseq/oligo_10x_multiome/oligo_astro_samples_all.txt"
# 
# # samps = c("astro_P0015_ini_S0028","astro_P0015_rec1_S0029",
# #           'astro_P0016_ini_S0030', 'astro_P0016_rec1_S0031')
# 
# sample.info = read.table(samp.file, header = F, colClasses = c("character","character","numeric","numeric"), 
#                          stringsAsFactors = F)
# sample.info = sample.info[order(sample.info$V2),]
# samples.newnames = sample.info[,2]
# 
# if(cancer!="all"){
#   samples.newnames = sort(samples.newnames[grepl(cancer,samples.newnames,fixed=T)])
# }else{
#   samples.newnames = sort(samples.newnames)
# }
# 
# # remove samples with no malignant cells
# samples.newnames = samples.newnames[!(samples.newnames %in% c("oligo_P0003_rec2_S0007","astro_P0014_ini_S0026",
#                                                               "oligo_P0022_rec1_S0043","oligo_P0030_ini_S0058"))]
# 
# #also remove samples which don't have both RNA and ATAC
# samples.newnames = samples.newnames[!(samples.newnames %in% c("oligo_P0011_ini_S0022","astro_P0017_ini_S0032",
#                                                               "astro_P0018_rec1_S0035","astro_P0025_ini_S0048"))]
# 
# #read in integrated object with cell names, samples, integrated clusters and cell types
# oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_all_",cancer,".rds"))
# 
# 
# #keep only cells in samples with malignant cells and good quality ATAC
# cell.samples = sapply(names(oligo.ss$orig.ident), function(x) {
#   unlist(strsplit(x,"-"))[3]
# })
# oligo.ss$sample = cell.samples
# cells.to.keep = names(oligo.ss$orig.ident)[which(cell.samples %in% samples.newnames)]
# oligo.ss = subset(oligo.ss, cells = cells.to.keep)
# 
# #in both cases remove cells which are mac/mal doublet or macs which have just ingested mal
# cells.remove = names(oligo.ss$orig.ident)[which(oligo.ss$rpca_clusters==13)]
# oligo.ss = subset(oligo.ss, cells = names(oligo.ss$orig.ident)[which(!(names(oligo.ss$orig.ident) %in% cells.remove))])
# 
# 
# #redo clustering to reduce number of malignant clusters
# oligo.ss <- FindClusters(oligo.ss, resolution = 0.25, cluster.name = "rpca_clusters_0_25")
# 
# cl.to.use = "rpca_clusters_0_25"
# md.col = which(colnames(oligo.ss@meta.data)==cl.to.use)
# 
# # points = Embeddings(oligo.ss, reduction = "umap.har")[,1:2]
# # mean(points[oligo.ss$rpca_clusters==19,1])
# 
# if(cancer=="astro"){
#   #add back in cluster on non-malignant astrocytes
#   levels(oligo.ss@meta.data[,md.col]) = c(levels(oligo.ss@meta.data[,md.col]),
#                                           max(as.numeric(levels(oligo.ss@meta.data[,md.col]))) + 1)
#   oligo.ss@meta.data[which(oligo.ss$rpca_clusters==19),md.col] = as.character(max(as.numeric(levels(oligo.ss@meta.data[,md.col]))))
# }
# 
# 
# #assign a cell type to each cluster
# cl.uni = unique(oligo.ss@meta.data[,md.col])
# cl.cell.types = sapply(cl.uni, function(x) {
#   names(sort(table(oligo.ss$cell.type.sr.imm[which(oligo.ss@meta.data[,md.col]==x)]),decreasing = T))[1]
# })
# 
# #assign clusters as malignant or not
# cl.malig.types = sapply(cl.uni, function(x) {
#   names(sort(table(oligo.ss$iCNV.malig.ident[which(oligo.ss@meta.data[,md.col]==x)]),decreasing = T))[1]
# })
# 
# rpca_cell_types = sapply(oligo.ss@meta.data[,md.col], function(x) {
#   cl.cell.types[which(cl.uni==x)]
# })
# rpca_malig_types = sapply(oligo.ss@meta.data[,md.col], function(x) {
#   cl.malig.types[which(cl.uni==x)]
# })
# oligo.ss$rpca_cell_types = rpca_cell_types
# oligo.ss$rpca_malig_types = rpca_malig_types
# oligo.ss$rpca_malig_cell_types = paste(oligo.ss$rpca_malig_types, oligo.ss$rpca_cell_types,
#                                        oligo.ss@meta.data[,md.col], sep="_")
# 
# 
# if(cancer=="oligo"){
#   df = data.frame(umaphar_1 = Embeddings(oligo.ss,"umap.har")[,1],
#                   umaphar_2 = Embeddings(oligo.ss,"umap.har")[,2], SNV = oligo.ss$IDH.mut,
#                   p1.ATAC = oligo.ss$p1.av.ATAC, q19.ATAC = oligo.ss$q19.av.ATAC)
# }else{
#   df = data.frame(umaphar_1 = Embeddings(oligo.ss,"umap.har")[,1],
#                   umaphar_2 = Embeddings(oligo.ss,"umap.har")[,2], SNV = oligo.ss$IDH.mut,
#                   q7.ATAC = oligo.ss$q7.av.ATAC)
# }
# 
# levs = unique(df$SNV)[!is.na(unique(df$SNV))]
# cols = gg_color_hue(length(levs))
# df$color = sapply(df$SNV, function(x) {
#   if(is.na(x)){return("gray70")}else{cols[which(levs==x)]}
# })
# #sort to put all NAs first, i.e., to plot all the other things last
# df = df[order(df$SNV, na.last = F),]
# 
# 
# #plot UMAP to check everything working correctly
# pt.sz = 0.002
# p1 = DimPlot(oligo.ss, group.by = "orig.ident", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p2 = DimPlot(oligo.ss, group.by = "nCNA", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p3 = DimPlot(oligo.ss, group.by = cl.to.use, reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p4 = DimPlot(oligo.ss, group.by = "rpca_cell_types", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p5 = DimPlot(oligo.ss, group.by = "rpca_malig_types", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p6 = DimPlot(oligo.ss, group.by = "rpca_malig_cell_types", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# p7 = DimPlot(oligo.ss, group.by = "rpca_clusters", reduction = 'umap.har', pt.size = pt.sz, shuffle=T, raster=F)
# if(cancer=="astro"){
#   p8 = FeaturePlot(oligo.ss, features = "q7.av.ATAC", reduction = 'umap.har', pt.size = pt.sz, raster=F) +
#     scale_colour_gradientn(colours = pal, na.value = "grey90")
# }else{
#   p8 = FeaturePlot(oligo.ss, features = "p1.av.ATAC", reduction = 'umap.har', pt.size = pt.sz, raster=F) +
#     scale_colour_gradientn(colours = pal, na.value = "grey90")
# }
# p9 = ggplot(df, aes(x = umaphar_1, y = umaphar_2, color = SNV)) + geom_point()
# print(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3))
# ggsave(paste0(data.dir,"/figures_scenicplus/umap_all_cells_",cancer,"_pre_cistopic.pdf"), width=20, height=12)
# 
# 
# #make dataframe with relevant meta data
# df.all = data.frame(cell_type = oligo.ss$rpca_malig_cell_types, sample = oligo.ss$sample)
# rownames(df.all) = names(oligo.ss$orig.ident)
# 
# 
# fwrite(df.all, file = paste0("scenicplus_",cancer,"/input_data/cell_annotations.tsv.gz"), 
#        sep = "\t", col.names = T, row.names = T)
# 
# #also write sample file
# write.table(samples.newnames, paste0("scenicplus_",cancer,"/input_data/sample_list.txt"),
#             col.names = F, row.names = F, quote = F)
# 



# df.list = list()
# for(i in 1:length(samps)){
#   
#   print(i)
#   
#   oligo = readRDS(paste0("analysis_",samps[i],"/oligo_filt_final_",samps[i],".rds"))
#   
#   cell.names = paste0(names(oligo$orig.ident),"-",oligo$orig.ident[1])
#   
#   cell.type = oligo$cell.type.sr.cl.imm.postmerge
#   cell.type[!is.na(oligo$traj.cell.type.new)] = oligo$traj.cell.type.new[!is.na(oligo$traj.cell.type.new)]
#   
#   df = data.frame(cell_type = cell.type, sample = oligo$orig.ident)
#   rownames(df) = cell.names
#   df.list[[i]] = df
#   
# }
# 
# df.all = do.call("rbind",df.list)
# 
# fwrite(df.all, file = "cell_annotations.tsv.gz", sep = "\t", col.names = T, row.names = T)
# 
# #also write sample file
# write.table(samps, paste0("scenic_plus_output/sample_list.txt"),
#             col.names = F, row.names = F, quote = F)
