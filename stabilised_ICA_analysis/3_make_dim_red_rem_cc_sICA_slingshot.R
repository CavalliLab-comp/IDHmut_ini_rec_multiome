##script to project cells into ICA space

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
data.dir = as.character(args[1])
samp.file = as.character(args[2])
cancer = as.character(args[3])

setwd(data.dir)

library(data.table)
library(gplots)
library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(pracma)
library(igraph)
library(RColorBrewer)
library(gage)
library(annotate)
library(clusterProfiler)
library(fgsea)
library(patchwork)

if(cancer=="astro"){
  cor.thresh = 0.45
}else{
  cor.thresh = 0.5
}
#cluster resolution for ICA_cell_states
res = 0.18


#directory with ICA input
fig.dir.ica = paste0("figures_sICA_",cancer,"_",cor.thresh)

#directory for STREAM output
fig.dir = paste0("figures_STREAM/",cancer,"_ss_rem_cc")
dir.create(fig.dir, recursive=T)


#ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#set to "malig" if want to run on malignant cells
fig.str = "malig"

if(cancer=="astro"){
  oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
}else{
  oligo.ss = readRDS(paste0("oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
}

oligo = SplitObject(oligo.ss, split.by = "orig.ident")
samps = sapply(oligo, function(x) unlist(strsplit(names(x$orig.ident)[1],"-"))[3])


#load in average of each dimension
av = "med"
df.mg.av = fread(paste0(data.dir,"/",fig.dir.ica,"/all_ICA_components_",av,"_",cancer,".tsv.gz"), data.table=F)
rownames(df.mg.av) = df.mg.av[,1]
df.mg.av = df.mg.av[,2:ncol(df.mg.av)]
nICAdims = ncol(df.mg.av)


#keep cell cycle dimensions for now, but they will be removed before running slingshot
if(cancer=="oligo"){
  #dim 8 is cell cycle
  cc.dim = c(5)
}else if(cancer=="astro"){
  #dim 2 is cell cycle
  cc.dim = c(2)
}else{
  cc.dim = c(3,8)
}


unmixing.mat = pinv(as.matrix(df.mg.av))
colnames(unmixing.mat) = rownames(df.mg.av)

fwrite(as.data.frame(unmixing.mat),
       file = paste0(fig.dir,"/unmixing_mat_rem_cc_",cancer,".tsv.gz"),
       sep = "\t", append = F, col.names = T, quote = F)


unmixing.mat = as.matrix(fread(paste0(fig.dir,"/unmixing_mat_rem_cc_",cancer,".tsv.gz"),
                               data.table = F))


#iterate through samples, select genes in intersect
#centre at 0, then project into new space


#list of lists
#one list per sample, and then list for each md within
data.proj = list()
for(i in 1:length(samps)){
  print(oligo[[i]]$orig.ident[1])
  
  data.filt = oligo[[i]]@assays$SCT@data[sapply(rownames(oligo[[i]]@assays$SCT@data),
                                                function(x) { x %in% rownames(df.mg.av) }),]
  #centre
  data.filt2 = t(apply(as.matrix(data.filt), 1, function(x) { x - mean(x) }))

  #project
  data.proj[[i]] = unmixing.mat %*% data.filt2

  #if i=1, write list of genes in intersect for universe in overenrichment analysis
  if(i==1){
    write.table(rownames(df.mg.av), row.names = F,col.names = F, quote=F,
                file = paste0(fig.dir,"/",cancer,"_samps_gene_universe.txt"))
  }
}



#all cells for both cancers
j = "all"
dir.create(paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j))

Ncell.samp = sapply(oligo, function(x) {
  length(x$orig.ident)
})
ss = lapply(seq(length(Ncell.samp)), function(x) {
  seq(Ncell.samp[x])
})

data.ss = do.call(cbind,lapply(seq(length(ss)), function(x) {
  data.proj[[x]][,ss[[x]]]
}))


#add row names for each dimension
rownames(data.ss) = colnames(df.mg.av)
#make input file for STREAM
fwrite(as.data.frame(data.ss[-cc.dim,]),
       file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,"/",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
       sep = "\t", row.names = T, col.names = T, quote = F)


#add ICA embedding and ica UMAP to full object
#this has full components including cell cycle
ica_assay <- CreateAssay5Object(data = as.matrix(data.ss))
oligo.ss[["ICA"]] <- ica_assay

#also put as DR so can directly use for clustering etc
ica.mat <- t(as.matrix(data.ss[-cc.dim,]))  # cells x components
#rename columns temporarily
colnames(ica.mat) = paste0("ICA_",1:ncol(ica.mat))
ica.dr <- CreateDimReducObject(
  embeddings = ica.mat, key = "ICA_"
)
oligo.ss[["ica"]] <- ica.dr

u = uwot::umap(ica.mat)
umap.dr <- CreateDimReducObject(embeddings = as.matrix(u[,1:2]), key = "UMAPICA_")
oligo.ss[["umap_ica"]] <- umap.dr


oligo.ss <- FindNeighbors(oligo.ss, dims = 1:ncol(ica.mat), reduction = 'ica', 
                          assay = "RNA", graph.name = c("ica_nn"))
oligo.ss <- FindClusters(oligo.ss, resolution = res, graph = "ica_nn", cluster.name = paste0("ica_snn_res.",res))
#find relevant meta data column
md.col.res = which(colnames(oligo.ss@meta.data)==paste0("ica_snn_res.",res))


#assign cells from tiny clusters to nearest
if(cancer=="oligo"){
  cells.cl5 = names(oligo.ss$orig.ident)[which(oligo.ss@meta.data[,md.col.res] %in% c(4,5,6))]
}else{
  cells.cl5 = names(oligo.ss$orig.ident)[which(oligo.ss@meta.data[,md.col.res] %in% c(4,5))]
}
knn = oligo.ss@graphs$ica_nn
#for each cell, find cluster of nearest neighbour and assign most frequent
new.cls = sapply(cells.cl5, function(x) {
  neighbs = colnames(knn)[knn[x, ] > 0]
  cl.neigh = oligo.ss@meta.data[which(rownames(oligo.ss@meta.data) %in% neighbs),md.col.res]
  names(sort(table(cl.neigh), decreasing = T))[1]
})
for(i in 1:length(cells.cl5)){
  oligo.ss@meta.data[which(names(oligo.ss$orig.ident)==cells.cl5[i]),md.col.res] = new.cls[i]
}
oligo.ss@meta.data[,md.col.res] = droplevels(oligo.ss@meta.data[,md.col.res])

p1 = DimPlot(oligo.ss, group.by = paste0("ica_snn_res.",res), reduction = "umap_ica", raster=F, pt.size = 0.001, shuffle=T)
p2 = DimPlot(oligo.ss, group.by = "ICA_score_states", reduction = "umap_ica", raster=F, pt.size = 0.001, shuffle=T)
print(p1 + p2)
ggsave(paste0("figures/umap_ica_projection_clusters_",cancer,".pdf"), width = 12)


p1 = FeaturePlot(oligo.ss, features = "OPCML", reduction = "umap_ica", raster=F, pt.size = 0.001)
p2 = FeaturePlot(oligo.ss, features = "KCNN3", reduction = "umap_ica", raster=F, pt.size = 0.001)
print(p1 + p2 + plot_layout(ncol = 2))
ggsave(paste0("figures/umap_ica_projection_features_",cancer,".pdf"), width = 12)



markers = FindAllMarkers(oligo.ss, group.by = paste0("ica_snn_res.",res), assay="SCT",slot="data")

#convert clusters to ICA cell states
oligo.ss$ICA_cell_states = sapply(oligo.ss@meta.data[,md.col.res], function(x) {
  if(cancer=="oligo"){
    if(x==0){
      return("OPC-like")
    }else if(x==1){
      return("intermediate")
    }else if(x==2){
      return("astro-like")
    }else if(x==3){
      return("neuro-like")
    }
  }else if(cancer=="astro"){
    if(x==0){
      return("intermediate")
    }else if(x==1){
      return("astro-like")
    }else if(x==2){
      return("OPC-like")
    }else if(x==3){
      return("neuro-like")
    }
  }
})

p1 = DimPlot(oligo.ss, group.by = "ICA_cell_states", reduction = "umap_ica", raster=F, pt.size = 0.001, shuffle=T)
p2 = DimPlot(oligo.ss, group.by = "ICA_score_states", reduction = "umap_ica", raster=F, pt.size = 0.001, shuffle=T)
print(p1 + p2)
ggsave(paste0("figures/umap_ica_projection_clusters_",cancer,".pdf"), width = 12)

if(cancer=="astro"){
  saveRDS(oligo.ss, paste0("oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))
}else{
  saveRDS(oligo.ss, paste0("oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))
}

#add ICA_cell_states to df so can be added to separate objects
df = fread(paste0("df_ICA_score_states_",cancer,".tsv.gz"), data.table = F, header = T)
rownames(df) = df[,1]
df[,1] = NULL
df$ICA_cell_states = oligo.ss$ICA_cell_states
fwrite(df, paste0("df_ICA_score_states_",cancer,".tsv.gz"), col.names = T, row.names = T, quote = F)


#also get relevant metadata 
if(cancer=="astro"){
  md.names = c("cyc.status","CytoTRACE","icnv.score","ICA_cell_states","ICA_score_states",
               "atlas.label","gest_week","cell.type.sr.imm","ini_rec","grade","treatment")
}else{
  md.names = c("cyc.status","CytoTRACE","icnv.score","atlas.label","ICA_cell_states","ICA_score_states",
               "gest_week","cell.type.sr.imm","ini_rec","grade","treatment","p1.av","p1.av.ATAC","q19.av.ATAC")
}
df = data.frame(oligo.ss@meta.data[,which(colnames(oligo.ss@meta.data) %in% md.names)])
df$cell = rownames(df)
df$sample = sapply(df$cell, function(x) {unlist(strsplit(x,"-"))[3]})
df = cbind(df, Embeddings(oligo.ss,"umap_ica"))
colnames(df)[(ncol(df)-1):ncol(df)] = c("UMAP1","UMAP2")
fwrite(df, file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,"/metadata_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
       sep = "\t", col.names=T, row.names=F, quote = F)



#make separate metadata file with ICA scores
md.cols = which(grepl("^ICA_.+_1$",colnames(oligo.ss@meta.data)))
fwrite(as.data.frame(oligo.ss@meta.data[,md.cols]), file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
                                                              "/metadata_ICA_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
       sep = "\t", col.names=T, row.names=T, quote = F)


#write full expression data
data.full.ss = list()
for(i in 1:length(samps)){
  data.filt = oligo[[i]]@assays$SCT@data[sapply(rownames(oligo[[i]]@assays$SCT@data),
                                                function(x) {x %in% rownames(df.mg.av)}),]
  data.full.ss[[i]] = data.filt[,ss[[i]]]
}


data.full.ss = do.call(cbind, data.full.ss)
#detect variable genes
oligo.temp = CreateSeuratObject(data.full.ss)
oligo.temp = FindVariableFeatures(oligo.temp, selection.method = "vst", nfeatures = 5000)
to_print = data.full.ss[rownames(data.full.ss) %in% VariableFeatures(oligo.temp),]

fwrite(as.data.frame(as.matrix(to_print)),
       file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,"/",cancer,"_samps_ss_full_data_rem_cc_",j,".tsv.gz"),
       sep = "\t", row.names = T, col.names = T, quote = F)







#####
##old code
#####

# #subsample cells and plot
# Nsamp = 20000
# Nrep = 1
# Ncell.samp = sapply(oligo, function(x) {
#   length(x$orig.ident)
# })
# N.per.samp = rmultinom(1,Nsamp,Ncell.samp)
# 
# 
# for(j in 1:Nrep){
#   dir.create(paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j))
# 
#   ss = lapply(seq(length(oligo)), function(x) {
#     sort(sample(Ncell.samp[x], N.per.samp[x,], replace=F))
#   })
# 
#   data.ss = do.call(cbind,lapply(seq(length(oligo)), function(x) {
#     data.proj[[x]][,ss[[x]]]
#   }))
#   
#   #make metadata list, one element per md
#   md.list.annot = list()
#   for(i in 1:length(md.names)){
#     md.list.annot[[i]] = unlist(lapply(seq(length(oligo)), function(x) {
#       md.list.samp[[x]][[i]][ss[[x]]]
#     }))
#   }
#   names(md.list.annot) = md.names
#   
#   cols = rep("green", Nsamp)
#   cols[unlist(md.list.annot[which(names(md.list.annot)=="tirosh.cell.type.new")])=="astro_new"] = "red"
#   cols[unlist(md.list.annot[which(names(md.list.annot)=="tirosh.cell.type.new")])=="undifferentiated"] = "blue"
#   png(paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#              "/pairwise_ica_components_ss_tirosh_",j,"_",cancer,".png"),
#       width = 1000, height = 1000)
#   pairs(t(data.ss), col = cols, pch = 20, cex = 0.5)
#   dev.off()
# 
#   cols = rep("cyan", Nsamp)
#   cols[unlist(md.list.annot[which(names(md.list.annot)=="cyc.status")])=="cycling"] = "red"
#   png(paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#              "/pairwise_ica_components_ss_cyc_status_",j,"_",cancer,".png"),
#       width = 1000, height = 1000)
#   pairs(t(data.ss), col = cols, pch = 20, cex = 0.5)
#   dev.off()
# 
#   pal = rainbow(length(ss))
#   cols = unlist(lapply(seq(length(ss)), function(x) {
#     rep(pal[x],length(ss[[x]]))
#   }))
#   #subsample again so we can see more clearly if samples are mixed
#   ss2 = sample(ncol(data.ss),ncol(data.ss)/5,replace=F)
#   data.ss2 = data.ss[,ss2]
#   cols2 = cols[ss2]
#   png(paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#              "/pairwise_ica_components_ss_samp_",j,"_",cancer,".png"),
#       width = 1000, height = 1000)
#   pairs(t(data.ss2), col = cols2, pch = 20, cex = 0.5)
#   dev.off()
# 
# 
#   #add row names for each dimension
#   rownames(data.ss) = colnames(df.mg.av)
# 
#   #make input file for STREAM
#   fwrite(as.data.frame(data.ss),
#          file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,"/",
#                        cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
#          sep = "\t", row.names = T, col.names = T, quote = F)
# 
# 
#   #make metadata file
#   df = data.frame(md.list.annot)
#   df$cell = colnames(data.ss)
#   df$sample = sapply(df$cell, function(x) {unlist(strsplit(x,"-"))[3]})
#   fwrite(df, file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#                                           "/metadata_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
#          sep = "\t", col.names=T, row.names=F, quote = F)
#   
#   #make separate metadata file with ICA scores
#   md.cnames = paste0("ICA_",1:nICAdims,"_score_med")
#   rel.ICA = lapply(seq(length(ss)), function(x) {
#     #ICA column names
#     inds = sapply(md.cnames, function(y) {which(colnames(oligo[[x]]@meta.data)==y)})
#     return(oligo[[x]]@meta.data[ss[[x]],inds])
#   })
#   fwrite(as.data.frame(do.call("rbind",rel.ICA)), file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#                                           "/metadata_ICA_",cancer,"_samps_ss_rem_cc_",j,".tsv.gz"),
#          sep = "\t", col.names=T, row.names=T, quote = F)
# 
# 
#   #write full data
#   data.full.ss = list()
#   for(i in 1:length(samps)){
# 
#     data.filt = oligo[[i]]@assays$SCT@data[sapply(rownames(oligo[[i]]@assays$SCT@data),
#                                                   function(x) {x %in% rownames(df.mg.av)}),]
#     data.full.ss[[i]] = data.filt[,ss[[i]]]
#   }
# 
#   data.full.ss = do.call(cbind, data.full.ss)
#   #detect variable genes
#   oligo.ss = CreateSeuratObject(data.full.ss)
#   oligo.ss = FindVariableFeatures(oligo.ss, selection.method = "vst", nfeatures = 5000)
#   to_print = data.full.ss[rownames(data.full.ss) %in% VariableFeatures(oligo.ss),]
# 
#   fwrite(as.data.frame(as.matrix(to_print)),
#          file = paste0(fig.dir,"/",cancer,"_ss_rem_cc_",j,
#                        "/",cancer,"_samps_ss_full_data_rem_cc_",j,".tsv.gz"),
#          sep = "\t", row.names = T, col.names = T, quote = F)
# }



# #do version with large number of cells
# Ncell.samp = sapply(oligo, function(x) {
#   length(x$orig.ident)
# })
# Nsamp = 65000
# N.per.samp = rmultinom(1,Nsamp,Ncell.samp)



# #if astro, use all cells, if oligo sample as above just more cells
# if(cancer=="astro"){
#   ss = lapply(seq(length(Ncell.samp)), function(x) {
#     seq(Ncell.samp[x])
#   })
# }else if(cancer=="oligo"){
#   ss = lapply(seq(length(Ncell.samp)), function(x) {
#     sort(sample(Ncell.samp[x], N.per.samp[x,], replace=F))
#   })
# }
