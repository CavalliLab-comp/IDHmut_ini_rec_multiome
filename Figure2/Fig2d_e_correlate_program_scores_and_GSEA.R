

library(data.table)
library(gplots)
library(ggplot2)
library(dplyr)
library(Seurat)
library(Signac)
library(igraph)
library(RColorBrewer)
library(gage)
library(annotate)
library(clusterProfiler)
library(fgsea)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)


data.dir = "/data/users/alyne/gliomult/svn/analyse/results"
setwd(data.dir)


#read in group membership
df.mem.list = list()
df.annot.list = list()
df.mg.av.list = list()
count = 1
for(cancer in c("astro","oligo")){
  if(cancer=="astro"){
    cor.thresh = 0.45
  }else{
    cor.thresh = 0.5
  }
  fig.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
  fig.dir = paste0("figures_sICA_",cancer,"_",cor.thresh)
  df.mem.list[[count]] = read.table(paste0(fig.dir,"/sub_graph_members_split_",cancer,".tsv"),
                      sep = "\t", header = T)
  df.annot.list[[count]] = as.data.frame(read_xlsx(paste0(fig.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
  colnames(df.annot.list[[count]]) = df.annot.list[[count]][1,]
  df.annot.list[[count]] = df.annot.list[[count]][-1,]
  df.mg.av.list[[count]] = fread(paste0(fig.dir,"/all_ICA_components_med_",cancer,".tsv.gz"),
         sep = "\t", header = T, data.table=F)
  rownames(df.mg.av.list[[count]]) = df.mg.av.list[[count]]$V1
  df.mg.av.list[[count]]$V1 = NULL
  count = count + 1
}

prog.order.list = list()
prog.order.list[[1]] = c("Reactive/MES-like","Astro−like1", "OPC−like1",
                         "OPC−like2","OPC−like3","Gliogenesis","Neuro−like",
                         "Cell cycle","Ribosome proteins","Signalling","Stress","Chromatin remodelling",
                         "lnRNA","RNA processing","chr2","chr3","chr4","chr5")
prog.order.list[[2]] = c("Reactive/MES-like","Astro−like1","Astro−like2","RG/astro synapse","OPC/neuro-like","OPC−like1",
                         "OPC−like2","OPC−like3","Synapse/neuron projection","Neural/developmental",
                         "Cell cycle","Ribosome proteins","MYC targets","Splicing","RNA/proteostatic regulation",
                         "lnRNA","RNA processing","chr1","chr2","chr3","chr4","chr5","chr7")


#use intersection of genes
all_rows <- names(table(unlist(lapply(df.mg.av.list, rownames))))[which(table(unlist(lapply(df.mg.av.list, rownames)))==2)]
# Create expanded matrices filled with 0
astro_exp <- matrix(0, nrow=length(all_rows), ncol=ncol(df.mg.av.list[[1]]),
                    dimnames=list(all_rows, colnames(df.mg.av.list[[1]])))
oligo_exp <- matrix(0, nrow=length(all_rows), ncol=ncol(df.mg.av.list[[2]]),
                    dimnames=list(all_rows, colnames(df.mg.av.list[[2]])))

# Fill existing values
astro_exp <- as.matrix(df.mg.av.list[[1]])[match(all_rows,rownames(df.mg.av.list[[1]])),]
oligo_exp <- as.matrix(df.mg.av.list[[2]])[match(all_rows,rownames(df.mg.av.list[[2]])),]
merged <- cbind(astro_exp, oligo_exp)

#sort dashes
prog.order.list = lapply(prog.order.list, function(x) {gsub("\u2212", "-", x)})
prog.orig.order.list = lapply(1:2, function(x) {
  colnames(df.annot.list[[x]])
})
prog.orig.order.list[[1]] = paste0("A-",prog.orig.order.list[[1]])
prog.orig.order.list[[2]] = paste0("O-",prog.orig.order.list[[2]])

colnames(merged) = c(prog.orig.order.list[[1]],prog.orig.order.list[[2]])


M = cor(merged, method="pearson")
col_fun <- colorRamp2(c(-0.45, 0, 0.45), c("blue","white","red"))
ht = Heatmap(as.matrix(M),
             name = "Correlation",
             col = col_fun, row_split = 4, column_split = 4,
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2")
pdf(paste0(fig.dir,"/corr_ICA_components_astro_oligo_pearson_intersect.pdf"), width = 12, height = 11, useDingbats = FALSE)
draw(ht)
dev.off()





#make combined figure showing enrichment of Nomura (GBM) and Johnson (IDHmut) NMF components
merged = merged[,!grepl("-chr",colnames(merged))]

#do GSEA of Nomura GBM signatures in IC components
GBM.genesets = as.data.frame(read_xlsx("Nomura_GBM_genesets.xlsx"))
IDH.genesets = as.data.frame(read_xlsx("Johnson_IDHmut_genesets.xlsx"))
colnames(GBM.genesets) = paste0("GBM_",colnames(GBM.genesets))
colnames(IDH.genesets) = paste0("IDH_",colnames(IDH.genesets))
glioma.genesets = cbind(GBM.genesets,IDH.genesets)


gsea.output = apply(as.matrix(merged), 2, function(x) {
  gsea.input = x
  names(gsea.input) = rownames(merged)
  fgsea(glioma.genesets, gsea.input, nproc=1, nPermSimple = 10000, scoreType = "pos")
})



gsea.pval = do.call("rbind",lapply(gsea.output, function(x) (-1)*log10(x$padj)))
colnames(gsea.pval) = gsea.output[[1]]$pathway
col_fun <- colorRamp2(c(0, 10), c("white","blue"))
ht = Heatmap(as.matrix(gsea.pval), name = "GSEA -log10(adj pval)",
             col = col_fun, #row_split = 3, column_split = 3,
             clustering_distance_rows = "euclidean",
             clustering_distance_columns = "euclidean",
             clustering_method_rows = "ward.D2",
             clustering_method_columns = "ward.D2",
             rect_gp = grid::gpar(col = "black", lwd = 0.5),
             column_names_max_height = grid::unit(12, "cm"))
pdf(paste0(fig.dir,"/GSEA_Nomura_Johnson_padj.pdf"), width = 9, height = 8, useDingbats = FALSE)
draw(ht)
dev.off()

