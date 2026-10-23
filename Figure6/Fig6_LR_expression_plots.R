library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Signac)
library(readxl)
library(tidyr)
library(ggpubr)
library(rstatix)
library(textshape)
library(compositions)
library(easyCODA)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(scCustomize)


### to change
setwd("/mnt/beegfs/home/alyne/stagein")
fig.dir = '/mnt/beegfs/home/alyne/stageout/figures_myeloid'
dir.create(fig.dir)

#read in full object
all.myeloid = readRDS("/mnt/beegfs/home/alyne/stagein/astro_oligo_data/astro_oligo_dataset_all.rds")


N.per.cat = 20
Ncell = 1000


#####
## First myeloid expression by subset
####


#LR expression to plot in myeloid cells
dot.genes.mac.final = c("SPP1","ITGB1","ITGAV","ITGA9","ITGB5","PDGFB","PDGFC","HLA-E","TLR4","ABCA1","CD74")


#combined dot plot with key LR
DotPlot(all.myeloid.filt, features = dot.genes.mac.final, group.by = "cluster_annotation",
        col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  scale_colour_viridis_c(option = "magma", direction=-1) + labs(y="",x="") + theme_bw() + scale_size(range = c(0, 6), limits = c(20, 90)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.text=element_text(size=12))
ggsave(paste0(fig.dir,"/dotplot_all_myeloid_cells_combined_scaled.pdf"), width = 10, height = 3.8)




######
## Then malignant cell expression by cell state
######


#load in Seurat objects with malignant cell states and gene expression

cancer = "oligo"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
setwd(data.dir)
cor.thresh = 0.5

#directory with sICA metaprograms
sica.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)

#set to "malig" if want to run on malignant cells
fig.str = "malig"
oligo = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,".rds"))

#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(sica.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.oligo = colnames(df.annot)



cancer = "astro"
data.dir = paste0("/mnt/beegfs/home/alyne/stagein/scenicplus_",cancer)
cor.thresh = 0.45

#directory with sICA metaprograms
sica.dir = paste0("/mnt/beegfs/home/alyne/stageout/figures_sICA_",cancer,"_",cor.thresh)
astro = readRDS(paste0(data.dir,"/oligo.ss_integrated_RPCA_",fig.str,"_",cancer,"_doub_rem.rds"))

#read in annotated sICA components
df.annot = as.data.frame(read_xlsx(paste0(sica.dir,"/top_gene_med_",cancer,"_annot.xlsx")))
colnames(df.annot) = df.annot[1,]
df.annot = df.annot[-1,]
prog.orig.order.astro = colnames(df.annot)


#find top metaprogram projection dimension for each cell
oligo$top_ICA_oligo = sapply(seq(length(colnames(oligo))), function(x) {prog.orig.order.oligo[which.max(oligo@assays$ICA@layers$data[,x])]})
astro$top_ICA_astro = sapply(seq(length(colnames(astro))), function(x) {prog.orig.order.astro[which.max(astro@assays$ICA@layers$data[,x])]})


#redefine cells as "Reactive/MES-like" if that's their highest dimension
oligo$ICA_mes = oligo$ICA_score_states
oligo$ICA_mes[which(oligo$top_ICA_oligo=="Reactive/MES-like")] = "Reactive/MES-like"

astro$ICA_mes = astro$ICA_score_states
astro$ICA_mes[which(astro$top_ICA_astro=="Reactive/MES-like")] = "Reactive/MES-like"


#define ligands/receptors to be plotted in malignant cells
dot.genes.malig = c("CD44","TNC","FN1","LAMA4","LTBP1","PDGFRA","KLRC2","ANGPTL2","RIMS2","MIF")


astro$ICA_mes = factor(astro$ICA_mes, levels = c("astro-like","Reactive/MES-like","OPC-like","neuro-like","intermediate"))
oligo$ICA_mes = factor(oligo$ICA_mes, levels = c("astro-like","Reactive/MES-like","OPC-like","intermediate"))

#set these to identities so they are used in DotPlot
Idents(astro) <- "ICA_mes"
Idents(oligo) <- "ICA_mes"


p1 = DotPlot(astro, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12))+#, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_viridis_c(option = "magma", direction=-1) + guides(colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
  labs(y = "",x = "",colour = "Average Expression",size = "Percent Expressed")
p2 = DotPlot(oligo, features = dot.genes.malig, group.by = "ICA_mes",
             col.min = -0.5, col.max = 1.5, dot.scale = 8, scale=T) +
  theme_bw() + scale_size(range = c(0, 6), limits = c(25, 100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_viridis_c(option = "magma", direction=-1) + theme(axis.text=element_text(size=12)) + 
  guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
  labs(y = "", x = "", colour = "Average Expression", size = "Percent Expressed")
p1 + p2 + plot_layout(ncol = 2)
ggsave(paste0(fig.dir,"/dotplot_all_malig_cells_mes_by_cancer_scaled.pdf"), width = 18, height = 3.8)









